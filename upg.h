#pragma once

// ============================================================================
// upg.h - multi-linked (combinatorial-map) loop, now with the topology surgery
// primitives that R1 / R2 / auto-tighten are built from.
//
// A crossing is ONE vertex bead with four slots; each slot points at the
// adjacent edge bead (neighborOfSlot) and records which other slot its strand
// continues to (throughSlot). All links are INDICES into LoopModel::beads.
//
// R1 / R2 need no face computation: a removable kink is a self-loop arc, and a
// removable bigon is two crossing-free arcs between the same two crossings -
// both read straight off the links. Removal is local relinking (CollapseRegion)
// plus a smooth rest-length collapse animation.
//
// This header declares the whole model; the implementation is split across:
//   upg_build.cpp    coordinate utils, grid topology map, BuildLoopFromPolygon
//   upg_surgery.cpp  link-rewiring primitives + dynamic resampling
//   upg_faces.cpp    arc / face queries and ClassifyFace (R1/R2/R3 detection)
//   upg_moves.cpp    the R1 / R2 / R3 moves
//   upg_render.cpp   raylib drawing
// (sigma_import.* and regions.* provide an alternative builder and the live
//  pin/region tracking; Sandbox_UPG.cpp is the physics + game-loop driver.)
// ============================================================================

#include <raylib.h>
#include <vector>
#include <string>
#include <cmath>
#include <set>
#include <map>
#include <random>
#include <functional>

// ============================================================================
// 1. CONSTANTS  (canvas-pixel units)
// ============================================================================

const int   SCREEN_WIDTH = 1200;
const int   SCREEN_HEIGHT = 1000;
// The board is GRID_N x GRID_N cells. This is the single knob for grid size:
// the render scale, the face-id machinery (GetFaceId / BuildTopologyMap), and
// the union-find size all derive their dimensions from it. The 3-cell slack
// keeps a 1.5-cell margin on each side of the knot, as the original 6/9 did.
const int   GRID_N = 6;
const float GRID_SIZE = SCREEN_WIDTH / (GRID_N + 3.0f);
const float KNOT_PIXEL_WIDTH = GRID_N * GRID_SIZE;
const float PADDING_X = (SCREEN_WIDTH - KNOT_PIXEL_WIDTH) / 2.0f;
const float PADDING_Y = (SCREEN_HEIGHT - KNOT_PIXEL_WIDTH) / 2.0f;

const float GEOM_EPS = 1e-6f;

// Physics
const float STEP_SIZE = 0.5f;
const float RADIUS = 0.25f * GRID_SIZE;
const float REST_LENGTH = 2.0f * RADIUS;
const float SPRING_CONSTANT = 1.0f;

// Keep arcs taut and evenly spaced (see GameState::Step).
//   K_SMOOTH: Laplacian smoothing force on each EDGE bead toward the midpoint of
//     its two neighbours. With the length springs this straightens folds and
//     distributes beads evenly; keep it well below SPRING_CONSTANT so it tidies
//     without collapsing arcs onto their chords.
//   EDGE_TENSION_SCALE: run the length-spring rest length a touch under
//     REST_LENGTH so the whole loop sits under mild tension and slack does not
//     accumulate. Directly-linked beads are exempt from the overlap projection, so
//     this doesn't fight it. 1.0 == no tension; lower == tighter (don't go so low
//     that short arcs are pulled below the overlap floor everywhere).
const float K_SMOOTH = 0.20f;
const float EDGE_TENSION_SCALE = 0.90f;

// Target even spacing for ResampleArcsEven: each arc is resampled (beads added OR
// removed) to round(arcLength / ARC_EVEN_SPACING) - 1 beads, so a too-sparse arc
// fills in and an over-packed one sheds beads, both landing at even spacing. Set
// to the tension equilibrium (EDGE_TENSION_SCALE * REST_LENGTH) so the resampled
// count matches what the springs want and the two don't fight. Larger = sparser.
const float ARC_EVEN_SPACING = REST_LENGTH * EDGE_TENSION_SCALE;

const float COLLISION_STIFFNESS = 1.0f;
const float INTEGRATION_STEP = 0.15f;
const float CENTER_SPEED = 0.15f;
const float PIN_BUFFER_RADIUS = REST_LENGTH * 1.5f;

const int   UPDATES_PER_FRAME = 10;

// Self-avoidance: two strands repel everywhere EXCEPT within this radius of a
// shared crossing, where they are allowed to coincide (that is the crossing).
const float CROSS_GUARD = REST_LENGTH * 1.0f;

// Dynamic resampling. Self-avoidance is bead-to-bead, so a wide gap between two
// consecutive edge beads is a hole another strand can slip through. Whenever an
// edge link stretches past this length we insert a fresh edge bead at its
// midpoint, keeping the strand densely enough sampled that collision holds.
const float EDGE_SPLIT_LENGTH = REST_LENGTH * 1.0f;

// Dynamic simplification (the inverse of the split above): an edge bead is
// redundant and gets spliced out when its strand doubles back too sharply at it,
// or when it sits in a pile of three mutually-touching edge beads. Only acts on
// consecutive strand beads, so strands merely crossing/piling are never fused.
const float EDGE_TOUCH_LENGTH = REST_LENGTH;   // disks overlap within this centre distance
const float MIN_ARC_ANGLE_DEG = 45.0f;         // interior angle sharper than this => drop the apex bead

// Minimum edge beads guaranteed on EVERY crossing-to-crossing arc. EnsureMinArcBeads
// tops up any arc below this each settled frame, so the diagram never develops a
// 0- or 1-bead arc between two crossings. A floor of 2 both keeps region tracing /
// rendering well-formed and (crucially) guarantees the R3 converging walk always has
// room: an empty or one-bead triangle arc used to collapse both walking T-junctions
// onto one bead and dereference beads[-1]. Kept at 2 (not higher) so it does not
// fight SimplifyDenseArcs -- a 2-bead arc has a crossing on one side of every bead,
// so neither the hairpin nor the pile rule can remove it, hence no add/remove churn.
const int   MIN_ARC_BEADS = 3;

// Absolute minimum an arc may be thinned to, even when its two crossings are so
// close that MIN_ARC_BEADS beads can't fit without overlap. The adaptive cap in
// ArcMinTarget lets a short arc drop to this rather than being padded past what it
// can hold; 2 keeps the arc well-formed for region tracing and the R3 walk.
const int   ARC_HARD_FLOOR = 2;

// Face-specific bead floors (each >= MIN_ARC_BEADS). A monogon (an R1 kink's self-
// loop) and each arc of a bigon (an R2 lens) are padded to these so the loop / lens
// reads as a rounded face and its collapse animates smoothly instead of snapping
// from a 1-2 bead sliver.
const int   MONOGON_MIN_BEADS = 6;
const int   BIGON_MIN_BEADS = 6;

// Overlap thinning: two CONSECUTIVE edge beads on the same arc closer than this are
// redundant and the excess is spliced out -- but never below the arc's own floor,
// so thinning can't fight the padding above. Well under REST_LENGTH, so normally
// spaced beads (tension holds them near 0.9 * REST_LENGTH) are kept; this only
// catches genuine piles/overlaps. Also sets how small a gap the fill pass can
// close: the fill threshold below is 2.2x this, so a lower value lets smaller gaps
// be filled (at the cost of gentler pile removal).
const float BEAD_OVERLAP_DIST = REST_LENGTH * 0.6f;

// Gap filling: when the length springs stretch a link between two edge beads past
// this, ResampleStretchedArcs inserts a bead at its midpoint (NormalizeArcBeads
// only pads to a bead-count floor, so it never fills a gap on a long, sparse arc).
// Kept above 2 * BEAD_OVERLAP_DIST (with a margin) so each midpoint half lands
// farther than the thinning distance and is not immediately removed again -- so
// filling and thinning can't oscillate. Lower BOTH this and BEAD_OVERLAP_DIST
// together if you want smaller gaps filled.
const float ARC_MAX_LINK = 2.2f * BEAD_OVERLAP_DIST;

// Hard non-overlap. The self-avoidance force above is a soft penalty that the
// displacement clamp can cap, so deep overlaps still leak through. After every
// integration substep a position-projection pass directly separates any pair of
// bead disks that overlap, guaranteeing none do -- except the cases it must
// leave alone: directly-linked beads (which rest exactly touching), a crossing
// with its own slot-neighbours, and the two strands that genuinely cross at a
// crossing (identified topologically, not by a radius, so a third strand merely
// passing nearby is still separated). Min centre separation is one bead diameter;
// a few Gauss-Seidel sweeps per substep clear even tight clusters.
const float BEAD_MIN_SEP = REST_LENGTH;   // = 2 * RADIUS: disks just touch, never overlap
const int   OVERLAP_RESOLVE_ITERS = 5;    // sweeps per substep; more clears tight clusters faster

// ---- overlap projection: Jacobi + under-relaxation (anti-jitter) -----------
// The old projection wrote each correction in place (Gauss-Seidel), so the fix a
// bead received depended on bead index order and it never fully settled -- a
// residual, order-biased push that buzzed against the length springs. The Jacobi
// path accumulates every pair's correction, then applies the per-bead AVERAGE
// scaled by OVERLAP_RELAX < 1, so squeezed beads glide apart symmetrically
// instead of oscillating. Set USE_JACOBI_OVERLAP=false for the legacy path.
const bool  USE_JACOBI_OVERLAP = true;
const float OVERLAP_RELAX = 0.5f;   // under-relaxation (0<w<=1); lower = smoother, needs more iters

// ---- sleep gating (anti-jitter) --------------------------------------------
// A settled loop never reaches an exact fixed point, so integrating it forever
// leaves it micro-vibrating at the float-noise floor. When bead AND pin motion
// stay under SLEEP_EPSILON px for SLEEP_FRAMES consecutive frames with nothing
// animating, freeze the simulation until an interaction wakes it. SLEEP_EPSILON
// is a per-frame max-displacement threshold; well below one pixel is invisible.
const bool  SLEEP_ENABLE = true;
const float SLEEP_EPSILON = 0.08f;
const int   SLEEP_FRAMES = 12;

// Pins settle by repulsion balance: pushed away from every nearby strand bead
// and crossing, a pin sits at the open centre of its face and cannot drift past
// a bounding strand (moving toward a wall only increases that wall's push).
const float PIN_SETTLE_RANGE = GRID_SIZE * 1.0f;   // how far a pin "sees" walls
const float PIN_SETTLE_STIFF = 0.45f;
const float PIN_PIN_RANGE = GRID_SIZE * 1.1f;   // pins separate from each other
const float PIN_PIN_STIFF = 0.6f;
const float PIN_MAX_STEP = RADIUS * 0.4f;      // smooth, jitter-free motion

// Smooth removal: rest length scales 1 -> 0, then the region is spliced out.
const float TRANSITION_COLLAPSE_RATE = 0.045f;
const int   COLLAPSE_FRAMES = 26;

// R3 keeps all three crossings, so it never collapses to a point (that pile
// relaxes chaotically and tangles the loop). Instead the triangle shrinks only
// to a small non-degenerate size, the legs swap there (a small jump because the
// crossings are close but distinct), then it regrows gently from that floor.
const float R3_COLLAPSE_FLOOR = 0.5f;   // smallest triangle scale during R3
const float R3_REGROW_RATE = 0.02f;   // gentler than collapse so physics keeps up
const int   R3_SHRINK_FRAMES = 54;
const int   R3_SLIDE_FRAMES = 58;       // guided flip into the clean layout
const int   R3_GROW_FRAMES = 54;
const int   MORPH_FRAMES = 140;      // R3 morph: crossings slide to the arc midpoints (slow, so the merge reads clearly)

const int   R3_WALK_GAP = 5;       // frames between T-junction hops during the R3 walk (shared by walk + merge timing)

// Before an R3 walk, each of the three triangle arcs is padded up to this many
// edge beads (spread on its longest links). More beads => more, smaller hops =>
// a smoother converging glide, and it guarantees room for two distinct T-junctions
// to walk toward the midpoint without landing on the same bead. These extra beads
// are transient: SimplifyDenseArcs prunes the surplus once the move settles. Must
// be >= 2; larger is smoother but slower (each extra bead adds ~R3_WALK_GAP frames).
const int   R3_MIN_ARC_BEADS = 4;

// ----- R3 walk outward-routing experiments (try each independently) -----
// The animated R3 T-junction walk rewires topology but leaves bead POSITIONS to
// the physics relaxer, which interpolates each migrating stem along a straight
// chord -- and that chord can cut across the triangle interior and the other two
// arcs, where the bead-to-bead self-avoidance is sparse enough to be threaded.
// Both options bias the walk to the OUTER side of the triangle (the half-plane
// away from the centroid captured in r3Center) so strands wrap the region:
//   Option 1 (R3_NUDGE_OUTWARD): nudge each relocated junction / merged crossing
//            and the stem it drags outward, on every hop and at every merge.
//   Option 2 (R3_PREBOW_STEM):   drop an outward guide bead on each stem the
//            instant it is detached, so it starts pre-bent around the outside.
// Set one false to isolate the other; both true combines them.
const bool  R3_NUDGE_OUTWARD = false;    // Option 1
const bool  R3_PREBOW_STEM = false;    // Option 2
const float R3_OUTWARD_NUDGE = 1.0f * REST_LENGTH;   // per-hop / per-merge offset (Option 1)
const float R3_STEM_BOW = 1.0f * REST_LENGTH;   // guide-bead outward bow   (Option 2)

//   Option 3 (R3_NORMAL_STEM): instead of pushing radially from the centroid,
//            drop the stem PERPENDICULAR to the local bar (the segment through the
//            junction's two bar neighbours), on the outward side. Because it tracks
//            the local arc rather than a fixed centroid ray, it stays correct on
//            curved / non-convex arcs and on skinny faces where the centroid sits
//            almost on the arc (where the radial ray drifts sideways instead of
//            away). At the merge the junction becomes a crossing, so "perpendicular"
//            degrades to "straight through": each transverse leg is seated opposite
//            its through-partner (a clean X) and the crossing is offset along the
//            outerA--outerB chord's normal. No-op outside an R3 walk.
// Option 3 composes with the others; it can fully replace Option 1's per-hop nudge.
const bool  R3_NORMAL_STEM = false;                 // Option 3
const float R3_STEM_NORMAL_LEN = 1.0f * REST_LENGTH;   // per-hop / per-merge normal offset (Option 3)
// If the one-shot placement above is tugged back by the springs before the next
// hop, ApplyStemNormalForces adds a standing per-substep angular spring holding the
// stem on the bar normal. Start soft so it biases without fighting the edge springs.
const bool  R3_NORMAL_STEM_SPRING = true;          // Option 3b: continuous angular hold (Option-1 leg-hold rides on this)
const float R3_STEM_NORMAL_K = 0.15f;              // angular-spring stiffness (Option 3b)

// Option 1: past the stem bead, also hold the first R3_STEM_LEG_SPAN beads of the
// outer LEG on the outward side of the junction's arc, so the whole leg -- not
// just the stem -- is kept out of the trigon. One-sided (only a leg bead that has
// drifted inward is pushed) with per-hop falloff, so it biases the leg outward
// without straightening it. KeepOutOfTrigon stays as a hard backstop.
const int   R3_STEM_LEG_SPAN = 5;                  // leg beads past the stem to hold (0 = stem only)
const float R3_STEM_LEG_MARGIN = RADIUS;           // clearance the leg keeps outside its arc line

// ----- Smooth-transition toggles for the manual (pin-click) R-moves -----
// R1 + R2: route the click through the SAME shrink-to-zero collapse the auto-
// tightener uses (arc rest lengths ramp 1->0 via the transitions map, then the
// region is spliced by CollapseRegion) instead of the bead-at-a-time eat (R1)
// or the single-frame surgery (R2). Set false for the legacy paths.
// R3: ease each T-junction hop -- the promoted junction starts at the OLD
// junction's position and springs forward one bead, so the branch point slides
// instead of teleporting. Independent of the outward-routing options above.
const bool  SMOOTH_R1R2_VIA_COLLAPSE = true;
const bool  R3_EASE_HOPS = false;

// R1 monogon collapse eats the self-loop one edge bead at a time; this is the
// frame gap between successive bead removals (the loop visibly retracts into
// the crossing, then the crossing is welded out on the final tick).
const int   R1_COLLAPSE_GAP = 5;

// Auto-tighten
const float TENSION_SETTLED = 0.6f;

// Stability / centering: a gentle force pulls the loop's centroid to the canvas
// centre (pure translation, no shape distortion), and a clamp caps how far any
// bead can move in one substep so stiff collisions can't explode.
const float CENTERING_STRENGTH = 0.05f;
const float MAX_DISPLACEMENT = RADIUS * 0.3f;

inline bool IsGridAligned(float v) {
    return std::fabs(v - std::round(v)) < 1e-4f;
}

// ============================================================================
// 2. CORE TYPES
// ============================================================================

class UnionFind {
public:
    std::vector<int> p;
    explicit UnionFind(int n) { p.resize(n); for (int i = 0; i < n; i++) p[i] = i; }
    int  find(int i) { return p[i] == i ? i : (p[i] = find(p[i])); }
    void unite(int i, int j) { p[find(i)] = p[find(j)]; }
};

// A smooth rest-length transition: scale multiplies REST_LENGTH so the physics
// animates a collapse (rate < 0) or growth (rate > 0).
struct BeadTransition {
    float scale = 1.0f;
    float rate = 0.0f;
    float minScale = 0.0f;   // collapse floor; R1/R2 use 0 (full collapse), R3 a small value
};

struct Bead {
    int     id = -1;
    Vector2 pos = { 0.0f, 0.0f };
    Vector2 force = { 0.0f, 0.0f };
    bool    isVertex = false;
    bool    active = true;
    bool    driven = false;

    int edge = -1;
    int neighbors[2] = { -1, -1 };

    int crossing = -1;
    int neighborOfSlot[4] = { -1, -1, -1, -1 };
    int throughSlot[4] = { -1, -1, -1, -1 };
};

struct Pin {
    int     id = -1;
    Vector2 pos = { 0.0f, 0.0f };
    bool    active = true;
    int     cell = -1;
};

// Result of walking a strand from a crossing slot to the next crossing.
struct ArcWalk {
    int              endVertex = -1;
    std::vector<int> beads;
};

// An arc leaving a crossing at a known slot.
struct SlotArc {
    int slot;
    int endVertex;
    std::vector<int> beads;
};

// A face classified for a click: type 1/2/3 = R1/R2/R3 (0 = none), the
// bordering crossings, the beads to collapse (R1/R2), and the three bounding
// arcs (R3), so the move uses the exact arcs of the clicked face.
struct FaceClass {
    int type = 0;
    std::vector<int> crossings;
    std::vector<int> doomed;
    std::vector<std::vector<int>> triArcs;   // R3: {arc(A,B), arc(B,C), arc(C,A)}
};

// One outer-leg swap of an R3 move: detach bead bu from crossing U and bead bv
// from crossing V, then cross-attach (bu->V, bv->U).
struct R3Swap { int U, bu, V, bv; };

// A planned R3: the three triangle arcs to shrink/regrow, the two leg swaps to
// apply at the triple point, and the triangle centre.
struct R3Plan {
    bool             ok = false;
    std::vector<int> arcBeads;
    std::vector<R3Swap> swaps;
    std::vector<int> cross;                  // [X, Y, Z]
    std::vector<std::vector<int>> arcs;      // [ab(X->Y), bc(Y->Z), ca(Z->X)]
    Vector2          center = { 0.0f, 0.0f };
};

struct LoopModel {
    std::vector<Bead> beads;
    std::vector<int>  crossingBead;
    std::vector<Pin>  pins;

    std::set<std::string> walls;
    UnionFind             uf{ GRID_N * GRID_N + 1 };

    std::map<int, std::vector<int>> regionBeads;
    std::map<int, std::set<int>>    regionCrossings;

    std::map<int, BeadTransition>   transitions;

    std::mt19937 rng{ std::random_device{}() };

    // Transient state for the animated R3 T-junction walk: the triangle centroid
    // captured at move start, and a flag true only while the walk's scheduled
    // callbacks are still in flight. The outward-routing options above read these
    // to bow each migrating junction / stem around the OUTSIDE of the triangle.
    // Inert (no geometric change) whenever r3Active is false.
    Vector2 r3Center{ 0.0f, 0.0f };
    bool    r3Active = false;
    // Boundary of the trigon face being flipped: corners + arc beads, in cycle
    // order (V0, arc0, V1, arc1, V2, arc2), captured when the walk starts.
    // KeepOutOfTrigon uses it to hold the outer strands out of the region.
    std::vector<int> r3Ring;
    // Debug: beads ejected by the most recent KeepOutOfTrigon pass, for the F3
    // overlay. Refilled every substep while r3Active; read only by DrawR3Debug.
    std::vector<int> r3Ejected;

    int Through(int vertexBead, int fromBead) const {
        const Bead& v = beads[vertexBead];
        for (int i = 0; i < 4; ++i)
            if (v.neighborOfSlot[i] == fromBead)
                return v.neighborOfSlot[v.throughSlot[i]];
        return -1;
    }
    bool IsVertexIdx(int i) const { return i >= 0 && beads[i].isVertex; }

    // topology surgery (implemented in upg_surgery.cpp)
    ArcWalk WalkFromSlot(int crossingBeadIdx, int slot) const;
    std::vector<int> FindSelfLoop(int crossingBeadIdx) const;
    std::vector<ArcWalk> ArcsBetween(int A, int B) const;
    void ReplaceLink(int x, int oldN, int newN);
    void CollapseRegion(const std::vector<int>& doomed);

    // dynamic resampling (implemented in upg_surgery.cpp): keep edge beads dense enough
    // that the bead-to-bead self-avoidance can't be threaded through a gap.
    int  InsertEdgeBead(int a, int b);        // new bead at the midpoint of link a-b
    int  SplitLink(int x, int y);             // new edge bead on link x-y (y may be a crossing)
    int  SpliceBeadBetween(int u, int v);     // new edge bead on link u-v (EITHER may be a crossing)
    int  EnsureMinArcBeads(int minBeads);     // top every crossing-to-crossing arc up to minBeads; returns beads added
    int  ArcMinTarget(int crossingBeadIdx, int slot, int otherCrossing) const; // face-aware bead floor for one arc
    int  NormalizeArcBeads(float overlapDist);// per arc: pad up to the face floor, then thin overlaps; returns net edits
    int  ResampleStretchedArcs(float maxLen); // split every over-long edge link; returns count
    int  ResampleArcsEven(float spacing);     // per arc: add/remove to an even, length-based bead count; returns net edits
    int  RemoveEdgeBead(int b);               // splice an edge bead out, relinking its two neighbors
    int  SimplifyDenseArcs(float touchLen, float minAngleDeg); // drop kinked/piled beads; returns count

    // Crossing-overlap guard support. For each active crossing, the edge beads on
    // its four incident arcs that lie within CROSS_GUARD of it, tagged by which of
    // the crossing's two strands (0/1) they belong to. Two beads may legitimately
    // overlap only where some crossing has them on its two DIFFERENT strands --
    // i.e. they are the two strands actually crossing there.
    struct CrossingStrands { int crossBead; std::map<int, int> strandOf; };
    std::vector<CrossingStrands> BuildCrossingStrandMap() const;

    // arc / face queries (implemented in upg_faces.cpp)
    std::vector<SlotArc> ArcsFrom(int crossingBeadIdx) const;
    int  SlotOf(int vertexBead, int bead) const;
    bool ValidBigon(int A, int B, const SlotArc& i, const SlotArc& j) const;
    FaceClass ClassifyFace(Vector2 p) const;   // which R-move the point's face wants

    void Monogon_Collapse_R1(int V, const FaceClass& face,
        std::function<void(int, std::function<void()>)> schedule = nullptr);
    void DeleteCrossingAndInsertEdgeBeads(int c, const std::set<int>& lensBeads);
    void Bigon_Swap_R2(int V1, int V2, const FaceClass& face,
        std::function<void(int, std::function<void()>)> schedule);


    // R3 (triangle flip): plan the swaps from the clicked face, apply them.
    R3Plan PlanR3(const FaceClass& face) const;
    void   ApplyR3(const std::vector<R3Swap>& swaps);

    // R3 Micro-Surgery Tools
    void RoleSwapCrossingAndEdge(int V, int E);
    void SlideOuterBeadIntoArc(int crossingVertexIdx, int slotTowardsArc);
    void SlideCrossingOverBead(int c, int slot);
    int  FindSlotTowardsCrossing(int crossingA, int crossingB) const;
    void StartR3MicroSurgery(int c0, int c1, int c2);
    void CreateTConnection(int V, int O, int A);
    void TraverseTConnection(int T, int A);
    void MergeTConnections(int T1, int T2);   // fuse two adjacent T-junctions back into one crossing
    void ApplyStemNormalForces(Vector2 center, float k);  // Option 3b: hold each T-junction stem on the bar normal (per substep)
    void KeepOutOfTrigon(float margin);                   // push non-boundary beads out of the R3 trigon face (per substep)
    void RebuildCrossingList();               // crossingBead = every active vertex (call after R3 surgery)
    void Perform_R3_Walk_S3(int V0, int V1, int V2, FaceClass& f, std::function<void(int, std::function<void()>)> schedule);
    void Perform_R3_Walk_S1(int V, int S1, int S2, std::vector<int>& arc_b1, std::vector<int>& arc_b2,
        std::function<void(int, std::function<void()>)> schedule = nullptr);

    // NOT a Reidemeister move. The symmetric "split / contract / expand" morph:
    // swaps the outer legs on ALL THREE triangle arcs (R3 swaps only two), which
    // turns the three monogons into three bigons -- i.e. morphs the diagram into
    // the trefoil shadow. Kept for the animation; it changes the knot type.
    R3Plan PlanFakeR3(const FaceClass& face) const;
};



// ============================================================================
// 3. BUILD / UTILITY / RENDER  (build in upg_build.cpp, render in upg_render.cpp)
// ============================================================================

Vector2 GridCoordToCanvasCoord(float x, float y);
int     GetFaceId(float cx, float cy, UnionFind& uf);
void    BuildTopologyMap(const std::vector<Vector2>& gridPoints,
    std::set<std::string>& walls, UnionFind& uf);

LoopModel BuildLoopFromPolygon(const std::vector<Vector2>& gridPoints,
    const std::vector<Vector2>& faceCoords,
    const std::vector<std::pair<int, int>>& pinCells);

void Draw_Debug_Screen(const LoopModel& m);
void DrawR3Debug(const LoopModel& m);   // F3 overlay: trigon outline, centroid, ejected beads
void DrawWireSpline(const LoopModel& model);
void DrawPins(const LoopModel& model);
void DrawInterface(float tension, bool autoRunning);