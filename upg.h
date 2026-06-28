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
const float SPRING_CONSTANT = 2.0f;
const float COLLISION_STIFFNESS = 1.5f;
const float INTEGRATION_STEP = 0.15f;
const float CENTER_SPEED = 0.15f;
const float PIN_BUFFER_RADIUS = REST_LENGTH * 1.5f;

const int   UPDATES_PER_FRAME = 5;

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
const int   OVERLAP_RESOLVE_ITERS = 1;

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

const int   R3_WALK_GAP = 20;       // frames between T-junction hops during the R3 walk (shared by walk + merge timing)

// R1 monogon collapse eats the self-loop one edge bead at a time; this is the
// frame gap between successive bead removals (the loop visibly retracts into
// the crossing, then the crossing is welded out on the final tick).
const int   R1_COLLAPSE_GAP = 20;

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

    int Through(int vertexBead, int fromBead) const {
        const Bead& v = beads[vertexBead];
        for (int i = 0; i < 4; ++i)
            if (v.neighborOfSlot[i] == fromBead)
                return v.neighborOfSlot[v.throughSlot[i]];
        return -1;
    }
    bool IsVertexIdx(int i) const { return i >= 0 && beads[i].isVertex; }

    // topology surgery (implemented in upg.cpp)
    ArcWalk WalkFromSlot(int crossingBeadIdx, int slot) const;
    std::vector<int> FindSelfLoop(int crossingBeadIdx) const;
    std::vector<ArcWalk> ArcsBetween(int A, int B) const;
    void ReplaceLink(int x, int oldN, int newN);
    void CollapseRegion(const std::vector<int>& doomed);

    // dynamic resampling (implemented in upg.cpp): keep edge beads dense enough
    // that the bead-to-bead self-avoidance can't be threaded through a gap.
    int  InsertEdgeBead(int a, int b);        // new bead at the midpoint of link a-b
    int  SplitLink(int x, int y);             // new edge bead on link x-y (y may be a crossing)
    int  ResampleStretchedArcs(float maxLen); // split every over-long edge link; returns count
    int  RemoveEdgeBead(int b);               // splice an edge bead out, relinking its two neighbors
    int  SimplifyDenseArcs(float touchLen, float minAngleDeg); // drop kinked/piled beads; returns count

    // Crossing-overlap guard support. For each active crossing, the edge beads on
    // its four incident arcs that lie within CROSS_GUARD of it, tagged by which of
    // the crossing's two strands (0/1) they belong to. Two beads may legitimately
    // overlap only where some crossing has them on its two DIFFERENT strands --
    // i.e. they are the two strands actually crossing there.
    struct CrossingStrands { int crossBead; std::map<int, int> strandOf; };
    std::vector<CrossingStrands> BuildCrossingStrandMap() const;

    // arc / face queries (implemented in upg.cpp)
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
// 3. BUILD / UTILITY / RENDER  (implemented in upg.cpp)
// ============================================================================

Vector2 GridCoordToCanvasCoord(float x, float y);
int     GetFaceId(float cx, float cy, UnionFind& uf);
void    BuildTopologyMap(const std::vector<Vector2>& gridPoints,
    std::set<std::string>& walls, UnionFind& uf);

LoopModel BuildLoopFromPolygon(const std::vector<Vector2>& gridPoints,
    const std::vector<Vector2>& faceCoords,
    const std::vector<std::pair<int, int>>& pinCells);

void Draw_Debug_Screen(const LoopModel& m);
void DrawWireSpline(const LoopModel& model);
void DrawPins(const LoopModel& model);
void DrawInterface(float tension, bool autoRunning);