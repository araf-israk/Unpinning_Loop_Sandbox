// ============================================================================
// upg_moves.cpp - the Reidemeister moves themselves.
//
// Each move is planned from a ClassifyFace result (see upg_faces.cpp) and run as
// local link surgery (see upg_surgery.cpp), optionally spread over frames via a
// `schedule` callback so it animates.
//
//   R1  Monogon_Collapse_R1
//          Eat a self-loop kink one edge bead at a time, then weld the
//          crossing's two outer legs straight through.
//
//   R2  DeleteCrossingAndInsertEdgeBeads / Bigon_Swap_R2
//          Open a bigon by cross-joining each lens leg to the OTHER strand's
//          outer leg - the smoothing that opens the lens while keeping the loop
//          connected (joining the two lens legs to each other would pinch off a
//          stray component; that was the original R2 bug).
//
//   R3  PlanR3 / ApplyR3
//          The pure leg-swap form of the triangle flip: swap the outer legs of
//          the two arcs meeting the stationary crossing.
//       Perform_R3_Walk_S3 / Perform_R3_Walk_S1, CreateTConnection,
//       TraverseTConnection, MergeTConnections, RebuildCrossingList
//          The animated "T-junction walk" form actually used on the live click
//          path: split each corner crossing into two T-junctions, walk them to
//          the arc midpoints, and merge each kissing pair into a new crossing.
//       PlanFakeR3
//          NOT a Reidemeister move - the symmetric 3-arc swap that morphs the
//          three monogons into the trefoil shadow (changes the knot type; kept
//          only for the animation).
//
// EXPERIMENTAL / NOT WIRED IN: RoleSwapCrossingAndEdge, SlideOuterBeadIntoArc,
// SlideCrossingOverBead and StartR3MicroSurgery are leftover R3 micro-surgery
// scaffolding. DispatchAtPin uses Perform_R3_Walk_S3, not these. In particular
// RoleSwapCrossingAndEdge still has hard-coded bead indices (3, 29, 30, 31) from
// a debugging session and is NOT correct in general; it is preserved verbatim to
// avoid changing behaviour and must be rewritten before being used.
// ============================================================================

#include "upg.h"

#include <algorithm>
#include <cstdio>

namespace {
    // Seat bead `b` so it leaves `junction` PERPENDICULAR to the bar (the segment
    // barA--barB), on the side away from `center`, at distance `delta`. Tracking the
    // local bar orientation (rather than a fixed centroid ray) keeps the stem off the
    // arc even where the arc curves or the centroid sits nearly on it. Used by the
    // walk and the merge (Option 3, R3_NORMAL_STEM).
    void PlaceBeadNormalToBar(Bead& b, Vector2 junction, Vector2 barA, Vector2 barB,
        Vector2 center, float delta) {
        float tx = barB.x - barA.x, ty = barB.y - barA.y;   // bar tangent
        float L = std::hypot(tx, ty);
        if (L <= GEOM_EPS) return;
        float nx = -ty / L, ny = tx / L;                    // unit normal to the bar
        float ox = junction.x - center.x, oy = junction.y - center.y;
        if (nx * ox + ny * oy < 0.0f) { nx = -nx; ny = -ny; }  // pick the outward side
        b.pos.x = junction.x + nx * delta;
        b.pos.y = junction.y + ny * delta;
    }

    // Even-odd point-in-polygon test.
    bool PtInPolyLocal(Vector2 p, const std::vector<Vector2>& poly) {
        bool in = false; size_t n = poly.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++)
            if (((poly[i].y > p.y) != (poly[j].y > p.y)) &&
                (p.x < (poly[j].x - poly[i].x) * (p.y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x))
                in = !in;
        return in;
    }

    // Closest point on segment a-b to p; writes the squared distance to d2out.
    Vector2 ClosestOnSeg(Vector2 p, Vector2 a, Vector2 b, float& d2out) {
        float vx = b.x - a.x, vy = b.y - a.y, wx = p.x - a.x, wy = p.y - a.y;
        float c1 = vx * wx + vy * wy, c2 = vx * vx + vy * vy;
        float t = (c2 > GEOM_EPS) ? c1 / c2 : 0.0f;
        if (t < 0.0f) t = 0.0f; if (t > 1.0f) t = 1.0f;
        Vector2 q = { a.x + t * vx, a.y + t * vy };
        float dx = p.x - q.x, dy = p.y - q.y; d2out = dx * dx + dy * dy; return q;
    }
} // namespace



// R1 monogon removal, animated. The clicked face is a self-loop hanging off one
// crossing V: V leaves at slot sIn, runs through the loop's edge beads, and
// returns to V at slot sOut. sIn and sOut are adjacent slots (the loop sits in
// one corner), so they belong to V's two DIFFERENT strands; the through-partner
// of each is therefore one of V's two outside legs. Removing the kink deletes V
// and the loop and welds those two outside legs straight through.
//
// Rather than collapse the whole region at once, we eat the loop one bead at a
// time: every R1_COLLAPSE_GAP frames the front edge bead is spliced out and V's
// slot advances onto the next, so the loop visibly retracts into the crossing.
// The final tick removes V together with the last loop bead and welds the legs.
// With no scheduler the whole collapse runs instantly.
void LoopModel::Monogon_Collapse_R1(int V, const FaceClass& face,
    std::function<void(int, std::function<void()>)> schedule) {
    if (face.type != 1) return;
    if (V < 0 || !beads[V].active || !beads[V].isVertex) return;

    // Use the exact loop the classifier identified for the clicked face. For R1
    // it stored the ordered loop beads followed by the crossing, so strip V back
    // out. (A crossing can carry two self-loops -- a figure-eight has a monogon
    // on each side -- so re-searching with FindSelfLoop could pick the wrong one;
    // fall back to it only when the face carries no bead list.)
    std::vector<int> loop;
    for (int b : face.doomed) if (b != V) loop.push_back(b);
    if (loop.empty()) loop = FindSelfLoop(V);
    if (loop.empty()) return;

    // The two slots the loop occupies on V (kept distinct so a single-bead loop,
    // whose one bead links V twice, still resolves to two slots).
    int sIn = -1, sOut = -1;
    for (int s = 0; s < 4; ++s)
        if (beads[V].neighborOfSlot[s] == loop.front()) { sIn = s; break; }
    for (int s = 0; s < 4; ++s)
        if (s != sIn && beads[V].neighborOfSlot[s] == loop.back()) { sOut = s; break; }
    if (sIn < 0 || sOut < 0) return;

    // The two outside legs are the through-partners of the loop slots. Capture
    // them now: the eat only ever touches slot sIn, so these stay valid until the
    // final weld. (Defensive: if the loop slots were through-partners the diagram
    // is malformed for R1, so bail rather than weld garbage.)
    int outSlotIn = beads[V].throughSlot[sIn];
    int outSlotOut = beads[V].throughSlot[sOut];
    if (outSlotIn == sOut || outSlotOut == sIn) return;
    int outsideIn = beads[V].neighborOfSlot[outSlotIn];
    int outsideOut = beads[V].neighborOfSlot[outSlotOut];

    // The terminal surgery: retire V and the last loop bead, weld the legs, and
    // refresh the crossing list. Captured by value so it survives as a callback.
    auto finish = [this, V, outsideIn, outsideOut](int last) {
        if (outsideIn >= 0) ReplaceLink(outsideIn, V, outsideOut);
        if (outsideOut >= 0) ReplaceLink(outsideOut, V, outsideIn);
        beads[V].active = false;
        beads[V].isVertex = false;
        beads[V].crossing = -1;
        for (int s = 0; s < 4; ++s) { beads[V].neighborOfSlot[s] = -1; beads[V].throughSlot[s] = -1; }
        if (last >= 0) {
            beads[last].active = false;
            beads[last].neighbors[0] = beads[last].neighbors[1] = -1;
        }
        RebuildCrossingList();
        };

    // Splice the front loop bead out, joining its two live neighbours directly.
    // Unlike RemoveEdgeBead this carries no parallel-edge guard: during the eat V
    // stays linked to the loop's tail via slot sOut, so the general guard would
    // wrongly refuse the last beads. Here each spliced bead's neighbours are known
    // to be exactly V (advancing on slot sIn) and the next loop bead.
    auto splice = [this](int b) {
        int n0 = beads[b].neighbors[0], n1 = beads[b].neighbors[1];
        if (n0 >= 0) ReplaceLink(n0, b, n1);
        if (n1 >= 0) ReplaceLink(n1, b, n0);
        beads[b].active = false;
        beads[b].neighbors[0] = beads[b].neighbors[1] = -1;
        };

    const int k = (int)loop.size();
    int last = loop.back();

    if (!schedule) {
        // Instant: splice every bead but the last, then weld V + the last out.
        for (int i = 0; i + 1 < k; ++i) splice(loop[i]);
        finish(last);
        return;
    }

    // Eat one edge bead every R1_COLLAPSE_GAP frames. When splice(loop[i]) fires,
    // V's slot sIn has advanced onto loop[i], so its neighbours are V and loop[i+1].
    for (int i = 0; i + 1 < k; ++i) {
        int b = loop[i];
        schedule((i + 1) * R1_COLLAPSE_GAP, [splice, b] { splice(b); });
    }
    // Final tick: V now reaches `last` from both loop slots; weld it out.
    schedule(k * R1_COLLAPSE_GAP, [finish, last] { finish(last); });
}

// Smooth crossing `c` out of the clicked bigon. `lensBeads` is the set of edge
// beads on the two arcs bounding that bigon, so the two slots of `c` whose
// adjacent bead is in the set are the lens (inner) legs; the other two are the
// external legs, and they are exactly the through-partners of the lens legs
// (ValidBigon guarantees the two lens arcs are different strands here).
//
// Two planar smoothings exist. Joining the two lens legs to EACH OTHER pinches
// the lens off into a stray closed loop (R2 then changes the component count --
// the original bug). Joining each lens leg to its own strand's external leg
// leaves the two new arcs geometrically crossing. The correct R2 smoothing is
// the third option: cross-join each lens leg to the OTHER strand's external leg,
// which opens the lens while keeping the loop connected.
void LoopModel::DeleteCrossingAndInsertEdgeBeads(int c, const std::set<int>& lensBeads) {
    if (c < 0 || c >= (int)beads.size() || !beads[c].isVertex || !beads[c].active) return;

    // The two lens legs: slots whose incident arc bead belongs to this bigon.
    int La = -1, Lb = -1;
    for (int s = 0; s < 4; ++s) {
        int n = beads[c].neighborOfSlot[s];
        if (n >= 0 && lensBeads.count(n)) { (La < 0 ? La : Lb) = s; }
    }
    if (La < 0 || Lb < 0) return;                  // not a clean bigon corner
    if (beads[c].throughSlot[La] == Lb) return;    // lens arcs share a strand: not R2

    // External legs = through-partners of the lens legs.
    int Ea = beads[c].throughSlot[La];
    int Eb = beads[c].throughSlot[Lb];
    int legLa = beads[c].neighborOfSlot[La];
    int legLb = beads[c].neighborOfSlot[Lb];
    int legEa = beads[c].neighborOfSlot[Ea];
    int legEb = beads[c].neighborOfSlot[Eb];

    Vector2 cPos = beads[c].pos;
    beads[c].active = false;

    // Cross-join: lens leg La -> strand-B external (Eb); lens leg Lb -> strand-A
    // external (Ea). New edge beads sit at the old crossing position.
    auto weld = [&](int p, int q) {
        if (p < 0 || q < 0) return;
        Bead eb;
        eb.id = (int)beads.size();
        eb.isVertex = false;
        eb.active = true;
        eb.pos = cPos;
        eb.edge = (!beads[p].isVertex) ? beads[p].edge : beads[q].edge;
        eb.neighbors[0] = p;
        eb.neighbors[1] = q;
        int id = eb.id;
        beads.push_back(eb);
        ReplaceLink(p, c, id);
        ReplaceLink(q, c, id);
        };
    weld(legLa, legEb);
    weld(legLb, legEa);

    // Refresh the list of active crossing indices used by physics self-avoidance.
    RebuildCrossingList();
}

void LoopModel::Bigon_Swap_R2(int V1, int V2, const FaceClass& face,
    std::function<void(int, std::function<void()>)> schedule) {
    if (V1 < 0 || V1 >= (int)beads.size() || !beads[V1].isVertex || !beads[V1].active) return;
    if (V2 < 0 || V2 >= (int)beads.size() || !beads[V2].isVertex || !beads[V2].active) return;

    // The clicked face's two bounding arcs name exactly which legs are the lens,
    // so each crossing smooths toward the correct (component-preserving) side.
    std::set<int> lensBeads(face.doomed.begin(), face.doomed.end());

    DeleteCrossingAndInsertEdgeBeads(V1, lensBeads);
    DeleteCrossingAndInsertEdgeBeads(V2, lensBeads);
}


// ----------------------------------------------------------------------------
// R3 (triangle flip)
// ----------------------------------------------------------------------------

// Plan an R3 from the clicked triangle face. crossings[0] is taken as the
// stationary crossing X; the moving strand (Y-Z) slides across it. The two
// strands through X (to Y and to Z) flip order with X, which reduces to swapping
// outer legs: for arc X-Y swap X's and Y's outer legs, for arc X-Z swap X's and
// Z's. The third strand (Y-Z) is untouched. Using the exact arcs the classifier
// found (not an endpoint re-search) keeps this unambiguous after earlier moves.
R3Plan LoopModel::PlanR3(const FaceClass& face) const {
    R3Plan plan;
    if (face.type != 3 || face.crossings.size() != 3 || face.triArcs.size() != 3)
        return plan;

    int X = face.crossings[0];          // stationary crossing
    int Y = face.crossings[1];
    int Z = face.crossings[2];
    const std::vector<int>& ab = face.triArcs[0];   // X -> Y  (front~X, back~Y)
    const std::vector<int>& bc = face.triArcs[1];   // Y -> Z  (moving strand)
    const std::vector<int>& ca = face.triArcs[2];   // Z -> X  (front~Z, back~X)
    if (ab.empty() || bc.empty() || ca.empty()) return plan;

    // Outer leg of the strand that leaves crossing `c` through triangle-bead
    // `adj`: the through-partner slot's neighbour, i.e. the same strand
    // continuing out of the triangle.
    auto outer = [&](int c, int adj) -> int {
        int s = SlotOf(c, adj);
        if (s < 0) return -1;
        return beads[c].neighborOfSlot[beads[c].throughSlot[s]];
        };

    int bXY = outer(X, ab.front());     // X's outer leg on strand X-Y
    int bYX = outer(Y, ab.back());      // Y's outer leg on strand X-Y
    int bXZ = outer(X, ca.back());      // X's outer leg on strand X-Z
    int bZX = outer(Z, ca.front());     // Z's outer leg on strand X-Z

    if (bXY < 0 || bYX < 0 || bXZ < 0 || bZX < 0) return plan;
    if (bXY == bXZ || bYX == bZX || bXY == bYX || bXZ == bZX) return plan;

    plan.arcBeads = ab;
    plan.arcBeads.insert(plan.arcBeads.end(), bc.begin(), bc.end());
    plan.arcBeads.insert(plan.arcBeads.end(), ca.begin(), ca.end());
    plan.swaps = { { X, bXY, Y, bYX }, { X, bXZ, Z, bZX } };
    plan.cross = { X, Y, Z };
    plan.arcs = { ab, bc, ca };
    plan.center = { (beads[X].pos.x + beads[Y].pos.x + beads[Z].pos.x) / 3.0f,
                    (beads[X].pos.y + beads[Y].pos.y + beads[Z].pos.y) / 3.0f };
    plan.ok = true;
    return plan;
}

// Cross-attach the swapped outer legs. Self-inverse for the same swap list.
void LoopModel::ApplyR3(const std::vector<R3Swap>& swaps) {
    for (auto& s : swaps) {
        ReplaceLink(s.U, s.bu, s.bv);
        ReplaceLink(s.V, s.bv, s.bu);
        ReplaceLink(s.bu, s.U, s.V);
        ReplaceLink(s.bv, s.V, s.U);
    }
}

// The "split / contract / expand" morph from the sketch. This is NOT a valid
// Reidemeister move: where a real R3 swaps the outer legs of only the two arcs
// meeting one stationary crossing, this swaps the outer legs of ALL THREE arcs.
// That symmetric swap turns the three monogons into three bigons -- it morphs
// the unknot-with-three-kinks into the trefoil shadow, changing the knot type.
// It reuses the same R3Plan/ApplyR3 plumbing so the smooth animation applies.
R3Plan LoopModel::PlanFakeR3(const FaceClass& face) const {
    R3Plan plan;
    if (face.type != 3 || face.crossings.size() != 3 || face.triArcs.size() != 3)
        return plan;

    int X = face.crossings[0], Y = face.crossings[1], Z = face.crossings[2];
    const std::vector<int>& ab = face.triArcs[0];   // X -> Y
    const std::vector<int>& bc = face.triArcs[1];   // Y -> Z
    const std::vector<int>& ca = face.triArcs[2];   // Z -> X
    if (ab.empty() || bc.empty() || ca.empty()) return plan;

    auto outer = [&](int c, int adj) -> int {
        int s = SlotOf(c, adj);
        if (s < 0) return -1;
        return beads[c].neighborOfSlot[beads[c].throughSlot[s]];
        };

    int lXab = outer(X, ab.front()), lYab = outer(Y, ab.back());   // arc X-Y legs
    int lYbc = outer(Y, bc.front()), lZbc = outer(Z, bc.back());   // arc Y-Z legs
    int lZca = outer(Z, ca.front()), lXca = outer(X, ca.back());   // arc Z-X legs
    if (lXab < 0 || lYab < 0 || lYbc < 0 || lZbc < 0 || lZca < 0 || lXca < 0) return plan;

    plan.arcBeads = ab;
    plan.arcBeads.insert(plan.arcBeads.end(), bc.begin(), bc.end());
    plan.arcBeads.insert(plan.arcBeads.end(), ca.begin(), ca.end());
    plan.swaps = { { X, lXab, Y, lYab }, { Y, lYbc, Z, lZbc }, { Z, lZca, X, lXca } };
    plan.cross = { X, Y, Z };
    plan.arcs = { ab, bc, ca };
    plan.center = { (beads[X].pos.x + beads[Y].pos.x + beads[Z].pos.x) / 3.0f,
                    (beads[X].pos.y + beads[Y].pos.y + beads[Z].pos.y) / 3.0f };
    plan.ok = true;
    return plan;
}
// Transforms an outer edge bead (E) into a crossing by transferring the 
// transverse strand (slots 1 and 3) from the old crossing (V) onto E.
void LoopModel::RoleSwapCrossingAndEdge(int V, int E) {
    // Safety checks: V must be a crossing, E must be an edge bead
    if (V < 0 || E < 0 || !beads[V].isVertex || beads[E].isVertex) return;

    // 1. Identify the transverse strand beads (living on slots 1 and 3)
    int B1 = beads[V].neighborOfSlot[1];
    int B3 = beads[V].neighborOfSlot[3];

    // 2. Identify the direct strand configuration
    // Find which slot (0 or 2) E is currently connected to on V
    int slot_E = (beads[V].neighborOfSlot[0] == E) ? 0 :
        (beads[V].neighborOfSlot[2] == E) ? 2 : -1;

    if (slot_E == -1) return; // Error: E must be on the 0/2 strand 
    int slot_opposite = (slot_E == 0) ? 2 : 0;
    int A_opp = beads[V].neighborOfSlot[slot_opposite]; // The bead on the other side of V

    // Find E's outer neighbor (the one that isn't V)
    int X = (beads[E].neighbors[0] == V) ? beads[E].neighbors[1] : beads[E].neighbors[0];

    // --- EXECUTE THE SWAP ---

    // 3. Convert V into an Edge Bead
    beads[V].isVertex = false;
    beads[V].neighbors[0] = 29;
    beads[V].neighbors[1] = 3;

    // 4. Convert E into a Crossing
    beads[E].isVertex = true;

    // Connect the direct strand (slot 0 and 2)
    beads[E].neighborOfSlot[0] = 3;
    beads[E].neighborOfSlot[1] = 31;
    beads[E].neighborOfSlot[2] = -1;
    beads[E].neighborOfSlot[3] = -1;

    // Set the throughSlots exactly as requested (1 to 3, 2 to 0)
    beads[E].throughSlot[0] = 0;
    beads[E].throughSlot[2] = 1;
    beads[E].throughSlot[1] = -1;
    beads[E].throughSlot[3] = -1;

    // 5. Update the external transverse beads to point to the new crossing (E)
    ReplaceLink(30, V, E);
    ReplaceLink(B3, V, E);
}

// Detaches an outer bead (O) from a crossing (V) and connects it 
// to the side of an arc bead (A), converting A into a 3-way T-junction.
void LoopModel::CreateTConnection(int V, int O, int A) {
    // 1. Safety Checks
    if (V < 0 || O < 0 || A < 0) return;
    if (!beads[V].isVertex || beads[O].isVertex || beads[A].isVertex) return;

    // 2. Detach O from the original crossing V
    for (int s = 0; s < 4; s++) {
        if (beads[V].neighborOfSlot[s] == O) {
            beads[V].neighborOfSlot[s] = -1; // Sever the connection at V
            break;
        }
    }

    // 3. Detach V from O, and point O towards its new home (A)
    if (beads[O].neighbors[0] == V) {
        beads[O].neighbors[0] = A;
    }
    else if (beads[O].neighbors[1] == V) {
        beads[O].neighbors[1] = A;
    }

    // 4. Gather A's current straight-line neighbors
    int A_prev = beads[A].neighbors[0];
    int A_next = beads[A].neighbors[1];

    // 5. Promote A into a T-Junction (A crossing with a dead end)
    beads[A].isVertex = true;

    // Wire the straight arc through slots 0 and 2
    beads[A].neighborOfSlot[0] = A_prev;
    beads[A].neighborOfSlot[2] = A_next;

    // Plug the new side-strand (O) into slot 1
    beads[A].neighborOfSlot[1] = O;

    // Leave slot 3 dead (empty) to make it a strict T-connection
    beads[A].neighborOfSlot[3] = -1;

    // 6. Route the internal weave logic
    beads[A].throughSlot[0] = 2;
    beads[A].throughSlot[2] = 0;

    // Route the T-stem into the dead end to cap it off safely
    beads[A].throughSlot[1] = 3;
    beads[A].throughSlot[3] = 1;

    // (Note: You do not need to call ReplaceLink for A_prev or A_next, 
    // because their pointers are still correctly aimed at index 'A'!)
    //beads[V].isVertex = false;
}

//.
// Slide an existing T-junction one bead along its bar: the junction moves from
// T to the adjacent bar bead A, and its stem follows (it now branches at A). The
// old junction T becomes a plain bar bead again. A must be one of T's two bar
// neighbours (slot 0 or slot 2) and an ordinary (non-vertex) bead. This just
// relocates the branch point; it does not change the diagram.
void LoopModel::TraverseTConnection(int T, int A) {
    if (T < 0 || A < 0) return;
    if (!beads[T].active || !beads[A].active) return;      // a between-hop prune may have severed one
    if (!beads[T].isVertex || beads[A].isVertex) return;   // T = junction, A = plain bar bead

    // A must sit on T's bar (slot 0 or slot 2). slot 1 is the stem, slot 3 dead.
    int slotToA = (beads[T].neighborOfSlot[0] == A) ? 0 :
        (beads[T].neighborOfSlot[2] == A) ? 2 : -1;
    if (slotToA == -1) return;

    int stem = beads[T].neighborOfSlot[1];                       // branch strand bead
    int barOther = beads[T].neighborOfSlot[slotToA == 0 ? 2 : 0];    // T's other bar neighbour
    int A_far = (beads[A].neighbors[0] == T) ? beads[A].neighbors[1]
        : beads[A].neighbors[0];

    // Demote T back to a plain bar bead, now sitting between barOther and A.
    beads[T].isVertex = false;
    for (int s = 0; s < 4; ++s) { beads[T].neighborOfSlot[s] = -1; beads[T].throughSlot[s] = -1; }
    beads[T].neighbors[0] = barOther;
    beads[T].neighbors[1] = A;

    // Promote A into the T-junction: bar = T <-> A_far, stem = the moved branch.
    beads[A].isVertex = true;
    beads[A].neighborOfSlot[0] = T;
    beads[A].neighborOfSlot[2] = A_far;
    beads[A].neighborOfSlot[1] = stem;
    beads[A].neighborOfSlot[3] = -1;
    beads[A].throughSlot[0] = 2; beads[A].throughSlot[2] = 0;   // bar passes straight through
    beads[A].throughSlot[1] = 3; beads[A].throughSlot[3] = 1;   // stem capped into the dead slot

    // The stem strand now attaches to A instead of T.
    if (stem >= 0) ReplaceLink(stem, T, A);

    // Option 3: drop the stem on the local bar's outward NORMAL. A's bar now runs
    // T (slot 0, behind) -> A_far (slot 2, ahead); the stem leaves at ~90deg to it
    // on the outward side, following the arc as the junction walks it. r3Center is
    // used only to pick the outward sign, not the push direction.
    if (r3Active && R3_NORMAL_STEM && stem >= 0) {
        PlaceBeadNormalToBar(beads[stem], beads[A].pos, beads[T].pos, beads[A_far].pos,
            r3Center, R3_STEM_NORMAL_LEN);
    }
}


void LoopModel::MergeTConnections(int T1, int T2) {
    if (T1 < 0 || T2 < 0 || T1 == T2) return;
    if (!beads[T1].active || !beads[T2].active) return;
    if (!beads[T1].isVertex || !beads[T2].isVertex) return;
    // Both must be strict T-junctions: bar on slots 0/2, stem on slot 1, slot 3 dead.
    if (beads[T1].neighborOfSlot[3] != -1 || beads[T2].neighborOfSlot[3] != -1) return;

    // They must be adjacent: T2 sits on one of T1's bar slots (and vice-versa).
    int slotT1toT2 = (beads[T1].neighborOfSlot[0] == T2) ? 0 :
        (beads[T1].neighborOfSlot[2] == T2) ? 2 : -1;
    int slotT2toT1 = (beads[T2].neighborOfSlot[0] == T1) ? 0 :
        (beads[T2].neighborOfSlot[2] == T1) ? 2 : -1;
    if (slotT1toT2 == -1 || slotT2toT1 == -1) return;

    int outerA = beads[T1].neighborOfSlot[slotT1toT2 == 0 ? 2 : 0];   // T1's far bar leg
    int outerB = beads[T2].neighborOfSlot[slotT2toT1 == 0 ? 2 : 0];   // T2's far bar leg
    int s1 = beads[T1].neighborOfSlot[1];                         // T1's stem
    int s2 = beads[T2].neighborOfSlot[1];                         // T2's stem

    // Need four real, distinct legs to make a clean crossing.
    if (outerA < 0 || outerB < 0 || s1 < 0 || s2 < 0) return;
    if (outerA == outerB || s1 == s2) return;

    // Weave: strand 1 = outerA<->s2 (slots 0<->2), strand 2 = outerB<->s1 (slots 1<->3).
    beads[T1].isVertex = true;
    beads[T1].neighborOfSlot[0] = outerA;
    beads[T1].neighborOfSlot[1] = s2;
    beads[T1].neighborOfSlot[2] = outerB;
    beads[T1].neighborOfSlot[3] = s1;
    beads[T1].throughSlot[0] = 1; beads[T1].throughSlot[1] = 0;
    beads[T1].throughSlot[2] = 3; beads[T1].throughSlot[3] = 2;

    // outerA and s1 already reference T1; only T2's two outward legs need repointing.
    ReplaceLink(outerB, T2, T1);
    ReplaceLink(s2, T2, T1);

    // Seat the crossing midway between the old junctions, then retire T2.
    beads[T1].pos = { (beads[T1].pos.x + beads[T2].pos.x) * 0.5f,
                      (beads[T1].pos.y + beads[T2].pos.y) * 0.5f };

    // Option 3: at the fused crossing there is no stem/bar -- both strands pass
    // straight through -- so "perpendicular" becomes "straight". Offset the crossing
    // along the outerA--outerB chord's outward normal, then seat each transverse leg
    // opposite its through-partner (strand outerA--T1--s2, strand outerB--T1--s1) so
    // each strand runs straight through T1 as a clean X rather than kinking.
    if (r3Active && R3_NORMAL_STEM) {
        PlaceBeadNormalToBar(beads[T1], beads[T1].pos, beads[outerA].pos, beads[outerB].pos,
            r3Center, R3_STEM_NORMAL_LEN);
        auto seatOpposite = [&](int leg, int partner) {
            if (leg < 0 || partner < 0) return;
            float dx = beads[T1].pos.x - beads[partner].pos.x;
            float dy = beads[T1].pos.y - beads[partner].pos.y;
            float L = std::hypot(dx, dy);
            if (L <= GEOM_EPS) return;
            beads[leg].pos.x = beads[T1].pos.x + dx / L * REST_LENGTH;
            beads[leg].pos.y = beads[T1].pos.y + dy / L * REST_LENGTH;
            };
        seatOpposite(s2, outerA);   // strand outerA--T1--s2
        seatOpposite(s1, outerB);   // strand outerB--T1--s1
    }
    beads[T2].active = false;
    beads[T2].isVertex = false;
    beads[T2].crossing = -1;
    for (int s = 0; s < 4; ++s) { beads[T2].neighborOfSlot[s] = -1; beads[T2].throughSlot[s] = -1; }
    beads[T2].neighbors[0] = beads[T2].neighbors[1] = -1;
}

// Option 3b: a soft, length-preserving angular spring. At every live strict
// T-junction it pulls the stem toward the bar's outward normal, so the branch is
// held at ~90deg across all physics substeps instead of only being placed once per
// hop (which the edge springs can undo before the next hop). Call once per substep
// from the physics loop while r3Active; inert otherwise. `center` is the triangle
// centroid (used only to choose the outward side); `k` is the spring stiffness.
void LoopModel::ApplyStemNormalForces(Vector2 center, float k) {
    for (int i = 0; i < (int)beads.size(); ++i) {
        Bead& J = beads[i];
        if (!J.active || !J.isVertex) continue;
        if (J.neighborOfSlot[3] != -1) continue;          // strict T-junctions only (dead slot 3)
        int bar0 = J.neighborOfSlot[0], bar2 = J.neighborOfSlot[2], stem = J.neighborOfSlot[1];
        if (bar0 < 0 || bar2 < 0 || stem < 0) continue;

        float tx = beads[bar2].pos.x - beads[bar0].pos.x;  // bar tangent
        float ty = beads[bar2].pos.y - beads[bar0].pos.y;
        float L = std::hypot(tx, ty);
        if (L <= GEOM_EPS) continue;
        float nx = -ty / L, ny = tx / L;                   // unit normal to the bar
        float ox = J.pos.x - center.x, oy = J.pos.y - center.y;
        if (nx * ox + ny * oy < 0.0f) { nx = -nx; ny = -ny; }   // outward side

        float sx = beads[stem].pos.x - J.pos.x, sy = beads[stem].pos.y - J.pos.y;
        float sl = std::hypot(sx, sy);
        if (sl <= GEOM_EPS) continue;
        float tX = J.pos.x + nx * sl, tY = J.pos.y + ny * sl;   // same radius, on the normal
        beads[stem].force.x += (tX - beads[stem].pos.x) * k;
        beads[stem].force.y += (tY - beads[stem].pos.y) * k;

        // Option 1: walk the outer leg past the stem and hold each bead on the
        // outward side of the bar line through J. One-sided -- only a bead that has
        // drifted to the inward (centroid) side is pushed out -- with per-hop
        // falloff, so legs already outside keep their natural shape. This is what
        // stops the leg (not just the stem) from sweeping across the trigon.
        int prev = i, cur = stem;
        for (int h = 1; h <= R3_STEM_LEG_SPAN; ++h) {
            const Bead& C = beads[cur];
            if (C.isVertex) break;
            int nxt = (C.neighbors[0] == prev) ? C.neighbors[1] : C.neighbors[0];
            if (nxt < 0 || nxt >= (int)beads.size() || !beads[nxt].active || beads[nxt].isVertex) break;
            float d = (beads[nxt].pos.x - J.pos.x) * nx + (beads[nxt].pos.y - J.pos.y) * ny;  // + = outward of bar line
            if (d < R3_STEM_LEG_MARGIN) {
                float w = 1.0f - float(h) / float(R3_STEM_LEG_SPAN + 1);   // closest leg beads held strongest
                float push = (R3_STEM_LEG_MARGIN - d) * k * w;
                beads[nxt].force.x += nx * push;
                beads[nxt].force.y += ny * push;
            }
            prev = cur; cur = nxt;
        }
    }
}

// Hard containment for the animated R3 walk: the trigon face (boundary =
// r3Ring, the corners + arc beads captured at move start) must stay empty, so
// any active bead that is NOT on that boundary and falls inside the polygon is
// projected out past its nearest edge (+ margin), pushed radially away from the
// centroid r3Center. This is what holds the outer-leg stems -- the strands that
// would otherwise be dragged across the region and clip through the arcs --
// outside the triangle for the whole walk. Anchors (crossings, dragged beads)
// and beads in a collapse transition are left in place, matching the exemptions
// ResolveOverlaps uses, so the two positional passes never fight the same bead.
// The migrating junctions are themselves ring beads and so are exempt; only
// their branch strands get pushed out. Call once per physics substep while
// r3Active; inert otherwise.
void LoopModel::KeepOutOfTrigon(float margin) {
    r3Ejected.clear();
    if (!r3Active || r3Ring.size() < 3) return;

    std::vector<Vector2> poly; poly.reserve(r3Ring.size());
    std::vector<char> onRing(beads.size(), 0);
    for (int b : r3Ring)
        if (b >= 0 && b < (int)beads.size() && beads[b].active) { poly.push_back(beads[b].pos); onRing[b] = 1; }
    if (poly.size() < 3) return;

    const size_t n = poly.size();
    for (int i = 0; i < (int)beads.size(); ++i) {
        Bead& p = beads[i];
        if (!p.active || onRing[i]) continue;
        if (p.isVertex || p.driven || transitions.count(i)) continue;   // immovable / mid-collapse: skip (see ResolveOverlaps)
        if (!PtInPolyLocal(p.pos, poly)) continue;

        // Nearest boundary edge, then re-seat p just outside it.
        float best2 = 1e30f; Vector2 bestQ = p.pos;
        for (size_t a = 0, b = n - 1; a < n; b = a++) {
            float d2; Vector2 q = ClosestOnSeg(p.pos, poly[b], poly[a], d2);
            if (d2 < best2) { best2 = d2; bestQ = q; }
        }
        float ox = bestQ.x - r3Center.x, oy = bestQ.y - r3Center.y;   // outward = away from centroid
        float L = std::hypot(ox, oy);
        if (L <= GEOM_EPS) { ox = p.pos.x - r3Center.x; oy = p.pos.y - r3Center.y; L = std::hypot(ox, oy); }
        if (L <= GEOM_EPS) continue;
        p.pos.x = bestQ.x + ox / L * margin;
        p.pos.y = bestQ.y + oy / L * margin;
        r3Ejected.push_back(i);
    }
}

// crossingBead must list exactly the live crossings. The T-junction R3 demotes the
// three old corners to edge beads and promotes three new mid-arc crossings, so after
// that surgery the build-time list is stale (it feeds the physics crossing-guard and
// auto-tighten). Rederive it from the beads themselves.
void LoopModel::RebuildCrossingList() {
    crossingBead.clear();
    for (int i = 0; i < (int)beads.size(); ++i)
        if (beads[i].active && beads[i].isVertex) crossingBead.push_back(i);
}

// Moves the outer leg bead across the crossing vertex V into the internal arc
void LoopModel::SlideOuterBeadIntoArc(int V, int slot_in) {
    // 1. Safety check: Ensure V is a valid crossing vertex
    if (V < 0 || !beads[V].isVertex) return;

    // 2. Identify the slots and the involved beads
    int slot_out = beads[V].throughSlot[slot_in];
    int O = beads[V].neighborOfSlot[slot_out]; // The outer bead to be moved
    int I = beads[V].neighborOfSlot[slot_in];   // The first inner bead of the arc

    // Safety check: The outer neighbor must be a normal edge bead to slide it
    if (O < 0 || beads[O].isVertex) return;

    // Find the next bead further out on the external leg
    int Onext = (beads[O].neighbors[0] == V) ? beads[O].neighbors[1] : beads[O].neighbors[0];

    // Helper lambda to safely update a neighbor connection on any bead type (vertex or edge)
    auto ReplaceNeighbor = [&](int beadIdx, int oldNeighbor, int newNeighbor) {
        if (beadIdx < 0) return;
        if (beads[beadIdx].isVertex) {
            for (int s = 0; s < 4; ++s) {
                if (beads[beadIdx].neighborOfSlot[s] == oldNeighbor) {
                    beads[beadIdx].neighborOfSlot[s] = newNeighbor;
                    break;
                }
            }
        }
        else {
            if (beads[beadIdx].neighbors[0] == oldNeighbor) beads[beadIdx].neighbors[0] = newNeighbor;
            else if (beads[beadIdx].neighbors[1] == oldNeighbor) beads[beadIdx].neighbors[1] = newNeighbor;
        }
        };

    // 3. Perform the Topology Surgery

    // Update the outer-most anchor to point to the crossing V instead of O
    ReplaceNeighbor(Onext, O, V);

    // Update the inner-most anchor to point to O instead of V
    ReplaceNeighbor(I, V, O);

    // Update the crossing slots: outer slot now sees Onext, inner slot now sees O
    beads[V].neighborOfSlot[slot_out] = Onext;
    beads[V].neighborOfSlot[slot_in] = O;

    // Update the moved bead O to sit directly between V and I
    beads[O].neighbors[0] = V;
    beads[O].neighbors[1] = I;
}


void LoopModel::SlideCrossingOverBead(int c, int slot) {
    int b = beads[c].neighborOfSlot[slot];
    if (b < 0 || beads[b].isVertex) return;

    int throughSlot = beads[c].throughSlot[slot];
    int nextB = beads[c].neighborOfSlot[throughSlot];

    int prevB = (beads[b].neighbors[0] == c) ? beads[b].neighbors[1] : beads[b].neighbors[0];

    ReplaceLink(prevB, b, c);
    beads[c].neighborOfSlot[slot] = prevB;

    beads[c].neighborOfSlot[throughSlot] = b;
    beads[b].neighbors[0] = c;
    beads[b].neighbors[1] = nextB;

    ReplaceLink(nextB, c, b);
}

int LoopModel::FindSlotTowardsCrossing(int crossingA, int crossingB) const {
    // Check all 4 ports on crossing A
    for (int slot = 0; slot < 4; slot++) {
        int currentBead = beads[crossingA].neighborOfSlot[slot];
        int prevBead = crossingA;

        // Walk along the normal edge beads
        while (currentBead != -1 && !beads[currentBead].isVertex) {
            int nextBead = (beads[currentBead].neighbors[0] == prevBead) ?
                beads[currentBead].neighbors[1] :
                beads[currentBead].neighbors[0];
            prevBead = currentBead;
            currentBead = nextBead;
        }

        // If the vertex we hit at the end of this strand is crossingB, return slot
        if (currentBead == crossingB) {
            return slot;
        }
    }
    return -1;
}

void LoopModel::Perform_R3_Walk_S3(int V0, int V1, int V2, FaceClass& f, std::function<void(int, std::function<void()>)> schedule) {
    // Degenerate-arc guard, hardened. The converging-walk surgery needs each
    // triangle arc to carry two DISTINCT T-junctions (one per end) that meet and
    // merge. A one-bead arc collapses both ends onto the same bead (leaving an
    // unmerged junction whose dead slot is later dereferenced as beads[-1]); an
    // empty arc -- two corners linked directly -- has no bead to walk at all. So we
    // now PAD every arc up to R3_MIN_ARC_BEADS before touching it, splitting the
    // longest link each time (SpliceBeadBetween can lengthen even a 0-bead arc, so
    // empty arcs are repaired rather than aborting the whole move). Padding also
    // makes the walk smoother: more beads => more, smaller hops.
    if (f.triArcs.size() == 3) {
        const int Vk[3] = { V0, V1, V2 };
        // The arc's endpoint crossings: derive from the live links when the arc has
        // beads; fall back to the corner order (triArcs[k] joins Vk[k] -> Vk[k+1])
        // for an empty arc, whose endpoints can't be read from beads.
        auto adjacentCorner = [&](int bead, int except) -> int {
            for (int c : Vk) {
                if (c == except || c < 0 || !beads[c].isVertex) continue;
                for (int s = 0; s < 4; ++s) if (beads[c].neighborOfSlot[s] == bead) return c;
            }
            return -1;
            };
        auto padArc = [&](std::vector<int>& arc, int cfFallback, int cbFallback, int target) {
            int cf = arc.empty() ? cfFallback : adjacentCorner(arc.front(), -1);
            int cb = arc.empty() ? cbFallback : adjacentCorner(arc.back(), cf);
            if (cf < 0) cf = cfFallback;
            if (cb < 0) cb = cbFallback;
            int guard = 0;
            while ((int)arc.size() < target && guard++ < target + 8) {
                std::vector<int> nodes; nodes.reserve(arc.size() + 2);
                nodes.push_back(cf);
                for (int b : arc) nodes.push_back(b);
                nodes.push_back(cb);
                int bestI = -1; float bestLen2 = -1.0f;
                for (int i = 0; i + 1 < (int)nodes.size(); ++i) {
                    Vector2 p = beads[nodes[i]].pos, q = beads[nodes[i + 1]].pos;
                    float d2 = (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y);
                    if (d2 > bestLen2) { bestLen2 = d2; bestI = i; }
                }
                if (bestI < 0) break;
                int m = SpliceBeadBetween(nodes[bestI], nodes[bestI + 1]);
                if (m < 0) break;                       // not adjacent (bad endpoints): stop
                arc.insert(arc.begin() + bestI, m);     // link i lands at arc index i
            }
            };
        for (int k = 0; k < 3; ++k)
            padArc(f.triArcs[k], Vk[k], Vk[(k + 1) % 3], R3_MIN_ARC_BEADS);
    }
    for (auto& arc : f.triArcs)
        if ((int)arc.size() < 2) { printf("[R3] arc too short after padding; move skipped\n"); return; }

    // The arc's endpoint bead names the exact leg, so slot-of-endpoint is
    // unambiguous even when an earlier R3 left >1 path between two crossings.
    auto slotForArc = [&](int V, const std::vector<int>& arc) -> int {
        if (arc.empty()) return -1;
        int f0 = arc.front(), b0 = arc.back();
        for (int s = 0; s < 4; ++s) {
            int n = beads[V].neighborOfSlot[s];
            if (n >= 0 && (n == f0 || n == b0)) return s;
        }
        return -1;
        };
    // Capture the triangle centroid from the corner positions BEFORE the surgery
    // below demotes them, and mark the walk active for its (possibly scheduled)
    // duration: this arms the normal-stem option / angular spring and suppresses
    // arc re-normalisation mid-walk. r3Active is cleared once the last merge has
    // run (see the schedule tail below / the synchronous else-branch).
    r3Center = { (beads[V0].pos.x + beads[V1].pos.x + beads[V2].pos.x) / 3.0f,
                 (beads[V0].pos.y + beads[V1].pos.y + beads[V2].pos.y) / 3.0f };
    r3Active = true;
    // Boundary cycle V0 -> arc0 -> V1 -> arc1 -> V2 -> arc2 (triArcs are stored
    // c_k -> c_{k+1}, so this concatenation is a simple loop). The corners demote
    // to bar beads and the junctions walk along the arcs, but every index here
    // stays a live boundary bead, so KeepOutOfTrigon can rebuild the polygon from
    // current positions each substep.
    r3Ring.clear();
    r3Ring.push_back(V0);
    r3Ring.insert(r3Ring.end(), f.triArcs[0].begin(), f.triArcs[0].end());
    r3Ring.push_back(V1);
    r3Ring.insert(r3Ring.end(), f.triArcs[1].begin(), f.triArcs[1].end());
    r3Ring.push_back(V2);
    r3Ring.insert(r3Ring.end(), f.triArcs[2].begin(), f.triArcs[2].end());
    Perform_R3_Walk_S1(V0, slotForArc(V0, f.triArcs[0]), slotForArc(V0, f.triArcs[2]), f.triArcs[0], f.triArcs[2], schedule);
    Perform_R3_Walk_S1(V1, slotForArc(V1, f.triArcs[0]), slotForArc(V1, f.triArcs[1]), f.triArcs[0], f.triArcs[1], schedule);
    Perform_R3_Walk_S1(V2, slotForArc(V2, f.triArcs[2]), slotForArc(V2, f.triArcs[1]), f.triArcs[2], f.triArcs[1], schedule);

    // Each arc now carries two T-junctions walking to its midpoint. Once they have
    // met, weld each kissing pair into the move's new crossing. An arc of length n
    // is walked in n/2 hops (last hop at n/2 * gap), so merge it one gap later.
    // The merge callback re-locates the pair at run time (the junctions are the
    // arc's only vertices), capturing a COPY of the arc so it stays valid after the
    // caller's FaceClass is gone. After all three pairs fuse, the old corner
    // crossings are demoted and three new ones exist, so refresh crossingBead.
    auto mergeArc = [this](std::vector<int> arc) {
        int t1 = -1, t2 = -1;
        for (int b : arc)
            if (b >= 0 && beads[b].active && beads[b].isVertex) { (t1 < 0 ? t1 : t2) = b; }
        if (t1 >= 0 && t2 >= 0) MergeTConnections(t1, t2);
        };

    if (schedule) {
        int lastMerge = 0;
        for (int k = 0; k < 3; ++k) {
            std::vector<int> arc = f.triArcs[k];
            int when = ((int)arc.size() / 2 + 1) * R3_WALK_GAP;
            lastMerge = std::max(lastMerge, when);
            schedule(when, [mergeArc, arc] { mergeArc(arc); });
        }
        // Topology is final one gap after the last merge; keep r3Active (and so
        // the trigon containment) running a few extra gaps so anything still
        // trapped inside is expelled during the post-merge settle, not left
        // clipped when the flag drops.
        schedule(lastMerge + R3_WALK_GAP, [this] { RebuildCrossingList(); });
        schedule(lastMerge + 8 * R3_WALK_GAP, [this] { r3Active = false; r3Ring.clear(); });
    }
    else {
        mergeArc(f.triArcs[0]);
        mergeArc(f.triArcs[1]);
        mergeArc(f.triArcs[2]);
        RebuildCrossingList();
        r3Active = false;
        r3Ring.clear();
    }
}

void LoopModel::Perform_R3_Walk_S1(int V, int S1, int S2, std::vector<int>& arc_b1, std::vector<int>& arc_b2,
    std::function<void(int, std::function<void()>)> schedule) {
    // Guards: V must be a crossing and S1/S2 real slots (FindSlotTowardsCrossing
    // returns -1 when there's no clean arc, and neighborOfSlot[-1] is OOB). Also
    // need at least one bead on each arc to read its end.
    if (V < 0 || !beads[V].isVertex) return;
    if (S1 < 0 || S1 > 3 || S2 < 0 || S2 > 3) return;
    if (arc_b1.empty() || arc_b2.empty()) return;

    int Ba = beads[V].neighborOfSlot[S1];
    int Bb = beads[V].neighborOfSlot[S2];
    int Ya = beads[V].neighborOfSlot[beads[V].throughSlot[S1]];
    int Yb = beads[V].neighborOfSlot[beads[V].throughSlot[S2]];
    if (Ba < 0 || Bb < 0 || Ya < 0 || Yb < 0) return;

    // Smooth crossing V open into two T-junctions, one on each leg (Ba, Bb)...
    CreateTConnection(V, Yb, Ba);
    CreateTConnection(V, Ya, Bb);
    // ...then turn the leftover 2-valent V into a plain bar bead joining the two
    // new junctions. After the two CreateTConnection calls, Ba and Bb are the only
    // beads that still reference V (their bars run V<->arc), and V's links to Ya/Yb
    // were severed. So V's correct neighbours are exactly {Ba, Bb} - but only if Ba
    // and Bb really do point back at V. Verify that and repair the back-links so the
    // demoted bead is never left half-connected.
    beads[V].isVertex = false;
    beads[V].crossing = -1;
    for (int s = 0; s < 4; ++s) { beads[V].neighborOfSlot[s] = -1; beads[V].throughSlot[s] = -1; }
    beads[V].neighbors[0] = Ba;
    beads[V].neighbors[1] = Bb;
    // Guarantee the bond is mutual: each junction must list V on its bar (slot 0/2).
    auto ensureBarLink = [&](int T, int back) {
        if (T < 0 || !beads[T].isVertex) return;
        if (beads[T].neighborOfSlot[0] != back && beads[T].neighborOfSlot[2] != back) {
            // bar slot that isn't the stem(1)/dead(3) and is currently dangling -> point it at V
            if (beads[T].neighborOfSlot[0] < 0) beads[T].neighborOfSlot[0] = back;
            else if (beads[T].neighborOfSlot[2] < 0) beads[T].neighborOfSlot[2] = back;
        }
        };
    ensureBarLink(Ba, V);
    ensureBarLink(Bb, V);

    // The junction sits on the arc end adjacent to V. The triArcs are stored
    // c0->cN, so if Ba/Bb is the arc's FRONT we walk forward (index up), else the
    // arc runs the other way and we walk backward (index down). Either way the
    // step count is bounded by the arc length - the old hard-coded 3 ran past the
    // end of a 3-bead arc. Walk j and j+-1, never reading outside [0, size).
    // Walk the junction along the arc. If `schedule` was supplied, spread the
    // hops over time - one TraverseTConnection every `gap` frames - so the move
    // animates; otherwise do them all at once (original behaviour). Each hop
    // captures only the two bead indices + `this`, so it stays valid after this
    // function (and the caller's arc vector) has returned.
    const int gap = R3_WALK_GAP;   // frames between hops
    auto walkAlong = [&](std::vector<int>& arc, int startBead) {
        const int n = (int)arc.size();              // FULL arc length (do NOT halve the index bound)
        if (n < 2) return;                          // nothing to slide across
        const int steps = n / 2;                    // walk only as far as the arc midpoint
        std::vector<std::pair<int, int>> hops;      // ordered (from,to), junction end first
        if (arc.front() == startBead)               // junction at front -> walk up toward middle
            for (int j = 0; j < steps && j + 1 < n; ++j) hops.push_back({ arc[j], arc[j + 1] });
        else                                        // junction at back -> walk down toward middle
            for (int j = n - 1; j > n - 1 - steps && j - 1 >= 0; --j) hops.push_back({ arc[j], arc[j - 1] });

        for (size_t k = 0; k < hops.size(); ++k) {
            int from = hops[k].first, to = hops[k].second;
            if (schedule)
                schedule((int)(k + 1) * gap, [this, from, to] { TraverseTConnection(from, to); });
            else
                TraverseTConnection(from, to);
        }
        };
    walkAlong(arc_b1, Ba);
    walkAlong(arc_b2, Bb);
}

void LoopModel::StartR3MicroSurgery(int c0, int c1, int c2) {
    // Find the slot that points from c0 inward towards c1
    int slotTowardsC1 = FindSlotTowardsCrossing(c0, c1);

    if (slotTowardsC1 != -1) {
        printf("Sliding crossing [0] %d towards crossing %d on slot %d\n", c0, c1, slotTowardsC1);
        //SlideOuterBeadIntoArc(c0, slotTowardsC1);
        RoleSwapCrossingAndEdge(c0, beads[c0].neighborOfSlot[1]);
    }
    else {
        printf("Error: c0 and c1 do not share a direct arc!\n");
    }
}