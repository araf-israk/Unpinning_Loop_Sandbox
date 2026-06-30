// ============================================================================
// upg_surgery.cpp - topology surgery primitives on the bead necklace.
//
// The low-level link-rewriting operations every Reidemeister move is built from,
// plus the dynamic resampling that keeps arcs densely enough sampled that the
// physics self-avoidance cannot be threaded through a gap:
//
//   * WalkFromSlot / FindSelfLoop / ArcsBetween
//                       - follow a strand from a crossing slot to the next
//                         crossing; detect a removable self-loop (R1 kink).
//   * ReplaceLink       - repoint one bead's reference oldN -> newN (handles
//                         both crossing slots and edge-bead neighbours).
//   * CollapseRegion    - splice a connected doomed region out, welding each
//                         entering strand to the one it continues to on the far
//                         side (the shared R1 + R2 removal primitive).
//   * InsertEdgeBead / SplitLink / ResampleStretchedArcs
//                       - add edge beads so no link is long enough for another
//                         strand to slip through.
//   * RemoveEdgeBead / SimplifyDenseArcs
//                       - the inverse: drop redundant beads at hairpins / piles.
//   * BuildCrossingStrandMap
//                       - tag beads near each crossing by which of its two
//                         strands they sit on, so the de-overlap pass exempts
//                         exactly the two strands that genuinely cross there.
// ============================================================================

#include "upg.h"

#include <algorithm>
#include <set>

// ----------------------------------------------------------------------------
// topology surgery
// ----------------------------------------------------------------------------

// Follow the arc leaving crossing `c` at `slot` until the next crossing.
ArcWalk LoopModel::WalkFromSlot(int c, int slot) const {
    ArcWalk w;
    int from = c, cur = beads[c].neighborOfSlot[slot];
    while (cur >= 0 && !beads[cur].isVertex) {
        w.beads.push_back(cur);
        int nxt = (beads[cur].neighbors[0] == from) ? beads[cur].neighbors[1]
            : beads[cur].neighbors[0];
        from = cur; cur = nxt;
    }
    w.endVertex = cur;   // crossing (or -1)
    return w;
}

// A self-loop arc at `c` (leaves c and returns to c with no other crossing) is
// an R1-removable kink. Returns its beads, or empty if none.
std::vector<int> LoopModel::FindSelfLoop(int c) const {
    for (int s = 0; s < 4; ++s) {
        if (beads[c].neighborOfSlot[s] < 0) continue;
        ArcWalk w = WalkFromSlot(c, s);
        if (w.endVertex == c) return w.beads;
    }
    return {};
}

// All crossing-free arcs running from A to B (found from A's four slots).
std::vector<ArcWalk> LoopModel::ArcsBetween(int A, int B) const {
    std::vector<ArcWalk> out;
    for (int s = 0; s < 4; ++s) {
        if (beads[A].neighborOfSlot[s] < 0) continue;
        ArcWalk w = WalkFromSlot(A, s);
        if (w.endVertex == B) out.push_back(w);
    }
    return out;
}

void LoopModel::ReplaceLink(int x, int oldN, int newN) {
    if (beads[x].isVertex) {
        for (int s = 0; s < 4; ++s) if (beads[x].neighborOfSlot[s] == oldN) beads[x].neighborOfSlot[s] = newN;
    }
    else {
        for (int k = 0; k < 2; ++k) if (beads[x].neighbors[k] == oldN) beads[x].neighbors[k] = newN;
    }
}

// Splice a connected doomed region out of the loop: each strand entering the
// region from outside is reconnected to the strand it continues to on the far
// side, then the doomed beads are deactivated. Works for both R1 (crossing +
// self-loop) and R2 (two crossings + the two lens arcs).
void LoopModel::CollapseRegion(const std::vector<int>& doomedV) {
    std::set<int> D(doomedV.begin(), doomedV.end());
    auto inD = [&](int x) { return x >= 0 && D.count(x) > 0; };

    auto neighborsOf = [&](int x) {
        std::vector<int> r;
        if (beads[x].isVertex) for (int s = 0; s < 4; ++s) r.push_back(beads[x].neighborOfSlot[s]);
        else { r.push_back(beads[x].neighbors[0]); r.push_back(beads[x].neighbors[1]); }
        return r;
        };

    std::set<int> done;
    for (int b = 0; b < (int)beads.size(); ++b) {
        if (inD(b) || !beads[b].active || done.count(b)) continue;
        for (int entry : neighborsOf(b)) {
            if (!inD(entry) || done.count(b)) continue;
            // walk through D to the strand's exit
            int from = b, cur = entry;
            while (inD(cur)) {
                int nxt = beads[cur].isVertex ? Through(cur, from)
                    : (beads[cur].neighbors[0] == from ? beads[cur].neighbors[1] : beads[cur].neighbors[0]);
                from = cur; cur = nxt;
                if (cur < 0) break;
            }
            if (cur < 0) continue;
            int exitBead = cur, exitFrom = from;        // exitBead is outside D
            ReplaceLink(b, entry, exitBead);
            ReplaceLink(exitBead, exitFrom, b);
            done.insert(b); done.insert(exitBead);
        }
    }
    for (int d : D) beads[d].active = false;
}

// Insert a fresh edge bead at the midpoint of the link between two adjacent,
// non-vertex edge beads a and b, rewiring so the strand reads a <-> new <-> b.
// Returns the new bead's index, or -1 if a/b are not a splittable edge pair.
int LoopModel::InsertEdgeBead(int a, int b) {
    if (a < 0 || b < 0 || a >= (int)beads.size() || b >= (int)beads.size()) return -1;
    const Bead& ba = beads[a];
    const Bead& bb = beads[b];
    if (ba.isVertex || bb.isVertex || !ba.active || !bb.active) return -1;
    bool aLinksB = (ba.neighbors[0] == b || ba.neighbors[1] == b);
    bool bLinksA = (bb.neighbors[0] == a || bb.neighbors[1] == a);
    if (!aLinksB || !bLinksA) return -1;

    // Build the new bead from a/b BEFORE push_back (which may reallocate and
    // dangle these references). It inherits a's edge id and sits at the midpoint.
    Bead nb;
    nb.id = (int)beads.size();
    nb.isVertex = false;
    nb.active = true;
    nb.edge = ba.edge;
    nb.pos = { (ba.pos.x + bb.pos.x) * 0.5f, (ba.pos.y + bb.pos.y) * 0.5f };
    nb.neighbors[0] = a;
    nb.neighbors[1] = b;

    int nbi = nb.id;
    beads.push_back(nb);                 // ba/bb references are now invalid

    // ReplaceLink re-fetches by index, so it is safe after the push_back.
    ReplaceLink(a, b, nbi);
    ReplaceLink(b, a, nbi);
    return nbi;
}

// Insert one edge bead on the link between edge bead x and its neighbour y, where
// y may be a crossing (repointed via its slot) or another edge bead (via its
// neighbours). Unlike InsertEdgeBead this allows a crossing endpoint, so it can
// lengthen a one-bead arc whose single bead touches two crossings. Returns the
// new bead index, or -1.
int LoopModel::SplitLink(int x, int y) {
    if (x < 0 || y < 0 || x >= (int)beads.size() || y >= (int)beads.size()) return -1;
    if (beads[x].isVertex || !beads[x].active || !beads[y].active) return -1;
    if (beads[x].neighbors[0] != y && beads[x].neighbors[1] != y) return -1;

    Bead nb;
    nb.id = (int)beads.size();
    nb.isVertex = false;
    nb.active = true;
    nb.edge = beads[x].edge;
    nb.pos = { (beads[x].pos.x + beads[y].pos.x) * 0.5f,
               (beads[x].pos.y + beads[y].pos.y) * 0.5f };
    nb.neighbors[0] = x;
    nb.neighbors[1] = y;
    int m = nb.id;
    beads.push_back(nb);                 // x/y references now invalid; use indices

    if (beads[x].neighbors[0] == y) beads[x].neighbors[0] = m;
    else                            beads[x].neighbors[1] = m;

    if (beads[y].isVertex) {
        for (int s = 0; s < 4; ++s)
            if (beads[y].neighborOfSlot[s] == x) { beads[y].neighborOfSlot[s] = m; break; }
    }
    else {
        if (beads[y].neighbors[0] == x) beads[y].neighbors[0] = m;
        else                            beads[y].neighbors[1] = m;
    }
    return m;
}

// One resampling pass: every link between two edge beads longer than maxLen gets
// a bead at its midpoint. Links are deduplicated by endpoint pair (orientation
// of neighbors[] is not globally consistent after surgery), and beads currently
// mid-transition (collapsing) are skipped so a split can't corrupt a removal.
// A long link splits once per call and the springs relax it over later frames.
int LoopModel::ResampleStretchedArcs(float maxLen) {
    const float maxLen2 = maxLen * maxLen;
    std::set<std::pair<int, int>> seen;
    std::vector<std::pair<int, int>> toSplit;

    const int n = (int)beads.size();
    for (int i = 0; i < n; ++i) {
        const Bead& b = beads[i];
        if (b.isVertex || !b.active || transitions.count(i)) continue;
        for (int k = 0; k < 2; ++k) {
            int j = b.neighbors[k];
            if (j < 0 || j >= n) continue;
            const Bead& bj = beads[j];
            if (bj.isVertex || !bj.active || transitions.count(j)) continue;
            auto key = std::minmax(i, j);
            if (!seen.insert({ key.first, key.second }).second) continue;
            float dx = bj.pos.x - b.pos.x, dy = bj.pos.y - b.pos.y;
            if (dx * dx + dy * dy > maxLen2) toSplit.push_back({ i, j });
        }
    }
    for (auto& [a, b] : toSplit) InsertEdgeBead(a, b);
    return (int)toSplit.size();
}

// Splice a redundant edge bead out: its two strand-neighbours are linked directly
// and it is deactivated. Only edge beads are removed; vertices are never touched.
// Refuses cases that would corrupt the strand (a lone neighbour, a 2-cycle, or a
// pair already linked, which would otherwise create a duplicate edge).
int LoopModel::RemoveEdgeBead(int b) {
    if (b < 0 || b >= (int)beads.size()) return -1;
    if (beads[b].isVertex || !beads[b].active) return -1;
    int a = beads[b].neighbors[0], c = beads[b].neighbors[1];
    if (a < 0 || c < 0 || a == c) return -1;
    if (!beads[a].active || !beads[c].active) return -1;

    auto linked = [&](int x, int y) {
        if (beads[x].isVertex) {
            for (int s = 0; s < 4; ++s) if (beads[x].neighborOfSlot[s] == y) return true;
            return false;
        }
        return beads[x].neighbors[0] == y || beads[x].neighbors[1] == y;
        };
    if (linked(a, c) || linked(c, a)) return -1;   // already adjacent: removing b would double the edge

    ReplaceLink(a, b, c);
    ReplaceLink(c, b, a);
    beads[b].active = false;
    return b;
}

// One simplification pass. For each edge bead b with strand-neighbours a and c,
// drop b when EITHER the interior angle a-b-c is sharper than minAngleDeg (the
// strand hairpins back on itself there), OR a, b, c are three edge beads all
// within touchLen of one another (a pile). Candidates are gathered first, then
// removed while skipping any region already merged this pass, so no two adjacent
// beads are spliced at once. Remaining density is cleaned over later frames.
int LoopModel::SimplifyDenseArcs(float touchLen, float minAngleDeg) {
    const float kPi = 3.14159265358979323846f;
    const float touch2 = touchLen * touchLen;
    const float cosThresh = std::cos(minAngleDeg * kPi / 180.0f); // angle<thresh <=> cos>cosThresh

    auto linked = [&](int x, int y) {
        if (beads[x].isVertex) {
            for (int s = 0; s < 4; ++s) if (beads[x].neighborOfSlot[s] == y) return true;
            return false;
        }
        return beads[x].neighbors[0] == y || beads[x].neighbors[1] == y;
        };

    std::vector<int> doomed;
    const int n = (int)beads.size();
    for (int b = 0; b < n; ++b) {
        const Bead& B = beads[b];
        if (B.isVertex || !B.active || transitions.count(b)) continue;
        int a = B.neighbors[0], c = B.neighbors[1];
        if (a < 0 || c < 0 || a == c) continue;
        if (!beads[a].active || !beads[c].active) continue;
        if (transitions.count(a) || transitions.count(c)) continue;
        if (linked(a, c) || linked(c, a)) continue;   // tiny loop / parallel edge: leave it

        Vector2 pa = beads[a].pos, pb = B.pos, pc = beads[c].pos;
        float ux = pa.x - pb.x, uy = pa.y - pb.y;     // apex -> a
        float vx = pc.x - pb.x, vy = pc.y - pb.y;     // apex -> c
        float lu = std::hypot(ux, uy), lv = std::hypot(vx, vy);

        bool sharp = false;
        if (lu > GEOM_EPS && lv > GEOM_EPS)
            sharp = (ux * vx + uy * vy) / (lu * lv) > cosThresh;

        bool pile = false;                            // three edge beads mutually touching
        if (!beads[a].isVertex && !beads[c].isVertex) {
            float dac = (pa.x - pc.x) * (pa.x - pc.x) + (pa.y - pc.y) * (pa.y - pc.y);
            pile = (lu * lu < touch2) && (lv * lv < touch2) && (dac < touch2);
        }

        if (sharp || pile) doomed.push_back(b);
    }

    std::set<int> touched;
    int removed = 0;
    for (int b : doomed) {
        if (!beads[b].active) continue;
        int a = beads[b].neighbors[0], c = beads[b].neighbors[1];
        if (a < 0 || c < 0) continue;
        if (touched.count(a) || touched.count(b) || touched.count(c)) continue;
        if (RemoveEdgeBead(b) >= 0) { touched.insert(a); touched.insert(b); touched.insert(c); ++removed; }
    }
    return removed;
}

// Tag the edge beads near each crossing by which of the crossing's two strands
// they sit on. From crossing C, the strand entering slot 0 leaves through
// throughSlot[0], so {0, throughSlot[0]} is one strand and the other two slots
// are the other. Each incident arc is walked outward from C and its beads are
// recorded until the arc leaves the CROSS_GUARD zone (or hits the next vertex).
// The result lets the de-overlap pass exempt exactly the two strands that cross
// at C -- and nothing else that merely drifts into the radius.
std::vector<LoopModel::CrossingStrands> LoopModel::BuildCrossingStrandMap() const {
    std::vector<CrossingStrands> out;
    const float G2 = CROSS_GUARD * CROSS_GUARD;
    for (int cb : crossingBead) {
        if (cb < 0 || cb >= (int)beads.size() || !beads[cb].active) continue;
        CrossingStrands cs; cs.crossBead = cb;
        Vector2 cpos = beads[cb].pos;
        int t0 = beads[cb].throughSlot[0];
        for (int s = 0; s < 4; ++s) {
            int strand = (s == 0 || s == t0) ? 0 : 1;
            int from = cb, cur = beads[cb].neighborOfSlot[s];
            int guard = 0;
            while (cur >= 0 && ++guard < 100000 && beads[cur].active && !beads[cur].isVertex) {
                float dx = beads[cur].pos.x - cpos.x, dy = beads[cur].pos.y - cpos.y;
                if (dx * dx + dy * dy > G2) break;       // past the crossing zone
                cs.strandOf[cur] = strand;
                int nxt = (beads[cur].neighbors[0] == from) ? beads[cur].neighbors[1]
                    : beads[cur].neighbors[0];
                from = cur; cur = nxt;
            }
        }
        out.push_back(std::move(cs));
    }
    return out;
}
