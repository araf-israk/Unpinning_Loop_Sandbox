// ============================================================================
// upg_faces.cpp - arc / face queries and Reidemeister classification.
//
// Reads the live bead geometry (so it stays correct after any move rewires the
// loop) to answer "what move does the face under this point want?":
//
//   * ArcsFrom / SlotOf / ValidBigon
//                   - enumerate the crossing-free arcs leaving a crossing and
//                     test whether two of them bound a genuine 2-gon (they must
//                     be different strands at BOTH crossings).
//   * ClassifyFace  - find the smallest face containing a point and classify it
//                     by how many crossings border it: 1 -> R1 (self-loop
//                     monogon), 2 -> R2 (bigon), 3 -> R3 (triangle). Reports the
//                     beads to collapse (R1/R2) or the three bounding arcs (R3)
//                     so the move acts on exactly the clicked face.
// ============================================================================

#include "upg.h"

// ----------------------------------------------------------------------------
// arc / face queries
// ----------------------------------------------------------------------------

std::vector<SlotArc> LoopModel::ArcsFrom(int c) const {
    std::vector<SlotArc> r;
    for (int s = 0; s < 4; ++s) {
        if (beads[c].neighborOfSlot[s] < 0) continue;
        ArcWalk w = WalkFromSlot(c, s);
        r.push_back({ s, w.endVertex, w.beads });
    }
    return r;
}

int LoopModel::SlotOf(int v, int bead) const {
    for (int s = 0; s < 4; ++s) if (beads[v].neighborOfSlot[s] == bead) return s;
    return -1;
}

// Two arcs bound a genuine 2-gon only if they are different strands at BOTH
// crossings (else they are one strand and the "bigon" is degenerate).
bool LoopModel::ValidBigon(int A, int B, const SlotArc& i, const SlotArc& j) const {
    if (i.slot == j.slot || i.beads.empty() || j.beads.empty()) return false;
    if (beads[A].throughSlot[i.slot] == j.slot) return false;
    int bi = SlotOf(B, i.beads.back()), bj = SlotOf(B, j.beads.back());
    if (bi < 0 || bj < 0 || beads[B].throughSlot[bi] == bj) return false;
    return true;
}

namespace {

    bool PointInPoly(Vector2 p, const std::vector<Vector2>& poly) {
        bool in = false;
        for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
            if (((poly[i].y > p.y) != (poly[j].y > p.y)) &&
                (p.x < (poly[j].x - poly[i].x) * (p.y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x))
                in = !in;
        return in;
    }

    float PolyArea(const std::vector<Vector2>& poly) {
        float a = 0.0f;
        for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
            a += (poly[j].x + poly[i].x) * (poly[j].y - poly[i].y);
        return std::fabs(a) * 0.5f;
    }

} // namespace

// Find the face the point p lies in, classify it by how many crossings border
// it (1->R1, 2->R2, 3->R3), and report the beads to collapse for R1/R2. Tests
// the live geometry, so it stays correct after moves rewire the loop.
FaceClass LoopModel::ClassifyFace(Vector2 p) const {
    FaceClass best;
    float bestArea = 1e30f;
    const int nc = (int)crossingBead.size();

    auto toPoly = [&](const std::vector<int>& path) {
        std::vector<Vector2> v; v.reserve(path.size());
        for (int b : path) v.push_back(beads[b].pos);
        return v;
        };
    auto consider = [&](int type, const std::vector<int>& path,
        std::vector<int> doomed, std::vector<int> crossings,
        std::vector<std::vector<int>> triArcs = {}) {
            if (path.size() < 3) return;
            std::vector<Vector2> poly = toPoly(path);
            if (!PointInPoly(p, poly)) return;
            float area = PolyArea(poly);
            if (area < bestArea) {
                bestArea = area; best.type = type;
                best.doomed = std::move(doomed); best.crossings = std::move(crossings);
                best.triArcs = std::move(triArcs);
            }
        };

    // R1: a self-loop arc bounds a 1-crossing face
    for (int a = 0; a < nc; ++a) {
        int A = crossingBead[a];
        if (!beads[A].active) continue;
        std::vector<int> loop = FindSelfLoop(A);
        if (loop.empty()) continue;
        std::vector<int> path = { A };
        path.insert(path.end(), loop.begin(), loop.end());
        std::vector<int> doomed = loop; doomed.push_back(A);
        consider(1, path, doomed, { A });
    }

    // R2: two different-strand arcs bound a 2-crossing face
    for (int a = 0; a < nc; ++a) {
        int A = crossingBead[a]; if (!beads[A].active) continue;
        std::vector<SlotArc> fromA = ArcsFrom(A);
        for (int b = a + 1; b < nc; ++b) {
            int B = crossingBead[b]; if (!beads[B].active) continue;
            std::vector<const SlotArc*> ab;
            for (auto& s : fromA) if (s.endVertex == B) ab.push_back(&s);
            for (size_t x = 0; x < ab.size(); ++x)
                for (size_t y = x + 1; y < ab.size(); ++y) {
                    if (!ValidBigon(A, B, *ab[x], *ab[y])) continue;
                    std::vector<int> path = { A };
                    path.insert(path.end(), ab[x]->beads.begin(), ab[x]->beads.end());
                    path.push_back(B);
                    for (auto it = ab[y]->beads.rbegin(); it != ab[y]->beads.rend(); ++it) path.push_back(*it);
                    std::vector<int> doomed = { A, B };
                    doomed.insert(doomed.end(), ab[x]->beads.begin(), ab[x]->beads.end());
                    doomed.insert(doomed.end(), ab[y]->beads.begin(), ab[y]->beads.end());
                    consider(2, path, doomed, { A, B });
                }
        }
    }

    // R3: three crossings pairwise joined by crossing-free arcs (triangle)
    for (int a = 0; a < nc; ++a) {
        int A = crossingBead[a]; if (!beads[A].active) continue;
        for (int b = a + 1; b < nc; ++b) {
            int B = crossingBead[b]; if (!beads[B].active) continue;
            for (int c = b + 1; c < nc; ++c) {
                int C = crossingBead[c]; if (!beads[C].active) continue;
                auto shortest = [&](int U, int V, std::vector<int>& out) {
                    size_t bestn = (size_t)-1; bool ok = false;
                    for (auto& s : ArcsFrom(U)) if (s.endVertex == V && s.beads.size() < bestn) {
                        out = s.beads; bestn = s.beads.size(); ok = true;
                    }
                    return ok;
                    };
                std::vector<int> ab, bc, ca;
                if (!shortest(A, B, ab) || !shortest(B, C, bc) || !shortest(C, A, ca)) continue;
                std::vector<int> path = { A };
                path.insert(path.end(), ab.begin(), ab.end()); path.push_back(B);
                path.insert(path.end(), bc.begin(), bc.end()); path.push_back(C);
                path.insert(path.end(), ca.begin(), ca.end());
                consider(3, path, {}, { A, B, C }, { ab, bc, ca });
            }
        }
    }
    return best;
}
