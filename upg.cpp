// ============================================================================
// upg.cpp - coordinate utils + region map (build only), polygon -> necklace
// builder, traversal-based render, and the topology-surgery primitives.
// ============================================================================

#include "upg.h"

#include <algorithm>
#include <unordered_map>

// ----------------------------------------------------------------------------
// utilities
// ----------------------------------------------------------------------------

Vector2 GridCoordToCanvasCoord(float x, float y) {
    return { PADDING_X + x * GRID_SIZE, PADDING_Y + y * GRID_SIZE };
}

int GetFaceId(float cx, float cy, UnionFind& uf) {
    int i = (int)std::floor(cx);
    int j = (int)std::floor(cy);
    if (i < 0 || i >= GRID_N || j < 0 || j >= GRID_N) return uf.find(GRID_N * GRID_N);
    return uf.find(j * GRID_N + i);
}

void BuildTopologyMap(const std::vector<Vector2>& pts,
    std::set<std::string>& walls, UnionFind& uf) {
    for (size_t k = 0; k < pts.size(); k++) {
        Vector2 s = pts[k], e = pts[(k + 1) % pts.size()];
        if (s.x == e.x) {
            for (int y = (int)std::min(s.y, e.y); y < (int)std::max(s.y, e.y); y++)
                walls.insert("V," + std::to_string((int)s.x) + "," + std::to_string(y));
        }
        else {
            for (int x = (int)std::min(s.x, e.x); x < (int)std::max(s.x, e.x); x++)
                walls.insert("H," + std::to_string(x) + "," + std::to_string((int)s.y));
        }
    }
    auto noWall = [&](const std::string& k) { return !walls.contains(k); };
    const int OUT = GRID_N * GRID_N;
    const std::string LAST = std::to_string(GRID_N);
    for (int j = 0; j < GRID_N; j++) {
        for (int i = 0; i < GRID_N; i++) {
            int id = j * GRID_N + i;
            if (i < GRID_N - 1 && noWall("V," + std::to_string(i + 1) + "," + std::to_string(j))) uf.unite(id, id + 1);
            if (i == GRID_N - 1 && noWall("V," + LAST + "," + std::to_string(j)))        uf.unite(id, OUT);
            if (i == 0 && noWall("V,0," + std::to_string(j)))                            uf.unite(id, OUT);
            if (j < GRID_N - 1 && noWall("H," + std::to_string(i) + "," + std::to_string(j + 1))) uf.unite(id, id + GRID_N);
            if (j == GRID_N - 1 && noWall("H," + std::to_string(i) + "," + LAST))        uf.unite(id, OUT);
            if (j == 0 && noWall("H," + std::to_string(i) + ",0"))                       uf.unite(id, OUT);
        }
    }
}

// ----------------------------------------------------------------------------
// polygon -> necklace
// ----------------------------------------------------------------------------

namespace {

    std::string PtKey(Vector2 p) {
        return std::to_string((int)std::lround(p.x)) + "," + std::to_string((int)std::lround(p.y));
    }

    std::vector<Vector2> SamplePath(const std::vector<Vector2>& path) {
        std::vector<float> cum = { 0.0f };
        for (size_t i = 0; i + 1 < path.size(); ++i) {
            float len = std::abs(path[i + 1].x - path[i].x) + std::abs(path[i + 1].y - path[i].y);
            cum.push_back(cum.back() + len);
        }
        const float total = cum.back();

        std::vector<Vector2> out;
        size_t seg = 0;
        for (float d = STEP_SIZE; d < total - 1e-4f; d += STEP_SIZE) {
            while (seg + 1 < cum.size() && d > cum[seg + 1]) ++seg;
            const float f = (d - cum[seg]) / (cum[seg + 1] - cum[seg]);
            out.push_back({ path[seg].x + (path[seg + 1].x - path[seg].x) * f,
                            path[seg].y + (path[seg + 1].y - path[seg].y) * f });
        }
        if (out.empty() && path.size() >= 2)
            out.push_back({ (path.front().x + path.back().x) * 0.5f,
                            (path.front().y + path.back().y) * 0.5f });
        return out;
    }

    std::vector<Vector2> AdjacentCells(float gx, float gy) {
        if (!IsGridAligned(gx)) return { {gx, gy - 0.5f}, {gx, gy + 0.5f} };
        if (!IsGridAligned(gy)) return { {gx - 0.5f, gy}, {gx + 0.5f, gy} };
        return { {gx - 0.5f, gy - 0.5f}, {gx + 0.5f, gy - 0.5f},
                 {gx - 0.5f, gy + 0.5f}, {gx + 0.5f, gy + 0.5f} };
    }

} // namespace

LoopModel BuildLoopFromPolygon(const std::vector<Vector2>& pts,
    const std::vector<Vector2>& faceCoords,
    const std::vector<std::pair<int, int>>& pinCells) {

    LoopModel m;
    BuildTopologyMap(pts, m.walls, m.uf);

    const int M = (int)pts.size();

    std::unordered_map<std::string, int> count, crossingOf;
    for (auto& p : pts) count[PtKey(p)]++;
    std::vector<Vector2> crossingGrid;
    for (auto& p : pts) {
        std::string k = PtKey(p);
        if (count[k] >= 2 && !crossingOf.count(k)) {
            crossingOf[k] = (int)crossingGrid.size();
            crossingGrid.push_back(p);
        }
    }
    const int numCrossings = (int)crossingGrid.size();
    auto isCrossing = [&](Vector2 p) { return count[PtKey(p)] >= 2; };
    auto crossIdx = [&](Vector2 p) { return crossingOf[PtKey(p)]; };

    m.crossingBead.resize(numCrossings);
    for (int c = 0; c < numCrossings; ++c) {
        Bead b;
        b.id = (int)m.beads.size();
        b.isVertex = true;
        b.crossing = c;
        b.pos = GridCoordToCanvasCoord(crossingGrid[c].x, crossingGrid[c].y);
        m.crossingBead[c] = (int)m.beads.size();
        m.beads.push_back(b);
    }

    int s0 = 0; while (s0 < M && !isCrossing(pts[s0])) ++s0;
    struct Edge { int startC, endC; std::vector<Vector2> path; int firstBead = -1, lastBead = -1; };
    std::vector<Edge> edges;
    {
        Edge cur; cur.startC = crossIdx(pts[s0]); cur.path.push_back(pts[s0]);
        for (int step = 1; step <= M; ++step) {
            Vector2 p = pts[(s0 + step) % M];
            cur.path.push_back(p);
            if (isCrossing(p)) {
                cur.endC = crossIdx(p);
                edges.push_back(cur);
                Edge ne; ne.startC = crossIdx(p); ne.path.push_back(p); cur = ne;
            }
        }
    }

    std::map<int, std::set<int>> regionSet;
    for (size_t ei = 0; ei < edges.size(); ++ei) {
        Edge& e = edges[ei];
        std::vector<Vector2> grid = SamplePath(e.path);
        const int c = (int)grid.size();
        const int base = (int)m.beads.size();
        for (int i = 0; i < c; ++i) {
            Bead b;
            b.id = base + i;
            b.isVertex = false;
            b.edge = (int)ei;
            b.pos = GridCoordToCanvasCoord(grid[i].x, grid[i].y);
            b.neighbors[0] = (i == 0) ? m.crossingBead[e.startC] : base + i - 1;
            b.neighbors[1] = (i == c - 1) ? m.crossingBead[e.endC] : base + i + 1;
            m.beads.push_back(b);
            for (auto& cell : AdjacentCells(grid[i].x, grid[i].y))
                regionSet[GetFaceId(cell.x, cell.y, m.uf)].insert(base + i);
        }
        e.firstBead = base;
        e.lastBead = base + c - 1;
    }

    std::vector<int> visit(numCrossings, 0);
    const int n = (int)edges.size();
    for (int e = 0; e < n; ++e) {
        const Edge& ce = edges[e];
        const Edge& ne = edges[(e + 1) % n];
        const int c = ce.endC;
        const int vb = m.crossingBead[c];
        const int sIn = 2 * visit[c];
        const int sOut = sIn + 1;
        m.beads[vb].neighborOfSlot[sIn] = ce.lastBead;
        m.beads[vb].neighborOfSlot[sOut] = ne.firstBead;
        m.beads[vb].throughSlot[sIn] = sOut;
        m.beads[vb].throughSlot[sOut] = sIn;
        ++visit[c];
    }

    for (int c = 0; c < numCrossings; ++c) {
        const int vb = m.crossingBead[c];
        const float gx = crossingGrid[c].x, gy = crossingGrid[c].y;
        std::set<int> roots;
        for (float dx : { -0.5f, 0.5f })
            for (float dy : { -0.5f, 0.5f })
                roots.insert(GetFaceId(gx + dx, gy + dy, m.uf));
        for (int r : roots) { regionSet[r].insert(vb); m.regionCrossings[r].insert(vb); }
    }

    for (auto& [root, s] : regionSet)
        m.regionBeads[root].assign(s.begin(), s.end());

    std::set<int> claimed;
    for (auto& [cell, pinIdx] : pinCells) {
        int root = m.uf.find(cell);
        Pin pin;
        pin.id = (int)m.pins.size();
        pin.cell = cell;
        pin.pos = GridCoordToCanvasCoord(faceCoords[pinIdx].x, faceCoords[pinIdx].y);
        pin.active = claimed.insert(root).second;
        m.pins.push_back(pin);
    }

    return m;
}

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
    beads[T2].active = false;
    beads[T2].isVertex = false;
    beads[T2].crossing = -1;
    for (int s = 0; s < 4; ++s) { beads[T2].neighborOfSlot[s] = -1; beads[T2].throughSlot[s] = -1; }
    beads[T2].neighbors[0] = beads[T2].neighbors[1] = -1;
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
    // Degenerate-arc guard. The converging-walk surgery needs each triangle arc to
    // carry two DISTINCT T-junctions (one per end) that meet and merge. A one-bead
    // arc forces both ends onto the same bead, leaving an unmerged junction with a
    // dangling dead slot (later dereferenced as beads[-1]); an empty arc makes one
    // crossing's surgery early-return while the others proceed. So: abort cleanly
    // on an empty arc, and subdivide any arc up to two beads before touching it.
    for (auto& arc : f.triArcs)
        if (arc.empty()) { printf("[R3] degenerate empty arc; move skipped\n"); return; }
    for (auto& arc : f.triArcs)
        while ((int)arc.size() < 2) {
            int b = arc.back();
            int inward = (arc.size() >= 2) ? arc[arc.size() - 2] : -1;
            int out = (beads[b].neighbors[0] == inward) ? beads[b].neighbors[1]
                : beads[b].neighbors[0];
            int m = SplitLink(b, out);
            if (m < 0) break;
            arc.push_back(m);
        }
    for (auto& arc : f.triArcs)
        if ((int)arc.size() < 2) { printf("[R3] arc too short to subdivide; move skipped\n"); return; }

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
        schedule(lastMerge + R3_WALK_GAP, [this] { RebuildCrossingList(); });
    }
    else {
        mergeArc(f.triArcs[0]);
        mergeArc(f.triArcs[1]);
        mergeArc(f.triArcs[2]);
        RebuildCrossingList();
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

// ----------------------------------------------------------------------------
// rendering
// ----------------------------------------------------------------------------

void Draw_Debug_Screen(const LoopModel& m) {
    Bead last_b, cur_show_b;
    for (auto& b : m.beads) {
        if (!b.active) continue;
        Vector2 pos = b.pos;
        std::string id = std::to_string(b.id);

        DrawText(id.c_str(), (int)(pos.x - MeasureText(id.c_str(), 30) / 2.0f),
            (int)(pos.y - 5), 30, RED);

        if (b.driven == true) {
            cur_show_b = b;

            std::string drag_info[11];
            drag_info[0] = "Dragging Bead Info";
            drag_info[1] = "ID: " + std::to_string(cur_show_b.id);
            drag_info[2] = "Pos: " + std::to_string((int)cur_show_b.pos.x) +
                " | " + std::to_string((int)cur_show_b.pos.y);

            drag_info[3] = "Force: " + std::to_string((int)cur_show_b.force.x) +
                " | " + std::to_string((int)cur_show_b.force.y);

            drag_info[4] = "IsVertex: " + std::to_string(cur_show_b.isVertex);
            drag_info[5] = "Driven: " + std::to_string(cur_show_b.driven);
            drag_info[6] = "Edge: " + std::to_string(cur_show_b.edge);
            drag_info[7] = "Neighbors: {" + std::to_string(cur_show_b.neighbors[0]) +
                " | " + std::to_string(cur_show_b.neighbors[1]) + "}";
            drag_info[8] = "Crossing: " + std::to_string(cur_show_b.crossing);
            drag_info[9] = "NeighborOfSlot: {" + std::to_string(cur_show_b.neighborOfSlot[0]) +
                " | " + std::to_string(cur_show_b.neighborOfSlot[1]) +
                " | " + std::to_string(cur_show_b.neighborOfSlot[2]) +
                " | " + std::to_string(cur_show_b.neighborOfSlot[3]) + "}";

            drag_info[10] = "ThroughSlot: {" + std::to_string(cur_show_b.throughSlot[0]) +
                " | " + std::to_string(cur_show_b.throughSlot[1]) +
                " | " + std::to_string(cur_show_b.throughSlot[2]) +
                " | " + std::to_string(cur_show_b.throughSlot[3]) + "}";

            for (int k = 0; k < 11; k++) {
                DrawText(drag_info[k].c_str(), (int)50,
                    (int)SCREEN_HEIGHT - (40 * 10 - (k * 30)), 20, YELLOW);

            }
        }



    }
}

void DrawWireSpline(const LoopModel& m) {
    const int N = (int)m.beads.size();
    std::vector<char> seen(N, 0);                 // arc beads already placed in a rope
    std::set<long long> usedEnd;                  // open-rope endpoints (vertex, inNeighbor)
    auto endKey = [](int v, int n) { return ((long long)v << 32) | (unsigned)n; };

    auto drawRope = [&](const std::vector<Vector2>& pts, bool closed) {
        const int nn = (int)pts.size();
        if (nn < 2) return;
        for (auto& p : pts) DrawRing(p, RADIUS - 3.0f, RADIUS + 3.0f, 0, 360, 1, BLACK);
        if (nn == 2) { DrawLineEx(pts[0], pts[1], 3.0f, WHITE); return; }
        std::vector<Vector2> c; c.reserve(nn + 3);
        if (closed) {            // wrap the Catmull-Rom phantom points around the loop
            c.push_back(pts.back()); c.insert(c.end(), pts.begin(), pts.end());
            c.push_back(pts[0]); c.push_back(pts[1]);
        }
        else {                   // open strand: duplicate the endpoints, do NOT wrap
            c.push_back(pts.front()); c.insert(c.end(), pts.begin(), pts.end());
            c.push_back(pts.back());
        }
        DrawSplineCatmullRom(c.data(), (int)c.size(), 10.0f, WHITE);
        };

    // Walk a rope starting at vertex v0, stepping first to neighbor n0. Follows
    // Through() at vertices and neighbors at arc beads, stopping at a dead end
    // (Through -> -1, e.g. a T-junction stem or a severed leg) or on return to
    // the start. A 4-valent crossing is passed straight through; a 3-valent
    // T-junction passes its bar through and dead-ends its stem.
    auto walkRope = [&](int v0, int n0, bool& closed) {
        std::vector<Vector2> pts{ m.beads[v0].pos };
        int prev = v0, cur = n0, guard = 0;
        closed = false;
        while (cur >= 0 && ++guard < 200000) {
            pts.push_back(m.beads[cur].pos);
            int nx;
            if (m.beads[cur].isVertex) nx = m.Through(cur, prev);
            else { seen[cur] = 1; nx = (m.beads[cur].neighbors[0] == prev) ? m.beads[cur].neighbors[1] : m.beads[cur].neighbors[0]; }
            if (nx < 0) { if (m.beads[cur].isVertex) usedEnd.insert(endKey(cur, prev)); break; }
            if (cur == v0 && nx == n0) { closed = true; break; }
            prev = cur; cur = nx;
        }
        return pts;
        };

    // 1. Open ropes: a dead end is a vertex slot whose strand can't continue
    //    (Through -> -1). Each such slot begins one open rope.
    for (int v = 0; v < N; ++v) {
        if (!m.beads[v].active || !m.beads[v].isVertex) continue;
        for (int s = 0; s < 4; ++s) {
            int n = m.beads[v].neighborOfSlot[s];
            if (n < 0 || !m.beads[n].active) continue;
            if (m.Through(v, n) != -1) continue;
            if (usedEnd.count(endKey(v, n))) continue;
            usedEnd.insert(endKey(v, n));
            bool closed; auto pts = walkRope(v, n, closed);
            drawRope(pts, false);
        }
    }

    // 2. Closed ropes: any arc bead not yet covered lies on a pure closed loop.
    for (int start = 0; start < N; ++start) {
        if (m.beads[start].isVertex || !m.beads[start].active || seen[start]) continue;
        std::vector<Vector2> pts;
        int prev = m.beads[start].neighbors[0], cur = start, guard = 0;
        bool closed = true;
        do {
            pts.push_back(m.beads[cur].pos);
            if (!m.beads[cur].isVertex) seen[cur] = 1;
            int nx = m.beads[cur].isVertex ? m.Through(cur, prev)
                : (m.beads[cur].neighbors[0] == prev ? m.beads[cur].neighbors[1] : m.beads[cur].neighbors[0]);
            if (nx < 0) { closed = false; break; }
            prev = cur; cur = nx;
        } while (cur != start && ++guard < 100000);
        drawRope(pts, closed);
    }
}

void DrawPins(const LoopModel& m) {
    for (auto& p : m.pins) {
        if (!p.active) continue;
        DrawCircleV(p.pos, RADIUS, RED);
        DrawCircleLines((int)p.pos.x, (int)p.pos.y, RADIUS, RED);
        std::string lbl = std::to_string(p.id);
        DrawText(lbl.c_str(), (int)(p.pos.x - MeasureText(lbl.c_str(), 10) / 2.0f),
            (int)(p.pos.y - 5), 20, WHITE);
    }
}

void DrawInterface(float tension, bool autoRunning) {
    std::string t = "Total Tension: " + std::to_string(tension).substr(0, 5);
    DrawText(t.c_str(), SCREEN_WIDTH / 2 - MeasureText(t.c_str(), 20) / 2, 5, 20, RED);

    Color       btn = autoRunning ? RED : BLUE;
    const char* txt = autoRunning ? "Stop Auto-Tightening" : "Start Auto-Tightening";
    DrawRectangle(SCREEN_WIDTH / 2 - 110, 32, 220, 36, btn);
    DrawText(txt, SCREEN_WIDTH / 2 - MeasureText(txt, 16) / 2, 42, 16, WHITE);

    const char* hint = "Click pin: R1/R2 remove, or R3 flip (3-crossing face)  |  Drag: move bead  |  Button: auto-tighten";
    DrawText(hint, SCREEN_WIDTH / 2 - MeasureText(hint, 11) / 2, SCREEN_HEIGHT - 20, 11, GRAY);
}