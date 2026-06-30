// ============================================================================
// upg_build.cpp - construction & coordinate utilities.
//
// Turns geometry into a LoopModel "bead necklace":
//   * GridCoordToCanvasCoord  - grid-cell coords -> screen pixels.
//   * BuildTopologyMap / GetFaceId
//                             - rasterise the polygon's edges into grid "walls"
//                               and union-find the cells they separate, so each
//                               grid cell knows which planar face it belongs to.
//                               Used only at build time to seed the region maps
//                               and pins; the live game later finds faces
//                               geometrically (see upg_faces.cpp / regions.cpp).
//   * BuildLoopFromPolygon    - trace a rectilinear grid polygon (corners +
//                               doubled crossing vertices) into crossing beads
//                               and edge-bead chains, wiring each crossing's
//                               four slots and populating the region/pin maps.
//
// This is one of two model builders; sigma_import.cpp is the combinatorial-map
// alternative. See upg.h for the bead / slot data model.
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
