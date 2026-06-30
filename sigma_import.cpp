// ============================================================================
// sigma_import.cpp - implementation of BuildLoopFromSigma.
// ============================================================================

#include "sigma_import.h"
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <stdexcept>

// Number of interior edge beads for a grid-space segment of the given length,
// spaced at the physics rest length; at least one so every slot links a bead.
static int EdgeBeadCount(float gridLen) {
    const float restGrid = REST_LENGTH / GRID_SIZE;   // rest spacing in grid units
    int n = (int)std::lround(gridLen / restGrid) - 1;
    return n < 1 ? 1 : n;
}

LoopModel BuildLoopFromSigma(const std::vector<std::array<int, 4>>& sigma,
    const std::vector<Vector2>& crossingGrid) {

    LoopModel m;
    const int V = (int)sigma.size();
    assert((int)crossingGrid.size() == V && "one grid coord per crossing");

    // Locate every dart: dartLoc[d] = (vertex, slot 0..3). Also recover E and
    // assert the dart set is exactly +-1..+-E (standard epsilon).
    std::unordered_map<int, std::pair<int, int>> dartLoc;
    int maxAbs = 0;
    for (int v = 0; v < V; ++v)
        for (int s = 0; s < 4; ++s) {
            int d = sigma[v][s];
            assert(d != 0 && !dartLoc.count(d) && "dart appears twice / is zero");
            dartLoc[d] = { v, s };
            maxAbs = std::max(maxAbs, std::abs(d));
        }
    const int E = maxAbs;
    assert((int)dartLoc.size() == 2 * E && "dart set is not +-1..+-E");

    // Crossing beads, in sigma order. throughSlot[s] = (s+2)%4 = sigma^2.
    m.crossingBead.resize(V);
    for (int v = 0; v < V; ++v) {
        Bead b;
        b.id = (int)m.beads.size();
        b.isVertex = true;
        b.crossing = v;
        b.pos = GridCoordToCanvasCoord(crossingGrid[v].x, crossingGrid[v].y);
        for (int s = 0; s < 4; ++s) b.throughSlot[s] = (s + 2) % 4;
        m.crossingBead[v] = (int)m.beads.size();
        m.beads.push_back(b);
    }

    // One edge bead chain per edge k, from the tail vertex (holding -k) to the
    // head vertex (holding +k). Endpoints link to the crossing beads; the
    // incident slot at each end is wired to the nearest edge bead.
    for (int k = 1; k <= E; ++k) {
        auto [vH, sH] = dartLoc[+k];
        auto [vT, sT] = dartLoc[-k];
        Vector2 gT = crossingGrid[vT], gH = crossingGrid[vH];
        float gridLen = std::hypot(gH.x - gT.x, gH.y - gT.y);
        int cnt = EdgeBeadCount(gridLen);

        const int base = (int)m.beads.size();
        for (int i = 0; i < cnt; ++i) {
            float t = (cnt == 1) ? 0.5f : (float)(i + 1) / (float)(cnt + 1);
            Vector2 g = { gT.x + (gH.x - gT.x) * t, gT.y + (gH.y - gT.y) * t };
            Bead b;
            b.id = base + i;
            b.isVertex = false;
            b.edge = k;
            b.pos = GridCoordToCanvasCoord(g.x, g.y);
            b.neighbors[0] = (i == 0) ? m.crossingBead[vT] : base + i - 1;
            b.neighbors[1] = (i == cnt - 1) ? m.crossingBead[vH] : base + i + 1;
            m.beads.push_back(b);
        }
        m.beads[m.crossingBead[vT]].neighborOfSlot[sT] = base;            // chain front
        m.beads[m.crossingBead[vH]].neighborOfSlot[sH] = base + cnt - 1;  // chain back
    }

    return m;
}

std::vector<Vector2> RingLayout(int numCrossings, float radiusGrid, Vector2 centerGrid) {
    std::vector<Vector2> out(numCrossings);
    for (int v = 0; v < numCrossings; ++v) {
        float a = 2.0f * PI * (float)v / (float)numCrossings;
        out[v] = { centerGrid.x + radiusGrid * std::cos(a),
                   centerGrid.y + radiusGrid * std::sin(a) };
    }
    return out;
}

int SeedPinsFromCells(LoopModel& m) {
    m.pins.clear();
    const int outside = m.uf.find(GRID_N * GRID_N);

    // group every interior grid cell by the region (union-find root) it lies in
    std::map<int, std::vector<std::pair<int, int>>> cells;
    for (int j = 0; j < GRID_N; ++j)
        for (int i = 0; i < GRID_N; ++i) {
            int r = m.uf.find(j * GRID_N + i);
            if (r != outside) cells[r].push_back({ i, j });
        }

    // pin each region at the member cell nearest its cell-centroid (interior)
    for (auto& [root, cs] : cells) {
        float cx = 0.0f, cy = 0.0f;
        for (auto& [i, j] : cs) { cx += i + 0.5f; cy += j + 0.5f; }
        cx /= (float)cs.size(); cy /= (float)cs.size();
        std::pair<int, int> best = cs.front(); float bestD = 1e30f;
        for (auto& [i, j] : cs) {
            float dx = i + 0.5f - cx, dy = j + 0.5f - cy, d = dx * dx + dy * dy;
            if (d < bestD) { bestD = d; best = { i, j }; }
        }
        Pin pin;
        pin.id = (int)m.pins.size();
        pin.pos = GridCoordToCanvasCoord(best.first + 0.5f, best.second + 0.5f);
        pin.active = true;
        pin.cell = best.second * GRID_N + best.first;
        m.pins.push_back(pin);
    }
    return (int)m.pins.size();
}

int BuildFacesFromPhi(LoopModel& m, const std::vector<std::array<int, 4>>& sigma) {
    const int V = (int)sigma.size();

    // dart -> (vertex, slot), and edge label -> its bead chain.
    std::unordered_map<int, std::pair<int, int>> dartLoc;
    for (int v = 0; v < V; ++v)
        for (int s = 0; s < 4; ++s) dartLoc[sigma[v][s]] = { v, s };
    std::map<int, std::vector<int>> edgeBeads;
    for (const Bead& b : m.beads)
        if (!b.isVertex && b.edge >= 1) edgeBeads[b.edge].push_back(b.id);

    // phi(d) = sigma^{-1}(eps(d)) = (sigma*eps)^{-1}.
    auto phi = [&](int d) {
        auto [v, s] = dartLoc[-d];
        return sigma[v][(s + 3) % 4];
        };

    m.regionBeads.clear();
    m.regionCrossings.clear();
    std::map<int, std::set<int>> beadSet;

    std::set<int> seen;
    int face = 0;
    for (auto& [d0, loc] : dartLoc) {
        if (seen.count(d0)) continue;
        for (int d = d0; !seen.count(d); d = phi(d)) {
            seen.insert(d);
            int v = dartLoc[d].first;
            int cb = m.crossingBead[v];
            m.regionCrossings[face].insert(cb);
            beadSet[face].insert(cb);
            for (int eb : edgeBeads[std::abs(d)]) beadSet[face].insert(eb);
        }
        ++face;
    }
    for (auto& [f, s] : beadSet) m.regionBeads[f].assign(s.begin(), s.end());
    return face;
}

int SeedPinsPerFace(LoopModel& m, const std::vector<std::array<int, 4>>& sigma) {
    const int V = (int)sigma.size();
    std::unordered_map<int, std::pair<int, int>> dartLoc;
    for (int v = 0; v < V; ++v)
        for (int s = 0; s < 4; ++s) dartLoc[sigma[v][s]] = { v, s };
    auto phi = [&](int d) { auto [v, s] = dartLoc[-d]; return sigma[v][(s + 3) % 4]; };

    // ordered corner crossings of each face (one phi orbit = one face boundary)
    std::set<int> seen;
    std::vector<std::vector<int>> faceCorners;
    for (auto& [d0, l] : dartLoc) {
        if (seen.count(d0)) continue;
        std::vector<int> vs;
        for (int d = d0; !seen.count(d); d = phi(d)) {
            seen.insert(d);
            int vv = dartLoc[d].first;
            if (vs.empty() || vs.back() != vv) vs.push_back(vv);
        }
        if (vs.size() > 1 && vs.front() == vs.back()) vs.pop_back();
        faceCorners.push_back(vs);
    }

    auto cpos = [&](int v) { return m.beads[m.crossingBead[v]].pos; };
    auto area = [&](const std::vector<int>& f) {
        float a = 0.0f; size_t n = f.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++) {
            Vector2 p = cpos(f[i]), q = cpos(f[j]);
            a += (q.x + p.x) * (q.y - p.y);
        }
        return std::fabs(a) * 0.5f;
        };

    // The outer (unbounded) face is the one whose corner polygon has the largest
    // area -- it traces the convex hull. It gets no pin (its boundary-centroid
    // would land in the middle of the diagram, not in the exterior region).
    int outer = -1; float maxA = -1.0f;
    for (size_t f = 0; f < faceCorners.size(); ++f) {
        float A = area(faceCorners[f]);
        if (A > maxA) { maxA = A; outer = (int)f; }
    }

    // Pin each bounded face at its corner centroid. Tutte inner faces are convex,
    // so the centroid is guaranteed inside the region; SettlePins refines it.
    m.pins.clear();
    for (size_t f = 0; f < faceCorners.size(); ++f) {
        if ((int)f == outer || faceCorners[f].empty()) continue;
        Vector2 c = { 0.0f, 0.0f };
        for (int v : faceCorners[f]) { Vector2 p = cpos(v); c.x += p.x; c.y += p.y; }
        c.x /= (float)faceCorners[f].size(); c.y /= (float)faceCorners[f].size();
        Pin pin;
        pin.id = (int)m.pins.size();
        pin.pos = c;
        pin.active = true;
        pin.cell = (int)f;
        m.pins.push_back(pin);
    }
    return (int)m.pins.size();
}

std::vector<Vector2> TutteLayout(const std::vector<std::array<int, 4>>& sigma,
    Vector2 centerGrid, float radiusGrid) {
    const int V = (int)sigma.size();

    std::unordered_map<int, std::pair<int, int>> dartLoc;
    for (int v = 0; v < V; ++v)
        for (int s = 0; s < 4; ++s) dartLoc[sigma[v][s]] = { v, s };

    // neighbour vertices of each crossing (with multiplicity), across its darts
    std::vector<std::vector<int>> nbr(V);
    for (int v = 0; v < V; ++v)
        for (int s = 0; s < 4; ++s)
            nbr[v].push_back(dartLoc[-sigma[v][s]].first);

    auto phi = [&](int d) { auto [v, s] = dartLoc[-d]; return sigma[v][(s + 3) % 4]; };

    // outer face = the longest phi orbit; collect its distinct corner vertices
    std::set<int> seen; std::vector<int> outer; size_t bestLen = 0;
    for (auto& [d0, l] : dartLoc) {
        if (seen.count(d0)) continue;
        std::vector<int> orb, verts;
        for (int d = d0; !seen.count(d); d = phi(d)) {
            seen.insert(d); orb.push_back(d);
            int vv = dartLoc[d].first;
            if (std::find(verts.begin(), verts.end(), vv) == verts.end()) verts.push_back(vv);
        }
        if (orb.size() > bestLen) { bestLen = orb.size(); outer = verts; }
    }

    std::vector<Vector2> pos(V, { 0.0f, 0.0f });
    if (outer.size() < 3)                      // degenerate: fall back to a ring
        return RingLayout(V, radiusGrid, centerGrid);

    // pin the outer face to a regular convex polygon, relax the rest
    std::vector<char> fixed(V, 0);
    for (size_t i = 0; i < outer.size(); ++i) {
        float a = 2.0f * PI * (float)i / (float)outer.size();
        pos[outer[i]] = { std::cos(a), std::sin(a) };
        fixed[outer[i]] = 1;
    }
    for (int it = 0; it < 2000; ++it)
        for (int v = 0; v < V; ++v) {
            if (fixed[v]) continue;
            Vector2 s = { 0.0f, 0.0f };
            for (int u : nbr[v]) { s.x += pos[u].x; s.y += pos[u].y; }
            float k = (float)nbr[v].size();
            pos[v] = { s.x / k, s.y / k };
        }

    std::vector<Vector2> out(V);
    for (int v = 0; v < V; ++v)
        out[v] = { centerGrid.x + radiusGrid * pos[v].x,
                   centerGrid.y + radiusGrid * pos[v].y };
    return out;
}
// ============================================================================
// Orthogonal-layout builder (folded in from sigma_to_grid.cpp + upg_sigma_ortho.cpp).
//   sigma -> Multiloop (validate) -> ComputeGridLayout -> SigmaGridPolygon
//   -> BuildLoopFromPolygon -> LoopModel
// ============================================================================

namespace {
    inline bool SamePt(const GridPoint& a, const GridPoint& b) { return a.x == b.x && a.y == b.y; }
}

std::vector<GridPoint> SigmaGridPolygon(const Multiloop& loop,
                                        const GridLayout& layout,
                                        int startEdge) {
    if (loop.numStrands != 1)
        throw std::runtime_error("SigmaGridPolygon: diagram has multiple components; "
                                 "BuildLoopFromPolygon takes a single closed polygon");

    // 1. Edge order along the strand: walk rho (straight-ahead) once around.
    std::vector<int> edgeOrder;
    {
        int start = startEdge, h = start, guard = 0;
        do { edgeOrder.push_back(std::abs(h)); h = loop.rho.Apply(h); }
        while (h != start && ++guard < 4 * loop.numEdges + 8);
    }
    if ((int)edgeOrder.size() != loop.numEdges)
        throw std::runtime_error("SigmaGridPolygon: strand did not cover every edge once");

    auto endpoints = [&](int k, GridPoint& a, GridPoint& b) {
        const std::vector<GridPoint>& p = layout.edgePaths.at(k - 1);
        a = p.front(); b = p.back();
        };

    // 2. Starting pen = the crossing shared by the last and first edges.
    GridPoint penPos;
    {
        GridPoint a0, b0, aL, bL;
        endpoints(edgeOrder.front(), a0, b0);
        endpoints(edgeOrder.back(), aL, bL);
        if      (SamePt(a0, aL) || SamePt(a0, bL)) penPos = a0;
        else if (SamePt(b0, aL) || SamePt(b0, bL)) penPos = b0;
        else throw std::runtime_error("SigmaGridPolygon: first/last edges are not adjacent");
    }

    // 3. Stitch, orienting each edge so it starts at the pen (geometry decides
    //    direction, so a path stored either way is handled).
    std::vector<GridPoint> poly;
    auto append = [&](const GridPoint& q) {
        if (poly.empty() || !SamePt(poly.back(), q)) poly.push_back(q);
        };
    for (int k : edgeOrder) {
        const std::vector<GridPoint>& path = layout.edgePaths.at(k - 1);
        if (SamePt(path.front(), penPos)) {
            for (const GridPoint& q : path) append(q);
            penPos = path.back();
        }
        else if (SamePt(path.back(), penPos)) {
            for (auto it = path.rbegin(); it != path.rend(); ++it) append(*it);
            penPos = path.front();
        }
        else {
            throw std::runtime_error("SigmaGridPolygon: edge does not continue from pen "
                                     "(discontinuous strand or bad layout)");
        }
    }
    // 4. Drop the wrap-around duplicate (front repeated as back).
    if (poly.size() >= 2 && SamePt(poly.front(), poly.back())) poly.pop_back();
    return poly;
}

LoopModel BuildLoopFromSigmaOrthogonal(const Permutation& sigma, int exteriorFace) {
    const std::string err = ValidateMultiloopSigma(sigma);
    if (!err.empty()) throw std::runtime_error(err);

    Multiloop loop(sigma);

    // Default exterior = the face bounded by the most half-edges (the outer
    // boundary in a typical drawing); callers may override.
    if (exteriorFace < 0) {
        exteriorFace = 0;
        size_t best = 0;
        for (int f = 0; f < (int)loop.faces.size(); ++f)
            if (loop.faces[f].size() > best) { best = loop.faces[f].size(); exteriorFace = f; }
    }

    GridLayout layout = ComputeGridLayout(loop, exteriorFace);
    std::vector<GridPoint> gp = SigmaGridPolygon(loop, layout);

    // Centre the integer drawing on the GRID_N board and flip y (orthogonal is
    // y-up; UPG's canvas is y-down). Centring is cosmetic: the spring physics
    // relaxes to rest-length spacing and re-centres the centroid, so a layout
    // larger than GRID_N still builds correctly. Coords stay integers, so the
    // crossing detection in BuildLoopFromPolygon (which rounds to int) is exact.
    const float offX = (GRID_N - layout.width)  / 2.0f;
    const float offY = (GRID_N - layout.height) / 2.0f;

    std::vector<Vector2> pts;
    pts.reserve(gp.size());
    for (const GridPoint& p : gp)
        pts.push_back({ (float)p.x + offX, (float)(layout.height - p.y) + offY });

    // Empty pinCells: Init(LoopModel) seeds pins geometrically via ReconcilePins.
    return BuildLoopFromPolygon(pts, {}, {});
}

LoopModel BuildLoopFromSigmaString(const std::string& cycleNotation, std::string& errorOut) {
    errorOut.clear();
    try {
        Permutation sigma = Permutation::FromString(cycleNotation);
        if (sigma.MaxAbsLabel() == 0)
            throw std::runtime_error("no crossings parsed - enter cycles like "
                                     "(16,5,-1,-6)(12,1,-13,-2)...");
        return BuildLoopFromSigmaOrthogonal(sigma, -1);
    }
    catch (const std::exception& e) {
        errorOut = e.what();
        return LoopModel{};
    }
}
