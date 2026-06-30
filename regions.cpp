// ============================================================================
// regions.cpp - see regions.h.
// ============================================================================

#include "regions.h"
#include <unordered_map>
#include <cmath>
#include <algorithm>

// ---- small geometry helpers (local; upg_faces.cpp's are file-private) -------

static bool PtInPoly(Vector2 p, const std::vector<Vector2>& poly) {
    bool in = false; size_t n = poly.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++)
        if (((poly[i].y > p.y) != (poly[j].y > p.y)) &&
            (p.x < (poly[j].x - poly[i].x) * (p.y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x))
            in = !in;
    return in;
}

static float PolyArea(const std::vector<Vector2>& poly) {
    float a = 0.0f; size_t n = poly.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++)
        a += (poly[j].x + poly[i].x) * (poly[j].y - poly[i].y);
    return std::fabs(a) * 0.5f;
}

static float Dist2(Vector2 a, Vector2 b) {
    float dx = a.x - b.x, dy = a.y - b.y; return dx * dx + dy * dy;
}

static float SegDist2(Vector2 p, Vector2 a, Vector2 b) {
    float vx = b.x - a.x, vy = b.y - a.y, wx = p.x - a.x, wy = p.y - a.y;
    float c1 = vx * wx + vy * wy;
    if (c1 <= 0) return Dist2(p, a);
    float c2 = vx * vx + vy * vy;
    if (c2 <= c1) return Dist2(p, b);
    float t = c1 / c2;
    Vector2 proj = { a.x + t * vx, a.y + t * vy };
    return Dist2(p, proj);
}

// active neighbours of a bead (4 for a crossing, 2 for an edge bead)
static std::vector<int> ActiveNeighbors(const LoopModel& m, int u) {
    std::vector<int> ns;
    const Bead& b = m.beads[u];
    if (b.isVertex) {
        for (int s = 0; s < 4; ++s) {
            int v = b.neighborOfSlot[s];
            if (v >= 0 && v < (int)m.beads.size() && m.beads[v].active) ns.push_back(v);
        }
    }
    else {
        for (int k = 0; k < 2; ++k) {
            int v = b.neighbors[k];
            if (v >= 0 && v < (int)m.beads.size() && m.beads[v].active) ns.push_back(v);
        }
    }
    return ns;
}

// ---- region enumeration -----------------------------------------------------

std::vector<Region> EnumerateRegions(const LoopModel& m) {
    const int N = (int)m.beads.size();

    // Arriving at v from u, the next edge bounding the same face is the one
    // turning the least counter-clockwise from the back-edge (v->u). At an edge
    // bead this just continues; at a crossing it picks the correct sector.
    auto nextNode = [&](int u, int v) -> int {
        const Vector2& pv = m.beads[v].pos;
        float a0 = std::atan2(m.beads[u].pos.y - pv.y, m.beads[u].pos.x - pv.x);
        int best = -1; float bestDelta = 1e30f;
        for (int w : ActiveNeighbors(m, v)) {
            float aw = std::atan2(m.beads[w].pos.y - pv.y, m.beads[w].pos.x - pv.x);
            float d = aw - a0;
            while (d <= 1e-5f) d += 2.0f * PI;     // strictly > 0; the back-edge sorts last
            while (d > 2.0f * PI) d -= 2.0f * PI;
            if (d < bestDelta) { bestDelta = d; best = w; }
        }
        return best;
    };

    std::unordered_map<long long, bool> seen;
    auto key = [&](int u, int v) { return (long long)u * N + v; };

    std::vector<Region> faces;
    for (int u = 0; u < N; ++u) {
        if (!m.beads[u].active) continue;
        for (int v : ActiveNeighbors(m, u)) {
            if (seen.count(key(u, v))) continue;
            std::vector<int> boundary;
            int cu = u, cv = v, guard = 0;
            while (!seen.count(key(cu, cv)) && guard++ < 4 * N + 8) {
                seen[key(cu, cv)] = true;
                boundary.push_back(cu);
                int w = nextNode(cu, cv);
                if (w < 0) break;
                cu = cv; cv = w;
            }
            if (boundary.size() >= 2) { Region r; r.boundary = std::move(boundary); faces.push_back(std::move(r)); }
        }
    }

    int outer = -1; float maxA = -1.0f;
    for (size_t i = 0; i < faces.size(); ++i) {
        std::vector<Vector2> poly; poly.reserve(faces[i].boundary.size());
        for (int b : faces[i].boundary) poly.push_back(m.beads[b].pos);
        faces[i].area = PolyArea(poly);
        if (faces[i].area > maxA) { maxA = faces[i].area; outer = (int)i; }
    }
    if (outer >= 0) faces[outer].bounded = false;
    return faces;
}

// ---- pole of inaccessibility ------------------------------------------------

Vector2 PoleOfInaccessibility(const std::vector<Vector2>& poly) {
    float minx = 1e30f, miny = 1e30f, maxx = -1e30f, maxy = -1e30f;
    for (auto& p : poly) {
        minx = std::min(minx, p.x); maxx = std::max(maxx, p.x);
        miny = std::min(miny, p.y); maxy = std::max(maxy, p.y);
    }
    Vector2 best = { (minx + maxx) * 0.5f, (miny + maxy) * 0.5f };
    float bestD = -1.0f;
    float x0 = minx, y0 = miny, x1 = maxx, y1 = maxy;
    const int G = 12;
    for (int iter = 0; iter < 4; ++iter) {
        for (int i = 0; i <= G; ++i)
            for (int j = 0; j <= G; ++j) {
                Vector2 q = { x0 + (x1 - x0) * i / G, y0 + (y1 - y0) * j / G };
                if (!PtInPoly(q, poly)) continue;
                float d = 1e30f; size_t n = poly.size();
                for (size_t a = 0, b = n - 1; a < n; b = a++) d = std::min(d, SegDist2(q, poly[b], poly[a]));
                if (d > bestD) { bestD = d; best = q; }
            }
        float rx = (x1 - x0) / G, ry = (y1 - y0) / G;
        x0 = best.x - rx; x1 = best.x + rx; y0 = best.y - ry; y1 = best.y + ry;
    }
    return best;
}

// ---- reconcile --------------------------------------------------------------

int ReconcilePins(LoopModel& m, bool spawnMissing) {
    auto regions = EnumerateRegions(m);

    struct B { std::vector<Vector2> poly; Vector2 cent; };
    std::vector<B> bs;
    for (auto& rg : regions) {
        if (!rg.bounded || rg.boundary.size() < 3) continue;
        B x; Vector2 c = { 0.0f, 0.0f };
        x.poly.reserve(rg.boundary.size());
        for (int b : rg.boundary) { Vector2 p = m.beads[b].pos; x.poly.push_back(p); c.x += p.x; c.y += p.y; }
        c.x /= (float)rg.boundary.size(); c.y /= (float)rg.boundary.size();
        x.cent = c; bs.push_back(std::move(x));
    }
    const int R = (int)bs.size();
    if (R == 0) return 0;

    // bucket each active pin into the region that contains it
    std::vector<std::vector<int>> in(R);
    std::vector<int> homeless;
    for (int pi = 0; pi < (int)m.pins.size(); ++pi) {
        if (!m.pins[pi].active) continue;
        int found = -1;
        for (int ri = 0; ri < R; ++ri) if (PtInPoly(m.pins[pi].pos, bs[ri].poly)) { found = ri; break; }
        if (found >= 0) in[found].push_back(pi);
        else            homeless.push_back(pi);
    }
    // escapees: attach to the nearest region (they will be snapped inside below)
    for (int pi : homeless) {
        int best = -1; float bd = 1e30f;
        for (int ri = 0; ri < R; ++ri) { float d = Dist2(m.pins[pi].pos, bs[ri].cent); if (d < bd) { bd = d; best = ri; } }
        if (best >= 0) in[best].push_back(pi);
        else m.pins[pi].active = false;
    }

    // one pin per region: keep the most central, snap it inside, drop the rest,
    // spawn for empties.
    for (int ri = 0; ri < R; ++ri) {
        auto& pins = in[ri];
        if (pins.empty()) {
            if (spawnMissing) {
                Pin np;
                np.id = (int)m.pins.size();
                np.pos = PoleOfInaccessibility(bs[ri].poly);
                np.active = true;
                np.cell = ri;
                m.pins.push_back(np);
            }
            continue;
        }
        int keep = pins[0]; float bd = Dist2(m.pins[pins[0]].pos, bs[ri].cent);
        for (int pi : pins) { float d = Dist2(m.pins[pi].pos, bs[ri].cent); if (d < bd) { bd = d; keep = pi; } }
        if (!PtInPoly(m.pins[keep].pos, bs[ri].poly))
            m.pins[keep].pos = PoleOfInaccessibility(bs[ri].poly);
        m.pins[keep].cell = ri;
        for (int pi : pins) if (pi != keep) m.pins[pi].active = false;
    }

    int cnt = 0; for (auto& p : m.pins) if (p.active) ++cnt;
    return cnt;
}
