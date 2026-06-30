// ============================================================================
// upg_render.cpp - raylib rendering of the live diagram.
//
//   * DrawWireSpline    - trace every rope (open strands that dead-end at a
//                         T-junction stem, plus pure closed loops) and draw it
//                         as a Catmull-Rom spline with a ring at each bead.
//   * DrawPins          - the red region pins.
//   * DrawInterface     - tension readout, auto-tighten button, help text.
//   * Draw_Debug_Screen - per-bead id labels and a dump of the dragged bead's
//                         links (off by default in the render loop).
// ============================================================================

#include "upg.h"

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
        for (auto& p : pts) DrawRing(p, RADIUS - 3.0f, RADIUS + 3.0f, 0, 360, 1, DARKGRAY);
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
        DrawCircleV(p.pos, RADIUS*0.75, RED);
        //DrawCircleLines((int)p.pos.x, (int)p.pos.y, RADIUS, RED);
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