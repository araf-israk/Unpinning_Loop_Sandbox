// ============================================================================
// Sandbox_UPG.cpp - stateful simulation over the index-linked LoopModel.
//
//   * physics  : springs from explicit links; global self-avoidance with a
//                crossing-guard; centroid centring + displacement clamp so the
//                loop stays put and stable.
//   * dispatch : clicking a pin finds the FACE it sits in, counts that face's
//                bordering crossings, and classifies the move (1->R1, 2->R2,
//                3->R3) before acting.
//   * R1 / R2  : smooth - only the arc's edge beads shrink (the crossings keep
//                their outer springs), then the region is spliced out.
// ============================================================================

#include "upg.h"
#include "sigma_import.h"
#include "regions.h"

#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <functional>
#include <string>
#include <unordered_map>

static float Dist2(Vector2 a, Vector2 b) {
    float dx = a.x - b.x, dy = a.y - b.y; return dx * dx + dy * dy;
}

static bool Adjacent(const LoopModel& m, int i, int j) {
    const Bead& a = m.beads[i];
    if (!a.isVertex) { if (a.neighbors[0] == j || a.neighbors[1] == j) return true; }
    else for (int s = 0; s < 4; ++s) if (a.neighborOfSlot[s] == j) return true;
    const Bead& b = m.beads[j];
    if (!b.isVertex) { if (b.neighbors[0] == i || b.neighbors[1] == i) return true; }
    else for (int s = 0; s < 4; ++s) if (b.neighborOfSlot[s] == i) return true;
    return false;
}

struct CollapseJob { std::vector<int> doomed; int frames; Vector2 center; };
struct Candidate { std::vector<int> doomed; Vector2 center; const char* kind; };

struct GameState {
    LoopModel model;
    float     tension = 0.0f;
    int       dragging = -1;
    Vector2   dragOffset = { 0.0f, 0.0f };

    bool                     autoRun = false;
    std::vector<CollapseJob> jobs;
    std::set<int>            pending;

    // ---- frame-delayed callbacks -----------------------------------------
    // A function scheduled to fire some number of frames (Update calls) later.
    // `interval` is the gap between fires; `repeats` is how many fires remain
    // (1 = one-shot, -1 = forever).
    struct DelayedTask {
        int                   framesLeft;
        int                   interval;
        int                   repeats;
        std::function<void()> fn;
    };
    std::vector<DelayedTask> delayedTasks;

    // Run `fn` once, `delay` frames from now. e.g. tickDelay(5, [&]{ ... });
    void tickDelay(int delay, std::function<void()> fn) {
        if (delay < 1) delay = 1;
        delayedTasks.push_back({ delay, delay, 1, std::move(fn) });
    }

    // Run `fn` every `interval` frames, `times` times (-1 = forever). The first
    // fire happens `interval` frames from now. e.g. tickEvery(5, [&]{ ... }, 4);
    void tickEvery(int interval, std::function<void()> fn, int times = -1) {
        if (interval < 1) interval = 1;
        delayedTasks.push_back({ interval, interval, times, std::move(fn) });
    }

    // Advance every scheduled task by one frame and fire the ones that are due.
    // Callbacks are collected first, then run, so a callback may safely schedule
    // more tasks (they start counting next frame) without invalidating anything.
    void TickDelays() {
        std::vector<std::function<void()>> due;
        for (size_t i = 0; i < delayedTasks.size();) {
            DelayedTask& t = delayedTasks[i];
            if (--t.framesLeft <= 0) {
                due.push_back(t.fn);
                if (t.repeats > 0) --t.repeats;
                if (t.repeats == 0) { delayedTasks.erase(delayedTasks.begin() + i); continue; }
                t.framesLeft = t.interval;   // re-arm a repeating task
            }
            ++i;
        }
        for (auto& fn : due) fn();
    }

    void Init(const std::vector<Vector2>& gridPoints,
        const std::vector<Vector2>& faceCoords,
        const std::vector<std::pair<int, int>>& pinCells) {
        model = BuildLoopFromPolygon(gridPoints, faceCoords, pinCells);
    }

    // Take a model built elsewhere (e.g. BuildLoopFromSigma + BuildFacesFromPhi
    // + SeedPinsPerFace).
    void Init(LoopModel m) {
        model = std::move(m);
        ReconcilePins(model, true);   // one pin per bounded region, topologically
    }

    // ---- physics ----------------------------------------------------------

    void AddSpring(int a, int b) {
        if (a < 0 || b < 0 || !model.beads[a].active || !model.beads[b].active) return;
        Vector2& pa = model.beads[a].pos;
        Vector2& pb = model.beads[b].pos;
        float dx = pb.x - pa.x, dy = pb.y - pa.y;
        float d = std::hypot(dx, dy);
        if (d < GEOM_EPS) return;
        float rest = REST_LENGTH;
        if (auto it = model.transitions.find(a); it != model.transitions.end()) rest *= it->second.scale;
        if (auto it = model.transitions.find(b); it != model.transitions.end()) rest *= it->second.scale;
        float disp = d - rest;
        tension += std::abs(disp);
        float mag = disp * SPRING_CONSTANT;
        float fx = (dx / d) * mag, fy = (dy / d) * mag;
        model.beads[a].force.x += fx; model.beads[a].force.y += fy;
        model.beads[b].force.x -= fx; model.beads[b].force.y -= fy;
    }

    void Step() {
        tension = 0.0f;
        for (auto& b : model.beads) b.force = { 0.0f, 0.0f };

        // (R3 move machinery removed - to be rewritten. Collisions always on.)
        const bool collisionsOff = false;

        // 1. springs (orientation-independent: add each edge link exactly once)
        for (int i = 0; i < (int)model.beads.size(); ++i) {
            Bead& b = model.beads[i];
            if (b.isVertex || !b.active) continue;
            for (int k = 0; k < 2; ++k) {
                int nb = b.neighbors[k];
                if (nb >= 0 && (model.beads[nb].isVertex || nb > i)) AddSpring(i, nb);
            }
        }
        // 1b. vertex<->vertex springs: when two crossings/T-junctions link
        // directly (no arc bead between them, e.g. a T-junction's bar meeting a
        // crossing), the arc-side pass above can't add that spring. Add it here.
        for (int i = 0; i < (int)model.beads.size(); ++i) {
            Bead& b = model.beads[i];
            if (!b.isVertex || !b.active) continue;
            for (int s = 0; s < 4; ++s) {
                int n = b.neighborOfSlot[s];
                if (n > i && model.beads[n].isVertex && model.beads[n].active) AddSpring(i, n);
            }
        }

        // 2. (collision is handled by the hard projection in step 6)

        // 3. pins push outward
        if (!collisionsOff) {
            for (auto& pin : model.pins) {
                if (!pin.active) continue;
                for (auto& b : model.beads) {
                    if (!b.active) continue;
                    float dx = b.pos.x - pin.pos.x, dy = b.pos.y - pin.pos.y;
                    float d = std::hypot(dx, dy);
                    if (d > GEOM_EPS && d < PIN_BUFFER_RADIUS) {
                        float mag = (PIN_BUFFER_RADIUS - d) * COLLISION_STIFFNESS;
                        b.force.x += (dx / d) * mag;
                        b.force.y += (dy / d) * mag;
                    }
                }
            }
        }

        // 4. centring: translate the whole loop toward the canvas centre so it
        //    neither drifts nor wanders. Uniform force => no shape distortion.
        Vector2 cen = { 0.0f, 0.0f }; int n = 0;
        for (auto& b : model.beads) if (b.active) { cen.x += b.pos.x; cen.y += b.pos.y; ++n; }
        if (n) {
            cen.x /= n; cen.y /= n;
            Vector2 push = { (SCREEN_WIDTH * 0.5f - cen.x) * CENTERING_STRENGTH,
                             (SCREEN_HEIGHT * 0.5f - cen.y) * CENTERING_STRENGTH };
            for (auto& b : model.beads) if (b.active) { b.force.x += push.x; b.force.y += push.y; }
        }

        // 5. integrate, clamping per-substep displacement for stability
        const float maxD2 = MAX_DISPLACEMENT * MAX_DISPLACEMENT;
        for (auto& b : model.beads) {
            if (!b.active || b.driven) continue;
            float dx = b.force.x * INTEGRATION_STEP, dy = b.force.y * INTEGRATION_STEP;
            float m2 = dx * dx + dy * dy;
            if (m2 > maxD2) { float s = MAX_DISPLACEMENT / std::sqrt(m2); dx *= s; dy *= s; }
            b.pos.x += dx; b.pos.y += dy;
        }

        // 6. hard non-overlap projection (runs last so the substep ends separated)
        if (!collisionsOff) ResolveOverlaps(OVERLAP_RESOLVE_ITERS);
    }

    // Hard non-overlap: push every pair of nearby beads apart until their disks
    // just touch (BEAD_MIN_SEP). A uniform grid (cell = BEAD_MIN_SEP) makes it
    // near-linear -- each bead only tests the 3x3 block of cells around it.
    // Exempt: directly-linked beads (they rest touching) and any two beads both
    // sitting within CROSS_GUARD of the same crossing (that is the crossing,
    // where two strands are meant to coincide). A crossing or dragged bead is an
    // immovable anchor, so its partner takes the whole correction.
    void ResolveOverlaps(int iters) {
        const float minSep = BEAD_MIN_SEP;
        const float minSep2 = minSep * minSep;
        const float guard2 = CROSS_GUARD * CROSS_GUARD;
        const float cell = minSep;

        std::vector<int> crossings;
        for (int c : model.crossingBead)
            if (c >= 0 && model.beads[c].active) crossings.push_back(c);
        auto sharesCrossing = [&](int i, int j) -> bool {
            for (int c : crossings) {
                Vector2 cp = model.beads[c].pos;
                if (Dist2(model.beads[i].pos, cp) <= guard2 &&
                    Dist2(model.beads[j].pos, cp) <= guard2) return true;
            }
            return false;
            };
        auto key = [](int ix, int iy) {
            return ((int64_t)ix << 32) ^ (int64_t)(uint32_t)iy;
            };

        for (int it = 0; it < iters; ++it) {
            std::unordered_map<int64_t, std::vector<int>> grid;
            grid.reserve(model.beads.size() * 2);
            for (int i = 0; i < (int)model.beads.size(); ++i) {
                if (!model.beads[i].active || model.transitions.count(i)) continue;
                int ix = (int)std::floor(model.beads[i].pos.x / cell);
                int iy = (int)std::floor(model.beads[i].pos.y / cell);
                grid[key(ix, iy)].push_back(i);
            }

            bool moved = false;
            for (int i = 0; i < (int)model.beads.size(); ++i) {
                Bead& bi = model.beads[i];
                if (!bi.active || model.transitions.count(i)) continue;
                int ix = (int)std::floor(bi.pos.x / cell);
                int iy = (int)std::floor(bi.pos.y / cell);
                for (int gx = ix - 1; gx <= ix + 1; ++gx)
                    for (int gy = iy - 1; gy <= iy + 1; ++gy) {
                        auto it = grid.find(key(gx, gy));
                        if (it == grid.end()) continue;
                        for (int j : it->second) {
                            if (j <= i) continue;
                            Bead& bj = model.beads[j];
                            float dx = bj.pos.x - bi.pos.x, dy = bj.pos.y - bi.pos.y;
                            float d2 = dx * dx + dy * dy;
                            if (d2 >= minSep2) continue;
                            if (Adjacent(model, i, j)) continue;
                            if (sharesCrossing(i, j)) continue;

                            float d = std::sqrt(d2);
                            if (d < GEOM_EPS) {
                                std::uniform_real_distribution<float> jit(-1.0f, 1.0f);
                                dx = jit(model.rng); dy = jit(model.rng); d = std::hypot(dx, dy);
                                if (d < GEOM_EPS) { dx = 1.0f; dy = 0.0f; d = 1.0f; }
                            }
                            float pen = minSep - d, ux = dx / d, uy = dy / d;
                            bool ai = bi.isVertex || bi.driven;
                            bool aj = bj.isVertex || bj.driven;
                            float wi = ai ? 0.0f : (aj ? 1.0f : 0.5f);
                            float wj = aj ? 0.0f : (ai ? 1.0f : 0.5f);
                            bi.pos.x -= ux * pen * wi; bi.pos.y -= uy * pen * wi;
                            bj.pos.x += ux * pen * wj; bj.pos.y += uy * pen * wj;
                            moved = true;
                        }
                    }
            }
            if (!moved) break;
        }
    }

    // Settle each pin at the open centre of its face by repulsion balance.
    // It is pushed away from every nearby strand bead AND crossing, so inside a
    // closed face the pushes cancel at the centre and any drift toward a wall is
    // resisted by that wall -- the pin stays in its region without chasing a
    // (jumpy) centroid. Pins also separate from one another. Overdamped + a
    // step clamp makes the motion smooth.
    void SettlePins() {
        const float R2s = PIN_SETTLE_RANGE * PIN_SETTLE_RANGE;
        const float PPs = PIN_PIN_RANGE * PIN_PIN_RANGE;
        const float maxStep2 = PIN_MAX_STEP * PIN_MAX_STEP;

        for (auto& pin : model.pins) {
            if (!pin.active) continue;
            Vector2 f = { 0.0f, 0.0f };

            // pushed off every nearby strand bead and crossing
            for (auto& b : model.beads) {
                if (!b.active) continue;
                float dx = pin.pos.x - b.pos.x, dy = pin.pos.y - b.pos.y;
                float d2 = dx * dx + dy * dy;
                if (d2 >= R2s) continue;
                float d = std::sqrt(d2);
                if (d < GEOM_EPS) {
                    std::uniform_real_distribution<float> jit(-0.5f, 0.5f);
                    dx = jit(model.rng); dy = jit(model.rng); d = std::hypot(dx, dy);
                }
                float mag = (PIN_SETTLE_RANGE - d) * PIN_SETTLE_STIFF;
                f.x += (dx / d) * mag; f.y += (dy / d) * mag;
            }

            // keep pins apart
            for (auto& q : model.pins) {
                if (&q == &pin || !q.active) continue;
                float dx = pin.pos.x - q.pos.x, dy = pin.pos.y - q.pos.y;
                float d2 = dx * dx + dy * dy;
                if (d2 >= PPs || d2 < GEOM_EPS) continue;
                float d = std::sqrt(d2);
                float mag = (PIN_PIN_RANGE - d) * PIN_PIN_STIFF;
                f.x += (dx / d) * mag; f.y += (dy / d) * mag;
            }

            float dx = f.x * INTEGRATION_STEP, dy = f.y * INTEGRATION_STEP;
            float m2 = dx * dx + dy * dy;
            if (m2 > maxStep2) { float s = PIN_MAX_STEP / std::sqrt(m2); dx *= s; dy *= s; }
            pin.pos.x += dx; pin.pos.y += dy;
        }
    }

    // ---- smooth removal ---------------------------------------------------

    Vector2 Centroid(const std::vector<int>& ids) const {
        Vector2 c = { 0.0f, 0.0f }; int n = 0;
        for (int id : ids) { c.x += model.beads[id].pos.x; c.y += model.beads[id].pos.y; ++n; }
        if (n) { c.x /= n; c.y /= n; }
        return c;
    }

    void StartCollapse(const std::vector<int>& doomed) {
        // Only the ARC edge beads shrink. Crossings are left untouched so their
        // outer springs keep full length and the rest of the loop is undistorted
        // while the kink / lens closes.
        for (int b : doomed) {
            if (!model.beads[b].isVertex)
                model.transitions[b] = { 1.0f, -TRANSITION_COLLAPSE_RATE };
            pending.insert(b);
        }
        jobs.push_back({ doomed, COLLAPSE_FRAMES, Centroid(doomed) });
    }

    void DeactivateNearestPin(Vector2 center) {
        int best = -1; float bestD = GRID_SIZE * 1.3f * GRID_SIZE * 1.3f;
        for (auto& pin : model.pins) {
            if (!pin.active) continue;
            float d = Dist2(pin.pos, center);
            if (d < bestD) { bestD = d; best = pin.id; }
        }
        if (best >= 0) model.pins[best].active = false;
    }

    void TickTransitions() {
        for (auto it = model.transitions.begin(); it != model.transitions.end();) {
            it->second.scale += it->second.rate;
            if (it->second.scale < it->second.minScale) it->second.scale = it->second.minScale;
            if (it->second.rate > 0.0f && it->second.scale >= 1.0f) { it = model.transitions.erase(it); continue; }
            ++it;
        }
        for (auto jit = jobs.begin(); jit != jobs.end();) {
            if (--jit->frames <= 0) {
                model.CollapseRegion(jit->doomed);
                for (int d : jit->doomed) { model.transitions.erase(d); pending.erase(d); }
                DeactivateNearestPin(jit->center);
                jit = jobs.erase(jit);
            }
            else ++jit;
        }
    }

    // ---- pin-click dispatch ----------------------------------------------

    // Click a pin: find the face it sits in, report its bordering crossings,
    // classify (1/2/3), and run the corresponding move smoothly.
    void DispatchAtPin(Vector2 p) {
        // Don't launch a new move while one is still animating (a scheduled R1
        // eat or R3 walk in flight, or a collapse job running): a second move
        // landing on the same shrinking region would corrupt the links.
        if (!delayedTasks.empty() || !jobs.empty()) return;

        FaceClass f = model.ClassifyFace(p);
        if (f.type == 0) { printf("[click] face has no R1/R2/R3 move\n"); fflush(stdout); return; }
        printf("[click] face borders %zu crossing(s) -> R%d\n", f.crossings.size(), f.type);
        fflush(stdout);
        if (f.type == 1 && !OverlapsPending(f.doomed)) {
            // Remove this face's red pin FIRST, then eat the monogon one edge
            // bead every R1_COLLAPSE_GAP frames and weld the crossing out on the
            // final tick.
            DeactivateNearestPin(Centroid(f.doomed));
            model.Monogon_Collapse_R1(f.crossings[0], f,
                [this](int d, std::function<void()> fn) { tickDelay(d, std::move(fn)); });
        }
        else if (f.type == 2 && !OverlapsPending(f.doomed)) {
            // Drop the bigon's pin, perform the strand-preserving swap, then let
            // the remaining pins re-settle into the new face layout (SettlePins
            // runs every quiescent frame; R2 is instant so it resumes at once).
            DeactivateNearestPin(Centroid(f.doomed));
            model.Bigon_Swap_R2(f.crossings[0], f.crossings[1], f,
                [this](int d, std::function<void()> fn) { tickDelay(d, std::move(fn)); });
        }
        else if (f.type == 3) {

            for (int k = 0; k < f.triArcs.size(); k++) {
                for (int i = 0; i < f.triArcs[k].size(); i++) {
                    printf(" %d |", f.triArcs[k][i]);
                }
                printf("\n");
            }
            int c0 = f.crossings[0]; int c1 = f.crossings[1]; int c2 = f.crossings[2];

            model.Perform_R3_Walk_S3(c0, c1, c2, f,
                [this](int d, std::function<void()> fn) { tickDelay(d, std::move(fn)); });


            printf("c0 -> %d | c1 -> %d | c2 -> %d\n", c0, c1, c2);

            printf("slot c0 -> c1 = %d | c1 -> c0 = %d\n", model.FindSlotTowardsCrossing(c0, c1),
                model.FindSlotTowardsCrossing(c1, c0));

            printf("slot c1 -> c2 = %d | c2 -> c1 = %d\n", model.FindSlotTowardsCrossing(c1, c2),
                model.FindSlotTowardsCrossing(c2, c1));

            printf("slot c2 -> c0 = %d | c0 -> c2 = %d\n", model.FindSlotTowardsCrossing(c2, c0),
                model.FindSlotTowardsCrossing(c0, c2));

            fflush(stdout);
        }
    }

    // ---- auto-tighten (R1/R2 only) ---------------------------------------

    bool OverlapsPending(const std::vector<int>& ids) const {
        for (int id : ids) if (pending.count(id)) return true;
        return false;
    }

    std::vector<Candidate> CollectCandidates() const {
        std::vector<Candidate> out;
        const int nc = (int)model.crossingBead.size();

        for (int vb : model.crossingBead) {
            if (!model.beads[vb].active) continue;
            std::vector<int> loop = model.FindSelfLoop(vb);
            if (loop.empty()) continue;
            std::vector<int> doomed = loop; doomed.push_back(vb);
            if (!OverlapsPending(doomed)) out.push_back({ doomed, Centroid(doomed), "R1" });
        }
        for (int a = 0; a < nc; ++a) {
            int A = model.crossingBead[a]; if (!model.beads[A].active) continue;
            std::vector<SlotArc> fromA = model.ArcsFrom(A);
            for (int b = a + 1; b < nc; ++b) {
                int B = model.crossingBead[b]; if (!model.beads[B].active) continue;
                std::vector<const SlotArc*> ab;
                for (auto& s : fromA) if (s.endVertex == B) ab.push_back(&s);
                // smallest valid lens
                int bx = -1, by = -1; size_t bestSz = (size_t)-1;
                for (size_t x = 0; x < ab.size(); ++x)
                    for (size_t y = x + 1; y < ab.size(); ++y)
                        if (model.ValidBigon(A, B, *ab[x], *ab[y])) {
                            size_t sz = ab[x]->beads.size() + ab[y]->beads.size();
                            if (sz < bestSz) { bestSz = sz; bx = (int)x; by = (int)y; }
                        }
                if (bx < 0) continue;
                std::vector<int> doomed = { A, B };
                doomed.insert(doomed.end(), ab[bx]->beads.begin(), ab[bx]->beads.end());
                doomed.insert(doomed.end(), ab[by]->beads.begin(), ab[by]->beads.end());
                if (!OverlapsPending(doomed)) out.push_back({ doomed, Centroid(doomed), "R2" });
            }
        }
        return out;
    }

    int CountCrossings() const {
        int n = 0; for (int vb : model.crossingBead) if (model.beads[vb].active) ++n; return n;
    }

    void AutoTick() {
        if (!autoRun || !jobs.empty()) return;
        std::vector<Candidate> cands = CollectCandidates();
        if (cands.empty()) { autoRun = false; return; }
        StartCollapse(cands.front().doomed);
        printf("[auto] %s removed, %d crossings left\n", cands.front().kind, CountCrossings() - 1);
        fflush(stdout);
    }

    // ---- input ------------------------------------------------------------

    void HandleInput(Vector2 mouse) {
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            bool onBtn = mouse.x > SCREEN_WIDTH / 2 - 110 && mouse.x < SCREEN_WIDTH / 2 + 110 &&
                mouse.y > 32 && mouse.y < 68;
            if (onBtn) { autoRun = !autoRun; return; }

            for (auto& pin : model.pins)
                if (pin.active && CheckCollisionPointCircle(mouse, pin.pos, RADIUS)) {
                    DispatchAtPin(pin.pos); return;
                }

            for (int i = 0; i < (int)model.beads.size(); ++i) {
                if (!model.beads[i].active) continue;
                if (CheckCollisionPointCircle(mouse, model.beads[i].pos, RADIUS)) {
                    dragging = i;
                    model.beads[i].driven = true;
                    dragOffset = { model.beads[i].pos.x - mouse.x, model.beads[i].pos.y - mouse.y };
                    break;
                }
            }
        }
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && dragging >= 0)
            model.beads[dragging].pos = { mouse.x + dragOffset.x, mouse.y + dragOffset.y };
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT) && dragging >= 0) {
            model.beads[dragging].driven = false;
            dragging = -1;
        }
    }

    void Update(Vector2 mouse) {
        HandleInput(mouse);
        for (int s = 0; s < UPDATES_PER_FRAME; ++s) Step();
        TickTransitions();
        TickDelays();
        AutoTick();
        // Re-settle pins only on a clean loop, not while a move sweeps strands
        // through their regions.
        if (jobs.empty()) {
            // Bead-adding is removed; only redundant kinked/piled beads are pruned.
            if (pending.empty() && delayedTasks.empty()) {
                model.SimplifyDenseArcs(EDGE_TOUCH_LENGTH, MIN_ARC_ANGLE_DEG);
                // Recompute regions from the live diagram and keep exactly one pin
                // per bounded region, each inside it. Runs only on a fully clean,
                // settled loop so the traced faces are well-formed.
                ReconcilePins(model, true);
            }
            SettlePins();
        }
    }

    void Render() const {
        BeginDrawing();
        ClearBackground(BLACK);
        DrawWireSpline(model);
        DrawInterface(tension, autoRun);
        DrawPins(model);
        DrawWireSpline(model);
        //Draw_Debug_Screen(model);
        EndDrawing();
    }
};

// ============================================================================
// START MENU
//
// A single text field for the vertex permutation sigma (cycle notation, e.g.
// "(16,5,-1,-6)(12,1,-13,-2)...") plus a Start button. Start parses the string,
// builds the model via the orthogonal-layout path (BuildLoopFromSigmaString),
// and on success hands the LoopModel to the game; on failure it shows the
// validator's message and stays on the menu. Enter also starts; Ctrl+V pastes.
// ============================================================================

struct StartMenu {
    // Pre-filled with the 10_1_18 example so the user can just press Start.
    std::string text =
        "(16,5,-1,-6)(12,1,-13,-2)(7,2,-8,-3)(11,6,-12,-7)"
        "(4,9,-5,-10)(15,10,-16,-11)(8,13,-9,-14)(3,14,-4,-15)";
    std::string error;
    float       backspaceTimer = 0.0f;

    // Centred widgets.
    Rectangle InputBox()  const { return { 150.0f, 470.0f, (float)SCREEN_WIDTH - 300.0f, 56.0f }; }
    Rectangle StartBtn()  const { return { SCREEN_WIDTH / 2.0f - 100.0f, 560.0f, 200.0f, 56.0f }; }

    // Returns true when Start (button or Enter) is pressed this frame.
    bool Update(Vector2 mouse) {
        // --- text entry ---
        int c;
        while ((c = GetCharPressed()) > 0)
            if (c >= 32 && c < 127) text.push_back((char)c);   // printable ASCII only

        // backspace with auto-repeat
        if (IsKeyPressed(KEY_BACKSPACE)) {
            if (!text.empty()) text.pop_back();
            backspaceTimer = 0.4f;
        }
        else if (IsKeyDown(KEY_BACKSPACE)) {
            backspaceTimer -= GetFrameTime();
            if (backspaceTimer <= 0.0f) { if (!text.empty()) text.pop_back(); backspaceTimer = 0.04f; }
        }

        // paste (Ctrl/Cmd + V)
        if ((IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_RIGHT_CONTROL)) && IsKeyPressed(KEY_V)) {
            const char* clip = GetClipboardText();
            if (clip) text += clip;
        }

        // --- start triggers ---
        bool start = IsKeyPressed(KEY_ENTER) || IsKeyPressed(KEY_KP_ENTER);
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && CheckCollisionPointRec(mouse, StartBtn()))
            start = true;
        return start;
    }

    void Draw(Vector2 mouse) const {
        BeginDrawing();
        ClearBackground(BLACK);

        const char* title = "THE UNPINNING GAME";
        DrawText(title, SCREEN_WIDTH / 2 - MeasureText(title, 48) / 2, 250, 48, WHITE);
        const char* sub = "Enter the vertex permutation sigma, then press Start.";
        DrawText(sub, SCREEN_WIDTH / 2 - MeasureText(sub, 20) / 2, 330, 20, GRAY);
        const char* ex = "e.g.  (16,5,-1,-6)(12,1,-13,-2)(7,2,-8,-3)(11,6,-12,-7) ...";
        DrawText(ex, SCREEN_WIDTH / 2 - MeasureText(ex, 16) / 2, 360, 16, DARKGRAY);

        // input box
        Rectangle box = InputBox();
        DrawRectangleRec(box, Fade(WHITE, 0.06f));
        DrawRectangleLinesEx(box, 2.0f, Fade(WHITE, 0.6f));

        const int   fs = 20, pad = 12;
        const int   innerW = (int)box.width - 2 * pad;
        const char* shown = text.empty()
            ? "(type or paste cycles here)"
            : text.c_str();
        const Color shownColor = text.empty() ? DARKGRAY : WHITE;
        const int   tw = MeasureText(text.c_str(), fs);
        const int   scroll = (tw > innerW) ? (tw - innerW) : 0;   // keep the tail visible

        BeginScissorMode((int)box.x + pad, (int)box.y, innerW, (int)box.height);
        DrawText(shown, (int)box.x + pad - scroll, (int)box.y + ((int)box.height - fs) / 2, fs, shownColor);
        EndScissorMode();

        // blinking caret at the end of the (scrolled) text
        if (((int)(GetTime() * 2.0)) % 2 == 0) {
            int caretX = (int)box.x + pad - scroll + (text.empty() ? 0 : tw);
            DrawRectangle(caretX, (int)box.y + 10, 2, (int)box.height - 20, WHITE);
        }

        // start button (hover highlight)
        Rectangle btn = StartBtn();
        bool hover = CheckCollisionPointRec(mouse, btn);
        DrawRectangleRec(btn, hover ? BLUE : Fade(BLUE, 0.7f));
        const char* bl = "START";
        DrawText(bl, (int)(btn.x + btn.width / 2) - MeasureText(bl, 24) / 2,
                 (int)(btn.y + btn.height / 2) - 12, 24, WHITE);

        // error message
        if (!error.empty()) {
            std::string msg = "Could not start: " + error;
            DrawText(msg.c_str(), SCREEN_WIDTH / 2 - MeasureText(msg.c_str(), 18) / 2, 650, 18, RED);
        }

        const char* foot = "Enter = Start    |    Ctrl+V = paste    |    Esc (in game) = back to menu";
        DrawText(foot, SCREEN_WIDTH / 2 - MeasureText(foot, 14) / 2, SCREEN_HEIGHT - 40, 14, GRAY);
        EndDrawing();
    }
};

// ============================================================================
// ENTRY POINT
// ============================================================================

#ifndef UPG_NO_MAIN
enum class AppState { Menu, Playing };

int main() {
    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "The Unpinning Game");
    SetTargetFPS(120);

    AppState  state = AppState::Menu;
    StartMenu menu;
    GameState gs;

    while (!WindowShouldClose()) {
        Vector2 mouse = GetMousePosition();

        if (state == AppState::Menu) {
            if (menu.Update(mouse)) {
                // Parse the sigma string and build the loop via the orthogonal
                // layout. On success, start; on failure, surface the message.
                std::string err;
                LoopModel m = BuildLoopFromSigmaString(menu.text, err);
                if (err.empty()) {
                    gs = GameState{};          // fresh game
                    gs.Init(std::move(m));
                    menu.error.clear();
                    state = AppState::Playing;
                }
                else {
                    menu.error = err;
                }
            }
            menu.Draw(mouse);
        }
        else { // Playing
            if (IsKeyPressed(KEY_ESCAPE)) { state = AppState::Menu; continue; }
            gs.Update(mouse);
            gs.Render();
        }
    }

    CloseWindow();
    return 0;
}
#endif