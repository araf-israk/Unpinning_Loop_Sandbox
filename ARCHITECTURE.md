# The Unpinning Game — code architecture

An interactive C++/raylib knot-diagram simulator. A knot diagram is held as a
**bead necklace**; a spring-mass physics loop relaxes it on screen, and clicking
a region's pin runs the Reidemeister move that region's face allows
(R1 / R2 / R3), simplifying the diagram toward the unknot.

This document describes the data model, the file layout (after splitting the old
monolithic `upg.cpp`), and how a click turns into a move.

---

## 1. The data model (`upg.h`)

Everything is one `std::vector<Bead>` (`LoopModel::beads`); all links are
**indices** into that vector, never pointers, so the array can grow without
invalidating the topology.

A `Bead` is one of two things:

- **Crossing (vertex) bead** — `isVertex == true`. Has four *slots*. For each
  slot `s`:
  - `neighborOfSlot[s]` = the index of the adjacent bead along that slot's dart.
  - `throughSlot[s]` = the slot the strand entering at `s` continues out of (an
    involution: the "straight through" partner). This is the only thing that
    encodes the over/under weave combinatorially.
- **Edge bead** — `isVertex == false`. A point on an arc, with exactly two
  strand neighbours `neighbors[0..1]`.

`LoopModel::Through(vertexBead, fromBead)` reads `throughSlot` to find where a
strand continues when it passes straight through a crossing. Most traversal code
is built on this plus the edge-bead neighbour walk.

Other state: `crossingBead` (indices of the live crossings — kept in sync by
`RebuildCrossingList`), `pins` (one red pin per bounded region), the build-time
grid `walls` / union-find / `regionBeads` / `regionCrossings` maps, and
`transitions` (per-bead rest-length scale for the smooth collapse/grow
animation).

Why this representation: a removable **kink** is just a self-loop arc (one
crossing, leaves and returns to itself), and a removable **bigon** is two
crossing-free arcs between the same two crossings — both are read straight off
the links with no separate face computation needed.

---

## 2. File layout

The old `upg.cpp` did six unrelated jobs in 1,500 lines. It is now split into
five focused translation units. **No logic changed** in the split — functions
were moved verbatim; only file-level comments and per-file `#include`s were
added. Replace `upg.cpp` in your build with these five files.

| File | Responsibility |
|------|----------------|
| `upg.h` | The whole model: constants, `Bead`/`Pin`/`LoopModel`, free-function decls. |
| `upg_build.cpp` | `GridCoordToCanvasCoord`, grid topology map (`BuildTopologyMap`/`GetFaceId`), and `BuildLoopFromPolygon` (polygon → necklace). |
| `upg_surgery.cpp` | Link-rewiring primitives (`WalkFromSlot`, `ReplaceLink`, `CollapseRegion`, …) and dynamic resampling (`InsertEdgeBead`, `SplitLink`, `SimplifyDenseArcs`, `BuildCrossingStrandMap`). |
| `upg_faces.cpp` | Arc/face queries (`ArcsFrom`, `ValidBigon`) and `ClassifyFace` — the R1/R2/R3 detector. |
| `upg_moves.cpp` | The R1/R2/R3 moves and the R3 T-junction-walk machinery. |
| `upg_render.cpp` | raylib drawing (`DrawWireSpline`, `DrawPins`, `DrawInterface`, `Draw_Debug_Screen`). |
| `sigma_import.h/.cpp` | Alternative builder from a combinatorial-map permutation σ, plus `TutteLayout`, face enumeration (`BuildFacesFromPhi`), and pin seeding. |
| `regions.h/.cpp` | Live region tracking: `EnumerateRegions` (trace faces from current geometry), `PoleOfInaccessibility`, `ReconcilePins` (one pin per bounded region). |
| `Sandbox_UPG.cpp` | The driver: `GameState` physics step, collision/pin settling, frame-scheduled animation callbacks, click dispatch, auto-tighten, and `main()`. |

Dependency direction: every `.cpp` includes `upg.h`; `sigma_import` and
`regions` add their own headers; `Sandbox_UPG` includes all three headers.

---

## 3. Two ways to build a model

1. **From a grid polygon** — `BuildLoopFromPolygon` (`upg_build.cpp`). Trace a
   rectilinear curve whose doubled points mark crossings; this reproduces the
   LooPindex `annotated.svg` drawing exactly. This is the path `main()` uses
   (`LooPindex_10_1_18_GridPoints`).
2. **From a σ permutation** — `BuildLoopFromSigma` (`sigma_import.cpp`) with a
   `TutteLayout`. Topologically identical, but arcs are straight crossing-to-
   crossing rather than rectilinear, so it does not match the drawing. Kept as
   an alternative (commented out in `main()`).

Either way the result is the same `LoopModel`, and `GameState::Init` then calls
`ReconcilePins` to seed one pin per bounded region.

---

## 4. The physics loop (`Sandbox_UPG.cpp`)

`GameState::Step()` runs `UPDATES_PER_FRAME` times per frame:

1. **Springs** along every edge link (rest length scaled by any active
   transition), added once per link regardless of orientation; plus direct
   vertex↔vertex springs for T-junction bars with no bead between them.
2. **Pins push outward** on nearby beads.
3. **Centring** — a uniform translation toward the canvas centre (no shape
   distortion) so the loop neither drifts nor wanders.
4. **Integrate** with a per-substep displacement clamp for stability.
5. **`ResolveOverlaps`** — a hard position-projection pass on a uniform grid that
   separates any overlapping bead disks, *exempting* directly-linked beads and
   the two strands that genuinely cross at a crossing (identified topologically
   via the crossing-guard radius, not by distance alone).

`SettlePins` then relaxes each pin to the open centre of its face by repulsion
balance. On a fully settled, quiescent loop the driver also runs
`SimplifyDenseArcs` and re-runs `ReconcilePins`.

---

## 5. Click → move pipeline

```
click pin
  └─ GameState::DispatchAtPin(pos)
       └─ LoopModel::ClassifyFace(pos)        [upg_faces.cpp]
            • smallest face under the point, classified by #bordering crossings
            • returns FaceClass{ type 1/2/3, crossings, doomed beads, triArcs }
       ├─ type 1 → Monogon_Collapse_R1(...)   [upg_moves.cpp]  (eat the kink)
       ├─ type 2 → Bigon_Swap_R2(...)         [upg_moves.cpp]  (open the lens)
       └─ type 3 → Perform_R3_Walk_S3(...)    [upg_moves.cpp]  (triangle flip)
```

Moves are animated by handing them a `schedule(framesLater, fn)` callback
(`GameState::tickDelay`), so the surgery is spread over frames while the physics
keeps the strands taut. With no scheduler the same code runs instantly (used by
the headless tests and auto-tighten).

- **R1** removes a self-loop one edge bead at a time, then welds the crossing's
  two outer legs straight through (`CollapseRegion`-style splice).
- **R2** opens a bigon by cross-joining each lens leg to the *other* strand's
  outer leg — the only smoothing that opens the lens without pinching off a
  stray component (the historical R2 bug) or leaving the arcs crossing.
- **R3** uses the "T-junction walk": split each corner crossing into two
  T-junctions, walk them to their arc's midpoint
  (`CreateTConnection` → `TraverseTConnection`), then `MergeTConnections` fuses
  each kissing pair into the move's new crossing.

`auto-tighten` (the button) repeatedly applies the smallest available R1/R2 with
no scheduler until none remain.

---

## 6. Notes & known rough edges

- `PlanR3` / `ApplyR3` are a pure leg-swap form of R3, and `PlanFakeR3` is the
  symmetric 3-arc swap that morphs the diagram into the trefoil shadow — that
  one is **not** a Reidemeister move (it changes the knot type) and is kept only
  for the animation. The live click path uses `Perform_R3_Walk_S3`, not these.
- `RoleSwapCrossingAndEdge`, `SlideOuterBeadIntoArc`, `SlideCrossingOverBead` and
  `StartR3MicroSurgery` are **experimental R3 micro-surgery scaffolding that is
  not wired into the game.** In particular `RoleSwapCrossingAndEdge` still
  contains hard-coded bead indices (`3, 29, 30, 31`) from a debugging session
  and is not correct in general. They are preserved verbatim in `upg_moves.cpp`
  (with a header note) so behaviour is unchanged; rewrite before use.

---

## 7. Building

Compile the eight `.cpp` files against raylib with a C++20 compiler (the code
uses `std::set::contains` and structured bindings):

```sh
g++ -std=c++20 *.cpp -o upg $(pkg-config --cflags --libs raylib)
```

(Replace `upg.cpp` with the five `upg_*.cpp` files in any existing project/IDE
configuration.) The repository was verified to compile and link cleanly after
the split, with an identical set of defined symbols to the original single-file
build.
