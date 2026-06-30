# The Unpinning Game - build & run

## Start menu
A text field for the vertex permutation sigma (cycle notation) plus a **Start**
button. Type/paste e.g.

    (16,5,-1,-6)(12,1,-13,-2)(7,2,-8,-3)(11,6,-12,-7)(4,9,-5,-10)(15,10,-16,-11)(8,13,-9,-14)(3,14,-4,-15)

and press **Start** (or Enter). The string is parsed, validated, laid out with
the orthogonal grid engine, and the game begins. Invalid input stays on the menu
and shows the reason. **Esc** in-game returns to the menu. **Ctrl+V** pastes.
SnapPy PD-style all-positive bracketed codes are also accepted.

Pipeline (all in sigma_import.*): `BuildLoopFromSigmaString`
-> `Permutation::FromString` -> `ValidateMultiloopSigma` -> `ComputeGridLayout`
-> `SigmaGridPolygon` (single straight-through traversal, crossings doubled)
-> `BuildLoopFromPolygon`.

## What changed in this revision
The permutation / multiloop / orthogonal code from the other project has been
**absorbed into our project and made self-contained**:
- The external `geometry.h` dependency is gone. It was only used by two unused
  helpers (`GridFacePolygon`/`GridFaceArea`); those were removed, so
  `orthogonal.cpp` now compiles with no external geometry header.
- The dead `RunPermutationTests` declaration was removed.
- The bridge + string adapter were folded into `sigma_import.*` (its natural
  home). The former `sigma_to_grid.*`, `sigma_ortho.h`, and `upg_sigma_ortho.cpp`
  glue files no longer exist.
- `permutation.*`, `multiloop.*`, `orthogonal.*` carry our banners and are kept
  as separate modules (orthogonal is one ~900-line cohesive algorithm; inlining
  it into another file would undo the earlier split).

## Files
Core game: `upg.h`, `upg_build.cpp`, `upg_surgery.cpp`, `upg_faces.cpp`,
`upg_moves.cpp`, `upg_render.cpp`, `regions.h/.cpp`.
Sigma -> loop: `sigma_import.h/.cpp` (Tutte builder + orthogonal builder + menu
entry point), `permutation.h/.cpp`, `multiloop.h/.cpp`, `orthogonal.h/.cpp`.
Driver/menu: `Sandbox_UPG.cpp`.
No `geometry.h` and no separate bridge/adapter files are needed.

## Build (with raylib installed)
```sh
g++ -std=c++20 *.cpp -o upg $(pkg-config --cflags --libs raylib)
./upg
```
C++20 required (std::set::contains, structured bindings).

## Verified headless
Real orthogonal `ComputeGridLayout` (not a stub) builds:
trefoil 3_1 (3 crossings), figure-8 4_1 (4), cinquefoil 5_1 (5), and
10_1_18 (8, single component). Malformed / 3-valent / empty input is rejected
with a clear message shown on the menu.

## Notes
- Single-component (knot) input: `BuildLoopFromPolygon` takes one closed polygon.
- A drawing larger than `GRID_N` still builds; the physics re-centres it. Raise
  `GRID_N` for very large knots if you want them to start fully on-screen.
