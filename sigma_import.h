#pragma once

// ============================================================================
// sigma_import.h - build a LoopModel from a combinatorial-map permutation
// representation (the form LooPindex / Plantri emit), instead of from a grid
// polygon.
//
// Darts are signed ints +-1..+-E. Edge k has ends +k and -k; the edge
// involution epsilon = (k,-k) is implicit and assumed standard. The vertex
// permutation sigma is given as one 4-cycle per crossing, listing the four
// incident darts in rotational order. Sign convention: +k is the dart oriented
// INTO its vertex (edge k's head), -k is oriented OUT (a tail), so a coherently
// oriented loop reads (+,+,-,-) cyclically at every crossing.
//
// Slot mapping into the Bead struct: slot s of a crossing = sigma cycle entry s,
// so neighborOfSlot[s] is the nearest edge bead along that dart and
// throughSlot[s] = (s+2)%4 -- the rotational-opposite entry, i.e. sigma^2, the
// strand that passes straight through. (BuildLoopFromPolygon instead pairs
// through-slots as 0<->1 / 2<->3; both are valid because the topology code reads
// throughSlot, never slot-index geometry.)
//
// Positions: crossing grid coords are supplied by the caller; edge beads are
// laid on the straight segment between the two crossings a dart joins, spaced at
// the physics rest length. regionBeads / pins are NOT populated (those need the
// grid walls); the faces are recoverable as the cycles of phi=(sigma*eps)^-1.
// ============================================================================

#include "upg.h"
#include "multiloop.h"     // Multiloop, ValidateMultiloopSigma
#include "orthogonal.h"    // GridLayout, GridPoint, ComputeGridLayout
#include <array>
#include <string>

// Build a LoopModel from sigma. `sigma[v]` is the rotational 4-cycle of vertex
// v; `crossingGrid[v]` is that vertex's position in grid-cell coordinates.
LoopModel BuildLoopFromSigma(const std::vector<std::array<int, 4>>& sigma,
    const std::vector<Vector2>& crossingGrid);

// Quick ring layout for testing when no embedding is on hand. NOTE: a ring seed
// is not a planar embedding; the physics relaxer may settle it into a different
// shadow. Supply a real layout (e.g. TutteLayout) for fidelity.
std::vector<Vector2> RingLayout(int numCrossings,
    float radiusGrid = 2.0f, Vector2 centerGrid = { 3.0f, 3.0f });

// Planar straight-line layout from sigma alone (Tutte barycentric embedding):
// the largest face is pinned to a convex polygon and every other crossing is
// relaxed to its neighbours' centroid. Produces a crossing-free embedding whose
// edges meet only at crossings -- exactly what the geometric ClassifyFace wants.
// Self-contained: no PD code or external layout tool needed.
std::vector<Vector2> TutteLayout(const std::vector<std::array<int, 4>>& sigma,
    Vector2 centerGrid = { GRID_N / 2.0f, GRID_N / 2.0f },
    float radiusGrid = (GRID_N - 2) / 2.0f);

// Enumerate faces by walking phi = (sigma*eps)^-1 and fill m.regionBeads and
// m.regionCrossings, keyed by face index 0..F-1. Each face is one phi orbit;
// regionBeads gets that face's bordering edge beads plus its corner crossing
// beads, regionCrossings just the corner crossings. Returns the face count.
// NOTE: these maps are presently write-only in UPG (the live code finds faces
// geometrically); this is the combinatorial parity of the polygon builder's
// region maps. Run after BuildLoopFromSigma.
int BuildFacesFromPhi(LoopModel& m, const std::vector<std::array<int, 4>>& sigma);

// One active pin per BOUNDED face, placed at the face's corner centroid (which
// is interior because Tutte inner faces are convex). The outer/unbounded face
// (largest corner-polygon area) gets no pin. Self-contained from sigma; clears
// existing pins. SettlePins then relaxes each pin to its region's open centre.
int SeedPinsPerFace(LoopModel& m, const std::vector<std::array<int, 4>>& sigma);

// The 10^1_18 vertex permutation, ready to hand to BuildLoopFromSigma.
inline std::vector<std::array<int, 4>> Sigma_10_1_18() {
    return {
        {16, 5, -1, -6}, {12, 1, -13, -2}, {7, 2, -8, -3}, {11, 6, -12, -7},
        {4, 9, -5, -10}, {15, 10, -16, -11}, {8, 13, -9, -14}, {3, 14, -4, -15}
    };
}

// The 10^1_18 loop as drawn on LooPindex: a rectilinear grid polygon on the 6x6
// grid (the canonical annotated.svg layout). Crossings are inserted as doubled
// vertices, so this feeds straight into BuildLoopFromPolygon and reproduces the
// LooPindex drawing exactly (rectilinear arcs, same crossing positions).
inline std::vector<Vector2> LooPindex_10_1_18_GridPoints() {
    // The drawing is 5 grid units wide (coords 0..5); centre it on the GRID_N grid.
    const float off = (GRID_N - 5) / 2;
    std::vector<Vector2> p = {
        {1,2},{1,3},{1,5},{3,5},{3,4},{3,3},{3,2},{3,1},{3,0},{0,0},{0,3},{1,3},
        {2,3},{3,3},{4,3},{5,3},{5,1},{3,1},{2,1},{2,2},{2,3},{2,4},{3,4},{4,4},
        {4,3},{4,2},{3,2},{2,2}
    };
    for (auto& q : p) { q.x += off; q.y += off; }
    return p;
}

// One active pin per BOUNDED region of a polygon-built model, placed at the
// centre of a grid cell that genuinely belongs to that region (so it is always
// interior, even for non-convex regions). The unbounded outer region gets none.
// Clears existing pins. SettlePins then relaxes each to its region's open centre.
int SeedPinsFromCells(LoopModel& m);
// ============================================================================
// Orthogonal-layout builder (folded in from the former sigma_to_grid / sigma_ortho
// glue). This is the rectilinear alternative to BuildLoopFromSigma + TutteLayout:
// it lays sigma out with the bend-minimised orthogonal engine (orthogonal.h) and
// feeds the resulting integer rectilinear polygon through BuildLoopFromPolygon,
// so the rendered game matches a real orthogonal drawing.
// ============================================================================

// Walk the multiloop's single strand (rho order) and stitch the orthogonal edge
// polylines into ONE closed rectilinear traversal in which every crossing
// appears exactly twice (the "doubled vertex" form BuildLoopFromPolygon wants).
// Edge orientation is decided geometrically, so it is independent of the
// dart-sign convention. Throws if the diagram is not a single component.
std::vector<GridPoint> SigmaGridPolygon(const Multiloop& loop,
                                        const GridLayout& layout,
                                        int startEdge = 1);

// Validate sigma, lay it out orthogonally (exteriorFace = -1 picks the largest
// face as the infinite region), and build the necklace. Throws on failure.
LoopModel BuildLoopFromSigmaOrthogonal(const Permutation& sigma, int exteriorFace = -1);

// Menu convenience: parse cycle notation like "(16,5,-1,-6)(12,1,-13,-2)..." and
// build the model. On failure fills `errorOut` (human-readable) instead of
// throwing and returns an empty model, so the UI can show the message.
LoopModel BuildLoopFromSigmaString(const std::string& cycleNotation, std::string& errorOut);
