#pragma once

// ============================================================================
// regions.h - topological region tracking for the red pins.
//
// A region is a face of the live 4-valent planar diagram. EnumerateRegions
// traces every face directly from the current bead geometry (the angular
// "next edge around the face" rule), so it stays correct after any R-move --
// no grid cells, no stale build-time data.
//
// ReconcilePins enforces the game invariant in one pass: exactly one active pin
// per bounded region, each one inside its region. It keeps the most-central pin
// in each region, deactivates duplicates, pulls escapees back in, and (when
// asked) spawns a pin at the region's pole of inaccessibility for any bounded
// region that has none. Call it on a clean (settled) loop.
// ============================================================================

#include "upg.h"

struct Region {
    std::vector<int> boundary;   // ordered boundary bead indices (one face cycle)
    bool             bounded = true;   // false = the unbounded outer region
    float            area = 0.0f;
};

// Trace all faces of the live diagram. The largest-area face is the outer one.
std::vector<Region> EnumerateRegions(const LoopModel& m);

// Most interior point of a polygon (pole of inaccessibility): the point farthest
// from the boundary. Robust for non-convex regions where a centroid escapes.
Vector2 PoleOfInaccessibility(const std::vector<Vector2>& poly);

// Enforce one active pin per bounded region. With spawnMissing, empty regions get
// a fresh pin at their pole. Returns the active pin count.
int ReconcilePins(LoopModel& m, bool spawnMissing);
