#include "raylib.h"
#include <vector>    // Required for std::vector
#include <string>    // Required for std::string and std::to_string
#include <cmath>     // Required for std::floor, std::hypot, std::abs, std::fmod, and std::pow
#include <memory>    // Required for std::shared_ptr and std::make_shared
#include <algorithm> // Required for std::min, std::max, std::find_if, and std::sort
#include <set>       // Required for std::set
#include <map>       // Required for std::map
#include <random>    // Required for std::random_device and std::mt19937

// ============================================================================
// 1. GLOBAL ENGINE CONFIGURATIONS & GAME PARAMETERS
// ============================================================================

// Window resolution constraints
const int SCREEN_WIDTH = 1200;
const int SCREEN_HEIGHT = 1200;

// Grid sizing: Map uses a 6x6 tile system centered nicely within 9 virtual grid blocks
const float GRID_SIZE = SCREEN_WIDTH / 9.0f;
const float KNOT_PIXEL_WIDTH = 6.0f * GRID_SIZE;

// Margins needed to perfectly position the grid game board in the screen center
const float PADDING_X = (SCREEN_WIDTH - KNOT_PIXEL_WIDTH) / 2.0f;
const float PADDING_Y = (SCREEN_HEIGHT - KNOT_PIXEL_WIDTH) / 2.0f;

// Physics parameters
const float STEP_SIZE = 0.5f;              // Inter-bead segment length spacing along the wire
const float RADIUS = 0.25f * GRID_SIZE;    // Visual and physical collision boundary size of nodes
const float REST_LENGTH = 2.0f * RADIUS;   // Ideal spring distance constraint between connected beads
const float SPRING_CONSTANT = 0.5f;        // Elasticity coefficient (Hooke's Law pull strength)
const float COLLISION_STIFFNESS = 1.0f;    // Force penalty scalar applied during grid region overlaps
const float DAMPING = 0.3f;                // Kinetic energy loss per frame to prevent infinite oscillations

// Automation control rules
const int UPDATES_PER_FRAME = 10;          // Number of physics substeps per frame for structural stability
const int MAX_FLOW_FRAMES = 150;           // Time allowance window for wire to settle before timeout
const float TENSION_THRESHOLD = 1.0f;      // Average loop tension requirement below which culling triggers

// Track automation states
enum class TightenState { IDLE, REMOVING, FLOWING };

// ============================================================================
// 2. DATA STRUCTURES & DATA ARCHITECTURE
// ============================================================================

// Represents an individual point mass in the structural simulation loop
struct Bead {
    int id;               // Unique object identifier
    Vector2 pos;          // Absolute rendering coordinates on the screen canvas
    std::string type;     // Classification flag: "face" (pins), "edge" (wire), "vertex" (intersections)
    Vector2 force = { 0.0f, 0.0f }; // Accumulated direction force vector applied in physics pass
    bool active = true;   // State tracker; false signals the bead has been culled from calculations
};

// Snapshot cache used to store state parameters for deterministic rollback/backtracking
struct BackupState {
    std::vector<int> sequenceIDs;      // Ordered list of node IDs forming the continuous loop
    std::vector<Vector2> nodePositions;// Coordinates of every structural simulation node
    std::vector<bool> nodeActive;      // Status states of every node in the master registry
};

// Disjoint-Set Union (DSU) engine for structural map partitioning
class UnionFind {
public:
    std::vector<int> p;
    UnionFind(int size) {
        p.resize(size);
        for (int i = 0; i < size; i++) p[i] = i; // Every sub-cell starts as its own unique group root
    }
    // Path compression lookups to find the primary parent identifier
    int find(int i) {
        return p[i] == i ? i : (p[i] = find(p[i]));
    }
    // Unifies two separated floor tile index sets into a single shared region identity
    void unite(int i, int j) {
        p[find(i)] = p[find(j)];
    }
};

// ============================================================================
// 3. UTILITY & GEOMETRIC CALCULATORS
// ============================================================================

// Translates a virtual 2D board index layout cleanly into absolute window display pixels
Vector2 GridCoordToCanvasCoord(float x, float y) {
    return { PADDING_X + x * GRID_SIZE, PADDING_Y + y * GRID_SIZE };
}

// Maps coordinate coordinates directly to their respective Disjoint-Set structural region roots
int GetFaceId(float cx, float cy, UnionFind& uf) {
    int i = std::floor(cx);
    int j = std::floor(cy);
    // Boundary check: Out of bounds requests resolve automatically to the outer perimeter (Zone 36)
    if (i < 0 || i >= 6 || j < 0 || j >= 6) return uf.find(36);
    return uf.find(j * 6 + i);
}

// Computes spatial cosine orientation angles using the dot product formula
float GetDotProduct(int seqIdx, const std::vector<std::shared_ptr<Bead>>& sequence) {
    int nSize = sequence.size();
    auto curr = sequence[seqIdx];
    auto prev = sequence[(seqIdx - 1 + nSize) % nSize];
    auto next = sequence[(seqIdx + 1) % nSize];

    // Build incoming and outgoing displacement direction vectors
    Vector2 d1 = { prev->pos.x - curr->pos.x, prev->pos.y - curr->pos.y };
    Vector2 d2 = { next->pos.x - curr->pos.x, next->pos.y - curr->pos.y };

    float len1 = std::hypot(d1.x, d1.y);
    float len2 = std::hypot(d2.x, d2.y);
    if (len1 == 0.0f || len2 == 0.0f) return 0.0f;

    // Returns cos(theta); values close to -1.0 identify flat, highly straight paths
    return ((d1.x * d2.x) + (d1.y * d2.y)) / (len1 * len2);
}

// Validation rules checking if a wire point can be safely removed without breaking topology
bool CanRemove(int seqIdx, const std::vector<std::shared_ptr<Bead>>& sequence, const std::set<int>& bannedBeads) {
    int nSize = sequence.size();
    if (sequence[seqIdx]->type == "vertex") return false;        // Never remove structural corners
    if (bannedBeads.contains(sequence[seqIdx]->id)) return false; // Skip nodes flagged by a rollback

    // Count line segments remaining between the neighboring vertex points
    int count = 1;
    int i = (seqIdx - 1 + nSize) % nSize;
    while (sequence[i]->type != "vertex") { count++; i = (i - 1 + nSize) % nSize; }

    i = (seqIdx + 1) % nSize;
    while (sequence[i]->type != "vertex") { count++; i = (i + 1) % nSize; }

    // Enforces a minimum of 2 checking beads per edge line segment to maintain reliable collision parsing
    return count > 1;
}

// ============================================================================
// 4. MAP PARSING & INITIALIZATION INITIALIZERS
// ============================================================================

// Identifies and catalogs structural wall coordinates across the map layout grid
void BuildTopologyMap(const std::vector<Vector2>& gridPoints, std::set<std::string>& walls, UnionFind& uf) {
    for (size_t k = 0; k < gridPoints.size(); k++) {
        Vector2 start = gridPoints[k];
        Vector2 end = gridPoints[(k + 1) % gridPoints.size()];

        // Scan line paths to identify horizontal vs vertical boundaries
        if (start.x == end.x) {
            for (int y = std::min(start.y, end.y); y < std::max(start.y, end.y); y++)
                walls.insert("V," + std::to_string((int)start.x) + "," + std::to_string(y));
        }
        else {
            for (int x = std::min(start.x, end.x); x < std::max(start.x, end.x); x++)
                walls.insert("H," + std::to_string(x) + "," + std::to_string((int)start.y));
        }
    }

    // Connect adjacent tiles where no wall is present to form open spatial loops
    for (int j = 0; j < 6; j++) {
        for (int i = 0; i < 6; i++) {
            int cellId = j * 6 + i;
            if (i < 5 && !walls.contains("V," + std::to_string(i + 1) + "," + std::to_string(j))) uf.unite(cellId, cellId + 1);
            else if (i == 5 && !walls.contains("V,6," + std::to_string(j))) uf.unite(cellId, 36);
            if (i == 0 && !walls.contains("V,0," + std::to_string(j))) uf.unite(cellId, 36);
            if (j < 5 && !walls.contains("H," + std::to_string(i) + "," + std::to_string(j + 1))) uf.unite(cellId, cellId + 6);
            else if (j == 5 && !walls.contains("H," + std::to_string(i) + ",6")) uf.unite(cellId, 36);
            if (j == 0 && !walls.contains("H," + std::to_string(i) + ",0")) uf.unite(cellId, 36);
        }
    }
}

// Populates the master game collections with tracking nodes along the puzzle path
void PopulateSimulationLoop(const std::vector<Vector2>& gridPoints, std::vector<std::shared_ptr<Bead>>& sequence,
    std::vector<std::shared_ptr<Bead>>& nodes, std::map<int, std::set<std::shared_ptr<Bead>>>& faceMap,
    UnionFind& uf, int& idCounter) {
    std::map<std::string, std::shared_ptr<Bead>> nodeMap;

    for (size_t i = 0; i < gridPoints.size(); i++) {
        Vector2 start = gridPoints[i];
        Vector2 end = gridPoints[(i + 1) % gridPoints.size()];
        float dist = std::max(std::abs(end.x - start.x), std::abs(end.y - start.y));
        float xDir = (end.x == start.x) ? 0.0f : (end.x - start.x) / std::abs(end.x - start.x);
        float yDir = (end.y == start.y) ? 0.0f : (end.y - start.y) / std::abs(end.y - start.y);

        // Subdivide target segments into individual beads based on STEP_SIZE
        for (int j = 0; j < dist / STEP_SIZE; j++) {
            float gridX = start.x + j * STEP_SIZE * xDir;
            float gridY = start.y + j * STEP_SIZE * yDir;
            std::string key = std::to_string(gridX) + "," + std::to_string(gridY);

            std::shared_ptr<Bead> bead;
            if (!nodeMap.contains(key)) {
                Vector2 coords = GridCoordToCanvasCoord(gridX, gridY);
                bead = std::make_shared<Bead>(idCounter++, coords, "edge");
                nodeMap[key] = bead;
                nodes.push_back(bead);
            }
            else {
                bead = nodeMap[key];
                bead->type = "vertex"; // Upgraded to vertex if path crosses itself
            }

            // Identify surrounding tiles to register spatial proximity links
            std::vector<Vector2> adjacentCells;
            if (std::fmod(gridX, 1.0f) != 0.0f) adjacentCells = { {gridX, gridY - 0.5f}, {gridX, gridY + 0.5f} };
            else if (std::fmod(gridY, 1.0f) != 0.0f) adjacentCells = { {gridX - 0.5f, gridY}, {gridX + 0.5f, gridY} };
            else adjacentCells = { {gridX - 0.5f, gridY - 0.5f}, {gridX + 0.5f, gridY - 0.5f}, {gridX - 0.5f, gridY + 0.5f}, {gridX + 0.5f, gridY + 0.5f} };

            for (auto& cell : adjacentCells) {
                int faceId = GetFaceId(cell.x, cell.y, uf);
                faceMap[faceId].insert(bead);
            }
            sequence.push_back(bead);
        }
    }
}

// ============================================================================
// 5. RENDERING ENGINE FUNCTIONS
// ============================================================================

// Iterates through active loop nodes and generates a seamless Catmull-Rom spline curve
void DrawWireSpline(const std::vector<std::shared_ptr<Bead>>& sequence) {
    std::vector<Vector2> activePoints;
    for (const auto& bead : sequence) {
        if (bead->active) activePoints.push_back(bead->pos);
    }

    int n = activePoints.size();
    // Safety check: Fall back to a standard line if the loop becomes too simple for curves
    if (n < 3) {
        if (n == 2) DrawLineEx(activePoints[0], activePoints[1], 3.0f, BLACK);
        return;
    }

    // Catmull-Rom boundary trick: duplicate edge points to smoothly close the spline loop
    std::vector<Vector2> closedPoints;
    closedPoints.reserve(n + 3);

    closedPoints.push_back(activePoints.back());                                       // Last element goes to front
    closedPoints.insert(closedPoints.end(), activePoints.begin(), activePoints.end()); // Core path elements
    closedPoints.push_back(activePoints[0]);                                           // First element goes to end
    closedPoints.push_back(activePoints[1]);                                           // Second element goes to end
    for (const auto& bead : sequence) DrawRing(bead->pos, RADIUS - 3.0f, RADIUS + 3.0f, 0, 360, 1, DARKGRAY);
    DrawSplineCatmullRom(closedPoints.data(), closedPoints.size(), 10.0f, WHITE);
}

// Draws the primary status indicators, real-time tension readings, and automation action buttons
void DrawInterface(TightenState state, float currentTension) {
    // 1. Text Metrics Displays
    std::string tensionText = "Total Tension: " + std::to_string(currentTension).substr(0, 5);
    DrawText(tensionText.c_str(), SCREEN_WIDTH / 2 - MeasureText(tensionText.c_str(), 20) / 2, 5, 20, RED);

    std::string statusText = "Status: IDLE";
    if (state == TightenState::REMOVING) statusText = "Status: Analyzing Vectors...";
    if (state == TightenState::FLOWING) statusText = "Status: Flowing...";
    if (state == TightenState::IDLE && currentTension > 50.0f) statusText = "Status: Finished (Max tension reached)";
    DrawText(statusText.c_str(), SCREEN_WIDTH / 2 - MeasureText(statusText.c_str(), 14) / 2, 28, 14, DARKGRAY);

    // 2. Control Toggle Box 
    Color btnColor = (state == TightenState::IDLE) ? BLUE : RED;
    std::string btnText = (state == TightenState::IDLE) ? "Start Auto-Tightening" : "Stop Tightening";
    DrawRectangle(SCREEN_WIDTH / 2 - 100, 45, 200, 40, btnColor);
    DrawText(btnText.c_str(), SCREEN_WIDTH / 2 - MeasureText(btnText.c_str(), 14) / 2, 57, 14, WHITE);
}

// Renders the solid green anchor pins with centering offset labels
void DrawPins(const std::vector<std::shared_ptr<Bead>>& faceBeads) {
    for (const auto& fBead : faceBeads) {
        DrawCircleV(fBead->pos, RADIUS, RED);
        DrawCircleLines(fBead->pos.x, fBead->pos.y, RADIUS, RED);
        std::string label = std::to_string(fBead->id);
        DrawText(label.c_str(), fBead->pos.x - MeasureText(label.c_str(), 10) / 2.0f, fBead->pos.y - 5, 20, WHITE);
    }
}

// ============================================================================
// 6. MAIN PROGRAM ENTRY & COORDINATOR LOOP
// ============================================================================
int main() {
    // Instantiate game window framework
    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "The Unpinning Game (C++ Refactored)");
    SetTargetFPS(60);

    // Initialize random number generation tools
    std::random_device rd;
    std::mt19937 gen(rd());

    // Layout configuration coordinates
    std::vector<Vector2> gridPoints = {
        {2, 0}, {4, 0}, {4, 3}, {1, 3}, {1, 5},
        {6, 5}, {6, 2}, {0, 2}, {0, 6}, {3, 6},
        {3, 1}, {5, 1}, {5, 4}, {2, 4}
    };
    std::vector<Vector2> faceCoords = {
        {2.5f, 0.5f}, {3.5f, 1.5f}, {4.5f, 1.5f}, {0.5f, 2.5f}, {2.5f, 2.5f},
        {3.5f, 2.5f}, {1.5f, 3.5f}, {2.5f, 3.5f}, {3.5f, 3.5f}, {3.5f, 4.5f}
    };

    // Tracking containers
    std::vector<std::shared_ptr<Bead>> faceBeads;
    std::vector<std::shared_ptr<Bead>> nodes;
    int idCounter = 1;

    // Build functional green anchor pins
    for (size_t k = 0; k < faceCoords.size(); k++) {
        Vector2 coords = GridCoordToCanvasCoord(faceCoords[k].x, faceCoords[k].y);
        auto bead = std::make_shared<Bead>(idCounter++, coords, "face");
        faceBeads.push_back(bead);
        nodes.push_back(bead);
    }

    // Run structural topology analysis setups
    std::set<std::string> walls;
    UnionFind uf(37);
    BuildTopologyMap(gridPoints, walls, uf);

    // Build the simulation wire loop pathing nodes
    std::vector<std::shared_ptr<Bead>> sequence;
    std::map<int, std::set<std::shared_ptr<Bead>>> faceMap;
    PopulateSimulationLoop(gridPoints, sequence, nodes, faceMap, uf, idCounter);

    // Map specific pins to their corresponding sub-grid regions
    faceMap[8].insert(faceBeads[0]);   faceMap[9].insert(faceBeads[1]);   faceMap[10].insert(faceBeads[2]);
    faceMap[15].insert(faceBeads[5]);  faceMap[22].insert(faceBeads[8]);  faceMap[14].insert(faceBeads[4]);
    faceMap[20].insert(faceBeads[7]);  faceMap[32].insert(faceBeads[3]);  faceMap[26].insert(faceBeads[6]);
    faceMap[29].insert(faceBeads[9]);

    // Track state variables
    TightenState tightenState = TightenState::IDLE;
    int flowFrames = 0;
    int consecutiveFailures = 0;
    float currentTension = 0.0f;
    std::set<int> bannedBeads;
    std::shared_ptr<Bead> draggingBead = nullptr;
    Vector2 dragOffset = { 0.0f, 0.0f };

    BackupState backup;
    int lastRemovedBeadId = -1;

    // Backup state engine capture macro logic
    auto CreateBackup = [&]() {
        backup.sequenceIDs.clear();
        for (auto& b : sequence) backup.sequenceIDs.push_back(b->id);
        backup.nodePositions.clear();
        backup.nodeActive.clear();
        for (auto& n : nodes) {
            backup.nodePositions.push_back(n->pos);
            backup.nodeActive.push_back(n->active);
        }
        };

    // Rollback recovery macro engine logic
    auto RestoreBackup = [&]() {
        sequence.clear();
        for (int id : backup.sequenceIDs) {
            auto it = std::find_if(nodes.begin(), nodes.end(), [id](const auto& n) { return n->id == id; });
            if (it != nodes.end()) sequence.push_back(*it);
        }
        for (size_t i = 0; i < nodes.size(); i++) {
            nodes[i]->pos = backup.nodePositions[i];
            nodes[i]->active = backup.nodeActive[i];
            nodes[i]->force = { 0.0f, 0.0f };
        }
        };

    // Selects and tries to remove flat nodes from the wire loop
    auto TryRemoveBead = [&]() -> bool {
        std::vector<std::vector<int>> edges;
        std::vector<int> currentEdge;

        int startIdx = 0;
        for (size_t i = 0; i < sequence.size(); i++) {
            if (sequence[i]->type == "vertex") { startIdx = i; break; }
        }

        for (size_t i = 0; i <= sequence.size(); i++) {
            int idx = (startIdx + i) % sequence.size();
            auto bead = sequence[idx];

            if (i > 0 && bead->type == "vertex") {
                if (!currentEdge.empty()) edges.push_back(currentEdge);
                currentEdge.clear();
            }
            else if (bead->type != "vertex" && bead->active) {
                if (CanRemove(idx, sequence, bannedBeads)) currentEdge.push_back(idx);
            }
        }

        if (edges.empty()) return false;

        // Prioritize longer edge paths using a random geometric scale
        std::sort(edges.begin(), edges.end(), [](const auto& a, const auto& b) { return a.size() > b.size(); });
        std::vector<float> probs;
        float totalProb = 0.0f;
        for (size_t i = 0; i < edges.size(); i++) {
            float p = std::pow(0.5f, i + 1);
            probs.push_back(p);
            totalProb += p;
        }

        std::uniform_real_distribution<float> dis(0.0f, totalProb);
        float randVal = dis(gen);
        float cumulative = 0.0f;
        std::vector<int> selectedEdge = edges.back();

        for (size_t i = 0; i < edges.size(); i++) {
            cumulative += probs[i];
            if (randVal <= cumulative) {
                selectedEdge = edges[i];
                break;
            }
        }

        // Isolate and delete the straightest node found along the edge path
        float maxDot = -INFINITY;
        std::vector<int> candidates;
        for (int idx : selectedEdge) {
            float dot = GetDotProduct(idx, sequence);
            if (dot > maxDot) {
                maxDot = dot;
                candidates = { idx };
            }
        }

        std::uniform_int_distribution<int> indexDis(0, candidates.size() - 1);
        int pickIdx = candidates[indexDis(gen)];
        auto beadToRemove = sequence[pickIdx];

        CreateBackup();
        lastRemovedBeadId = beadToRemove->id;
        sequence.erase(sequence.begin() + pickIdx);
        beadToRemove->active = false;

        return true;
        };

    // Primary Game Run Cycle Loop
    while (!WindowShouldClose()) {
        Vector2 mousePos = GetMousePosition();

        // --------------------------------------------------------------------
        // A. USER POLL CONTROL INPUT PIPELINE
        // --------------------------------------------------------------------
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            Rectangle btnBounds = { SCREEN_WIDTH / 2.0f - 100, 45, 200, 40 };
            if (CheckCollisionPointRec(mousePos, btnBounds)) {
                // Clicked Interface Button Box
                if (tightenState == TightenState::IDLE) {
                    tightenState = TightenState::REMOVING;
                    consecutiveFailures = 0;
                    bannedBeads.clear();
                }
                else {
                    tightenState = TightenState::IDLE;
                }
            }
            else {
                // Search for valid node element intersections to begin tracking drag operations
                for (auto& n : nodes) {
                    if (!n->active) continue;
                    if (CheckCollisionPointCircle(mousePos, n->pos, RADIUS)) {
                        draggingBead = n;
                        dragOffset = { n->pos.x - mousePos.x, n->pos.y - mousePos.y };
                        break;
                    }
                }
            }
        }

        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && draggingBead != nullptr) {
            draggingBead->pos = { mousePos.x + dragOffset.x, mousePos.y + dragOffset.y };
        }
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) {
            draggingBead = nullptr;
        }

        // --------------------------------------------------------------------
        // B. CORE PHYSICS CONSTRAINT ENGINE SUBSTEPS
        // --------------------------------------------------------------------
        for (int step = 0; step < UPDATES_PER_FRAME; step++) {
            currentTension = 0.0f;
            for (auto& n : nodes) n->force = { 0.0f, 0.0f };

            // Hooke's Elastic Tension Pass
            for (size_t i = 0; i < sequence.size(); i++) {
                auto n1 = sequence[i];
                auto n2 = sequence[(i + 1) % sequence.size()];
                float dx = n2->pos.x - n1->pos.x;
                float dy = n2->pos.y - n1->pos.y;
                float dist = std::hypot(dx, dy);

                if (dist == 0.0f) continue;
                float displacement = dist - REST_LENGTH;
                currentTension += std::abs(displacement);

                float forceMag = displacement * SPRING_CONSTANT;
                Vector2 f = { (dx / dist) * forceMag, (dy / dist) * forceMag };
                n1->force.x += f.x; n1->force.y += f.y;
                n2->force.x -= f.x; n2->force.y -= f.y;
            }

            // Spatial Face-Proximity Boundary Collisions Pass
            for (auto& [faceId, itemSet] : faceMap) {
                std::vector<std::shared_ptr<Bead>> arr;
                for (auto& b : itemSet) if (b->active) arr.push_back(b);

                for (size_t i = 0; i < arr.size(); i++) {
                    for (size_t j = i + 1; j < arr.size(); j++) {
                        auto b1 = arr[i], b2 = arr[j];
                        float dx = b2->pos.x - b1->pos.x;
                        float dy = b2->pos.y - b1->pos.y;
                        float dist = std::hypot(dx, dy);

                        if (dist == 0.0f) { // Random jitter injection to resolve complete overlaps
                            std::uniform_real_distribution<float> d(-0.5f, 0.5f);
                            dx = d(gen); dy = d(gen); dist = std::hypot(dx, dy);
                        }
                        if (dist < REST_LENGTH) {
                            float overlap = REST_LENGTH - dist;
                            float forceMag = overlap * COLLISION_STIFFNESS;
                            Vector2 f = { (dx / dist) * forceMag, (dy / dist) * forceMag };
                            b1->force.x -= f.x; b1->force.y -= f.y;
                            b2->force.x += f.x; b2->force.y += f.y;
                        }
                    }
                }
            }

            // Position Integration Pass
            for (auto& n : nodes) {
                if (n->active && n != draggingBead && n->type != "face") {
                    n->pos.x += n->force.x * DAMPING;
                    n->pos.y += n->force.y * DAMPING;
                }
            }

            // ----------------------------------------------------------------
            // C. AUTOMATION SYSTEM STATE CONTROLLER
            // ----------------------------------------------------------------
            if (tightenState == TightenState::REMOVING) {
                if (TryRemoveBead()) {
                    tightenState = TightenState::FLOWING;
                    flowFrames = 0;
                }
                else {
                    tightenState = TightenState::IDLE;
                }
            }
            else if (tightenState == TightenState::FLOWING) {
                flowFrames++;
                // Check if loop has settled below our tension threshold
                if ((currentTension / sequence.size()) < TENSION_THRESHOLD) {
                    consecutiveFailures = 0;
                    bannedBeads.clear();
                    tightenState = TightenState::REMOVING;
                }
                else if (flowFrames >= MAX_FLOW_FRAMES) { // Settle timeout reached; trigger rollback
                    RestoreBackup();
                    bannedBeads.insert(lastRemovedBeadId);
                    consecutiveFailures++;
                    tightenState = (consecutiveFailures >= 1) ? TightenState::IDLE : TightenState::REMOVING;
                }
            }
        }

        // --------------------------------------------------------------------
        // D. HARDWARE RENDER SYSTEM DRAW PASS
        // --------------------------------------------------------------------
        BeginDrawing();
        ClearBackground(BLACK);

        DrawWireSpline(sequence);      // Draw smooth line curves
        DrawInterface(tightenState, currentTension); // Render information layers
        DrawPins(faceBeads);           // Render the interactive target pins

        EndDrawing();
    }

    // Unload engine assets and close application context links
    CloseWindow();
    return 0;
}