#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>
#include <math.h>
#include <stdio.h>

#define G 200.0f // Gravitational constant
#define C_SPEED 1000.0f // Speed of light approximation
#define MAX_BODIES 500
#define TRAIL_LENGTH 200
#define DISK_INNER_RADIUS 80.0f
#define DISK_OUTER_RADIUS 600.0f

typedef struct {
    Vector3 position;
    Vector3 velocity;
    float mass;
    float radius;
    Color color;
    Vector3 trail[TRAIL_LENGTH];
    int trailIndex;
    bool active;
} Body;

typedef struct {
    Vector3 position;
    Vector3 velocity;
} State;

typedef struct {
    Vector3 dPosition; // velocity
    Vector3 dVelocity; // acceleration
} Derivative;

// Extract state from bodies for RK4
void getStates(Body bodies[], State states[]) {
    for(int i=0; i<MAX_BODIES; i++) {
        states[i].position = bodies[i].position;
        states[i].velocity = bodies[i].velocity;
    }
}

// Calculate derivatives (accelerations) for a given state
void calculateDerivatives(State states[], Derivative derivs[], Body bodies[], bool enableDrag) {
    for (int i = 0; i < MAX_BODIES; i++) {
        if (!bodies[i].active) {
            derivs[i].dPosition = (Vector3){0,0,0};
            derivs[i].dVelocity = (Vector3){0,0,0};
            continue;
        }

        derivs[i].dPosition = states[i].velocity; // dr/dt = v
        Vector3 force = { 0.0f, 0.0f, 0.0f };

        // Drag Force (Gas/Dust) - Accretion Disk Model
        if (enableDrag && bodies[0].active) {
            // Assume Body 0 is the center (Sun)
            float distToCenter = Vector3Distance(states[i].position, states[0].position);
            
            if (distToCenter > DISK_INNER_RADIUS && distToCenter < DISK_OUTER_RADIUS) {
                // Density function: Higher density closer to center
                float normalizedDist = (distToCenter - DISK_INNER_RADIUS) / (DISK_OUTER_RADIUS - DISK_INNER_RADIUS);
                float density = 1.0f - normalizedDist;
                if (density < 0.0f) density = 0.0f;
                
                float dragCoeff = 0.002f * density; 
                
                // F_drag = -c * rho * v
                Vector3 drag = Vector3Scale(states[i].velocity, -dragCoeff);
                force = Vector3Add(force, drag);
            }
        }

        for (int j = 0; j < MAX_BODIES; j++) {
            if (i == j || !bodies[j].active) continue;

            Vector3 direction = Vector3Subtract(states[j].position, states[i].position);
            float distance = Vector3Length(direction);
            
            float minDist = bodies[i].radius + bodies[j].radius; 
            if (distance < minDist) distance = minDist; 

            float forceMagnitude = (G * bodies[i].mass * bodies[j].mass) / (distance * distance);

            // Relativistic correction (Precession)
            Vector3 relVel = Vector3Subtract(states[j].velocity, states[i].velocity);
            // h = |r x v|
            Vector3 hVec = Vector3CrossProduct(direction, relVel);
            float h = Vector3Length(hVec);
            
            float correction = (3.0f * h * h) / (C_SPEED * C_SPEED * distance * distance);
            forceMagnitude *= (1.0f + correction);

            Vector3 forceVec = Vector3Scale(Vector3Normalize(direction), forceMagnitude);
            force = Vector3Add(force, forceVec);
        }
        derivs[i].dVelocity = Vector3Scale(force, 1.0f / bodies[i].mass); // dv/dt = a
    }
}

// Runge-Kutta 4th Order Integration
void integrateRK4(Body bodies[], float dt, bool enableDrag) {
    State initialStates[MAX_BODIES];
    getStates(bodies, initialStates);

    Derivative k1[MAX_BODIES], k2[MAX_BODIES], k3[MAX_BODIES], k4[MAX_BODIES];
    State tempStates[MAX_BODIES];

    // k1
    calculateDerivatives(initialStates, k1, bodies, enableDrag);

    // k2
    for(int i=0; i<MAX_BODIES; i++) {
        tempStates[i].position = Vector3Add(initialStates[i].position, Vector3Scale(k1[i].dPosition, dt * 0.5f));
        tempStates[i].velocity = Vector3Add(initialStates[i].velocity, Vector3Scale(k1[i].dVelocity, dt * 0.5f));
    }
    calculateDerivatives(tempStates, k2, bodies, enableDrag);

    // k3
    for(int i=0; i<MAX_BODIES; i++) {
        tempStates[i].position = Vector3Add(initialStates[i].position, Vector3Scale(k2[i].dPosition, dt * 0.5f));
        tempStates[i].velocity = Vector3Add(initialStates[i].velocity, Vector3Scale(k2[i].dVelocity, dt * 0.5f));
    }
    calculateDerivatives(tempStates, k3, bodies, enableDrag);

    // k4
    for(int i=0; i<MAX_BODIES; i++) {
        tempStates[i].position = Vector3Add(initialStates[i].position, Vector3Scale(k3[i].dPosition, dt));
        tempStates[i].velocity = Vector3Add(initialStates[i].velocity, Vector3Scale(k3[i].dVelocity, dt));
    }
    calculateDerivatives(tempStates, k4, bodies, enableDrag);

    // Update
    for(int i=0; i<MAX_BODIES; i++) {
        if (!bodies[i].active) continue;
        if (i == 0) continue; // Sun is stationary
        
        Vector3 dPos = Vector3Scale(Vector3Add(Vector3Add(k1[i].dPosition, Vector3Scale(k2[i].dPosition, 2.0f)), Vector3Add(Vector3Scale(k3[i].dPosition, 2.0f), k4[i].dPosition)), dt / 6.0f);
        Vector3 dVel = Vector3Scale(Vector3Add(Vector3Add(k1[i].dVelocity, Vector3Scale(k2[i].dVelocity, 2.0f)), Vector3Add(Vector3Scale(k3[i].dVelocity, 2.0f), k4[i].dVelocity)), dt / 6.0f);

        bodies[i].position = Vector3Add(bodies[i].position, dPos);
        bodies[i].velocity = Vector3Add(bodies[i].velocity, dVel);
    }
}

void handleCollisions(Body bodies[]) {
    for (int i = 0; i < MAX_BODIES; i++) {
        if (!bodies[i].active) continue;
        for (int j = i + 1; j < MAX_BODIES; j++) {
            if (!bodies[j].active) continue;

            if (CheckCollisionSpheres(bodies[i].position, bodies[i].radius, bodies[j].position, bodies[j].radius)) {
                // Merge j into i (Fusion)
                Body *b1 = &bodies[i];
                Body *b2 = &bodies[j];

                Vector3 momentum1 = Vector3Scale(b1->velocity, b1->mass);
                Vector3 momentum2 = Vector3Scale(b2->velocity, b2->mass);
                Vector3 totalMomentum = Vector3Add(momentum1, momentum2);
                float totalMass = b1->mass + b2->mass;

                if (i == 0) {
                    // Sun absorbs body, but stays stationary
                    b1->mass = totalMass;
                    b1->radius = cbrtf(powf(b1->radius, 3.0f) + powf(b2->radius, 3.0f));
                } else {
                    b1->velocity = Vector3Scale(totalMomentum, 1.0f / totalMass);
                    
                    // Weighted position
                    b1->position = Vector3Scale(Vector3Add(Vector3Scale(b1->position, b1->mass), Vector3Scale(b2->position, b2->mass)), 1.0f/totalMass);

                    b1->radius = cbrtf(powf(b1->radius, 3.0f) + powf(b2->radius, 3.0f));
                    b1->mass = totalMass;
                }
                
                b2->active = false;
            }
        }
    }
}

void explodeBody(Body bodies[], int index) {
    bodies[index].active = false;
    
    int fragments = 8;
    float newMass = bodies[index].mass / fragments;
    float newRadius = bodies[index].radius / 2.0f;
    if (newRadius < 2.0f) newRadius = 2.0f;
    
    for (int k = 0; k < fragments; k++) {
        for (int j = 0; j < MAX_BODIES; j++) {
            if (!bodies[j].active) {
                bodies[j].active = true;
                bodies[j].mass = newMass;
                bodies[j].radius = newRadius;
                bodies[j].color = (Color){ 
                    (unsigned char)fminf(255, bodies[index].color.r + GetRandomValue(-20, 20)),
                    (unsigned char)fminf(255, bodies[index].color.g + GetRandomValue(-20, 20)),
                    (unsigned char)fminf(255, bodies[index].color.b + GetRandomValue(-20, 20)),
                    255 
                };
                
                Vector3 offset = { (float)GetRandomValue(-5, 5), (float)GetRandomValue(-5, 5), (float)GetRandomValue(-5, 5) };
                bodies[j].position = Vector3Add(bodies[index].position, offset);
                
                Vector3 velSpread = { (float)GetRandomValue(-20, 20) / 10.0f, (float)GetRandomValue(-20, 20) / 10.0f, (float)GetRandomValue(-20, 20) / 10.0f };
                bodies[j].velocity = Vector3Add(bodies[index].velocity, velSpread);
                
                bodies[j].trailIndex = 0;
                for(int t=0; t<TRAIL_LENGTH; t++) bodies[j].trail[t] = bodies[j].position;
                break;
            }
        }
    }
}

void checkRocheLimit(Body bodies[]) {
    if (!bodies[0].active) return;
    
    for (int i = 1; i < MAX_BODIES; i++) {
        if (!bodies[i].active) continue;
        
        float dist = Vector3Distance(bodies[i].position, bodies[0].position);
        
        if (bodies[i].mass < 0.1f) continue; 
        
        float massRatio = bodies[0].mass / bodies[i].mass;
        float rocheLimit = 1.26f * bodies[i].radius * cbrtf(massRatio);
        
        if (bodies[i].radius > 3.0f && dist < rocheLimit) {
            explodeBody(bodies, i);
        }
    }
}

// Lagrange points are strictly defined for the restricted 3-body problem in a plane.
// We will project to the orbital plane of the largest planet for visualization.
void drawLagrangePoints(Body bodies[]) {
    if (!bodies[0].active) return;

    int heaviestIndex = -1;
    float maxMass = 0.0f;

    for (int i = 1; i < MAX_BODIES; i++) {
        if (bodies[i].active && bodies[i].mass > maxMass) {
            maxMass = bodies[i].mass;
            heaviestIndex = i;
        }
    }

    if (heaviestIndex == -1) return;

    Body *m1 = &bodies[0];
    Body *m2 = &bodies[heaviestIndex];

    Vector3 r1 = m1->position;
    Vector3 r2 = m2->position;
    Vector3 R_vec = Vector3Subtract(r2, r1);
    float R = Vector3Length(R_vec);
    Vector3 u = Vector3Scale(R_vec, 1.0f / R);

    float massRatio = m2->mass / m1->mass; 
    float hillRadius = R * cbrtf(massRatio / 3.0f);

    Vector3 l1 = Vector3Subtract(r2, Vector3Scale(u, hillRadius));
    Vector3 l2 = Vector3Add(r2, Vector3Scale(u, hillRadius));

    float l3_dist = R * (1.0f + (5.0f/12.0f) * massRatio);
    Vector3 l3 = Vector3Subtract(r1, Vector3Scale(u, l3_dist));

    // For L4/L5 we need a vector perpendicular to R_vec in the orbital plane.
    // We can use the velocity of m2 relative to m1 to find the orbital plane normal.
    Vector3 vRel = Vector3Subtract(m2->velocity, m1->velocity);
    Vector3 orbitalNormal = Vector3Normalize(Vector3CrossProduct(R_vec, vRel));
    // If velocity is parallel to position (falling straight in), this fails, but that's rare for planets.
    // Perpendicular vector in plane:
    Vector3 perp = Vector3Normalize(Vector3CrossProduct(orbitalNormal, u));

    // L4/L5 form equilateral triangles.
    // Position is R/2 along u, and sqrt(3)/2 * R along perp.
    float h_tri = 0.8660254f * R; // sin(60) * R
    
    Vector3 l4 = Vector3Add(Vector3Add(r1, Vector3Scale(u, R * 0.5f)), Vector3Scale(perp, h_tri));
    Vector3 l5 = Vector3Add(Vector3Add(r1, Vector3Scale(u, R * 0.5f)), Vector3Scale(perp, -h_tri));

    Color lColor = VIOLET;
    float lRadius = 5.0f;
    
    DrawSphere(l1, lRadius, lColor); 
    DrawSphere(l2, lRadius, lColor); 
    DrawSphere(l3, lRadius, lColor); 
    DrawSphere(l4, lRadius, lColor); 
    DrawSphere(l5, lRadius, lColor); 
}

// Helper to create a body in a stable circular orbit around a parent
Body createOrbitingBody(Body parent, float orbitRadius, float angleDeg, float mass, float radius, Color color) {
    float theta = angleDeg * DEG2RAD;
    
    // Position offset
    float dx = orbitRadius * cosf(theta);
    float dz = orbitRadius * sinf(theta);
    
    Vector3 position = Vector3Add(parent.position, (Vector3){ dx, 0.0f, dz });
    
    // Orbital velocity (circular) v = sqrt(GM / r)
    float vMag = sqrtf((G * parent.mass) / orbitRadius);
    
    // Velocity vector (tangent to circle)
    // If pos is (cos, sin), vel is (-sin, cos) for counter-clockwise orbit
    float vx = -vMag * sinf(theta);
    float vz = vMag * cosf(theta);
    
    Vector3 velocity = Vector3Add(parent.velocity, (Vector3){ vx, 0.0f, vz });
    
    Body b = {0};
    b.position = position;
    b.velocity = velocity;
    b.mass = mass;
    b.radius = radius;
    b.color = color;
    b.active = true;
    b.trailIndex = 0;
    for(int i=0; i<TRAIL_LENGTH; i++) b.trail[i] = position;
    
    return b;
}

void initBodies(Body bodies[]) {
    for(int i=0; i<MAX_BODIES; i++) bodies[i].active = false;

    // Sun
    bodies[0] = (Body){
        .position = { 0.0f, 0.0f, 0.0f },
        .velocity = { 0.0f, 0.0f, 0.0f },
        .mass = 10000.0f,
        .radius = 40.0f,
        .color = YELLOW,
        .active = true
    };
    
    // Planet 1
    bodies[1] = createOrbitingBody(bodies[0], 200.0f, 0.0f, 10.0f, 10.0f, BLUE);
    
    // Planet 2
    bodies[2] = createOrbitingBody(bodies[0], 500.0f, 0.0f, 8.0f, 8.0f, RED);
    
    // Planet 3 (Orange)
    bodies[3] = createOrbitingBody(bodies[0], 850.0f, 0.0f, 500.0f, 20.0f, ORANGE);
    
    // Moon of Planet 3
    bodies[4] = createOrbitingBody(bodies[3], 45.0f, 0.0f, 1.0f, 4.0f, WHITE);

    for (int i = 0; i < MAX_BODIES; i++) {
        for (int j = 0; j < TRAIL_LENGTH; j++) {
            bodies[i].trail[j] = bodies[i].position;
        }
        bodies[i].trailIndex = 0;
    }
}

typedef enum {
    STATE_SIMULATION,
    STATE_MENU,
    STATE_SETTINGS
} AppState;

void drawAccretionDisk(Vector3 center) {
    float inner = DISK_INNER_RADIUS;
    float outer = DISK_OUTER_RADIUS;
    int slices = 60;
    Color c = Fade(BLUE, 0.15f);

    rlBegin(RL_TRIANGLES);
    rlColor4ub(c.r, c.g, c.b, c.a);

    for (int i = 0; i < slices; i++) {
        float angle1 = (float)i / slices * 2.0f * PI;
        float angle2 = (float)(i + 1) / slices * 2.0f * PI;

        float sin1 = sinf(angle1);
        float cos1 = cosf(angle1);
        float sin2 = sinf(angle2);
        float cos2 = cosf(angle2);

        // Quad formed by 2 triangles
        // V1(in,1) -> V2(out,1) -> V3(in,2)
        rlVertex3f(center.x + inner*cos1, center.y, center.z + inner*sin1);
        rlVertex3f(center.x + outer*cos1, center.y, center.z + outer*sin1);
        rlVertex3f(center.x + inner*cos2, center.y, center.z + inner*sin2);

        // V2(out,1) -> V4(out,2) -> V3(in,2)
        rlVertex3f(center.x + outer*cos1, center.y, center.z + outer*sin1);
        rlVertex3f(center.x + outer*cos2, center.y, center.z + outer*sin2);
        rlVertex3f(center.x + inner*cos2, center.y, center.z + inner*sin2);
    }
    rlEnd();
}

void drawOrbit(Body bodies[], int targetIndex, int centerIndex) {
    if (!bodies[centerIndex].active || !bodies[targetIndex].active || targetIndex == centerIndex) return;

    Body *center = &bodies[centerIndex];
    Body *planet = &bodies[targetIndex];

    Vector3 rVec = Vector3Subtract(planet->position, center->position);
    Vector3 vVec = Vector3Subtract(planet->velocity, center->velocity);
    
    float r = Vector3Length(rVec);
    float v = Vector3Length(vVec);
    float mu = G * (center->mass + planet->mass);

    // Specific angular momentum h = r x v
    Vector3 hVec = Vector3CrossProduct(rVec, vVec);
    float h = Vector3Length(hVec);
    if (h < 0.1f) return; 

    // Eccentricity vector e = (v x h) / mu - r / |r|
    Vector3 vxh = Vector3CrossProduct(vVec, hVec);
    Vector3 eVec = Vector3Subtract(Vector3Scale(vxh, 1.0f/mu), Vector3Scale(rVec, 1.0f/r));
    float e = Vector3Length(eVec);

    // Semi-major axis a = 1 / (2/r - v^2/mu)
    float energy = v*v/2.0f - mu/r;
    if (fabs(energy) < 0.0001f) return; 
    float a = -mu / (2.0f * energy);

    if (e >= 1.0f || a < 0) return; // Hyperbolic/Parabolic

    // Basis vectors for orbital plane
    Vector3 n = Vector3Normalize(hVec); // Normal to plane
    Vector3 p; // Periapsis direction
    if (e > 0.001f) {
        p = Vector3Normalize(eVec);
    } else {
        if (fabs(n.y) < 0.9f) p = Vector3Normalize(Vector3CrossProduct((Vector3){0,1,0}, n));
        else p = Vector3Normalize(Vector3CrossProduct((Vector3){1,0,0}, n));
    }
    Vector3 q = Vector3CrossProduct(n, p); 

    // Draw Orbit Path
    rlBegin(RL_LINES);
    rlColor4ub(255, 255, 255, 60); 

    int segments = 100;
    Vector3 prevPos = {0};
    bool first = true;

    for (int i = 0; i <= segments; i++) {
        float theta = (float)i / segments * 2.0f * PI;
        float radius = a * (1.0f - e*e) / (1.0f + e * cosf(theta));
        
        Vector3 posInPlane = Vector3Add(Vector3Scale(p, radius * cosf(theta)), Vector3Scale(q, radius * sinf(theta)));
        Vector3 worldPos = Vector3Add(center->position, posInPlane);

        if (!first) {
            rlVertex3f(prevPos.x, prevPos.y, prevPos.z);
            rlVertex3f(worldPos.x, worldPos.y, worldPos.z);
        }
        prevPos = worldPos;
        first = false;
    }
    rlEnd();

    // Draw Periapsis (Green) and Apoapsis (Red)
    float r_peri = a * (1.0f - e);
    float r_apo = a * (1.0f + e);
    
    Vector3 posPeri = Vector3Add(center->position, Vector3Scale(p, r_peri));
    Vector3 posApo = Vector3Add(center->position, Vector3Scale(p, -r_apo)); 

    DrawSphere(posPeri, 3.0f, GREEN);
    DrawSphere(posApo, 3.0f, RED);
}

int findParentBody(Body bodies[], int subjectIndex) {
    if (subjectIndex == 0) return -1; 
    
    int bestParent = 0; 
    float maxForce = -1.0f;
    
    for (int i = 0; i < MAX_BODIES; i++) {
        if (i == subjectIndex || !bodies[i].active) continue;
        
        float dist = Vector3Distance(bodies[i].position, bodies[subjectIndex].position);
        if (dist < 1.0f) continue;
        
        float force = bodies[i].mass / (dist * dist);
        
        if (force > maxForce) {
            maxForce = force;
            bestParent = i;
        }
    }
    return bestParent;
}

void drawOrbitEditor(Body bodies[], int targetIndex, int screenWidth, int screenHeight) {
    if (targetIndex == 0 || !bodies[targetIndex].active) return;
    
    int parentIndex = findParentBody(bodies, targetIndex);
    if (parentIndex == -1) return; 
    
    Body *b = &bodies[targetIndex];
    Body *p = &bodies[parentIndex];
    
    // Calculate relative state
    Vector3 relPos = Vector3Subtract(b->position, p->position);
    Vector3 relVel = Vector3Subtract(b->velocity, p->velocity);
    
    float dist = Vector3Length(relPos);
    float speed = Vector3Length(relVel);
    float mu = G * (p->mass + b->mass);

    // Calculate orbital elements
    Vector3 hVec = Vector3CrossProduct(relPos, relVel);
    float h = Vector3Length(hVec);
    Vector3 n = Vector3Scale(hVec, 1.0f/h);
    
    Vector3 vxh = Vector3CrossProduct(relVel, hVec);
    Vector3 eVec = Vector3Subtract(Vector3Scale(vxh, 1.0f/mu), Vector3Scale(relPos, 1.0f/dist));
    float e = Vector3Length(eVec);
    
    // Handle circular orbits (e ~ 0)
    Vector3 eDir;
    if (e < 1e-4f) {
        e = 0.0f;
        eDir = Vector3Normalize(relPos); 
    } else {
        eDir = Vector3Normalize(eVec);
    }
    
    Vector3 qDir = Vector3CrossProduct(n, eDir);
    
    float energy = speed*speed/2.0f - mu/dist;
    float a = -mu / (2.0f * energy);
    
    // True Anomaly nu
    float nu = atan2f(Vector3DotProduct(relPos, qDir), Vector3DotProduct(relPos, eDir));

    // Calculate Angle (Argument of Periapsis relative to reference)
    Vector3 worldUp = {0, 1, 0};
    if (fabsf(n.y) > 0.95f) worldUp = (Vector3){1, 0, 0};
    Vector3 u = Vector3Normalize(Vector3CrossProduct(worldUp, n));
    Vector3 v = Vector3CrossProduct(n, u);
    float angle = atan2f(Vector3DotProduct(eDir, v), Vector3DotProduct(eDir, u));
    if (angle < 0) angle += 2*PI;

    // Calculate Inclination
    Vector3 refNormal = {0, 1, 0};
    float inclination = acosf(Vector3DotProduct(n, refNormal));
    Vector3 nodeVec = Vector3CrossProduct(refNormal, n);
    if (Vector3Length(nodeVec) < 0.001f) nodeVec = (Vector3){1, 0, 0};
    nodeVec = Vector3Normalize(nodeVec);

    // Derived parameters
    float r_peri = a * (1.0f - e);
    float r_apo = a * (1.0f + e);
    float period = 2.0f * PI * sqrtf(powf(a, 3.0f) / mu);

    // UI Layout
    int uiWidth = 280;
    int uiHeight = 400; // Increased height
    int uiX = screenWidth - uiWidth - 10;
    int uiY = screenHeight - uiHeight - 10;
    Rectangle uiRect = { uiX, uiY, uiWidth, uiHeight };
    
    DrawRectangleRec(uiRect, Fade(BLACK, 0.8f));
    DrawRectangleLinesEx(uiRect, 1, DARKGRAY);
    
    int startX = uiX + 10;
    int startY = uiY + 10;
    
    DrawText("ORBIT EDITOR", startX, startY, 20, WHITE);
    startY += 30;
    
    // Info
    DrawText(TextFormat("Period: %.1f s", period), startX, startY, 10, YELLOW);
    startY += 20;

    // Helper to apply changes
    float new_a = a;
    float new_e = e;
    float new_angle = angle;
    float new_inclination = inclination;
    bool changed = false;

    // 1. Eccentricity
    DrawText(TextFormat("Eccentricity: %.3f", e), startX, startY, 10, LIGHTGRAY);
    Rectangle eRect = { startX, startY + 15, 260, 20 };
    DrawRectangleRec(eRect, DARKGRAY);
    DrawRectangle(eRect.x, eRect.y, e * eRect.width, eRect.height, ORANGE);
    if (CheckCollisionPointRec(GetMousePosition(), eRect) && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        float val = (GetMouseX() - eRect.x) / eRect.width;
        if (val < 0.0f) val = 0.0f;
        if (val > 0.95f) val = 0.95f;
        new_e = val;
        changed = true;
    }
    startY += 45;

    // 2. Semi-major Axis
    DrawText(TextFormat("Semi-major Axis: %.1f", a), startX, startY, 10, LIGHTGRAY);
    Rectangle aRect = { startX, startY + 15, 260, 20 };
    DrawRectangleRec(aRect, DARKGRAY);
    float normA = a / 1000.0f;
    if (normA > 1.0f) normA = 1.0f;
    DrawRectangle(aRect.x, aRect.y, normA * aRect.width, aRect.height, BLUE);
    if (CheckCollisionPointRec(GetMousePosition(), aRect) && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        float val = (GetMouseX() - aRect.x) / aRect.width;
        new_a = val * 1000.0f;
        if (new_a < 20.0f) new_a = 20.0f;
        changed = true;
    }
    startY += 45;

    // 3. Periapsis
    DrawText(TextFormat("Periapsis: %.1f", r_peri), startX, startY, 10, GREEN);
    Rectangle pRect = { startX, startY + 15, 260, 20 };
    DrawRectangleRec(pRect, DARKGRAY);
    float normP = r_peri / 1000.0f;
    if (normP > 1.0f) normP = 1.0f;
    DrawRectangle(pRect.x, pRect.y, normP * pRect.width, pRect.height, GREEN);
    if (CheckCollisionPointRec(GetMousePosition(), pRect) && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        float val = (GetMouseX() - pRect.x) / pRect.width;
        float new_p = val * 1000.0f;
        if (new_p < 10.0f) new_p = 10.0f;
        if (new_p >= r_apo) new_p = r_apo - 1.0f; // Clamp
        
        new_a = (new_p + r_apo) / 2.0f;
        new_e = (r_apo - new_p) / (r_apo + new_p);
        changed = true;
    }
    startY += 45;

    // 4. Apoapsis
    DrawText(TextFormat("Apoapsis: %.1f", r_apo), startX, startY, 10, RED);
    Rectangle apRect = { startX, startY + 15, 260, 20 };
    DrawRectangleRec(apRect, DARKGRAY);
    float normAp = r_apo / 1000.0f;
    if (normAp > 1.0f) normAp = 1.0f;
    DrawRectangle(apRect.x, apRect.y, normAp * apRect.width, apRect.height, RED);
    if (CheckCollisionPointRec(GetMousePosition(), apRect) && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        float val = (GetMouseX() - apRect.x) / apRect.width;
        float new_ap = val * 1000.0f;
        if (new_ap <= r_peri) new_ap = r_peri + 1.0f; // Clamp
        
        new_a = (r_peri + new_ap) / 2.0f;
        new_e = (new_ap - r_peri) / (new_ap + r_peri);
        changed = true;
    }
    startY += 45;

    // 5. Rotation (Arg. Periapsis)
    DrawText(TextFormat("Rotation: %.0f deg", angle * RAD2DEG), startX, startY, 10, PURPLE);
    Rectangle angRect = { startX, startY + 15, 260, 20 };
    DrawRectangleRec(angRect, DARKGRAY);
    float normAng = angle / (2*PI);
    DrawRectangle(angRect.x, angRect.y, normAng * angRect.width, angRect.height, PURPLE);
    if (CheckCollisionPointRec(GetMousePosition(), angRect) && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        float val = (GetMouseX() - angRect.x) / angRect.width;
        new_angle = val * 2 * PI;
        changed = true;
    }
    startY += 45;

    // 6. Inclination
    DrawText(TextFormat("Inclination: %.0f deg", inclination * RAD2DEG), startX, startY, 10, MAGENTA);
    Rectangle incRect = { startX, startY + 15, 260, 20 };
    DrawRectangleRec(incRect, DARKGRAY);
    float normInc = inclination / PI; // 0 to 180
    DrawRectangle(incRect.x, incRect.y, normInc * incRect.width, incRect.height, MAGENTA);
    if (CheckCollisionPointRec(GetMousePosition(), incRect) && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        float val = (GetMouseX() - incRect.x) / incRect.width;
        new_inclination = val * PI;
        changed = true;
    }
    startY += 45;

    if (changed) {
        // Reconstruct vectors
        Vector3 new_eDir = Vector3Add(Vector3Scale(u, cosf(new_angle)), Vector3Scale(v, sinf(new_angle)));
        Vector3 new_qDir = Vector3CrossProduct(n, new_eDir);
        
        float slr = new_a * (1.0f - new_e * new_e);
        float new_r = slr / (1.0f + new_e * cosf(nu));
        
        float v_radial = sqrtf(mu/slr) * new_e * sinf(nu);
        float v_tangential = sqrtf(mu/slr) * (1.0f + new_e * cosf(nu));
        
        Vector3 r_hat = Vector3Add(Vector3Scale(new_eDir, cosf(nu)), Vector3Scale(new_qDir, sinf(nu)));
        Vector3 t_hat = Vector3CrossProduct(n, r_hat);
        
        Vector3 pos = Vector3Scale(r_hat, new_r);
        Vector3 vel = Vector3Add(Vector3Scale(r_hat, v_radial), Vector3Scale(t_hat, v_tangential));

        // Apply Inclination Change
        if (fabsf(new_inclination - inclination) > 0.001f) {
            float deltaInc = new_inclination - inclination;
            Quaternion q = QuaternionFromAxisAngle(nodeVec, deltaInc);
            pos = Vector3RotateByQuaternion(pos, q);
            vel = Vector3RotateByQuaternion(vel, q);
        }

        b->position = Vector3Add(p->position, pos);
        b->velocity = Vector3Add(p->velocity, vel);
    }
}

int main(void)
{
    const int screenWidth = 1200;
    const int screenHeight = 800;

    InitWindow(screenWidth, screenHeight, "Solar System Simulation 3D - Solaray");
    SetTargetFPS(60);
    SetExitKey(KEY_NULL);

    // 3D Camera
    Camera3D camera = { 0 };
    camera.position = (Vector3){ 0.0f, 400.0f, 400.0f };
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    // Orbit Camera State
    float camDist = 600.0f;
    Vector2 camAngle = { 0.0f, 1.0f }; // Theta (XZ plane), Phi (Y axis)
    int cameraTargetIndex = 0; // Default to Sun

    float timeScale = 1.0f;
    bool enableDrag = false;
    bool showLagrange = false;
    bool enableRoche = true;
    bool relativeView = false;

    Body bodies[MAX_BODIES];
    initBodies(bodies);

    bool creationMode = false;
    bool isDragging = false;
    Vector3 dragStartPos = {0};
    float newBodyMass = 10.0f;
    float spawnHeight = 0.0f; // For 3D spawning
    float spawnAngle = 0.0f; // Vertical launch angle in degrees

    // Lighting Shader
    Shader lightShader = LoadShaderFromMemory(
        "#version 330\n"
        "in vec3 vertexPosition;\n"
        "in vec2 vertexTexCoord;\n"
        "in vec3 vertexNormal;\n"
        "in vec4 vertexColor;\n"
        "out vec3 fragPosition;\n"
        "out vec3 fragNormal;\n"
        "out vec4 fragColor;\n"
        "uniform mat4 mvp;\n"
        "uniform mat4 matModel;\n"
        "uniform mat4 matNormal;\n"
        "void main() {\n"
        "    fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));\n"
        "    fragNormal = normalize(vec3(matNormal * vec4(vertexNormal, 1.0)));\n"
        "    fragColor = vertexColor;\n"
        "    gl_Position = mvp * vec4(vertexPosition, 1.0);\n"
        "}\n",
        
        "#version 330\n"
        "in vec3 fragPosition;\n"
        "in vec3 fragNormal;\n"
        "in vec4 fragColor;\n"
        "out vec4 finalColor;\n"
        "uniform vec3 lightPos;\n"
        "uniform vec4 lightColor;\n"
        "uniform vec4 ambientColor;\n"
        "uniform vec4 objectColor;\n"
        "uniform vec3 viewPos;\n"
        "void main() {\n"
        "    vec3 normal = normalize(fragNormal);\n"
        "    vec3 lightDir = normalize(lightPos - fragPosition);\n"
        "    float diff = max(dot(normal, lightDir), 0.0);\n"
        "    vec3 diffuse = diff * lightColor.rgb;\n"
        "    vec3 ambient = ambientColor.rgb;\n"
        "    vec3 result = (ambient + diffuse) * objectColor.rgb;\n"
        "    finalColor = vec4(result, objectColor.a);\n"
        "}\n"
    );

    // Get shader locations
    int lightPosLoc = GetShaderLocation(lightShader, "lightPos");
    int viewPosLoc = GetShaderLocation(lightShader, "viewPos");
    int ambientLoc = GetShaderLocation(lightShader, "ambientColor");
    int lightColorLoc = GetShaderLocation(lightShader, "lightColor");
    int objectColorLoc = GetShaderLocation(lightShader, "objectColor");

    // Set shader constants
    float ambient[4] = { 0.1f, 0.1f, 0.1f, 1.0f };
    SetShaderValue(lightShader, ambientLoc, ambient, SHADER_UNIFORM_VEC4);
    float lightColor[4] = { 1.0f, 1.0f, 0.9f, 1.0f }; // Sun color
    SetShaderValue(lightShader, lightColorLoc, lightColor, SHADER_UNIFORM_VEC4);

    // Create Sphere Model for Planets
    Mesh sphereMesh = GenMeshSphere(1.0f, 16, 16);
    Model sphereModel = LoadModelFromMesh(sphereMesh);
    sphereModel.materials[0].shader = lightShader;

    AppState currentState = STATE_SIMULATION;
    bool shouldExit = false;

    while (!shouldExit && !WindowShouldClose()) {
        if (IsKeyPressed(KEY_ESCAPE)) {
            if (currentState == STATE_SIMULATION) currentState = STATE_MENU;
            else if (currentState == STATE_MENU) currentState = STATE_SIMULATION;
            else if (currentState == STATE_SETTINGS) currentState = STATE_MENU;
        }

        if (currentState == STATE_SIMULATION) {
            if (IsKeyPressed(KEY_N)) {
                creationMode = !creationMode;
                spawnHeight = 0.0f;
            }

            // Camera Controls
            if (!creationMode) {
                if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                    Vector2 delta = GetMouseDelta();
                    camAngle.x -= delta.x * 0.005f;
                    camAngle.y -= delta.y * 0.005f;
                    if (camAngle.y < 0.01f) camAngle.y = 0.01f;
                    if (camAngle.y > PI - 0.01f) camAngle.y = PI - 0.01f;
                }
                camDist -= GetMouseWheelMove() * 50.0f;
                if (camDist < 50.0f) camDist = 50.0f;
            }

            // Update Camera Position
            if (!bodies[cameraTargetIndex].active) cameraTargetIndex = 0;
            Vector3 targetPos = bodies[cameraTargetIndex].position;
            camera.target = targetPos;
            
            camera.position.x = targetPos.x + camDist * sinf(camAngle.y) * cosf(camAngle.x);
            camera.position.y = targetPos.y + camDist * cosf(camAngle.y);
            camera.position.z = targetPos.z + camDist * sinf(camAngle.y) * sinf(camAngle.x);

            // Update Shader Uniforms
            SetShaderValue(lightShader, viewPosLoc, &camera.position, SHADER_UNIFORM_VEC3);
            if (bodies[0].active) {
                SetShaderValue(lightShader, lightPosLoc, &bodies[0].position, SHADER_UNIFORM_VEC3);
            }

            if (creationMode) {
                // Check if mouse is over UI (simple rect check)
                int uiWidth = 280;
                int uiHeight = 200;
                int uiX = screenWidth - uiWidth - 10;
                int uiY = 10;
                Rectangle uiRect = { uiX, uiY, uiWidth, uiHeight };
                bool mouseOverUI = CheckCollisionPointRec(GetMousePosition(), uiRect);

                if (!mouseOverUI) {
                    Ray ray = GetMouseRay(GetMousePosition(), camera);
                    
                    if (fabs(ray.direction.y) > 0.001f) {
                        float t = (spawnHeight - ray.position.y) / ray.direction.y;
                        if (t >= 0) {
                            Vector3 mouseWorldPos = Vector3Add(ray.position, Vector3Scale(ray.direction, t));
                            
                            if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
                                isDragging = true;
                                dragStartPos = mouseWorldPos;
                            }

                            if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT) && isDragging) {
                                isDragging = false;
                                for (int i = 0; i < MAX_BODIES; i++) {
                                    if (!bodies[i].active) {
                                        bodies[i].active = true;
                                        bodies[i].position = dragStartPos;
                                        
                                        Vector3 dragVec = Vector3Subtract(mouseWorldPos, dragStartPos);
                                        float speed = Vector3Length(dragVec);
                                        Vector3 dirXZ = Vector3Normalize(dragVec);
                                        
                                        // Calculate 3D velocity based on angle
                                        float angleRad = spawnAngle * DEG2RAD;
                                        float vy = speed * sinf(angleRad);
                                        float vxz = speed * cosf(angleRad);
                                        
                                        bodies[i].velocity = (Vector3){ vxz * dirXZ.x, vy, vxz * dirXZ.z };
                                        
                                        bodies[i].mass = newBodyMass;
                                        bodies[i].radius = sqrtf(newBodyMass) * 3.0f;
                                        if (bodies[i].radius < 5.0f) bodies[i].radius = 5.0f;
                                        bodies[i].color = (Color){ GetRandomValue(100, 255), GetRandomValue(100, 255), GetRandomValue(100, 255), 255 };
                                        
                                        for(int t=0; t<TRAIL_LENGTH; t++) bodies[i].trail[t] = bodies[i].position;
                                        bodies[i].trailIndex = 0;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (IsKeyPressed(KEY_RIGHT)) timeScale *= 2.0f;
            if (IsKeyPressed(KEY_LEFT)) timeScale *= 0.5f;
            if (IsKeyPressed(KEY_SPACE)) timeScale = (timeScale == 0.0f) ? 1.0f : 0.0f;

            float dt = GetFrameTime() * timeScale;
            
            int subSteps = (int)ceilf(fabs(timeScale)); 
            if (subSteps < 1) subSteps = 1;
            float subDt = dt / subSteps;

            for (int step = 0; step < subSteps; step++) {
                integrateRK4(bodies, subDt, enableDrag);
                handleCollisions(bodies);
                if (enableRoche) checkRocheLimit(bodies);
            }

            if (timeScale != 0.0f) {
                for (int i = 0; i < MAX_BODIES; i++) {
                    if (!bodies[i].active) continue;
                    bodies[i].trail[bodies[i].trailIndex] = bodies[i].position;
                    bodies[i].trailIndex = (bodies[i].trailIndex + 1) % TRAIL_LENGTH;
                }
            }
        }

        BeginDrawing();
        ClearBackground(BLACK);
        
        BeginMode3D(camera);
        
        // Custom Grid
        rlSetClipPlanes(1.0f, 10000.0f); // Increase draw distance
        
        int slices = 100;
        float spacing = 50.0f;
        float halfSize = slices * spacing / 2.0f;
        Color gridColor = Fade(LIGHTGRAY, 0.2f); // 20% visible
        float gridY = -1.0f; // Slightly below 0 to avoid z-fighting with accretion disk

        rlBegin(RL_LINES);
        rlColor4ub(gridColor.r, gridColor.g, gridColor.b, gridColor.a);
        for (int i = -slices/2; i <= slices/2; i++) {
            // Lines parallel to Z axis
            rlVertex3f(i * spacing, gridY, -halfSize);
            rlVertex3f(i * spacing, gridY, halfSize);

            // Lines parallel to X axis
            rlVertex3f(-halfSize, gridY, i * spacing);
            rlVertex3f(halfSize, gridY, i * spacing);
        }
        rlEnd();

        if (enableDrag && bodies[0].active) {
            drawAccretionDisk(bodies[0].position);
        }

        if (cameraTargetIndex != 0 && bodies[cameraTargetIndex].active) {
            if (relativeView) {
                int parent = findParentBody(bodies, cameraTargetIndex);
                if (parent != -1) drawOrbit(bodies, cameraTargetIndex, parent);
            } else {
                drawOrbit(bodies, cameraTargetIndex, 0);
            }
        }

        for (int i = 0; i < MAX_BODIES; i++) {
            if (!bodies[i].active) continue;
            for (int j = 0; j < TRAIL_LENGTH - 1; j++) {
                int idx = (bodies[i].trailIndex + j) % TRAIL_LENGTH;
                int nextIdx = (idx + 1) % TRAIL_LENGTH;
                if (Vector3DistanceSqr(bodies[i].trail[idx], bodies[i].trail[nextIdx]) > 0.001f)
                    DrawLine3D(bodies[i].trail[idx], bodies[i].trail[nextIdx], Fade(bodies[i].color, 0.4f));
            }
        }
        
        for (int i = 0; i < MAX_BODIES; i++) {
            if (!bodies[i].active) continue;
            
            if (i == 0) {
                // Sun is unlit
                DrawSphere(bodies[i].position, bodies[i].radius, bodies[i].color);
            } else {
                // Planets are lit
                Vector3 scale = { bodies[i].radius, bodies[i].radius, bodies[i].radius };
                float col[4] = { bodies[i].color.r/255.0f, bodies[i].color.g/255.0f, bodies[i].color.b/255.0f, bodies[i].color.a/255.0f };
                SetShaderValue(lightShader, objectColorLoc, col, SHADER_UNIFORM_VEC4);
                DrawModelEx(sphereModel, bodies[i].position, (Vector3){0,1,0}, 0.0f, scale, WHITE);
            }
        }

        if (creationMode && currentState == STATE_SIMULATION) {
            // 3D Visualization of Drag
            int uiWidth = 280;
            int uiHeight = 200;
            int uiX = screenWidth - uiWidth - 10;
            int uiY = 10;
            Rectangle uiRect = { uiX, uiY, uiWidth, uiHeight };

            Ray ray = GetMouseRay(GetMousePosition(), camera);
            if (fabs(ray.direction.y) > 0.001f) {
                float t = (spawnHeight - ray.position.y) / ray.direction.y;
                if (t >= 0) {
                    Vector3 mouseWorldPos = Vector3Add(ray.position, Vector3Scale(ray.direction, t));
                    
                    if (!CheckCollisionPointRec(GetMousePosition(), uiRect)) {
                        if (isDragging) {
                            // Calculate end point based on angle
                            Vector3 dragVec = Vector3Subtract(mouseWorldPos, dragStartPos);
                            float speed = Vector3Length(dragVec);
                            Vector3 dirXZ = Vector3Normalize(dragVec);
                            float angleRad = spawnAngle * DEG2RAD;
                            float vy = speed * sinf(angleRad);
                            float vxz = speed * cosf(angleRad);
                            Vector3 velocity = { vxz * dirXZ.x, vy, vxz * dirXZ.z };
                            
                            // Draw trajectory preview (simple line)
                            Vector3 endPos = Vector3Add(dragStartPos, velocity); // Just a visual vector
                            DrawLine3D(dragStartPos, endPos, RED);
                            DrawSphere(dragStartPos, 5.0f, GREEN);
                            
                            // Draw height reference
                            DrawLine3D(dragStartPos, (Vector3){dragStartPos.x, 0, dragStartPos.z}, Fade(GREEN, 0.3f));
                        } else {
                            DrawSphereWires(mouseWorldPos, sqrtf(newBodyMass)*3.0f, 8, 8, Fade(GRAY, 0.5f));
                        }
                    }
                }
            }
        }

        if (showLagrange) {
            drawLagrangePoints(bodies);
        }
        
        EndMode3D();
        
        // UI and 2D Overlays
        if (cameraTargetIndex != 0 && currentState == STATE_SIMULATION) {
            drawOrbitEditor(bodies, cameraTargetIndex, screenWidth, screenHeight);
        }

        if (creationMode && currentState == STATE_SIMULATION) {
            int uiWidth = 280;
            int uiHeight = 200;
            int uiX = screenWidth - uiWidth - 10;
            int uiY = 10;
            Rectangle uiRect = { uiX, uiY, uiWidth, uiHeight };
            
            DrawRectangleRec(uiRect, Fade(BLACK, 0.8f));
            DrawRectangleLinesEx(uiRect, 1, DARKGRAY);
            
            int startY = uiY + 10;
            int startX = uiX + 10;
            
            DrawText("PLANET CREATOR", startX, startY, 20, WHITE);
            startY += 30;

            // Mass
            DrawText(TextFormat("Mass: %.1f", newBodyMass), startX, startY, 10, LIGHTGRAY);
            Rectangle massRect = { startX, startY + 15, 260, 20 };
            DrawRectangleRec(massRect, DARKGRAY);
            DrawRectangle(massRect.x, massRect.y, (newBodyMass / 100.0f) * massRect.width, massRect.height, BLUE);
            if (CheckCollisionPointRec(GetMousePosition(), massRect) && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                float val = (GetMouseX() - massRect.x) / massRect.width;
                if (val < 0) val = 0; 
                if (val > 1) val = 1;
                newBodyMass = val * 100.0f;
                if (newBodyMass < 1.0f) newBodyMass = 1.0f;
            }
            startY += 45;

            // Height
            DrawText(TextFormat("Height Offset: %.1f", spawnHeight), startX, startY, 10, LIGHTGRAY);
            Rectangle heightRect = { startX, startY + 15, 260, 20 };
            DrawRectangleRec(heightRect, DARKGRAY);
            float normHeight = (spawnHeight + 200.0f) / 400.0f;
            DrawRectangle(heightRect.x, heightRect.y, normHeight * heightRect.width, heightRect.height, GREEN);
            if (CheckCollisionPointRec(GetMousePosition(), heightRect) && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                float val = (GetMouseX() - heightRect.x) / heightRect.width;
                if (val < 0) val = 0; 
                if (val > 1) val = 1;
                spawnHeight = (val * 400.0f) - 200.0f;
            }
            startY += 45;

            // Angle
            DrawText(TextFormat("Launch Angle: %.1f deg", spawnAngle), startX, startY, 10, LIGHTGRAY);
            Rectangle angleRect = { startX, startY + 15, 260, 20 };
            DrawRectangleRec(angleRect, DARKGRAY);
            float normAngle = (spawnAngle + 90.0f) / 180.0f;
            DrawRectangle(angleRect.x, angleRect.y, normAngle * angleRect.width, angleRect.height, RED);
            if (CheckCollisionPointRec(GetMousePosition(), angleRect) && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                float val = (GetMouseX() - angleRect.x) / angleRect.width;
                if (val < 0) val = 0; 
                if (val > 1) val = 1;
                spawnAngle = (val * 180.0f) - 90.0f;
            }
            
            // Instructions
            startY += 40;
            DrawText("Click & Drag in space to launch", startX, startY, 10, GRAY);
        }

        // ... (Keep existing UI logic mostly, but remove 2D specific hover for now or adapt it)
        // Adapting hover to 3D:
        if (currentState == STATE_SIMULATION) {
            Ray ray = GetMouseRay(GetMousePosition(), camera);
            int hitIndex = -1;
            float minHitDist = 1e9f;

            for (int i = 0; i < MAX_BODIES; i++) {
                if (!bodies[i].active) continue;
                RayCollision collision = GetRayCollisionSphere(ray, bodies[i].position, bodies[i].radius);
                if (collision.hit) {
                    if (collision.distance < minHitDist) {
                        minHitDist = collision.distance;
                        hitIndex = i;
                    }
                }
            }

            if (hitIndex != -1) {
                if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && !creationMode) {
                    cameraTargetIndex = hitIndex;
                }

                Body *b = &bodies[hitIndex];
                float speed = Vector3Length(b->velocity);
                float distToSun = Vector3Distance(b->position, bodies[0].position);
                float kineticE = 0.5f * b->mass * speed * speed;
                
                char infoText[512];
                sprintf(infoText, "Mass: %.1f\nSpeed: %.1f\nDist to Sun: %.1f\nKinetic E: %.1e\nPos: (%.0f, %.0f, %.0f)", 
                        b->mass, speed, distToSun, kineticE, b->position.x, b->position.y, b->position.z);
                
                Vector2 screenPos = GetWorldToScreen(b->position, camera);
                DrawRectangle(screenPos.x + 20, screenPos.y - 60, 220, 110, Fade(DARKGRAY, 0.9f));
                DrawRectangleLines(screenPos.x + 20, screenPos.y - 60, 220, 110, WHITE);
                DrawText(infoText, screenPos.x + 25, screenPos.y - 55, 10, WHITE);
            }
        }

        // ... (Rest of UI)
        DrawFPS(10, 10);
        DrawText(TextFormat("Time Scale: %.2fx", timeScale), 10, 30, 20, WHITE);
        DrawText("RK4 | N-Body | Collisions | Relativistic Precession", 10, 50, 20, GREEN);
        
        int statusY = 70;
        int x = 10;
        
        DrawText("Roche Limit:", x, statusY, 20, LIGHTGRAY);
        x += MeasureText("Roche Limit: ", 20);
        DrawText(enableRoche ? "ON" : "OFF", x, statusY, 20, enableRoche ? GREEN : RED);
        x += MeasureText("ON ", 20) + 10;
        
        DrawText("| Accretion:", x, statusY, 20, LIGHTGRAY);
        x += MeasureText("| Accretion: ", 20);
        DrawText(enableDrag ? "ON" : "OFF", x, statusY, 20, enableDrag ? GREEN : RED);
        x += MeasureText("ON ", 20) + 10;

        DrawText("| Lagrange:", x, statusY, 20, LIGHTGRAY);
        x += MeasureText("| Lagrange: ", 20);
        DrawText(showLagrange ? "ON" : "OFF", x, statusY, 20, showLagrange ? GREEN : RED);
        x += MeasureText("ON ", 20) + 10;

        DrawText("| Create (N):", x, statusY, 20, LIGHTGRAY);
        x += MeasureText("| Create (N): ", 20);
        DrawText(creationMode ? "ON" : "OFF", x, statusY, 20, creationMode ? GREEN : RED);
        
        x += MeasureText("ON ", 20) + 20;
        Rectangle btnRel = { x, statusY - 2, 140, 24 };
        bool hoverRel = CheckCollisionPointRec(GetMousePosition(), btnRel);
        DrawRectangleRec(btnRel, hoverRel ? GRAY : DARKGRAY);
        DrawText("Rel. View", btnRel.x + 10, btnRel.y + 2, 20, relativeView ? GREEN : WHITE);
        if (hoverRel && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) relativeView = !relativeView;
        
        if (creationMode && currentState == STATE_SIMULATION) {
            DrawText(TextFormat("Creation Mode: Click & Drag to launch. Scroll: Mass %.1f", newBodyMass), 10, 100, 20, YELLOW);
        }

        // Menu and Settings (Same as before)
        if (currentState == STATE_MENU) {
            DrawRectangle(0, 0, screenWidth, screenHeight, Fade(BLACK, 0.7f));
            int menuX = screenWidth / 2 - 100;
            int menuY = screenHeight / 2 - 100;
            
            DrawText("PAUSED", screenWidth/2 - MeasureText("PAUSED", 40)/2, menuY - 60, 40, WHITE);

            Vector2 mouse = GetMousePosition();
            
            Rectangle btnResume = { menuX, menuY, 200, 40 };
            bool hoverResume = CheckCollisionPointRec(mouse, btnResume);
            DrawRectangleRec(btnResume, hoverResume ? GRAY : DARKGRAY);
            DrawText("Resume", btnResume.x + 20, btnResume.y + 10, 20, WHITE);
            if (hoverResume && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) currentState = STATE_SIMULATION;

            Rectangle btnReset = { menuX, menuY + 50, 200, 40 };
            bool hoverReset = CheckCollisionPointRec(mouse, btnReset);
            DrawRectangleRec(btnReset, hoverReset ? GRAY : DARKGRAY);
            DrawText("Reset", btnReset.x + 20, btnReset.y + 10, 20, WHITE);
            if (hoverReset && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
                initBodies(bodies);
                currentState = STATE_SIMULATION;
            }

            Rectangle btnSettings = { menuX, menuY + 100, 200, 40 };
            bool hoverSettings = CheckCollisionPointRec(mouse, btnSettings);
            DrawRectangleRec(btnSettings, hoverSettings ? GRAY : DARKGRAY);
            DrawText("Settings", btnSettings.x + 20, btnSettings.y + 10, 20, WHITE);
            if (hoverSettings && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) currentState = STATE_SETTINGS;

            Rectangle btnQuit = { menuX, menuY + 150, 200, 40 };
            bool hoverQuit = CheckCollisionPointRec(mouse, btnQuit);
            DrawRectangleRec(btnQuit, hoverQuit ? RED : MAROON);
            DrawText("Quit", btnQuit.x + 20, btnQuit.y + 10, 20, WHITE);
            if (hoverQuit && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) shouldExit = true;
        }
        else if (currentState == STATE_SETTINGS) {
            DrawRectangle(0, 0, screenWidth, screenHeight, Fade(BLACK, 0.8f));
            int menuX = screenWidth / 2 - 150;
            int menuY = screenHeight / 2 - 100;
            
            DrawText("SETTINGS", screenWidth/2 - MeasureText("SETTINGS", 40)/2, menuY - 60, 40, WHITE);
            Vector2 mouse = GetMousePosition();

            Rectangle btnRoche = { menuX, menuY, 300, 40 };
            bool hoverRoche = CheckCollisionPointRec(mouse, btnRoche);
            DrawRectangleRec(btnRoche, hoverRoche ? GRAY : DARKGRAY);
            DrawText(TextFormat("Roche Limit: %s", enableRoche ? "ON" : "OFF"), btnRoche.x + 20, btnRoche.y + 10, 20, enableRoche ? GREEN : RED);
            if (hoverRoche && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) enableRoche = !enableRoche;

            Rectangle btnDrag = { menuX, menuY + 50, 300, 40 };
            bool hoverDrag = CheckCollisionPointRec(mouse, btnDrag);
            DrawRectangleRec(btnDrag, hoverDrag ? GRAY : DARKGRAY);
            DrawText(TextFormat("Accretion Drag: %s", enableDrag ? "ON" : "OFF"), btnDrag.x + 20, btnDrag.y + 10, 20, enableDrag ? GREEN : RED);
            if (hoverDrag && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) enableDrag = !enableDrag;

            Rectangle btnLag = { menuX, menuY + 100, 300, 40 };
            bool hoverLag = CheckCollisionPointRec(mouse, btnLag);
            DrawRectangleRec(btnLag, hoverLag ? GRAY : DARKGRAY);
            DrawText(TextFormat("Lagrange Points: %s", showLagrange ? "ON" : "OFF"), btnLag.x + 20, btnLag.y + 10, 20, showLagrange ? GREEN : RED);
            if (hoverLag && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) showLagrange = !showLagrange;

            Rectangle btnBack = { menuX, menuY + 160, 300, 40 };
            bool hoverBack = CheckCollisionPointRec(mouse, btnBack);
            DrawRectangleRec(btnBack, hoverBack ? GRAY : DARKGRAY);
            DrawText("Back", btnBack.x + 120, btnBack.y + 10, 20, WHITE);
            if (hoverBack && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) currentState = STATE_MENU;
        }

        EndDrawing();
    }
    CloseWindow();
    return 0;
}