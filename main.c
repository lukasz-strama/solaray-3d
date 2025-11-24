#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define RAYGUI_IMPLEMENTATION
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <raygui.h>
#pragma GCC diagnostic pop

#define G_REAL 6.67430e-11
#define M_SUN 1.989e30
#define M_EARTH 5.972e24
#define AU 1.496e11
#define C_SPEED_REAL 299792458.0

// Use Real Constants directly
#define G G_REAL
#define C_SPEED C_SPEED_REAL

#define MAX_BODIES 500
#define TRAIL_LENGTH 200
// Disk radii in meters (approximate for solar system scale)
#define DISK_INNER_RADIUS (0.3 * AU)
#define DISK_OUTER_RADIUS (50.0 * AU)

typedef struct {
    double x;
    double y;
    double z;
} DVector3;

// Double precision math helpers
DVector3 DVector3Add(DVector3 v1, DVector3 v2) { return (DVector3){ v1.x + v2.x, v1.y + v2.y, v1.z + v2.z }; }
DVector3 DVector3Subtract(DVector3 v1, DVector3 v2) { return (DVector3){ v1.x - v2.x, v1.y - v2.y, v1.z - v2.z }; }
DVector3 DVector3Scale(DVector3 v, double scale) { return (DVector3){ v.x * scale, v.y * scale, v.z * scale }; }
double DVector3Length(DVector3 v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }
double DVector3Distance(DVector3 v1, DVector3 v2) { return DVector3Length(DVector3Subtract(v1, v2)); }
double DVector3DotProduct(DVector3 v1, DVector3 v2) { return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z; }
DVector3 DVector3CrossProduct(DVector3 v1, DVector3 v2) { return (DVector3){ v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x }; }
DVector3 DVector3Normalize(DVector3 v) { double len = DVector3Length(v); if (len == 0) return (DVector3){0}; return DVector3Scale(v, 1.0/len); }
Vector3 ToVector3(DVector3 v) { return (Vector3){ (float)v.x, (float)v.y, (float)v.z }; }

typedef struct {
    DVector3 position;
    DVector3 velocity;
    DVector3 acceleration; // For Verlet integration
    double mass;
    double radius;
    Color color;
    DVector3 trail[TRAIL_LENGTH];
    int trailIndex;
    bool active;
} Body;

typedef struct {
    DVector3 position;
    DVector3 velocity;
} State;

typedef struct {
    DVector3 dPosition; // velocity
    DVector3 dVelocity; // acceleration
} Derivative;

typedef enum {
    STATE_SIMULATION,
    STATE_MENU,
    STATE_SETTINGS
} AppState;

typedef enum {
    INTEGRATOR_RK4,
    INTEGRATOR_VERLET
} IntegratorType;

typedef struct {
    Camera3D renderCamera; // Used for rendering only
    DVector3 camPos;       // Real world position
    DVector3 camTarget;    // Real world target position
    double camDist;        // Real world distance from target
    Vector2 camAngle;
    int cameraTargetIndex;

    float timeScale; // Multiplier for real time (1.0 = 1 sec/sec, but we usually want faster)
    bool enableDrag;
    bool showLagrange;
    bool enableRoche;
    bool relativeView;
    IntegratorType integrator;

    bool creationMode;
    bool isDragging;
    Vector3 dragStartPos; // Visual drag start (relative)
    double newBodyMass;   // Real mass
    double spawnHeight;   // Real height
    float spawnAngle;

    AppState currentState;
    bool shouldExit;

    Shader lightShader;
    Model sphereModel;
    int lightPosLoc;
    int viewPosLoc;
    int objectColorLoc;
} SimulationState;

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
            derivs[i].dPosition = (DVector3){0,0,0};
            derivs[i].dVelocity = (DVector3){0,0,0};
            continue;
        }

        derivs[i].dPosition = states[i].velocity; // dr/dt = v
        DVector3 force = { 0.0, 0.0, 0.0 };

        // Drag Force (Gas/Dust) - Accretion Disk Model
        if (enableDrag && bodies[0].active) {
            // Assume Body 0 is the center (Sun)
            double distToCenter = DVector3Distance(states[i].position, states[0].position);

            if (distToCenter > DISK_INNER_RADIUS && distToCenter < DISK_OUTER_RADIUS) {
                // Density function: Higher density closer to center
                double normalizedDist = (distToCenter - DISK_INNER_RADIUS) / (DISK_OUTER_RADIUS - DISK_INNER_RADIUS);
                double density = 1.0 - normalizedDist;
                if (density < 0.0) density = 0.0;

                double dragCoeff = 0.002 * density;

                // F_drag = -c * rho * v
                DVector3 drag = DVector3Scale(states[i].velocity, -dragCoeff);
                force = DVector3Add(force, drag);
            }
        }

        for (int j = 0; j < MAX_BODIES; j++) {
            if (i == j || !bodies[j].active) continue;

            DVector3 direction = DVector3Subtract(states[j].position, states[i].position);
            double distance = DVector3Length(direction);

            double minDist = bodies[i].radius + bodies[j].radius;
            if (distance < minDist) distance = minDist;

            double forceMagnitude = (G * bodies[i].mass * bodies[j].mass) / (distance * distance);

            // Relativistic correction (Precession)
            DVector3 relVel = DVector3Subtract(states[j].velocity, states[i].velocity);
            // h = |r x v|
            DVector3 hVec = DVector3CrossProduct(direction, relVel);
            double h = DVector3Length(hVec);

            double correction = (3.0 * h * h) / (C_SPEED * C_SPEED * distance * distance);
            forceMagnitude *= (1.0 + correction);

            DVector3 forceVec = DVector3Scale(DVector3Normalize(direction), forceMagnitude);
            force = DVector3Add(force, forceVec);
        }
        derivs[i].dVelocity = DVector3Scale(force, 1.0 / bodies[i].mass); // dv/dt = a
    }
}

// Runge-Kutta 4th Order Integration
void integrateRK4(Body bodies[], double dt, bool enableDrag) {
    State initialStates[MAX_BODIES];
    getStates(bodies, initialStates);

    Derivative k1[MAX_BODIES], k2[MAX_BODIES], k3[MAX_BODIES], k4[MAX_BODIES];
    State tempStates[MAX_BODIES];

    // k1
    calculateDerivatives(initialStates, k1, bodies, enableDrag);

    // k2
    for(int i=0; i<MAX_BODIES; i++) {
        tempStates[i].position = DVector3Add(initialStates[i].position, DVector3Scale(k1[i].dPosition, dt * 0.5));
        tempStates[i].velocity = DVector3Add(initialStates[i].velocity, DVector3Scale(k1[i].dVelocity, dt * 0.5));
    }
    calculateDerivatives(tempStates, k2, bodies, enableDrag);

    // k3
    for(int i=0; i<MAX_BODIES; i++) {
        tempStates[i].position = DVector3Add(initialStates[i].position, DVector3Scale(k2[i].dPosition, dt * 0.5));
        tempStates[i].velocity = DVector3Add(initialStates[i].velocity, DVector3Scale(k2[i].dVelocity, dt * 0.5));
    }
    calculateDerivatives(tempStates, k3, bodies, enableDrag);

    // k4
    for(int i=0; i<MAX_BODIES; i++) {
        tempStates[i].position = DVector3Add(initialStates[i].position, DVector3Scale(k3[i].dPosition, dt));
        tempStates[i].velocity = DVector3Add(initialStates[i].velocity, DVector3Scale(k3[i].dVelocity, dt));
    }
    calculateDerivatives(tempStates, k4, bodies, enableDrag);

    // Update
    for(int i=0; i<MAX_BODIES; i++) {
        if (!bodies[i].active) continue;
        if (i == 0) continue; // Sun is stationary

        DVector3 dPos = DVector3Scale(DVector3Add(DVector3Add(k1[i].dPosition, DVector3Scale(k2[i].dPosition, 2.0)), DVector3Add(DVector3Scale(k3[i].dPosition, 2.0), k4[i].dPosition)), dt / 6.0);
        DVector3 dVel = DVector3Scale(DVector3Add(DVector3Add(k1[i].dVelocity, DVector3Scale(k2[i].dVelocity, 2.0)), DVector3Add(DVector3Scale(k3[i].dVelocity, 2.0), k4[i].dVelocity)), dt / 6.0);

        bodies[i].position = DVector3Add(bodies[i].position, dPos);
        bodies[i].velocity = DVector3Add(bodies[i].velocity, dVel);
    }
}

// Calculate accelerations for Verlet integration
void calculateAccelerations(Body bodies[], bool enableDrag) {
    for (int i = 0; i < MAX_BODIES; i++) {
        if (!bodies[i].active) {
            bodies[i].acceleration = (DVector3){0,0,0};
            continue;
        }

        DVector3 force = { 0.0, 0.0, 0.0 };

        // Drag Force (Gas/Dust) - Accretion Disk Model
        if (enableDrag && bodies[0].active) {
            double distToCenter = DVector3Distance(bodies[i].position, bodies[0].position);

            if (distToCenter > DISK_INNER_RADIUS && distToCenter < DISK_OUTER_RADIUS) {
                double normalizedDist = (distToCenter - DISK_INNER_RADIUS) / (DISK_OUTER_RADIUS - DISK_INNER_RADIUS);
                double density = 1.0 - normalizedDist;
                if (density < 0.0) density = 0.0;

                double dragCoeff = 0.002 * density;
                DVector3 drag = DVector3Scale(bodies[i].velocity, -dragCoeff);
                force = DVector3Add(force, drag);
            }
        }

        for (int j = 0; j < MAX_BODIES; j++) {
            if (i == j || !bodies[j].active) continue;

            DVector3 direction = DVector3Subtract(bodies[j].position, bodies[i].position);
            double distance = DVector3Length(direction);

            double minDist = bodies[i].radius + bodies[j].radius;
            if (distance < minDist) distance = minDist;

            double forceMagnitude = (G * bodies[i].mass * bodies[j].mass) / (distance * distance);

            // Relativistic correction (Precession)
            DVector3 relVel = DVector3Subtract(bodies[j].velocity, bodies[i].velocity);
            DVector3 hVec = DVector3CrossProduct(direction, relVel);
            double h = DVector3Length(hVec);

            double correction = (3.0 * h * h) / (C_SPEED * C_SPEED * distance * distance);
            forceMagnitude *= (1.0 + correction);

            DVector3 forceVec = DVector3Scale(DVector3Normalize(direction), forceMagnitude);
            force = DVector3Add(force, forceVec);
        }
        bodies[i].acceleration = DVector3Scale(force, 1.0 / bodies[i].mass);
    }
}

// Velocity Verlet Integration (Symplectic, faster than RK4)
void integrateVerlet(Body bodies[], double dt, bool enableDrag) {
    // 1. First half-kick (Update Velocity using current Acceleration)
    for (int i = 0; i < MAX_BODIES; i++) {
        if (!bodies[i].active) continue;
        bodies[i].velocity = DVector3Add(bodies[i].velocity, DVector3Scale(bodies[i].acceleration, dt * 0.5));
    }

    // 2. Drift (Update Position using new Velocity)
    for (int i = 0; i < MAX_BODIES; i++) {
        if (!bodies[i].active) continue;
        bodies[i].position = DVector3Add(bodies[i].position, DVector3Scale(bodies[i].velocity, dt));
    }

    // 3. Recalculate Forces (Accelerations)
    calculateAccelerations(bodies, enableDrag);

    // 4. Second half-kick (Update Velocity using NEW Acceleration)
    for (int i = 0; i < MAX_BODIES; i++) {
        if (!bodies[i].active) continue;
        bodies[i].velocity = DVector3Add(bodies[i].velocity, DVector3Scale(bodies[i].acceleration, dt * 0.5));
    }
}

void handleCollisions(Body bodies[]) {
    for (int i = 0; i < MAX_BODIES; i++) {
        if (!bodies[i].active) continue;
        for (int j = i + 1; j < MAX_BODIES; j++) {
            if (!bodies[j].active) continue;

            if (DVector3Distance(bodies[i].position, bodies[j].position) < (bodies[i].radius + bodies[j].radius)) {
                // Merge j into i (Fusion)
                Body *b1 = &bodies[i];
                Body *b2 = &bodies[j];

                DVector3 momentum1 = DVector3Scale(b1->velocity, b1->mass);
                DVector3 momentum2 = DVector3Scale(b2->velocity, b2->mass);
                DVector3 totalMomentum = DVector3Add(momentum1, momentum2);
                double totalMass = b1->mass + b2->mass;

                if (i == 0) {
                    // Sun absorbs body, but stays stationary
                    b1->mass = totalMass;
                    b1->radius = cbrt(pow(b1->radius, 3.0) + pow(b2->radius, 3.0));
                } else {
                    b1->velocity = DVector3Scale(totalMomentum, 1.0 / totalMass);

                    // Weighted position
                    b1->position = DVector3Scale(DVector3Add(DVector3Scale(b1->position, b1->mass), DVector3Scale(b2->position, b2->mass)), 1.0/totalMass);

                    b1->radius = cbrt(pow(b1->radius, 3.0) + pow(b2->radius, 3.0));
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
    double newMass = bodies[index].mass / fragments;
    double newRadius = bodies[index].radius / 2.0;
    if (newRadius < 10000.0) newRadius = 10000.0; // Min fragment size 10km

    for (int k = 0; k < fragments; k++) {
        for (int j = 0; j < MAX_BODIES; j++) {
            if (!bodies[j].active) {
                bodies[j].active = true;
                bodies[j].mass = newMass;
                bodies[j].radius = newRadius;
                bodies[j].color = (Color){
                    (unsigned char)fmin(255, bodies[index].color.r + GetRandomValue(-20, 20)),
                    (unsigned char)fmin(255, bodies[index].color.g + GetRandomValue(-20, 20)),
                    (unsigned char)fmin(255, bodies[index].color.b + GetRandomValue(-20, 20)),
                    255
                };

                DVector3 offset = { (double)GetRandomValue(-5, 5) * newRadius, (double)GetRandomValue(-5, 5) * newRadius, (double)GetRandomValue(-5, 5) * newRadius };
                bodies[j].position = DVector3Add(bodies[index].position, offset);

                DVector3 velSpread = { (double)GetRandomValue(-20, 20) * 100.0, (double)GetRandomValue(-20, 20) * 100.0, (double)GetRandomValue(-20, 20) * 100.0 };
                bodies[j].velocity = DVector3Add(bodies[index].velocity, velSpread);

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

        double dist = DVector3Distance(bodies[i].position, bodies[0].position);

        if (bodies[i].mass < 1e20) continue; // Ignore small asteroids

        double massRatio = bodies[0].mass / bodies[i].mass;
        double rocheLimit = 1.26 * bodies[i].radius * cbrt(massRatio);

        if (bodies[i].radius > 100000.0 && dist < rocheLimit) {
            explodeBody(bodies, i);
        }
    }
}

// Lagrange points are strictly defined for the restricted 3-body problem in a plane.
// We will project to the orbital plane of the largest planet for visualization.
void drawLagrangePoints(Body bodies[], DVector3 camPos) {
    if (!bodies[0].active) return;

    int heaviestIndex = -1;
    double maxMass = 0.0;

    for (int i = 1; i < MAX_BODIES; i++) {
        if (bodies[i].active && bodies[i].mass > maxMass) {
            maxMass = bodies[i].mass;
            heaviestIndex = i;
        }
    }

    if (heaviestIndex == -1) return;

    Body *m1 = &bodies[0];
    Body *m2 = &bodies[heaviestIndex];

    DVector3 r1 = m1->position;
    DVector3 r2 = m2->position;
    DVector3 R_vec = DVector3Subtract(r2, r1);
    double R = DVector3Length(R_vec);
    DVector3 u = DVector3Scale(R_vec, 1.0 / R);

    double massRatio = m2->mass / m1->mass;
    double hillRadius = R * cbrt(massRatio / 3.0);

    DVector3 l1 = DVector3Subtract(r2, DVector3Scale(u, hillRadius));
    DVector3 l2 = DVector3Add(r2, DVector3Scale(u, hillRadius));

    double l3_dist = R * (1.0 + (5.0/12.0) * massRatio);
    DVector3 l3 = DVector3Subtract(r1, DVector3Scale(u, l3_dist));

    // For L4/L5 we need a vector perpendicular to R_vec in the orbital plane.
    // We can use the velocity of m2 relative to m1 to find the orbital plane normal.
    DVector3 vRel = DVector3Subtract(m2->velocity, m1->velocity);
    DVector3 orbitalNormal = DVector3Normalize(DVector3CrossProduct(R_vec, vRel));
    // If velocity is parallel to position (falling straight in), this fails, but that's rare for planets.
    // Perpendicular vector in plane:
    DVector3 perp = DVector3Normalize(DVector3CrossProduct(orbitalNormal, u));

    // L4/L5 form equilateral triangles.
    // Position is R/2 along u, and sqrt(3)/2 * R along perp.
    double h_tri = 0.8660254 * R; // sin(60) * R

    DVector3 l4 = DVector3Add(DVector3Add(r1, DVector3Scale(u, R * 0.5)), DVector3Scale(perp, h_tri));
    DVector3 l5 = DVector3Add(DVector3Add(r1, DVector3Scale(u, R * 0.5)), DVector3Scale(perp, -h_tri));

    Color lColor = VIOLET;
    float lRadius = 5.0f; // This should probably be scaled too, but let's keep it fixed for now or use visual size logic

    DrawSphere(ToVector3(DVector3Subtract(l1, camPos)), lRadius, lColor);
    DrawSphere(ToVector3(DVector3Subtract(l2, camPos)), lRadius, lColor);
    DrawSphere(ToVector3(DVector3Subtract(l3, camPos)), lRadius, lColor);
    DrawSphere(ToVector3(DVector3Subtract(l4, camPos)), lRadius, lColor);
    DrawSphere(ToVector3(DVector3Subtract(l5, camPos)), lRadius, lColor);
}

// Helper to create a body in a stable circular orbit around a parent
Body createOrbitingBody(Body parent, double orbitRadius, double angleDeg, double mass, double radius, Color color) {
    double theta = angleDeg * DEG2RAD;

    // Position offset
    double dx = orbitRadius * cos(theta);
    double dz = orbitRadius * sin(theta);

    DVector3 position = DVector3Add(parent.position, (DVector3){ dx, 0.0, dz });

    // Orbital velocity (circular) v = sqrt(GM / r)
    double vMag = sqrt((G * parent.mass) / orbitRadius);

    // Velocity vector (tangent to circle)
    // If pos is (cos, sin), vel is (-sin, cos) for counter-clockwise orbit
    double vx = -vMag * sin(theta);
    double vz = vMag * cos(theta);

    DVector3 velocity = DVector3Add(parent.velocity, (DVector3){ vx, 0.0, vz });

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
        .position = { 0.0, 0.0, 0.0 },
        .velocity = { 0.0, 0.0, 0.0 },
        .mass = M_SUN,
        .radius = 6.9634e8, // Sun Radius
        .color = YELLOW,
        .active = true
    };

    // Planet 1 (Mercury-like)
    bodies[1] = createOrbitingBody(bodies[0], 0.39 * AU, 0.0, 3.301e23, 2.4397e6, BLUE);

    // Planet 2 (Earth-like)
    bodies[2] = createOrbitingBody(bodies[0], 1.0 * AU, 0.0, M_EARTH, 6.371e6, RED);

    // Planet 3 (Jupiter-like)
    bodies[3] = createOrbitingBody(bodies[0], 5.2 * AU, 0.0, 1.898e27, 6.9911e7, ORANGE);

    // Moon of Planet 3 (Europa-like)
    bodies[4] = createOrbitingBody(bodies[3], 6.71e8, 0.0, 4.8e22, 1.56e6, WHITE);

    for (int i = 0; i < MAX_BODIES; i++) {
        for (int j = 0; j < TRAIL_LENGTH; j++) {
            bodies[i].trail[j] = bodies[i].position;
        }
        bodies[i].trailIndex = 0;
    }
}

void drawAccretionDisk(DVector3 center, DVector3 camPos) {
    float inner = DISK_INNER_RADIUS;
    float outer = DISK_OUTER_RADIUS;
    int slices = 60;
    Color c = Fade(BLUE, 0.15f);

    DVector3 relCenter = DVector3Subtract(center, camPos);

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
        rlVertex3f(relCenter.x + inner*cos1, relCenter.y, relCenter.z + inner*sin1);
        rlVertex3f(relCenter.x + outer*cos1, relCenter.y, relCenter.z + outer*sin1);
        rlVertex3f(relCenter.x + inner*cos2, relCenter.y, relCenter.z + inner*sin2);

        // V2(out,1) -> V4(out,2) -> V3(in,2)
        rlVertex3f(relCenter.x + outer*cos1, relCenter.y, relCenter.z + outer*sin1);
        rlVertex3f(relCenter.x + outer*cos2, relCenter.y, relCenter.z + outer*sin2);
        rlVertex3f(relCenter.x + inner*cos2, relCenter.y, relCenter.z + inner*sin2);
    }
    rlEnd();
}

void drawOrbit(Body bodies[], int targetIndex, int centerIndex, DVector3 camPos) {
    if (!bodies[centerIndex].active || !bodies[targetIndex].active || targetIndex == centerIndex) return;

    Body *center = &bodies[centerIndex];
    Body *planet = &bodies[targetIndex];

    DVector3 rVec = DVector3Subtract(planet->position, center->position);
    DVector3 vVec = DVector3Subtract(planet->velocity, center->velocity);

    double r = DVector3Length(rVec);
    double v = DVector3Length(vVec);
    double mu = G * (center->mass + planet->mass);

    // Specific angular momentum h = r x v
    DVector3 hVec = DVector3CrossProduct(rVec, vVec);
    double h = DVector3Length(hVec);
    if (h < 0.1) return;

    // Eccentricity vector e = (v x h) / mu - r / |r|
    DVector3 vxh = DVector3CrossProduct(vVec, hVec);
    DVector3 eVec = DVector3Subtract(DVector3Scale(vxh, 1.0/mu), DVector3Scale(rVec, 1.0/r));
    double e = DVector3Length(eVec);

    // Semi-major axis a = 1 / (2/r - v^2/mu)
    double energy = v*v/2.0 - mu/r;
    if (fabs(energy) < 0.0001) return;
    double a = -mu / (2.0 * energy);

    if (e >= 1.0 || a < 0) return; // Hyperbolic/Parabolic

    // Basis vectors for orbital plane
    DVector3 n = DVector3Normalize(hVec); // Normal to plane
    DVector3 p; // Periapsis direction
    if (e > 0.001) {
        p = DVector3Normalize(eVec);
    } else {
        if (fabs(n.y) < 0.9) p = DVector3Normalize(DVector3CrossProduct((DVector3){0,1,0}, n));
        else p = DVector3Normalize(DVector3CrossProduct((DVector3){1,0,0}, n));
    }
    DVector3 q = DVector3CrossProduct(n, p);

    // Draw Orbit Path
    rlBegin(RL_LINES);
    rlColor4ub(255, 255, 255, 60);

    int segments = 360;
    DVector3 prevPos = {0};
    bool first = true;

    for (int i = 0; i <= segments; i++) {
        double theta = (double)i / segments * 2.0 * PI;
        double radius = a * (1.0 - e*e) / (1.0 + e * cos(theta));

        DVector3 posInPlane = DVector3Add(DVector3Scale(p, radius * cos(theta)), DVector3Scale(q, radius * sin(theta)));
        DVector3 worldPos = DVector3Add(center->position, posInPlane);
        DVector3 relPos = DVector3Subtract(worldPos, camPos);

        if (!first) {
            rlVertex3f(prevPos.x, prevPos.y, prevPos.z);
            rlVertex3f(relPos.x, relPos.y, relPos.z);
        }
        prevPos = relPos;
        first = false;
    }
    rlEnd();

    // Draw Periapsis (Green) and Apoapsis (Red)
    double r_peri = a * (1.0 - e);
    double r_apo = a * (1.0 + e);

    DVector3 posPeri = DVector3Add(center->position, DVector3Scale(p, r_peri));
    DVector3 posApo = DVector3Add(center->position, DVector3Scale(p, -r_apo));

    DrawSphere(ToVector3(DVector3Subtract(posPeri, camPos)), 3.0f, GREEN);
    DrawSphere(ToVector3(DVector3Subtract(posApo, camPos)), 3.0f, RED);
}

int findParentBody(Body bodies[], int subjectIndex) {
    if (subjectIndex == 0) return -1;

    int bestParent = 0;
    double maxForce = -1.0;

    for (int i = 0; i < MAX_BODIES; i++) {
        if (i == subjectIndex || !bodies[i].active) continue;

        double dist = DVector3Distance(bodies[i].position, bodies[subjectIndex].position);
        if (dist < 1.0) continue;

        double force = bodies[i].mass / (dist * dist);

        if (force > maxForce) {
            maxForce = force;
            bestParent = i;
        }
    }
    return bestParent;
}

void drawOrbitEditor(Body bodies[], int targetIndex, int screenWidth, int screenHeight, SimulationState *state) {
    if (targetIndex == 0 || !bodies[targetIndex].active) return;

    int parentIndex = findParentBody(bodies, targetIndex);
    if (parentIndex == -1) return;

    Body *b = &bodies[targetIndex];
    Body *p = &bodies[parentIndex];

    // Calculate relative state
    DVector3 relPos = DVector3Subtract(b->position, p->position);
    DVector3 relVel = DVector3Subtract(b->velocity, p->velocity);

    double dist = DVector3Length(relPos);
    double speed = DVector3Length(relVel);
    double mu = G * (p->mass + b->mass);

    // Calculate orbital elements
    DVector3 hVec = DVector3CrossProduct(relPos, relVel);
    double h = DVector3Length(hVec);
    DVector3 n = DVector3Scale(hVec, 1.0/h);

    DVector3 vxh = DVector3CrossProduct(relVel, hVec);
    DVector3 eVec = DVector3Subtract(DVector3Scale(vxh, 1.0/mu), DVector3Scale(relPos, 1.0/dist));
    double e = DVector3Length(eVec);

    // Handle circular orbits (e ~ 0)
    DVector3 eDir;
    if (e < 1e-4) {
        e = 0.0;
        eDir = DVector3Normalize(relPos);
    } else {
        eDir = DVector3Normalize(eVec);
    }

    DVector3 qDir = DVector3CrossProduct(n, eDir);

    double energy = speed*speed/2.0 - mu/dist;
    double a = -mu / (2.0 * energy);

    // True Anomaly nu
    double nu = atan2(DVector3DotProduct(relPos, qDir), DVector3DotProduct(relPos, eDir));

    // Calculate Angle (Argument of Periapsis relative to reference)
    DVector3 worldUp = {0, 1, 0};
    if (fabs(n.y) > 0.95) worldUp = (DVector3){1, 0, 0};
    DVector3 u = DVector3Normalize(DVector3CrossProduct(worldUp, n));
    DVector3 v = DVector3CrossProduct(n, u);
    double angle = atan2(DVector3DotProduct(eDir, v), DVector3DotProduct(eDir, u));
    if (angle < 0) angle += 2*PI;

    // Calculate Inclination
    DVector3 refNormal = {0, 1, 0};
    double inclination = acos(DVector3DotProduct(n, refNormal));
    DVector3 nodeVec = DVector3CrossProduct(refNormal, n);
    if (DVector3Length(nodeVec) < 0.001) nodeVec = (DVector3){1, 0, 0};
    nodeVec = DVector3Normalize(nodeVec);

    // Derived parameters
    double r_peri = a * (1.0 - e);
    double r_apo = a * (1.0 + e);
    double period = 2.0 * PI * sqrt(pow(a, 3.0) / mu);

    // UI Layout
    int uiWidth = 280;
    int uiHeight = 420;
    int uiX = screenWidth - uiWidth - 10;
    int uiY = screenHeight - uiHeight - 10;
    Rectangle uiRect = { uiX, uiY, uiWidth, uiHeight };

    if (GuiWindowBox(uiRect, "ORBIT EDITOR")) {
        state->cameraTargetIndex = 0;
    }

    int startX = uiX + 10;
    int startY = uiY + 30;
    int controlWidth = 140;
    int labelWidth = 110;

    // Static state for UI
    static int unitSelection = 2; // 0: m, 1: km, 2: AU
    static bool unitEditMode = false;
    
    static bool editEcc = false;
    static bool editSemi = false;
    static bool editPeri = false;
    static bool editApo = false;
    static bool editRot = false;
    static bool editInc = false;

    static char textEcc[64] = "";
    static char textSemi[64] = "";
    static char textPeri[64] = "";
    static char textApo[64] = "";
    static char textRot[64] = "";
    static char textInc[64] = "";

    static int lastTarget = -1;
    if (lastTarget != targetIndex) {
        lastTarget = targetIndex;
        editEcc = false; editSemi = false; editPeri = false; 
        editApo = false; editRot = false; editInc = false;
        unitEditMode = false;
    }

    double unitScale = 1.0;
    const char* unitLabel = "m";
    if (unitSelection == 1) { unitScale = 1000.0; unitLabel = "km"; }
    else if (unitSelection == 2) { unitScale = AU; unitLabel = "AU"; }

    // Info
    GuiLabel((Rectangle){startX, startY, 260, 20}, TextFormat("Period: %.1f s", period));
    startY += 30;

    // Unit Selector
    Rectangle dropdownRect = {startX + 60, startY, 120, 20};
    GuiLabel((Rectangle){startX, startY, 50, 20}, "Units:");
    startY += 40;

    // Helper to apply changes
    double new_a = a;
    double new_e = e;
    double new_angle = angle;
    double new_inclination = inclination;
    bool changed = false;

    // 1. Eccentricity
    GuiLabel((Rectangle){startX, startY, labelWidth, 20}, "Eccentricity:");
    if (!editEcc) snprintf(textEcc, 64, "%.5f", e);
    if (GuiTextBox((Rectangle){startX + labelWidth, startY, controlWidth, 20}, textEcc, 64, editEcc)) {
        editEcc = !editEcc;
        if (!editEcc) {
            new_e = atof(textEcc);
            if (new_e < 0) new_e = 0;
            if (new_e >= 0.99) new_e = 0.99;
            changed = true;
        }
    }
    startY += 30;

    // 2. Semi-major Axis
    GuiLabel((Rectangle){startX, startY, labelWidth, 20}, TextFormat("Semi-major (%s):", unitLabel));
    if (!editSemi) snprintf(textSemi, 64, "%.4f", a / unitScale);
    if (GuiTextBox((Rectangle){startX + labelWidth, startY, controlWidth, 20}, textSemi, 64, editSemi)) {
        editSemi = !editSemi;
        if (!editSemi) {
            new_a = atof(textSemi) * unitScale;
            if (new_a < 1000.0) new_a = 1000.0;
            changed = true;
        }
    }
    startY += 30;

    // 3. Periapsis
    GuiLabel((Rectangle){startX, startY, labelWidth, 20}, TextFormat("Periapsis (%s):", unitLabel));
    if (!editPeri) snprintf(textPeri, 64, "%.4f", r_peri / unitScale);
    if (GuiTextBox((Rectangle){startX + labelWidth, startY, controlWidth, 20}, textPeri, 64, editPeri)) {
        editPeri = !editPeri;
        if (!editPeri) {
            double new_p = atof(textPeri) * unitScale;
            if (new_p >= r_apo) new_p = r_apo - 1000.0;
            new_a = (new_p + r_apo) / 2.0;
            new_e = (r_apo - new_p) / (r_apo + new_p);
            changed = true;
        }
    }
    startY += 30;

    // 4. Apoapsis
    GuiLabel((Rectangle){startX, startY, labelWidth, 20}, TextFormat("Apoapsis (%s):", unitLabel));
    if (!editApo) snprintf(textApo, 64, "%.4f", r_apo / unitScale);
    if (GuiTextBox((Rectangle){startX + labelWidth, startY, controlWidth, 20}, textApo, 64, editApo)) {
        editApo = !editApo;
        if (!editApo) {
            double new_ap = atof(textApo) * unitScale;
            if (new_ap <= r_peri) new_ap = r_peri + 1000.0;
            new_a = (r_peri + new_ap) / 2.0;
            new_e = (new_ap - r_peri) / (new_ap + r_peri);
            changed = true;
        }
    }
    startY += 30;

    // 5. Rotation
    GuiLabel((Rectangle){startX, startY, labelWidth, 20}, "Rotation (deg):");
    if (!editRot) snprintf(textRot, 64, "%.2f", angle * RAD2DEG);
    if (GuiTextBox((Rectangle){startX + labelWidth, startY, controlWidth, 20}, textRot, 64, editRot)) {
        editRot = !editRot;
        if (!editRot) {
            new_angle = atof(textRot) * DEG2RAD;
            changed = true;
        }
    }
    startY += 30;

    // 6. Inclination
    GuiLabel((Rectangle){startX, startY, labelWidth, 20}, "Inclination (deg):");
    if (!editInc) snprintf(textInc, 64, "%.2f", inclination * RAD2DEG);
    if (GuiTextBox((Rectangle){startX + labelWidth, startY, controlWidth, 20}, textInc, 64, editInc)) {
        editInc = !editInc;
        if (!editInc) {
            new_inclination = atof(textInc) * DEG2RAD;
            changed = true;
        }
    }
    startY += 30;

    // Draw Dropdown last to be on top
    if (GuiDropdownBox(dropdownRect, "Meters;Kilometers;AU", &unitSelection, unitEditMode)) {
        unitEditMode = !unitEditMode;
    }

    if (changed) {
        // Reconstruct vectors
        DVector3 new_eDir = DVector3Add(DVector3Scale(u, cos(new_angle)), DVector3Scale(v, sin(new_angle)));
        DVector3 new_qDir = DVector3CrossProduct(n, new_eDir);

        double slr = new_a * (1.0 - new_e * new_e);
        double new_r = slr / (1.0 + new_e * cos(nu));

        double v_radial = sqrt(mu/slr) * new_e * sin(nu);
        double v_tangential = sqrt(mu/slr) * (1.0 + new_e * cos(nu));

        DVector3 r_hat = DVector3Add(DVector3Scale(new_eDir, cos(nu)), DVector3Scale(new_qDir, sin(nu)));
        DVector3 t_hat = DVector3CrossProduct(n, r_hat);

        DVector3 pos = DVector3Scale(r_hat, new_r);
        DVector3 vel = DVector3Add(DVector3Scale(r_hat, v_radial), DVector3Scale(t_hat, v_tangential));

        // Apply Inclination Change
        if (fabs(new_inclination - inclination) > 0.001) {
            double deltaInc = new_inclination - inclination;
            Quaternion q = QuaternionFromAxisAngle(ToVector3(nodeVec), (float)deltaInc);
            Vector3 posF = Vector3RotateByQuaternion(ToVector3(pos), q);
            Vector3 velF = Vector3RotateByQuaternion(ToVector3(vel), q);
            pos = (DVector3){posF.x, posF.y, posF.z};
            vel = (DVector3){velF.x, velF.y, velF.z};
        }

        b->position = DVector3Add(p->position, pos);
        b->velocity = DVector3Add(p->velocity, vel);
    }
}

// Initialize simulation state, camera, and shaders
void InitSimulation(SimulationState *state) {
    state->camDist = 2.0 * AU;
    state->camPos = (DVector3){0, 2.0 * AU, 2.0 * AU};
    state->camAngle = (Vector2){ 0.0f, 1.0f };
    state->cameraTargetIndex = 0;

    state->timeScale = 10000.0f; // Start with faster time
    state->enableDrag = false;
    state->showLagrange = false;
    state->enableRoche = true;
    state->relativeView = false;
    state->integrator = INTEGRATOR_RK4;

    state->creationMode = false;
    state->isDragging = false;
    state->dragStartPos = (Vector3){0};
    state->newBodyMass = M_EARTH;
    state->spawnHeight = 0.0;
    state->spawnAngle = 0.0f;

    state->currentState = STATE_SIMULATION;
    state->shouldExit = false;

    // 3D Camera (Raylib camera used for rendering relative to camPos)
    state->renderCamera = (Camera3D){ 0 };
    state->renderCamera.position = (Vector3){ 0.0f, 0.0f, 0.0f };
    state->renderCamera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    state->renderCamera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    state->renderCamera.fovy = 45.0f;
    state->renderCamera.projection = CAMERA_PERSPECTIVE;

    // Lighting Shader
    state->lightShader = LoadShaderFromMemory(
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

    state->lightPosLoc = GetShaderLocation(state->lightShader, "lightPos");
    state->viewPosLoc = GetShaderLocation(state->lightShader, "viewPos");
    int ambientLoc = GetShaderLocation(state->lightShader, "ambientColor");
    int lightColorLoc = GetShaderLocation(state->lightShader, "lightColor");
    state->objectColorLoc = GetShaderLocation(state->lightShader, "objectColor");

    float ambient[4] = { 0.1f, 0.1f, 0.1f, 1.0f };
    SetShaderValue(state->lightShader, ambientLoc, ambient, SHADER_UNIFORM_VEC4);
    float lightColor[4] = { 1.0f, 1.0f, 0.9f, 1.0f };
    SetShaderValue(state->lightShader, lightColorLoc, lightColor, SHADER_UNIFORM_VEC4);

    Mesh sphereMesh = GenMeshSphere(1.0f, 16, 16);
    state->sphereModel = LoadModelFromMesh(sphereMesh);
    state->sphereModel.materials[0].shader = state->lightShader;
}

// Handle user input for camera, creation mode, and time scale
void HandleInput(SimulationState *state, Body bodies[], int screenWidth, int screenHeight) {
    (void)screenHeight; // Unused parameter
    if (IsKeyPressed(KEY_ESCAPE)) {
        if (state->currentState == STATE_SIMULATION) state->currentState = STATE_MENU;
        else if (state->currentState == STATE_MENU) state->currentState = STATE_SIMULATION;
        else if (state->currentState == STATE_SETTINGS) state->currentState = STATE_MENU;
    }

    if (state->currentState == STATE_SIMULATION) {
        if (IsKeyPressed(KEY_N)) {
            state->creationMode = !state->creationMode;
            state->spawnHeight = 0.0;
        }

        if (!state->creationMode) {
            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                Vector2 delta = GetMouseDelta();
                state->camAngle.x -= delta.x * 0.005f;
                state->camAngle.y -= delta.y * 0.005f;
                if (state->camAngle.y < 0.01f) state->camAngle.y = 0.01f;
                if (state->camAngle.y > PI - 0.01f) state->camAngle.y = PI - 0.01f;
            }
            
            float wheel = GetMouseWheelMove();
            if (wheel != 0) {
                state->camDist *= (1.0f - wheel * 0.1f);
                if (state->camDist < 1000.0) state->camDist = 1000.0; // Min distance 1km
            }
        }

        if (state->creationMode) {
            int uiWidth = 240;
            int uiHeight = 200;
            int uiX = screenWidth - uiWidth - 10;
            int uiY = 10;
            Rectangle uiRect = { uiX, uiY, uiWidth, uiHeight };
            bool mouseOverUI = CheckCollisionPointRec(GetMousePosition(), uiRect);

            if (!mouseOverUI) {
                // Ray from camera (which is at 0,0,0 in render space)
                Ray ray = GetMouseRay(GetMousePosition(), state->renderCamera);
                
                // We need to intersect with plane y = spawnHeight
                // Ray origin in world space is state->camPos
                // Ray dir is ray.direction
                
                // P = O + D*t
                // P.y = spawnHeight
                // O.y + D.y*t = spawnHeight
                // t = (spawnHeight - O.y) / D.y
                
                if (fabs(ray.direction.y) > 0.001f) {
                    double t = (state->spawnHeight - state->camPos.y) / ray.direction.y;
                    if (t >= 0) {
                        DVector3 mouseWorldPos = {
                            state->camPos.x + ray.direction.x * t,
                            state->camPos.y + ray.direction.y * t,
                            state->camPos.z + ray.direction.z * t
                        };

                        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
                            state->isDragging = true;
                            // Store relative drag start for visualization
                            state->dragStartPos = (Vector3){ (float)(mouseWorldPos.x - state->camPos.x), (float)(mouseWorldPos.y - state->camPos.y), (float)(mouseWorldPos.z - state->camPos.z) };
                        }

                        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT) && state->isDragging) {
                            state->isDragging = false;
                            for (int i = 0; i < MAX_BODIES; i++) {
                                if (!bodies[i].active) {
                                    bodies[i].active = true;
                                    bodies[i].position = mouseWorldPos;

                                    // Calculate velocity based on drag distance
                                    // We need to reconstruct drag start in world space
                                    DVector3 dragStartWorld = {
                                        state->camPos.x + state->dragStartPos.x,
                                        state->camPos.y + state->dragStartPos.y,
                                        state->camPos.z + state->dragStartPos.z
                                    };
                                    
                                    DVector3 dragVec = DVector3Subtract(mouseWorldPos, dragStartWorld);
                                    double dist = DVector3Length(dragVec);
                                    // Scale velocity: 1 screen unit -> huge velocity?
                                    // Let's say drag distance is proportional to orbital velocity at this distance
                                    double speed = dist / state->camDist * 30000.0; // Arbitrary scaling

                                    DVector3 dir = DVector3Normalize(dragVec);
                                    
                                    // Simple launch in direction
                                    bodies[i].velocity = DVector3Scale(dir, speed);

                                    bodies[i].mass = state->newBodyMass;
                                    bodies[i].radius = cbrt(state->newBodyMass / M_EARTH) * 6.371e6;
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

        if (IsKeyPressed(KEY_RIGHT)) state->timeScale *= 2.0f;
        if (IsKeyPressed(KEY_LEFT)) state->timeScale *= 0.5f;
        if (IsKeyPressed(KEY_SPACE)) state->timeScale = (state->timeScale == 0.0f) ? 10000.0f : 0.0f;
    }
}

// Update camera position and run physics simulation (RK4, collisions)
void UpdatePhysics(SimulationState *state, Body bodies[]) {
    if (state->currentState != STATE_SIMULATION) return;

    float dt = GetFrameTime() * state->timeScale;

    // Dynamic sub-stepping for stability
    // Target subDt around 50.0 seconds for good stability with moons
    double targetSubDt = 50.0;
    int subSteps = (int)ceil(dt / targetSubDt);
    if (subSteps < 10) subSteps = 10; // Minimum steps
    if (subSteps > 1000) subSteps = 1000; // Cap to prevent freeze
    double subDt = (double)dt / subSteps;

    for (int step = 0; step < subSteps; step++) {
        if (state->integrator == INTEGRATOR_RK4) {
            integrateRK4(bodies, subDt, state->enableDrag);
        } else {
            integrateVerlet(bodies, subDt, state->enableDrag);
        }
        
        handleCollisions(bodies);
        if (state->enableRoche) checkRocheLimit(bodies);
    }

    if (state->timeScale != 0.0f) {
        for (int i = 0; i < MAX_BODIES; i++) {
            if (!bodies[i].active) continue;
            bodies[i].trail[bodies[i].trailIndex] = bodies[i].position;
            bodies[i].trailIndex = (bodies[i].trailIndex + 1) % TRAIL_LENGTH;
        }
    }

    // Update Camera Position AFTER integration to prevent jitter/lag
    if (!bodies[state->cameraTargetIndex].active) state->cameraTargetIndex = 0;
    DVector3 targetPos = bodies[state->cameraTargetIndex].position;
    state->camTarget = targetPos;

    state->camPos.x = targetPos.x + state->camDist * sin(state->camAngle.y) * cos(state->camAngle.x);
    state->camPos.y = targetPos.y + state->camDist * cos(state->camAngle.y);
    state->camPos.z = targetPos.z + state->camDist * sin(state->camAngle.y) * sin(state->camAngle.x);

    // Update Shader Uniforms
    // We render relative to camera, so view pos is always 0,0,0
    Vector3 viewPos = {0,0,0};
    SetShaderValue(state->lightShader, state->viewPosLoc, &viewPos, SHADER_UNIFORM_VEC3);
    
    if (bodies[0].active) {
        DVector3 sunRel = DVector3Subtract(bodies[0].position, state->camPos);
        Vector3 sunPos = ToVector3(sunRel);
        SetShaderValue(state->lightShader, state->lightPosLoc, &sunPos, SHADER_UNIFORM_VEC3);
    }
}

// Render the 3D world (bodies, orbits, grid)
void Draw3DScene(SimulationState *state, Body bodies[]) {
    // Setup render camera at 0,0,0 looking at relative target
    DVector3 targetRel = DVector3Subtract(state->camTarget, state->camPos);
    state->renderCamera.target = ToVector3(targetRel);
    state->renderCamera.position = (Vector3){0,0,0};
    
    BeginMode3D(state->renderCamera);

    // Custom Grid centered on target
    rlSetClipPlanes(1.0f, 1e15f); // Increase far plane

    // Adaptive Grid
    double gridSpacing = pow(10.0, floor(log10(state->camDist / 2.0)));
    int slices = 100;
    double spacing = gridSpacing;
    double halfSize = slices * spacing / 2.0;
    Color gridColor = Fade(LIGHTGRAY, 0.2f);
    
    // Grid is on XZ plane of the target
    DVector3 gridCenter = state->camTarget;
    DVector3 gridRel = DVector3Subtract(gridCenter, state->camPos);

    rlBegin(RL_LINES);
    rlColor4ub(gridColor.r, gridColor.g, gridColor.b, gridColor.a);
    for (int i = -slices/2; i <= slices/2; i++) {
        // Lines parallel to Z
        DVector3 p1 = { gridRel.x + i * spacing, gridRel.y, gridRel.z - halfSize };
        DVector3 p2 = { gridRel.x + i * spacing, gridRel.y, gridRel.z + halfSize };
        rlVertex3f((float)p1.x, (float)p1.y, (float)p1.z);
        rlVertex3f((float)p2.x, (float)p2.y, (float)p2.z);

        // Lines parallel to X
        DVector3 p3 = { gridRel.x - halfSize, gridRel.y, gridRel.z + i * spacing };
        DVector3 p4 = { gridRel.x + halfSize, gridRel.y, gridRel.z + i * spacing };
        rlVertex3f((float)p3.x, (float)p3.y, (float)p3.z);
        rlVertex3f((float)p4.x, (float)p4.y, (float)p4.z);
    }
    rlEnd();

    if (state->enableDrag && bodies[0].active) {
        drawAccretionDisk(bodies[0].position, state->camPos);
    }

    if (state->cameraTargetIndex != 0 && bodies[state->cameraTargetIndex].active) {
        if (state->relativeView) {
            int parent = findParentBody(bodies, state->cameraTargetIndex);
            if (parent != -1) drawOrbit(bodies, state->cameraTargetIndex, parent, state->camPos);
        } else {
            drawOrbit(bodies, state->cameraTargetIndex, 0, state->camPos);
        }
    }

    for (int i = 0; i < MAX_BODIES; i++) {
        if (!bodies[i].active) continue;
        for (int j = 0; j < TRAIL_LENGTH - 1; j++) {
            int idx = (bodies[i].trailIndex + j) % TRAIL_LENGTH;
            int nextIdx = (idx + 1) % TRAIL_LENGTH;
            
            DVector3 p1d = DVector3Subtract(bodies[i].trail[idx], state->camPos);
            DVector3 p2d = DVector3Subtract(bodies[i].trail[nextIdx], state->camPos);
            
            Vector3 p1 = ToVector3(p1d);
            Vector3 p2 = ToVector3(p2d);
            
            if (Vector3DistanceSqr(p1, p2) > 0.001f)
                DrawLine3D(p1, p2, Fade(bodies[i].color, 0.4f));
        }
    }

    for (int i = 0; i < MAX_BODIES; i++) {
        if (!bodies[i].active) continue;

        double distToCam = DVector3Distance(bodies[i].position, state->camPos);
        
        // Visual Size Logic
        double minVisualSize = distToCam * 0.002; // tan(angle) approx
        double drawRadius = bodies[i].radius;
        if (drawRadius < minVisualSize) drawRadius = minVisualSize;

        DVector3 relPos = DVector3Subtract(bodies[i].position, state->camPos);
        Vector3 pos = ToVector3(relPos);

        // LOD: If object is very small on screen, draw a point
        // Approximate screen size in pixels (assuming 45 deg FOV and 800px height)
        // size_px = (diameter / dist) * (screen_height / (2 * tan(fov/2)))
        // tan(22.5) ~ 0.414. 800 / 0.828 ~ 966.
        double screenDiameter = (bodies[i].radius * 2.0 / distToCam) * 1000.0;

        if (i != 0 && screenDiameter < 3.0) {
            DrawPoint3D(pos, bodies[i].color);
        } else {
            if (i == 0) {
                DrawSphere(pos, (float)drawRadius, bodies[i].color);
            } else {
                Vector3 scale = { (float)drawRadius, (float)drawRadius, (float)drawRadius };
                float col[4] = { bodies[i].color.r/255.0f, bodies[i].color.g/255.0f, bodies[i].color.b/255.0f, bodies[i].color.a/255.0f };
                SetShaderValue(state->lightShader, state->objectColorLoc, col, SHADER_UNIFORM_VEC4);
                DrawModelEx(state->sphereModel, pos, (Vector3){0,1,0}, 0.0f, scale, WHITE);
            }
        }
    }

    if (state->creationMode && state->currentState == STATE_SIMULATION) {
        int uiWidth = 280;
        int uiHeight = 200;
        int uiX = GetScreenWidth() - uiWidth - 10;
        int uiY = 10;
        Rectangle uiRect = { uiX, uiY, uiWidth, uiHeight };

        Ray ray = GetMouseRay(GetMousePosition(), state->renderCamera);
        if (fabs(ray.direction.y) > 0.001f) {
            double t = (state->spawnHeight - state->camPos.y) / ray.direction.y;
            if (t >= 0) {
                DVector3 mouseWorldPos = {
                    state->camPos.x + ray.direction.x * t,
                    state->camPos.y + ray.direction.y * t,
                    state->camPos.z + ray.direction.z * t
                };
                
                DVector3 mouseRelPos = DVector3Subtract(mouseWorldPos, state->camPos);
                Vector3 mouseRel = ToVector3(mouseRelPos);

                if (!CheckCollisionPointRec(GetMousePosition(), uiRect)) {
                    if (state->isDragging) {
                        // dragStartPos is relative to camPos at start of drag
                        // But camPos might have moved? No, we assume camera is static during drag for simplicity or we should have stored world pos.
                        // In HandleInput we stored dragStartPos as relative to camPos.
                        // Let's assume we want to draw from dragStartPos (relative) to current mouse (relative)
                        
                        // Wait, dragStartPos in HandleInput was stored as relative to camPos.
                        // If camPos changes, this is wrong. But camPos only changes if we move camera.
                        // If we are dragging, we probably aren't moving camera (unless we allow both).
                        
                        Vector3 endPos = mouseRel;
                        DrawLine3D(state->dragStartPos, endPos, RED);
                        DrawSphere(state->dragStartPos, 5.0f, GREEN); // This 5.0f is tiny in real scale...
                        // We need visual size for this too?
                        // Or just use screen space drawing?
                        // Let's just draw a line.
                        
                    } else {
                        DrawSphereWires(mouseRel, 10.0f, 8, 8, Fade(GRAY, 0.5f)); // 10.0f is tiny
                    }
                }
            }
        }
    }

    if (state->showLagrange) {
        drawLagrangePoints(bodies, state->camPos);
    }

    EndMode3D();
}

// Render the 2D user interface (orbit editor, creation menu, settings)
void DrawUI(SimulationState *state, Body bodies[], int screenWidth, int screenHeight) {
    if (state->cameraTargetIndex != 0 && state->currentState == STATE_SIMULATION) {
        drawOrbitEditor(bodies, state->cameraTargetIndex, screenWidth, screenHeight, state);
    }

    if (state->creationMode && state->currentState == STATE_SIMULATION) {
        int uiWidth = 280;
        int uiHeight = 210;
        int uiX = screenWidth - uiWidth - 10;
        int uiY = 10;
        Rectangle uiRect = { uiX, uiY, uiWidth, uiHeight };

        if (GuiWindowBox(uiRect, "PLANET CREATOR")) {
            state->creationMode = false;
        }

        int startX = uiX + 10;
        int startY = uiY + 30;

        // Mass
        GuiLabel((Rectangle){startX, startY, 260, 20}, TextFormat("Mass: %.1e kg", state->newBodyMass));
        startY += 20;
        // Logarithmic slider for mass?
        // Just a simple slider for now, maybe scaling M_EARTH
        float tempMass = (float)(state->newBodyMass / M_EARTH);
        GuiSlider((Rectangle){startX, startY, 200, 20}, NULL, NULL, &tempMass, 0.1f, 1000.0f);
        state->newBodyMass = (double)tempMass * M_EARTH;
        startY += 30;

        // Height Offset
        GuiLabel((Rectangle){startX, startY, 260, 20}, TextFormat("Height Offset: %.1f AU", state->spawnHeight / AU));
        startY += 20;
        float tempHeight = (float)(state->spawnHeight / AU);
        GuiSlider((Rectangle){startX, startY, 200, 20}, NULL, NULL, &tempHeight, -2.0f, 2.0f);
        state->spawnHeight = (double)tempHeight * AU;
        startY += 30;

        // Launch Angle
        GuiLabel((Rectangle){startX, startY, 260, 20}, TextFormat("Launch Angle: %.1f deg", state->spawnAngle));
        startY += 20;
        float tempAngle = state->spawnAngle;
        GuiSlider((Rectangle){startX, startY, 200, 20}, NULL, NULL, &tempAngle, -90.0f, 90.0f);
        state->spawnAngle = tempAngle;
        startY += 30;

        GuiLabel((Rectangle){startX, startY, 260, 20}, "Click & Drag in space to launch");
    }

    if (state->currentState == STATE_SIMULATION) {
        // Draw Planet Indicators
        for (int i = 0; i < MAX_BODIES; i++) {
            if (!bodies[i].active) continue;
            
            DVector3 relPos = DVector3Subtract(bodies[i].position, state->camPos);
            Vector3 center = ToVector3(relPos);
            
            // Check if in front of camera
            Vector3 camForward = Vector3Normalize(Vector3Subtract(state->renderCamera.target, state->renderCamera.position));
            Vector3 toBody = Vector3Normalize(center);
            
            if (Vector3DotProduct(toBody, camForward) > 0.0f) {
                Vector2 screenPos = GetWorldToScreen(center, state->renderCamera);
                
                // Calculate visual size to decide if we need a marker
                double distToCam = DVector3Length(relPos);
                double screenDiameter = (bodies[i].radius * 2.0 / distToCam) * 1000.0; // Approx pixels
                
                if (screenDiameter < 10.0) { // If smaller than 10 pixels
                    DrawCircleLines((int)screenPos.x, (int)screenPos.y, 10.0f, Fade(bodies[i].color, 0.5f));
                }
            }
        }

        Ray ray = GetMouseRay(GetMousePosition(), state->renderCamera);
        int hitIndex = -1;
        float minHitDist = 1e9f;
        
        // Proximity selection variables
        int closestProximityIndex = -1;
        float minScreenDist = 30.0f; // 30 pixels radius tolerance for "near miss" clicks
        Vector2 mousePos = GetMousePosition();

        for (int i = 0; i < MAX_BODIES; i++) {
            if (!bodies[i].active) continue;
            
            DVector3 relPos = DVector3Subtract(bodies[i].position, state->camPos);
            Vector3 center = ToVector3(relPos);
            
            double distToCam = DVector3Length(relPos);
            double minVisualSize = distToCam * 0.002;
            double hitRadius = fmax(bodies[i].radius, minVisualSize);

            // 1. Exact Raycast (Priority)
            RayCollision collision = GetRayCollisionSphere(ray, center, (float)hitRadius);
            if (collision.hit) {
                if (collision.distance < minHitDist) {
                    minHitDist = collision.distance;
                    hitIndex = i;
                }
            }

            // 2. Screen Space Proximity (Fallback for small/distant objects)
            // Check if object is in front of camera
            Vector3 toBody = Vector3Normalize(center);
            Vector3 camForward = Vector3Normalize(Vector3Subtract(state->renderCamera.target, state->renderCamera.position));
            
            if (Vector3DotProduct(toBody, camForward) > 0.0f) {
                Vector2 screenPos = GetWorldToScreen(center, state->renderCamera);
                float screenDist = Vector2Distance(mousePos, screenPos);
                
                if (screenDist < minScreenDist) {
                    minScreenDist = screenDist;
                    closestProximityIndex = i;
                }
            }
        }

        // If no direct hit, use the closest proximity hit
        if (hitIndex == -1) {
            hitIndex = closestProximityIndex;
        }

        if (hitIndex != -1) {
            if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && !state->creationMode) {
                state->cameraTargetIndex = hitIndex;
            }

            Body *b = &bodies[hitIndex];
            double speed = DVector3Length(b->velocity);
            double distToSun = DVector3Distance(b->position, bodies[0].position);
            double kineticE = 0.5 * b->mass * speed * speed;

            char infoText[512];
            sprintf(infoText, "Mass: %.2e kg\nSpeed: %.1f km/s\nDist to Sun: %.2f AU\nKinetic E: %.1e J",
                    b->mass, speed / 1000.0, distToSun / AU, kineticE);

            // We need screen pos of the body
            DVector3 relPos = DVector3Subtract(b->position, state->camPos);
            Vector2 screenPos = GetWorldToScreen(ToVector3(relPos), state->renderCamera);
            
            DrawRectangle(screenPos.x + 20, screenPos.y - 60, 220, 110, Fade(DARKGRAY, 0.9f));
            DrawRectangleLines(screenPos.x + 20, screenPos.y - 60, 220, 110, WHITE);
            DrawText(infoText, screenPos.x + 25, screenPos.y - 55, 10, WHITE);
        }
    }

    DrawFPS(10, 10);
    DrawText(TextFormat("Time Scale: %.0fx", state->timeScale), 10, 30, 20, WHITE);
    DrawText(state->integrator == INTEGRATOR_RK4 ? "Integrator: RK4 (Accurate)" : "Integrator: Verlet (Fast)", 10, 50, 20, state->integrator == INTEGRATOR_RK4 ? GREEN : SKYBLUE);

    int statusY = 70;
    int x = 10;

    DrawText("Roche Limit:", x, statusY, 20, LIGHTGRAY);
    x += MeasureText("Roche Limit: ", 20);
    DrawText(state->enableRoche ? "ON" : "OFF", x, statusY, 20, state->enableRoche ? GREEN : RED);
    x += MeasureText("ON ", 20) + 10;

    DrawText("| Accretion:", x, statusY, 20, LIGHTGRAY);
    x += MeasureText("| Accretion: ", 20);
    DrawText(state->enableDrag ? "ON" : "OFF", x, statusY, 20, state->enableDrag ? GREEN : RED);
    x += MeasureText("ON ", 20) + 10;

    DrawText("| Lagrange:", x, statusY, 20, LIGHTGRAY);
    x += MeasureText("| Lagrange: ", 20);
    DrawText(state->showLagrange ? "ON" : "OFF", x, statusY, 20, state->showLagrange ? GREEN : RED);
    x += MeasureText("ON ", 20) + 10;

    DrawText("| Create (N):", x, statusY, 20, LIGHTGRAY);
    x += MeasureText("| Create (N): ", 20);
    DrawText(state->creationMode ? "ON" : "OFF", x, statusY, 20, state->creationMode ? GREEN : RED);

    x += MeasureText("ON ", 20) + 20;
    
    if (GuiButton((Rectangle){x, statusY - 2, 100, 24}, state->relativeView ? "Rel. View: ON" : "Rel. View: OFF")) {
        state->relativeView = !state->relativeView;
    }

    if (state->creationMode && state->currentState == STATE_SIMULATION) {
        DrawText(TextFormat("Creation Mode: Click & Drag to launch. Scroll: Mass %.1e", state->newBodyMass), 10, 100, 20, YELLOW);
    }

    if (state->currentState == STATE_MENU) {
        DrawRectangle(0, 0, screenWidth, screenHeight, Fade(BLACK, 0.7f));
        int menuX = screenWidth / 2 - 100;
        int menuY = screenHeight / 2 - 100;

        DrawText("PAUSED", screenWidth/2 - MeasureText("PAUSED", 40)/2, menuY - 60, 40, WHITE);

        if (GuiButton((Rectangle){menuX, menuY, 200, 40}, "Resume")) state->currentState = STATE_SIMULATION;
        if (GuiButton((Rectangle){menuX, menuY + 50, 200, 40}, "Reset")) {
            initBodies(bodies);
            state->currentState = STATE_SIMULATION;
        }
        if (GuiButton((Rectangle){menuX, menuY + 100, 200, 40}, "Settings")) state->currentState = STATE_SETTINGS;
        if (GuiButton((Rectangle){menuX, menuY + 150, 200, 40}, "Quit")) state->shouldExit = true;
    }
    else if (state->currentState == STATE_SETTINGS) {
        DrawRectangle(0, 0, screenWidth, screenHeight, Fade(BLACK, 0.8f));
        int menuX = screenWidth / 2 - 150;
        int menuY = screenHeight / 2 - 100;

        DrawText("SETTINGS", screenWidth/2 - MeasureText("SETTINGS", 40)/2, menuY - 60, 40, WHITE);

        GuiToggle((Rectangle){menuX, menuY, 300, 40}, "Roche Limit", &state->enableRoche);
        GuiToggle((Rectangle){menuX, menuY + 50, 300, 40}, "Accretion Drag", &state->enableDrag);
        GuiToggle((Rectangle){menuX, menuY + 100, 300, 40}, "Lagrange Points", &state->showLagrange);
        
        bool isVerlet = (state->integrator == INTEGRATOR_VERLET);
        bool prevVerlet = isVerlet;
        GuiToggle((Rectangle){menuX, menuY + 150, 300, 40}, "Use Verlet Integrator", &isVerlet);
        
        if (isVerlet != prevVerlet) {
            state->integrator = isVerlet ? INTEGRATOR_VERLET : INTEGRATOR_RK4;
            // Recalculate accelerations if switching to Verlet
            if (state->integrator == INTEGRATOR_VERLET) {
                calculateAccelerations(bodies, state->enableDrag);
            }
        }

        if (GuiButton((Rectangle){menuX, menuY + 210, 300, 40}, "Back")) state->currentState = STATE_MENU;
    }
}

int main(void)
{
    const int screenWidth = 1200;
    const int screenHeight = 800;

    InitWindow(screenWidth, screenHeight, "Solar System Simulation 3D - Solaray");
    SetTargetFPS(60);
    SetExitKey(KEY_NULL);

    SimulationState state;
    Body bodies[MAX_BODIES];
    initBodies(bodies);
    InitSimulation(&state);

    // Initialize accelerations for Verlet
    calculateAccelerations(bodies, state.enableDrag);

    while (!state.shouldExit && !WindowShouldClose()) {
        HandleInput(&state, bodies, screenWidth, screenHeight);
        UpdatePhysics(&state, bodies);

        BeginDrawing();
        ClearBackground(BLACK);

        Draw3DScene(&state, bodies);
        DrawUI(&state, bodies, screenWidth, screenHeight);

        EndDrawing();
    }
    CloseWindow();
    return 0;
}
