# Solar System Simulation 3D (N-Body RK4)

Educational project designed to explore numerical integration methods and astrophysical concepts in a fully interactive 3D environment. It is not a scientifically accurate simulator of the real solar system, but rather a sandbox for visualizing mathematical models of gravity, orbital mechanics, and celestial phenomena. Made also to learn `raylib`. 

Based on [2D version](https://github.com/lukasz-strama/solaray).

![Program Demo](demo.gif)

## Mathematical Models

### 1. Numerical Integration
The simulation uses the **Runge-Kutta 4th Order (RK4)** method for solving the differential equations of motion. This provides significantly higher stability and accuracy compared to Euler or Verlet integration, especially for eccentric orbits.

Given state $y$ (position and velocity) and derivative function $f(t, y)$:

$$
\begin{aligned}
    k_1 &= f(t_n, y_n) \\
    k_2 &= f\left(t_n + \frac{h}{2}, y_n + h \frac{k_1}{2}\right) \\
    k_3 &= f\left(t_n + \frac{h}{2}, y_n + h \frac{k_2}{2}\right) \\
    k_4 &= f\left(t_n + h, y_n + h k_3\right) \\
    y_{n+1} &= y_n + \frac{h}{6}\left(k_1 + 2k_2 + 2k_3 + k_4\right)
\end{aligned}
$$

### 2. Gravity & General Relativity
The base interaction is Newtonian gravity, augmented with a Post-Newtonian correction to simulate the **precession of the perihelion** (similar to Mercury's orbit).

$$ F_{total} = F_{Newton} \left( 1 + \frac{3h^2}{c^2 r^2} \right) $$

Where:
*   $F_{Newton} = G \frac{m_1 m_2}{r^2}$
*   $h = |\vec{r} \times \vec{v}|$ (Specific angular momentum)
*   $c$ is the speed of light (scaled for simulation visibility)

### 3. Collisions (Volume Conservation)
Collisions are perfectly inelastic. When two bodies merge, momentum is conserved, and the new radius is calculated assuming **Volume Conservation** (constant density).

$$ r_{new} = \sqrt[3]{r_1^3 + r_2^3} $$

### 4. Roche Limit (Tidal Disruption)
Bodies approaching a massive primary too closely will disintegrate due to tidal forces. The simulation uses the **Rigid Body** approximation for the Roche Limit:

$$ d_{limit} \approx 1.26 R_{primary} \sqrt[3]{\frac{M_{primary}}{M_{satellite}}} $$

If $d < d_{limit}$, the body explodes into smaller fragments.

### 5. Accretion Disk Drag
A simplified drag model simulates gas/dust density near the central body.
*   **Density Profile**: Linear gradient decreasing from $R_{in}$ to $R_{out}$.
*   **Drag Force**:

$$
\vec{F}_{drag} = -C_{drag} \cdot \rho(r) \cdot \vec{v}
$$

## Assumptions & Simplifications

*   **3D Space**: The simulation runs in a full 3D coordinate system. Orbits can have inclination and bodies can move freely on the Y-axis.
*   **Scaled Constants**: $G$ and $c$ (speed of light) are arbitrary game units chosen for visual clarity, not SI units.
*   **Perfect Spheres**: All bodies are treated as point masses for gravity and spheres for collisions.
*   **Lagrange Points**: Calculated analytically for the Sun and the most massive planet, then projected onto the orbital plane for visualization.

## Current Status

### Working Correctly
*   **N-Body Gravity**: All bodies attract all other bodies $O(N^2)$.
*   **3D Orbital Mechanics**: Support for inclined orbits and 3D velocity vectors.
*   **Orbital Stability**: RK4 maintains stable orbits for long durations.
*   **Precession**: Relativistic apsidal precession is observable.
*   **Creation Mode**: User can inject new bodies with specific mass, height offset, and launch angle.
*   **Orbit Editor**: Real-time adjustment of orbital parameters (Eccentricity, Semi-Major Axis, Inclination).

### Limitations / Known Issues
*   **Scale**: Distances and masses are not to scale with the real solar system.
*   **Fragmentation Physics**: When the Roche limit is breached, bodies spawn random fragments. This is a visual approximation.
*   **Performance**: The $O(N^2)$ complexity limits the simulation to approximately 500 bodies before frame rate drops on average hardware.
*   **No Atmospheric Effects**: Drag is a simple linear model; no complex fluid dynamics are simulated.

### TO-DO
*   **Improved Collision Handling**: More sophisticated fragmentation and debris field simulation.
*   **Enhanced Visuals**: Add textures, lighting effects, and particle systems.
*   **User Interface**: More intuitive controls and information display.
*   **Scale Constants**: Implement a more realistic scaling system for distances and masses.
*   **Loading Presets**: Predefined solar system configurations imitating real systems.

## Controls

*   **Mouse Wheel**: Zoom Camera
*   **Right Mouse Drag**: Rotate Camera
*   **Middle Mouse Drag**: Pan Camera
*   **Left Click (on Body)**: Select/Focus Camera on Body
*   **Arrow Keys**: Change Time Scale
*   **Space**: Pause/Resume Simulation
*   **N**: Toggle Creation Mode
    *   **Left Click & Drag**: Launch new body
    *   **UI Sliders**: Adjust Mass, Height Offset, Launch Angle
*   **ESC**: Open Menu

## Build

Requires `raylib` and a C compiler.

```bash
./nob
./main
```
