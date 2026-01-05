# 2D Molecular Dynamics Simulation (Lennard–Jones Fluid)

This repository contains a **2D molecular dynamics (MD) simulation** of interacting particles using a **Lennard–Jones (12–6) potential** with periodic boundary conditions.

The project includes:
- A **C++ implementation** used for numerical experiments and data generation
- A **Python implementation** that reproduces the same physics and provides a **real-time animation with energy and temperature plots**

The simulation models a classical Lennard–Jones fluid and demonstrates thermostatting, energy conservation, velocity distributions, and particle trajectories.

---

## Physical Model

### Particles
- Identical particles with mass \( m = 1 \)
- Motion restricted to two spatial dimensions
- Periodic boundary conditions in a square simulation box

### Interaction Potential
Particles interact via a **shifted Lennard–Jones potential**:

\[
U(r) = 4\left(\frac{1}{r^{12}} - \frac{1}{r^6}\right)
\]

The potential and its derivative are shifted at the cutoff radius \( r_c \) such that:
- \( U(r_c) = 0 \)
- \( \frac{dU}{dr}(r_c) = 0 \)

This guarantees continuity of both energy and force.

---

## Numerical Method

### Time Integration
- Velocity Verlet–type integration scheme
- Includes a **Nosé–Hoover–style thermostat**
- A thermostat variable dynamically adjusts the kinetic energy to control temperature

### Thermostat Parameters
- `Q`: thermostat mass (controls coupling strength)
- `T`: target temperature (constant or time-dependent)

---

## Repository Structure


---

## C++ Implementation

### Features
- Square-lattice initialization
- Maxwell–Boltzmann velocity distribution
- Removal of center-of-mass velocity
- Periodic boundary conditions (minimum image convention)
- Lennard–Jones forces with cutoff shifting
- Energy tracking (kinetic, potential, total)
- Velocity histograms
- Single-particle trajectories
- Temperature ramps
- Output to text files for analysis

### Compilation

Requirements:
- C++11 or newer
- Standard library only

Compile with:
```bash
g++ main.cpp functions.cpp -O2 -o md2d
