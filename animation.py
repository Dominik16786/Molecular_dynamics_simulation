import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, gridspec

N = 201               
L = 25.0              
dt = 0.002           
rcut = 2.7            
Tref = 119.0
T_target_default = 300.0 / Tref   
mass = 1.0
Q = 1e20              
nsteps = 4000       
frame_interval = 10   
np.random.seed(12345)


def init_positions(N, L):
    """Place N particles on (rough) square lattice inside box L x L."""
    k = int(np.ceil(np.sqrt(N)))
    delta = L / k
    x = np.zeros(N)
    y = np.zeros(N)
    for i in range(N):
        x[i] = (i % k + 0.5) * delta
        y[i] = (i // k + 0.5) * delta
    return x, y

def init_velocities(N, T):
    """Draw velocities from normal(0, sqrt(T)), then remove center-of-mass motion."""
    sigma = np.sqrt(T)
    vx = np.random.normal(0.0, sigma, size=N)
    vy = np.random.normal(0.0, sigma, size=N)
    vx -= vx.mean()
    vy -= vy.mean()
    return vx, vy

def compute_forces(x, y, L, rcut, mass):
    """
    Compute pairwise forces using LJ 12-6, with the same cutoff shifting
    (value & slope shifted) as in your C++ compute_forces.
    Returns (fx, fy, Epot)
    """
    N = x.size
    fx = np.zeros(N)
    fy = np.zeros(N)
    Epot = 0.0
    rcut2 = rcut * rcut

    inv_rcut = 1.0 / rcut
    inv_rcut6 = inv_rcut**6
    inv_rcut12 = inv_rcut6 * inv_rcut6
    Ucut = 4.0 * (inv_rcut12 - inv_rcut6)
    dUcut = -48.0 * inv_rcut12 * inv_rcut + 24.0 * inv_rcut6 * inv_rcut

    for i in range(N):
        xi = x[i]
        yi = y[i]
        for j in range(i+1, N):
            dx = xi - x[j]
            dy = yi - y[j]
            # periodic boundaries (minimum image)
            if dx >  L/2: dx -= L
            if dx < -L/2: dx += L
            if dy >  L/2: dy -= L
            if dy < -L/2: dy += L
            r2 = dx*dx + dy*dy
            if r2 < rcut2 and r2 > 1e-12:
                r = np.sqrt(r2)
                inv_r = 1.0 / r
                inv_r6 = inv_r**6
                inv_r12 = inv_r6 * inv_r6
                U = 4.0 * (inv_r12 - inv_r6)
                # derivative of U wrt r (as in C++ code)
                dU = -48.0 * inv_r12 * inv_r + 24.0 * inv_r6 * inv_r
                Ushift = U - Ucut - dUcut * (r - rcut)
                # force magnitude from shifted potential: F = -(dU - dUcut) / r
                F = -(dU - dUcut) / r
                fx[i] += F * dx
                fy[i] += F * dy
                fx[j] -= F * dx
                fy[j] -= F * dy
                Epot += Ushift

    # accelerations would be fx/mass etc outside
    return fx, fy, Epot

def kinetic_energy(vx, vy, mass):
    """Compute total kinetic energy 0.5 m v^2 per particle (mass uniform)."""
    return 0.5 * mass * np.sum(vx*vx + vy*vy)


class MDSystem:
    def __init__(self, N, L, dt, rcut, mass=1.0, xi_init=0.0):
        self.N = N
        self.L = L
        self.dt = dt
        self.rcut = rcut
        self.mass = mass
        self.x, self.y = init_positions(N, L)
        self.vx, self.vy = init_velocities(N, T_target_default)
        self.ax = np.zeros(N)
        self.ay = np.zeros(N)
        self.fx = np.zeros(N)
        self.fy = np.zeros(N)
        # initial forces and accelerations:
        self.fx, self.fy, _ = compute_forces(self.x, self.y, self.L, self.rcut, self.mass)
        self.ax = self.fx / self.mass
        self.ay = self.fy / self.mass
        self.xi = xi_init   # thermostat variable (static in C++ file-scope)
        self.nd = 2 * N

    def step(self, Q, T_target):
        """Perform one time step using the same algorithm as your C++ time_it."""
        dt = self.dt
        # half-step predictor for velocities including thermostat damping term
        vxp = self.vx + 0.5 * dt * (self.ax - self.xi * self.vx)
        vyp = self.vy + 0.5 * dt * (self.ay - self.xi * self.vy)

        # update positions
        self.x += vxp * dt
        self.y += vyp * dt
        # periodic boundaries
        self.x = np.mod(self.x, self.L)
        self.y = np.mod(self.y, self.L)

        # recompute forces at new positions
        self.fx, self.fy, Epot = compute_forces(self.x, self.y, self.L, self.rcut, self.mass)
        self.ax = self.fx / self.mass
        self.ay = self.fy / self.mass

        # iterative thermostat update (the loop of 5 iterations in your C++)
        for k in range(5):
            Ek = kinetic_energy(self.vx, self.vy, self.mass)
            # xi update (same formula)
            self.xi += dt * (2.0 * Ek - self.nd * T_target) / Q
            # update velocities using current xi and computed accelerations
            denom = 1.0 + 0.5 * dt * self.xi
            # avoid division by zero in pathological xi
            if denom == 0.0:
                denom = 1e-12
            self.vx = (vxp + 0.5 * dt * self.ax) / denom
            self.vy = (vyp + 0.5 * dt * self.ay) / denom

        Ek = kinetic_energy(self.vx, self.vy, self.mass)
        return Ek, Epot


sim = MDSystem(N=N, L=L, dt=dt, rcut=rcut, mass=mass, xi_init=0.0)
sim.vx, sim.vy = init_velocities(N, T_target_default)
sim.fx, sim.fy, Epot0 = compute_forces(sim.x, sim.y, sim.L, sim.rcut, sim.mass)
sim.ax = sim.fx / sim.mass
sim.ay = sim.fy / sim.mass

# arrays to collect time-series
nframes = nsteps // frame_interval
Ekin_arr = np.zeros(nframes)
Epot_arr = np.zeros(nframes)
Etot_arr = np.zeros(nframes)
Tinst_arr = np.zeros(nframes)
times = np.zeros(nframes)


fig = plt.figure(figsize=(10, 5))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.1])

ax_anim = fig.add_subplot(gs[0])
ax_anim.set_xlim(0, L)
ax_anim.set_ylim(0, L)
ax_anim.set_aspect('equal')
ax_anim.set_title('Particle positions')

scat = ax_anim.scatter(sim.x, sim.y, s=20)

ax_e = fig.add_subplot(gs[1])
line_kin, = ax_e.plot([], [], label='E_kin')
line_pot, = ax_e.plot([], [], label='E_pot')
line_tot, = ax_e.plot([], [], label='E_tot')
ax_e.set_xlim(0, nframes * dt * frame_interval)
ax_e.set_ylim(-5.0, 50.0)
ax_e.set_xlabel('time')
ax_e.set_ylabel('Energy')
ax_e.legend(loc='upper left')

# Temperature axis (shared x)
ax_T = ax_e.twinx()
line_T, = ax_T.plot([], [], 'k--', label='T_inst')
ax_T.set_ylabel('T_inst')

# text on animation showing step/time
time_text = ax_anim.text(0.02, 0.98, '', transform=ax_anim.transAxes, va='top')

def init_anim():
    scat.set_offsets(np.c_[sim.x, sim.y])
    line_kin.set_data([], [])
    line_pot.set_data([], [])
    line_tot.set_data([], [])
    line_T.set_data([], [])
    time_text.set_text('')
    return scat, line_kin, line_pot, line_tot, line_T, time_text

frame_index = 0
def update(frame):
    global frame_index
    for k in range(frame_interval):
        Ek, Epot = sim.step(Q=Q, T_target=T_target_default)
    t = (frame * frame_interval + 1) * dt

    Ekin_arr[frame] = Ek
    Epot_arr[frame] = Epot
    Etot_arr[frame] = Ek + Epot
    Tinst_arr[frame] = Ek / N
    times[frame] = t

    scat.set_offsets(np.c_[sim.x, sim.y])

    times_plot = times[:frame+1]
    line_kin.set_data(times_plot, Ekin_arr[:frame+1])
    line_pot.set_data(times_plot, Epot_arr[:frame+1])
    line_tot.set_data(times_plot, Etot_arr[:frame+1])
    line_T.set_data(times_plot, Tinst_arr[:frame+1])
    all_e = np.concatenate([Ekin_arr[:frame+1], Epot_arr[:frame+1], Etot_arr[:frame+1]])
    ymin, ymax = np.min(all_e) - 0.1*np.abs(np.min(all_e)+1e-12), np.max(all_e) + 0.1*np.abs(np.max(all_e)+1e-12)
    if ymax <= ymin:
        ymin, ymax = -1, 1
    ax_e.set_ylim(ymin, ymax)

    Tmax = max(1.0, np.max(Tinst_arr[:frame+1]) * 1.2)
    ax_T.set_ylim(0, Tmax)

    time_text.set_text(f't = {t:.3f}')
    frame_index += 1
    return scat, line_kin, line_pot, line_tot, line_T, time_text

anim = animation.FuncAnimation(fig, update, init_func=init_anim,
                               frames=nframes, interval=30, blit=False)

plt.tight_layout()
plt.show()
