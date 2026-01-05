#include "particle.h"
#include <vector>
#include <cmath>
#include <random>
#include<fstream>

void init_positions(std::vector<Particle>& p, double L){
    int N = p.size();
    int k = (int)std::ceil(std::sqrt((double)N));
    double delta = L / k;

    for (int i = 0; i < N; ++i) {
        p[i].x = (i % k + 0.5) * delta;
        p[i].y = (i / k + 0.5) * delta;
    }
}

void init_velocities(std::vector<Particle>& p, double T){
    static std::mt19937 gen(12345);
    std::normal_distribution<double> normal(0.0, 1.0);

    int N = p.size();
    double sigma = std::sqrt(T);

    for (int i = 0; i < N; i++){
        p[i].vx = sigma * normal(gen);
        p[i].vy = sigma * normal(gen);
    }

    double vxcm = 0.0, vycm = 0.0;
    for (int i = 0; i < N; i++){
        vxcm += p[i].vx;
        vycm += p[i].vy;
    }
    vxcm /= N;
    vycm /= N;

    for (int i = 0; i < N; i++){
        p[i].vx -= vxcm;
        p[i].vy -= vycm;
    }
}

double compute_forces(std::vector<Particle>& p, double L, double rcut){
    int N = p.size();
    double Epot = 0.0;
    double rcut2 = rcut * rcut;

    for (auto& pi : p){
        pi.fx = pi.fy = 0.0;
    }

    double inv_rcut = 1.0 / rcut;
    double inv_rcut6 = std::pow(inv_rcut, 6);
    double inv_rcut12 = inv_rcut6 * inv_rcut6;
    double Ucut = 4.0 * (inv_rcut12 - inv_rcut6);
    double dUcut = -48.0 * inv_rcut12 * inv_rcut + 24.0 * inv_rcut6 * inv_rcut;

    for (int i = 0; i < N; i++){
        for (int j = i + 1; j < N; j++){
            double dx = p[i].x - p[j].x;
            double dy = p[i].y - p[j].y;

            if (dx >  L/2) dx -= L;
            if (dx < -L/2) dx += L;
            if (dy >  L/2) dy -= L;
            if (dy < -L/2) dy += L;

            double r2 = dx*dx + dy*dy;
            if (r2 < rcut2 && r2 > 0.0){
                double r = std::sqrt(r2);
                double inv_r = 1.0 / r;
                double inv_r6 = std::pow(inv_r, 6);
                double inv_r12 = inv_r6 * inv_r6;

                double U = 4.0 * (inv_r12 - inv_r6);
                double dU = -48.0 * inv_r12 * inv_r + 24.0 * inv_r6 * inv_r;

                double Ushift = U - Ucut - dUcut * (r - rcut);
                double F = -(dU - dUcut) / r;

                p[i].fx += F * dx;
                p[i].fy += F * dy;
                p[j].fx -= F * dx;
                p[j].fy -= F * dy;

                Epot += Ushift;
            }
        }
    }

    for (auto& pi : p){
        pi.ax = pi.fx / pi.mass;
        pi.ay = pi.fy / pi.mass;
    }

    return Epot;
}

double kinetic_energy(const std::vector<Particle>& p){
    double Ek = 0.0;
    for (const auto& pi : p){
        Ek += 0.5 * (pi.vx*pi.vx + pi.vy*pi.vy);
    }
    return Ek;
}

void time_it(std::vector<Particle>& p, double dt, double L, double rcut, double Q, double T_target){
    static double xi = 0.0;
    int N = p.size();
    int nd = 2 * N;

    for (auto& pi : p) {
        pi.vxp = pi.vx + 0.5 * dt * (pi.ax - xi * pi.vx);
        pi.vyp = pi.vy + 0.5 * dt * (pi.ay - xi * pi.vy);
    }

    for (auto& pi : p) {
        pi.x += pi.vxp * dt;
        pi.y += pi.vyp * dt;
        if (pi.x < 0) pi.x += L;
        if (pi.x >= L) pi.x -= L;
        if (pi.y < 0) pi.y += L;
        if (pi.y >= L) pi.y -= L;
    }

    compute_forces(p, L, rcut);

    for (int k = 0; k < 5; k++){
        double Ek = kinetic_energy(p);
        xi += dt * (2.0 * Ek - nd * T_target) / Q;

        for (auto& pi : p) {
            pi.vx = (pi.vxp + 0.5 * dt * pi.ax) / (1.0 + 0.5 * dt * xi);
            pi.vy = (pi.vyp + 0.5 * dt * pi.ay) / (1.0 + 0.5 * dt * xi);
        }
    }
}
void velocity_distribution(std::vector<double>& vhist,const std::vector<Particle>& particles, double dv){
    int nv = vhist.size();
    for (const auto& p : particles){
        double v = std::sqrt(p.vx * p.vx + p.vy * p.vy);
        int k = static_cast<int>(std::floor(v / dv));
        if (k >= 0 && k < nv) {
            vhist[k] += 1.0;
        }
    }
}

void to_file(std::string filename, std::vector<Particle>& particles){
    std::ofstream out(filename);
    for (size_t i = 0; i < particles.size(); i++){
        out << i << " "
            << particles[i].x  << " "
            << particles[i].y  << " "
            << particles[i].vx << " "
            << particles[i].vy << " "
            << particles[i].ax << " "
            << particles[i].ay << " "
            << particles[i].fx << " "
            << particles[i].fy << "\n";
    }
    out.close();
}
