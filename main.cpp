#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "particle.h"

using namespace std;

double kinetic_energy(const vector<Particle>& p);
double compute_forces(vector<Particle>& p, double L, double rcut);
void init_positions(vector<Particle>& p, double L);
void init_velocities(vector<Particle>& p, double T);
void time_it(vector<Particle>& particles, double dt, double L, double rcut, double Q, double T_target);
void velocity_distribution(vector<double>& vhist, const vector<Particle>& particles, double dv);
void to_file(string filename, vector<Particle>& particles);

int main()
{
    const int N = 201;
    const double L = 25.0;
    const double dt = 0.002;
    const double rcut = 2.7;
    const double Tref = 119.0;
    const double T = 300.0 / Tref;

    vector<Particle> particles(N);
    for (auto& p : particles) p.mass = 1.0;

    init_positions(particles, L);
    init_velocities(particles, T);
    compute_forces(particles, L, rcut);
    to_file("init_positions.txt", particles);

    {
        vector<Particle> p = particles;
        int Nit = 10000;
        double Q = 1e20;

        ofstream energy("energy_task2a.txt");
        for (int it = 0; it < Nit; it++) {
            time_it(p, dt, L, rcut, Q, T);
            if (it % 10 == 0) {
                double Ekin = kinetic_energy(p);
                double Epot = compute_forces(p, L, rcut);
                energy << it * dt << " "
                       << Ekin << " "
                       << Epot << " "
                       << Ekin + Epot << "\n";
            }
        }
        energy.close();
    }

    {
        vector<Particle> p = particles;
        int Nit = 50000;
        double Q = 1e20;

        ofstream traj("trajectory_task2b.txt");
        int l = N / 2;

        for (int it = 0; it < Nit; it++) {
            time_it(p, dt, L, rcut, Q, T);
            if (it % 10 == 0)
                traj << p[l].x << " " << p[l].y << "\n";
        }
        traj.close();
    }

    {
        int Nit = 10000;
        vector<double> Qlist = {1.0, 0.1, 0.01};

        for (size_t idx = 0; idx < Qlist.size(); idx++) {
            vector<Particle> p = particles;
            double Q = Qlist[idx];

            ofstream energy("energy_task3_Q" + to_string(idx + 1) + ".txt");
            ofstream temp("temp_task3_Q" + to_string(idx + 1) + ".txt");

            for (int it = 0; it < Nit; it++) {
                time_it(p, dt, L, rcut, Q, T);
                if (it % 10 == 0) {
                    double Ekin = kinetic_energy(p);
                    double Epot = compute_forces(p, L, rcut);
                    double Tinst = Ekin / N;

                    energy << it * dt << " "
                           << Ekin << " "
                           << Epot << " "
                           << Ekin + Epot << "\n";

                    temp << it * dt << " " << Tinst << "\n";
                }
            }
            energy.close();
            temp.close();
        }
    }

    {
        vector<int> Nit_list = {1000, 10000, 100000};
        double Q = 10.0;

        int nv = 200;
        double vc = sqrt(2.0 * T);
        double vmax = 4.0 * vc;
        double dv = vmax / nv;

        for (size_t idx = 0; idx < Nit_list.size(); idx++) {
            vector<Particle> p = particles;
            vector<double> vhist(nv, 0.0);
            int Nit = Nit_list[idx];

            for (int it = 0; it < Nit; it++) {
                time_it(p, dt, L, rcut, Q, T);
                velocity_distribution(vhist, p, dv);
            }

            for (int k = 0; k < nv; k++)
                vhist[k] /= (Nit * N * dv);

            ofstream vhfile("vhist_task4_" + to_string(idx + 1) + ".txt");
            for (int k = 0; k < nv; k++) {
                double vmid = (k + 0.5) * dv;
                vhfile << vmid << " " << vhist[k] << "\n";
            }
            vhfile.close();
        }
    }

    {
        const double Q = 1.0;
        const double Tmax = 300.0 / Tref;
        const double Tmin = 5.0 / Tref;
        const int nT = 500;

        vector<int> Nit_list = {10000, 100000};

        for (int case_id = 0; case_id < 2; ++case_id) {

            int Nit = Nit_list[case_id];
            int Kit = Nit / nT;
            double dT = (Tmax - Tmin) / nT;

            vector<Particle> p = particles;
            ofstream temp_file("temp_task5_case" + to_string(case_id + 1) + ".txt");

            for (int it = 0; it < Nit; ++it) {

                double T_target = Tmax - dT * floor(double(it) / Kit);
                if (T_target < Tmin) T_target = Tmin;

                time_it(p, dt, L, rcut, Q, T_target);

                double Ekin = kinetic_energy(p);
                double Tinst = Ekin / N;

                temp_file << it * dt << " " << Tinst << "\n";
            }

            temp_file.close();
            to_file("positions_task5_case" + to_string(case_id + 1) + ".txt", p);
        }
    }

    return 0;
}
