#pragma once

class Particle {
public:
    double mass;
    double x, y;
    double vx, vy;
    double vxp, vyp;
    double ax, ay;
    double fx, fy;

    Particle() {
        mass = 1.0;
        x = y = 0.0;
        vx = vy = 0.0;
        vxp = vyp = 0.0;
        ax = ay = 0.0;
        fx = fy = 0.0;
    }
};
