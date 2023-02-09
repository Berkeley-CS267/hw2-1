#include "common.h"
#include <cmath>
#include <set>
#include <unordered_map>
#include <iostream>

#define EPS  0.001

std::unordered_map<int, std::set<int>> bins;
double bin_size = 3*cutoff;
int lda;

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}

int calculate_bin_number(double x, double y, double size, double bin_size, int lda){
    double quotient;
    quotient = x/bin_size;
    int column_index =(int) quotient;

    quotient = y/bin_size;
    int row_index=(int)quotient;
    //   std::cout << column_index + lda* row_index << "\n";
    return column_index + lda* row_index;

}


void init_simulation(particle_t* parts, int num_parts, double size) {
    double quotient= size/bin_size;
    lda = (int) ceil(quotient);
    int index;
    for (int i = 0; i < num_parts; ++i){
        index = calculate_bin_number(parts[i].x,parts[i].y, size, bin_size,lda);
        bins[index].insert(i);
    }
    // You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here
}

bool check_boundary(int row , int column){
    if ((row < 0) or (row>=lda)){
        return false;
    }
    if ((column < 0) or (column >=lda)){
        return false;
    }
    return true;
}
void apply_force_bin(int row, int column, int p, particle_t* parts, int lda){
    if (not check_boundary(row, column)){
        return;
    }
    for (auto it = bins[column+row*lda].begin(); it != bins[column+row*lda].end(); ++it){
        apply_force(parts[p], parts[*it]);
    }

}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Compute Forces
    // for (int i = 0; i < num_parts; ++i) {
    //     parts[i].ax = parts[i].ay = 0;
    //     for (int j = 0; j < num_parts; ++j) {
    //         apply_force(parts[i], parts[j]);
    //     }
    // }

    for(int i = 0; i < lda; ++i){
        for (int j = 0; j < lda; ++j){
            for (auto it = bins[i+j*lda].begin(); it != bins[i+j*lda].end(); ++it){
                parts[*it].ax = parts[*it].ay = 0;
                apply_force_bin(j-1,i-1, *it, parts, lda);
                apply_force_bin(j-1,i, *it, parts, lda);
                apply_force_bin(j-1,i+1, *it, parts, lda);
                apply_force_bin(j,i-1, *it, parts, lda);
                apply_force_bin(j,i, *it, parts, lda);
                apply_force_bin(j,i+1, *it, parts, lda);
                apply_force_bin(j+1,i-1, *it, parts, lda);
                apply_force_bin(j+1,i, *it, parts, lda);
                apply_force_bin(j+1,i+1, *it, parts, lda);
            }
        }
    }


    // Move Particles
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }

    for( int i=0; i< lda*lda;i++){
        bins[i].clear();
    }
    int index;
    for (int i = 0; i < num_parts; ++i){
        index = calculate_bin_number(parts[i].x,parts[i].y, size, bin_size, lda);
        bins[index].insert(i);
    }
}