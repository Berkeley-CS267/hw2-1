#include "common.h"
#include <cmath>
#include <set>
#include <unordered_map>
#include <iostream>
#include <bits/stdc++.h>

std::unordered_map<int, std::unordered_set<particle_t *>> bins;
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

int calculate_bin_number(double x, double y, double size, double bin_size, int lda){
    double quotient;
    quotient = x/bin_size;
    int column_index =(int) quotient;

    quotient = y/bin_size;
    int row_index=(int)quotient;
    //   std::cout << column_index + lda* row_index << "\n";
    return column_index + lda* row_index;

}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method

    int origin_bin = calculate_bin_number(p.x, p.y, size, bin_size, lda);
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

    int new_bin =calculate_bin_number(p.x, p.y, size, bin_size, lda);
    if(origin_bin == new_bin) return;

    bins[origin_bin].erase(&p);
    bins[new_bin].insert(&p);
}



void init_simulation(particle_t* parts, int num_parts, double size) {
    double quotient= size/bin_size;
    lda = (int) ceil(quotient);
    int index;
    const int space = ceil(1.5 * bin_size * bin_size * 1. / density);
    for(int i = 0; i< lda*lda; ++i){
        bins[i].reserve(space);
    }

    for (int i = 0; i < num_parts; ++i){
        index = calculate_bin_number(parts[i].x,parts[i].y, size, bin_size,lda);
        bins[index].insert(&parts[i]);
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

void apply_force_bin(int row, int column, particle_t* p, int lda){
    if (not check_boundary(row, column)){
        return;
    }
    for (auto it = bins[column+row*lda].begin(); it != bins[column+row*lda].end(); ++it){
        apply_force(*p, **it);
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
    for(int i = 0; i < num_parts; ++i){
        parts[i].ax = parts[i].ay = 0;
    }
    for(int i = 0; i < lda; ++i){
        for (int j = 0; j < lda; ++j){
            for (auto it = bins[i+j*lda].begin(); it != bins[i+j*lda].end(); ++it){
                apply_force_bin(j-1,i-1, *it, lda);
                apply_force_bin(j-1,i, *it, lda);
                apply_force_bin(j-1,i+1, *it, lda);
                apply_force_bin(j,i-1, *it, lda);
                apply_force_bin(j,i, *it,  lda);
                apply_force_bin(j,i+1, *it, lda);
                apply_force_bin(j+1,i-1, *it,  lda);
                apply_force_bin(j+1,i, *it, lda);
                apply_force_bin(j+1,i+1, *it,  lda);
            }
        }
    }
    // Move Particles
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }

}
