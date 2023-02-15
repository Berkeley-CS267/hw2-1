#include "common.h"
#include <cmath>
#include <set>
#include <unordered_map>
#include <iostream>
#include <bits/stdc++.h>

std::unordered_map<int, std::unordered_set<particle_t *>> bins;
std::unordered_map<int, int> particle_to_bin;
// std::unordered_map<int, std::tuple<double, double>> force_vectors;
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
    neighbor.ax -= coef * dx;
    neighbor.ay -= coef * dy;
    // return std::make_tuple(coef * dx, coef * dy);
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


void apply_force_self(particle_t& particle, particle_t& neighbor) {
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
    // return std::make_tuple(coef * dx, coef * dy);
}

// Integrate the ODE
void move(int particle_ind, particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method

    int origin_bin = particle_to_bin[particle_ind];
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

    particle_to_bin[particle_ind] = new_bin;
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
        particle_to_bin[i] = index;
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

void apply_force_bin(int row, int column, int row2, int column2, int lda){
    // Force every particle in each bin to interact
    if (!check_boundary(row, column) && !check_boundary(row2, column2)){
        return;
    }
    for (auto it2 = bins[column2+row2*lda].begin(); it2 != bins[column2+row2*lda].end(); ++it2){
        for (auto it = bins[column+row*lda].begin(); it != bins[column+row*lda].end(); ++it){
            // Interact particles
            apply_force(**it2, **it);
        }
    }
}

void apply_force_bin_self(int row, int column, int row2, int column2, int lda){
    // Force every particle in each bin to interact
    if (!check_boundary(row, column) && !check_boundary(row2, column2)){
        return;
    }
    for (auto it2 = bins[column2+row2*lda].begin(); it2 != bins[column2+row2*lda].end(); ++it2){
        for (auto it = bins[column+row*lda].begin(); it != bins[column+row*lda].end(); ++it){
            // Interact particles
            apply_force_self(**it2, **it);
        }
    }
}

/*
void apply_force_bin(std::unordered_set<particle_t*> *bin_1, std::unordered_set<particle_t*> *bin_2, bool same = false) {
    // We will calculate the outer_product of the two vectors
    // Matrix of dx,dy force corrections
    std::tuple<double, double> force_matrix[bin_1->size()][bin_2->size()];
    for (auto i = 0; i < bin_1->size(); ++i){
        for (auto j = 0; j < bin_2->size(); ++j){
            force_matrix[i][j] = apply_force(*bin_1->at(i), *bin_2->at(j));
        }
    }
    // Ok now reduce across for each particle
    for (auto i = 0; i < bin_1->size(); ++i){
        // For each particle in bin_1
        double dx_sum = 0.0;
        double dy_sum = 0.0;
        for (auto j = 0; j < bin_2->size(); ++j){
            auto item = force_matrix[i][j];
            dx_sum += std::get<0>(item);
            dy_sum += std::get<0>(item);
        }
        bin_1->at(i)->ax += dx_sum;
        bin_1->at(i)->ay += dy_sum;
    }
    if (!same) {
        // Now reduce down for each particle, but subtract bc wrong direction
        for (auto j = 0; j < bin_2->size(); ++j){
            // For each particle in bin_1
            double dx_sum = 0.0;
            double dy_sum = 0.0;
            for (auto i = 0; i < bin_2->size(); ++i){
                auto item = force_matrix[i][j];
                dx_sum -= std::get<0>(item);
                dy_sum -= std::get<0>(item);
            }
            bin_2->at(j)->ax += dx_sum;
            bin_2->at(j)->ay += dy_sum;
        }
    }
}
*/

void apply_force_bins(int row, int column, int lda){
    // Create a matrix of force for each bin to bin

    // Loop over 5 forward neighbor bins (3 below, 1 to left, and itself)
    apply_force_bin_self(row, column, row, column, lda);
    // right neighbor
    apply_force_bin(row, column, row, column + 1, lda);
    // bottom left neighbor
    apply_force_bin(row, column, row + 1, column - 1, lda);
    // bottom neighbor
    apply_force_bin(row, column, row + 1, column, lda);
    // bottom right neighbor
    apply_force_bin(row, column, row + 1, column + 1, lda);
}



void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Compute Forces
    for(int i = 0; i < num_parts; ++i){
        parts[i].ax = parts[i].ay = 0;
    }
    for(int i = 0; i < lda; ++i){
        for (int j = 0; j < lda; ++j){
            apply_force_bins(i, j, lda);
        }
    }
    // Move Particles
    for (int i = 0; i < num_parts; ++i) {
        move(i, parts[i], size);
    }

}