#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <cstring>
#include "common.h"


// =================
// Helper Functions
// =================

// I/O routines
void save(std::ofstream &fsave, particle_t *parts, int num_parts, double size) {
    static bool first = true;

    if (first) {
        fsave << num_parts << " " << size << std::endl;
        first = false;
    }

    for (int i = 0; i < num_parts; ++i) {
        fsave << parts[i].x << " " << parts[i].y << std::endl;
    }

    fsave << std::endl;
}


// Particle Initialization
void init_particles(particle_t *parts, int num_parts, double size, int part_seed) {
    srand48(part_seed);

    int sx = (int) ceil(sqrt((double) num_parts));
    int sy = (num_parts + sx - 1) / sx;

    int *shuffle = (int *) malloc(num_parts * sizeof(int));

    for (int i = 0; i < num_parts; ++i) {
        shuffle[i] = i;
    }

    for (int i = 0; i < num_parts; ++i) {
        // Make sure particles are not spatially sorted
        int j = lrand48() % (num_parts - i);
        int k = shuffle[j];
        shuffle[j] = shuffle[num_parts - i - 1];

        // Distribute particles evenly to ensure proper spacing
        parts[i].x = size * (1. + (k % sx)) / (1 + sx);
        parts[i].y = size * (1. + (k / sx)) / (1 + sy);

        // Assign random velocities within a bound
        parts[i].vx = drand48() * 2 - 1;
        parts[i].vy = drand48() * 2 - 1;
    }

    free(shuffle);
}


// Command Line Option Processing
int find_arg_idx(int argc, char **argv, const char *option) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], option) == 0) {
            return i;
        }
    }
    return -1;
}


int find_int_arg(int argc, char **argv, const char *option, int default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return std::stoi(argv[iplace + 1]);
    }

    return default_value;
}


char *find_string_option(int argc, char **argv, const char *option, char *default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return argv[iplace + 1];
    }

    return default_value;
}


// ==============
// Main Function
// ==============

int main(int argc, char **argv) {
    // Parse Args
    if (find_arg_idx(argc, argv, "-h") >= 0) {
        std::cout << "Options:" << std::endl;
        std::cout << "-h: see this help" << std::endl;
        std::cout << "-n <int>: set number of particles" << std::endl;
        std::cout << "-o <filename>: set the output file name" << std::endl;
        std::cout << "-s <int>: set particle initialization seed" << std::endl;
        return 0;
    }

    // Open Output File
    char *savename = find_string_option(argc, argv, "-o", nullptr);
    std::ofstream fsave(savename);

    // Initialize Particles
    int num_parts = find_int_arg(argc, argv, "-n", 1000);
    int part_seed = find_int_arg(argc, argv, "-s", 0);
    double size = sqrt(density * num_parts);

    particle_t *parts = new particle_t[num_parts];

    init_particles(parts, num_parts, size, part_seed);

    // Algorithm
    auto start_time = std::chrono::steady_clock::now();

#ifdef _OPENMP
#pragma omp parallel default(shared)
    {
#endif

        for (int step = 0; step < nsteps; ++step) {
            simulate_one_step(parts, num_parts, size);

            // Save state if necessary
#ifdef _OPENMP
#pragma omp master
#endif
            if (fsave.good() && (step % savefreq) == 0) {
                save(fsave, parts, num_parts, size);
            }
        }

#ifdef _OPENMP
    }
#endif

    auto end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();

    // Finalize
    std::cout << "Simulation Time = "
              << seconds
              << " seconds for "
              << num_parts
              << " particles."
              << std::endl;
    fsave.close();
    delete[] parts;
}
