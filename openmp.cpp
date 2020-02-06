#include "common.h"
#include <omp.h>

void init_simulation(particle_t* parts, int num_parts, double size, void **my_data) {
	// You can use this space to initialize data objects that you may need
	// This function will be called once before the algorithm begins
	// Do not do any particle simulation here
	// The my_data object will persist between calls to simulate_one_step
	// So you can create objects and store them there for later use
}

void simulate_one_step(particle_t* parts, int num_parts, double size, void *my_data) {
    // Write this function
}
