#include "main.h"

int main(int argc, char** argv)
{
    // Establish some physical parameters 
    const double epsilon = 1.0;
    const double temperature = 0.1 * epsilon; // Resonable temperature at which facets will be smooth
    const double mu = std::atof(argv[1]) * epsilon;
    const double radius = 20.0; // In units of lattice constant
    const double k_D = exp(-2.0 * epsilon / temperature); // A reasonable choice for the diffusion barrier

    // Simulation details
    int max_size = 150; // Size of lattice should just be set large enough so that there will be no interactions with the boundaries
    int trial_number = 0;
    std::stringstream filebase;
    filebase << "test" << std::fixed << std::setprecision(2) << temperature << "_" << mu << "_" << radius;
    
    // Create an instance of the nanocrystal object
    nanocrystal nc(sphere, radius, temperature, mu, epsilon, k_D, max_size, filebase.str(), "../TestFiles/", trial_number);
    
    // Run a simulation
    unsigned long long num_steps = 100; //5e6;
    int config_records = 5e2;
    int data_records = 1e4;
    int seed = 101; // Obviously get this from an external source for general use
    nc.kmc_run(num_steps, config_records, data_records, seed, false);
}   