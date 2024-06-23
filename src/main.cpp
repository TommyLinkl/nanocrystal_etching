#include "main.h"

int main(int argc, char** argv)
{
    const double epsilon = 1.0;
    const double temperature = 0.1 * epsilon;
    const double mu = -6.0 * epsilon;
    const double radius = 5.0; // In units of lattice constant
    const double k_D = exp(-3.0 * epsilon / temperature); // A reasonable choice for the diffusion barrier

    // Simulation details
    int max_size = 50; // Size of lattice (should just be set large enough so that there will be no interactions with the boundaries)
    int trial_number = 0;
    std::stringstream filebase;
    filebase << "test" << std::fixed << std::setprecision(2) << temperature << "_" << mu;
    blocked_crystal nanocrystal(sphere, radius, temperature, mu, epsilon, k_D, max_size, filebase.str(), "TestFiles/", trial_number);
}