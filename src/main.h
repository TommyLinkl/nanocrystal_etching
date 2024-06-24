//  Created by Sam Oaks-Leaf.
//  A header describing a nanocrystal class designed for etching in the BEP framework for surface diffusion.

#ifndef Etching_h
#define Etching_h

#include "helper.h"
#include <iostream>
#include <array>
#include <iomanip>
#include <random>
#include <functional>
#include <fstream>
#include <string>
#include <sstream>
#include <limits.h>
#include <algorithm>
#include <numeric>

class nanocrystal
{
    class atom
    {
        // Subclass that describes a site on the lattice
        public:
            int lattice_index, vector_index; // position in the lattice vector and in the active sites vector
            std::vector<int> poss_rate_indices; // vector storing inidces in active rates vector for processes ocurring at this site
            std::vector<int> neighbors;
            atom(int index_l, int index_v, nanocrystal&);
            atom();
    };
    private:
        // Some private constants that determine structure of lattice (for more general applications these matrices should be read from input file)
        
        const int zlat = 12;
        const int zlat2 = 54;
        // Three vectors set geometry of FCC lattice
        const std::vector<double> xtl_1_xyz = {0, 1, 1}; // Distance is measured in units of sqrt(2) * radius of metal atom or 1/sqrt(2) * lattice constant
        const std::vector<double> xtl_2_xyz = {1, 0, 1};
        const std::vector<double> xtl_3_xyz = {1, 1, 0};

        // A matrix gives neighbors relative to a single site
        const std::vector<std::vector<int> > order = { {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1},
        
        {1, -1, 0}, {-1, 1, 0}, {0, 1, -1}, {0, -1, 1}, {-1, 0, 1}, {1, 0, -1} };

        // And this one gives the neighbors of neighbors
        const std::vector<std::vector<int> > order2 = { {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1},
        
        {1, -1, 0}, {-1, 1, 0}, {0, 1, -1}, {0, -1, 1}, {-1, 0, 1}, {1, 0, -1},
        
        {2, 0, 0}, {-2, 0, 0}, {0, 2, 0}, {0, -2, 0}, {0, 0, 2}, {0, 0, -2},
        
        {1, 1, 0}, {-1, -1, 0}, {1, 0, 1}, {-1, 0, -1}, {0, 1, 1}, {0, -1, -1},

        {2, -2, 0}, {-2, 2, 0}, {2, 0, -2}, {-2, 0, 2}, {0, 2, -2}, {0, -2, 2},

        {2, -1, 0}, {-2, 1, 0}, {2, 0, -1}, {-2, 0, 1}, {0, 2, -1}, {0, -2, 1}, {-1, 2, 0}, {1, -2, 0}, {0, -1, 2}, {0, 1, -2}, {-1, 0, 2}, {1, 0, -2},
        
        {1, -1, 1}, {-1, 1, -1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, -1}, {-1, 1, 1},

        {2, -1, -1}, {-2, 1, 1}, {-1, 2, -1}, {1, -2, 1}, {-1, -1, 2}, {1, 1, -2} };

        const double max_time = INFINITY;
        const double sqrt_2 = sqrt(2.0);
    
    public:
        // Parameters that define interactions and rates in model
        double T, mu, k_D, epsilon; // Energy measured in units of eps_11, time in units of rate_1 (usually unless k_1 set to zero)

        // We track particle number, energy, and time elapsed throughout any simulation
        long double time, energy;
        int num_atoms;
        long long unsigned int count_steps = 0;
        double initial_size;

        // Lattice dimensions
        int max_size;
        int max_est;
        int max_size_sq;
        int N_tot;
        
        
        //*** Simulation Tools ***


        // Where to store data files
        std::string filebase;
        std::string direc;
        int trial_number;
        std::ofstream config_stream;
        std::ofstream data_stream;
        std::ofstream full_config_stream;

        // Random generator
        std::mt19937 contgen;
        std::uniform_real_distribution<double> contdis;

        // Containers for keeping track of currently active processes and associated rates
        std::array<int, 3> poss_rates;
        std::array<std::vector<double>, 3> partial_sums;
        std::array<std::vector<int>, 3 >sums_order;
        std::array<std::vector<long double>, 3> pre_computed_rates; // Pre-computed Rates of adsorption and desorption and diffusion
        std::array<std::vector<int>, 3> rate_counts;
        std::array<std::vector<std::vector<std::pair<int, int> > >, 3> active_rate_indices;
        std::array<std::vector<std::vector<int> >, 3> open_slots_active_rate_indices;
        std::array<std::vector<int>, 3> num_open_slots_active_rate_indices;
        
        std::vector<atom> active_sites;
        std::vector<int> open_slots_active_sites;
        int num_open_slots_active_sites;
        std::vector<int> vector_indices;

        // Containers for keeping track of the lattice composition itself
        std::vector<bool> sites;
        std::vector<char> neighbs;

        nanocrystal(std::function<bool(int, int, double)> f, double r, double temp, double chem_pot_1, double e_scale, double diff_quick, int r_max, std::string name, std::string dir, int trial_n);
        nanocrystal(std::string);

        std::vector<int> get_nearest_neighbors(int);
        std::vector<int> get_next_nearest_neighbors(int);
        void pre_compute_rates();
        void configure_surface();
        void add_site_to_new_surface(int, int);

        void add_active_site(int);
        void add_pre_computed_rates(int);
        int add_rate(int, int, int, int);
        
        void remove_active_site(int);
        void remove_pre_computed_rates(int);
        void remove_rate(int, int, int);

        void kmc_run(unsigned long long, int, int, int, bool);
        void prepare_run(bool, int);
        void kmc_step();
        void execute_move(int, int);
        void update_neighbors(int index, int new_val);

        void output_data();
        void print_config();
        void dump_full_config();

        bool brute_force_check();
        ~nanocrystal(){};
};

nanocrystal::nanocrystal(std::function<bool(int, int, double)> f, double r, double temp, double chem_pot, double e_scale, double diff_quick, int r_max, std::string name, std::string dir, int trial_n)
{
    trial_number = trial_n;
    filebase = name;
    direc = dir;
    T = temp; // Assumes temperature, chemical potential, and epsilon are all in the same units
    mu = chem_pot;
    epsilon = e_scale;
    
    k_D = diff_quick; // Assumes this kinetic factor is relative to the adsorption/desorption prefactor
    
    initial_size = r;
    max_size = r_max;

    max_size_sq = max_size * max_size;
    N_tot = max_size * max_size_sq;

    time = 0.0;
    count_steps = 0;

    num_atoms = 0;
    sites.reserve(N_tot);
    neighbs.reserve(N_tot);

    int lattice_index = 0;
    for (std::size_t i = 0; i < max_size; i++)
    {
        for (std::size_t j = 0; j < max_size; j++)
        {
            for (std::size_t k = 0; k < max_size; k++)
            {
                if (f(lattice_index, max_size, initial_size))
                {
                    sites.push_back(1);
                    num_atoms += 1;
                }
                else
                {
                    sites.push_back(0);
                }
                lattice_index++;
            }
        }
    }
    energy = 0.0;
    configure_surface();
}

nanocrystal::nanocrystal(std::string infilename) // For loading just a minimal config
{
    std::ifstream load_file(infilename, std::ios::in | std::ios::binary);
    if (load_file.is_open())
    {
        int double_size = sizeof(double);
        int ldouble_size = sizeof(long double);
        load_file.read(reinterpret_cast<char *>(&initial_size), double_size);
        load_file.read(reinterpret_cast<char *>(&T), double_size);
        load_file.read(reinterpret_cast<char *>(&mu), double_size);
        load_file.read(reinterpret_cast<char *>(&epsilon), double_size);
        load_file.read(reinterpret_cast<char *>(&k_D), double_size);
        load_file.read(reinterpret_cast<char *>(&time), ldouble_size);
    
        int int_size = sizeof(int);
        load_file.read(reinterpret_cast<char *>(&max_size), int_size);
        load_file.read(reinterpret_cast<char *>(&num_atoms), int_size);
        load_file.read(reinterpret_cast<char *>(&count_steps), sizeof(long long unsigned int));
        
        max_size_sq = max_size * max_size;
        N_tot = max_size * max_size_sq;
        
        sites.resize(N_tot);
        neighbs.reserve(N_tot);

        char temp_val;
        int temp_index;
        int char_size = sizeof(char);
        for (int i = 0; i < num_atoms; i++)
        {
            load_file.read(reinterpret_cast<char *>(&temp_val), char_size);
            load_file.read(reinterpret_cast<char *>(&temp_index), int_size);
            sites[temp_index] = temp_val;
        }
        load_file.close();
    }
    else
    {
        throw std::invalid_argument("Can't load dump file\n");
    }
    configure_surface();
}

nanocrystal::atom::atom(int index_l, int index_v, nanocrystal& bc)
{
    lattice_index = index_l;
    vector_index = index_v;
    poss_rate_indices.resize(bc.zlat + 1);
    neighbors = bc.get_next_nearest_neighbors(lattice_index);
}

nanocrystal::atom::atom()
{
    lattice_index = INT32_MAX;
    vector_index = 0;
}

std::vector<int> nanocrystal::get_nearest_neighbors(int lattice_index)
{
    std::array<int, 3> ijk = get_ijk(lattice_index, max_size);
    std::vector<int> neighbors(zlat);
    for (int l = 0; l < zlat; l++)
    {
        int ni = nnmod(ijk[0] + order[l][0], max_size);
        int nj = nnmod(ijk[1] + order[l][1], max_size);
        int nk = nnmod(ijk[2] + order[l][2], max_size);
        neighbors[l] = ni * max_size_sq + nj * max_size + nk;
    }
    return neighbors;
}

std::vector<int> nanocrystal::get_next_nearest_neighbors(int lattice_index)
{
    std::array<int, 3> ijk = get_ijk(lattice_index, max_size);
    std::vector<int> neighbors(zlat2);
    for (int l = 0; l < zlat2; l++)
    {
        int ni = nnmod(ijk[0] + order2[l][0], max_size);
        int nj = nnmod(ijk[1] + order2[l][1], max_size);
        int nk = nnmod(ijk[2] + order2[l][2], max_size);
        neighbors[l] = ni * max_size_sq + nj * max_size + nk;
    }
    return neighbors;
}

void nanocrystal::pre_compute_rates()
{
    // Pre computes all possible rates for adsorption, desorption and diffusion processes
    for (int i = 0; i < 3; i++)
    {
        poss_rates[i] = 0;
        pre_computed_rates[i].clear();
    }
    for (int num_neighb = 0; num_neighb < zlat; num_neighb++)
    {
        pre_computed_rates[0].push_back( exp( (epsilon * num_neighb + mu) / (2 * T) ) );
        poss_rates[0] += 1;
        pre_computed_rates[1].push_back( exp( -(epsilon * num_neighb + mu) / (2 * T) ) );
        poss_rates[1] += 1;
    }
    for (int diff_neighbs = -8; diff_neighbs <= 8; diff_neighbs++)
    {
        pre_computed_rates[2].push_back( k_D * exp( (diff_neighbs - 1) * epsilon / (2 * T) ) );
        poss_rates[2] += 1;
    }
}

void nanocrystal::configure_surface()
{
    // Determine rates for all processes
    pre_compute_rates();
    for (int i = 0; i < 3; i++)
    {
        partial_sums[i].resize(poss_rates[i]);
        sums_order[i].resize(poss_rates[i]);
        rate_counts[i].resize(poss_rates[i]);
        active_rate_indices[i].resize(poss_rates[i]);
        open_slots_active_rate_indices[i].resize(poss_rates[i]);
        num_open_slots_active_rate_indices[i].resize(poss_rates[i]);
        for (int j = 0; j < poss_rates[i]; j++)
        {
            active_rate_indices[i][j].push_back({0, zlat + 1});
        }
    }

    active_sites.push_back(atom());
    vector_indices.resize(N_tot);

    int lind = 0;
    for (int i = 0; i < max_size; i++ )
    {
        for (int j = 0; j < max_size; j++ )
        {
            for (int k = 0; k < max_size; k++ )
            {
                int val = sites[lind];
                add_site_to_new_surface(lind, val);
                lind++;
            }
        }
    }

    for (int i = 0; i < N_tot; i++)
    {
        if ((sites[i] > 0 && neighbs[i] < zlat) || (sites[i] == 0 && neighbs[i] > 0 ))
        {
            add_active_site(i);
        }
    }
}

void nanocrystal::add_site_to_new_surface(int index, int val)
{
    std::vector<int> current_neighbors = get_nearest_neighbors(index);
    int num_neighb = 0;
    for (int m = 0; m < zlat; m++)
    {
        int n_index = current_neighbors[m];
        num_neighb += sites[n_index];
    }

    neighbs.push_back(num_neighb);
    
    if (val == 1)
    {
        energy -= 0.5 * neighbs[index] * epsilon;
    }
}

void nanocrystal::add_pre_computed_rates(int vector_index)
{
    int lat_index = active_sites[vector_index].lattice_index;
    int val = sites[lat_index];
    int num_neighb = neighbs[lat_index];

    if (val == 0)
    {
        if (num_neighb > 0 && num_neighb < zlat)
        {
            active_sites[vector_index].poss_rate_indices[0] = add_rate(num_neighb, vector_index, 0, 0);
        }
    }
    else
    {
        active_sites[vector_index].poss_rate_indices[zlat] = add_rate(num_neighb, vector_index, zlat, 1);
        for (int l = 0; l < zlat; l++)
        {
            int n_index = active_sites[vector_index].neighbors[l];
            if (sites[n_index] == 0)
            {
                active_sites[vector_index].poss_rate_indices[l] = add_rate( (neighbs[n_index] - num_neighb + 8), vector_index, l, 2);
            }
            else
            {
                active_sites[vector_index].poss_rate_indices[l] = 0;
            }
        }
    }
}

int nanocrystal::add_rate(int rate_num, int vector_index, int rate_label, int rate_category)
{
    int new_index;
    int num_slots = num_open_slots_active_rate_indices[rate_category][rate_num];
    if (rate_category == 0 && !(rate_num > 0 && rate_num < zlat))
    {
        std::cout << "Pause to debug\n";
    }
    rate_counts[rate_category][rate_num]++;

    if (num_slots > 0)
    {
        new_index = open_slots_active_rate_indices[rate_category][rate_num][ (num_slots - 1) ];
        open_slots_active_rate_indices[rate_category][rate_num][ (num_slots-1) ] = 0;
        num_open_slots_active_rate_indices[rate_category][rate_num]--;
        active_rate_indices[rate_category][rate_num][new_index] = {vector_index, rate_label};
    }
    else
    {
        new_index = active_rate_indices[rate_category][rate_num].size();
        active_rate_indices[rate_category][rate_num].push_back({vector_index, rate_label});
    }
    return new_index;
}

void nanocrystal::add_active_site(int site_index)
{
    int new_atom_index;
    if (num_open_slots_active_sites > 0)
    {
        new_atom_index = open_slots_active_sites[(num_open_slots_active_sites - 1)];
        open_slots_active_sites[(num_open_slots_active_sites - 1)] = 0;
        num_open_slots_active_sites--;
        active_sites[new_atom_index].lattice_index = site_index;
        active_sites[new_atom_index].vector_index = new_atom_index;
        active_sites[new_atom_index].neighbors = get_next_nearest_neighbors(site_index);
    }
    else
    {
        new_atom_index = active_sites.size();
        active_sites.push_back(atom(site_index, new_atom_index, *this));
    }
    vector_indices[site_index] = new_atom_index;
    add_pre_computed_rates(new_atom_index);
}

void nanocrystal::execute_move(int vector_index, int rate_label)
{
    int lattice_index = active_sites[vector_index].lattice_index;
    int old_val = sites[lattice_index];
    int num_neighbs = neighbs[lattice_index];
    if (!old_val)
    {
        remove_pre_computed_rates(vector_index);
        energy -= num_neighbs * epsilon;
        num_atoms++;
        sites[lattice_index] = 1;
        update_neighbors(lattice_index, 1);
        add_pre_computed_rates(vector_index);
    }

    else if (rate_label == zlat)
    {
        remove_pre_computed_rates(vector_index);
        energy += num_neighbs * epsilon;
        num_atoms--;
        sites[lattice_index] = 0;
        update_neighbors(lattice_index, 0);
        add_pre_computed_rates(vector_index);
    }

    else
    {
        energy += num_neighbs * epsilon;
        int n_lattice_index = active_sites[vector_index].neighbors[rate_label];
        remove_pre_computed_rates(vector_index);
        sites[lattice_index] = 0;
        update_neighbors(lattice_index, 0);
        
        int n_vector_index = vector_indices[n_lattice_index];
        if (n_vector_index > 0) // If the site to which the atom diffuses has no neighbors, treat as a desorption
        {
            remove_pre_computed_rates(n_vector_index);
            sites[n_lattice_index] = 1;
            update_neighbors(n_lattice_index, 1);
            add_pre_computed_rates(n_vector_index);
            energy -= neighbs[n_lattice_index] * epsilon;
        }
    }
}

void nanocrystal::prepare_run(bool no_config, int seed)
{
    std::stringstream trial_prefix;
    trial_prefix << direc << "t" << trial_number;
   
    std::string full_config_filename = trial_prefix.str() + "full" + filebase + ".txt";
    full_config_stream.open( full_config_filename, std::ios::out | std::ios::binary );

    if (!no_config)
    {
        config_stream << std::fixed <<std::setprecision(16);
        std::string config_filename = trial_prefix.str() + "c" + filebase + ".txt";
        config_stream.open( config_filename, std::ofstream::out | std::ofstream::app );
    }

    std::string data_filename = trial_prefix.str() + "d" + filebase + ".txt";
    data_stream.open( data_filename, std::ofstream::out | std::ofstream::app );
    data_stream << std::scientific <<std::setprecision(16);

    if (!no_config)
    {
        config_stream <<  "T mu eps k_D max_size\n";
        config_stream << T << " " << mu << " " << epsilon << " " << k_D << " " << max_size << "\n";
        config_stream << "val index neighbs neighb_2 " << seed << "\n";
    }

    data_stream <<  "T mu eps k_D max_size\n";
    data_stream << T << " " << mu << " " << epsilon << " " << k_D  << " " << max_size << "\n";
    data_stream <<  "Steps t N E " << seed << "\n";
    
    int double_size = sizeof(double);
    full_config_stream.write(reinterpret_cast<const char*>(&initial_size), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&T), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&mu), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&epsilon), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&k_D), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&max_size), sizeof(int));
}

void nanocrystal::remove_rate(int rate_num, int vector_index, int rate_category)
{
    rate_counts[rate_category][rate_num]--;
    active_rate_indices[rate_category][rate_num][vector_index] = {0, zlat + 1};
    
    num_open_slots_active_rate_indices[rate_category][rate_num]++;
    if (num_open_slots_active_rate_indices[rate_category][rate_num] > open_slots_active_rate_indices[rate_category][rate_num].size())
    {
        open_slots_active_rate_indices[rate_category][rate_num].push_back(vector_index);
    }
    else
    {
        open_slots_active_rate_indices[rate_category][rate_num][num_open_slots_active_rate_indices[rate_category][rate_num] - 1] = vector_index;
    }
}

void nanocrystal::remove_pre_computed_rates(int vector_index)
{
    int lat_index = active_sites[vector_index].lattice_index;
    int num_neighb = neighbs[lat_index];
    int val = sites[lat_index];
    int current_index;

    if (val == 0)
    {
        // Remove adsorption rate
        current_index = active_sites[vector_index].poss_rate_indices[0];
        if (current_index > 0)
        {
            remove_rate(num_neighb, current_index, 0);
            active_sites[vector_index].poss_rate_indices[0] = 0;
        }
    }
    else
    {
        // Remove desorption rate
        current_index = active_sites[vector_index].poss_rate_indices[zlat];
        if (current_index > 0)
        {
            remove_rate(num_neighb, current_index, 1);
            active_sites[vector_index].poss_rate_indices[12] = 0;
        }
        // Remove all diffusion rates
        for (int l = 0; l < zlat; l++)
        {
            current_index = active_sites[vector_index].poss_rate_indices[l];
            if (current_index > 0)
            {
                int n_index = active_sites[vector_index].neighbors[l];
                remove_rate( (neighbs[n_index] - num_neighb + 8), current_index, 2);
                active_sites[vector_index].poss_rate_indices[l] = 0;
            }
        }
    }
}

void nanocrystal::remove_active_site(int vector_index)
{
    num_open_slots_active_sites++;
    if (num_open_slots_active_sites > open_slots_active_sites.size())
    {
        open_slots_active_sites.push_back(vector_index);
    }
    else
    {
        open_slots_active_sites[(num_open_slots_active_sites - 1)] = vector_index;
    }
    
    vector_indices[active_sites[vector_index].lattice_index] = 0;
    active_sites[vector_index].lattice_index = N_tot;
    active_sites[vector_index].vector_index = 0;
}

void nanocrystal::update_neighbors(int lattice_index, int new_val)
{
    int vector_index = vector_indices[lattice_index];
    for (int l = 0; l < zlat2; l++)
    {
        int n_lattice_index = active_sites[vector_index].neighbors[l];
        int n_vector_index = vector_indices[n_lattice_index];
        if (n_vector_index > 0)
        {
            remove_pre_computed_rates(n_vector_index);
        }
    }
    for (int l = 0; l < zlat; l++)
    {
        int n_lattice_index = active_sites[vector_index].neighbors[l];
        int n_vector_index = vector_indices[n_lattice_index];
        if (new_val)
        {
            neighbs[n_lattice_index] += 1;
        }
        else
        {
            neighbs[n_lattice_index] -= 1;
        }
    }
    for (int l = 0; l < zlat2; l++)
    {
        int n_lattice_index = active_sites[vector_index].neighbors[l];
        int n_vector_index = vector_indices[n_lattice_index];
        int n_val = sites[n_lattice_index];
        int num_neighb = neighbs[n_lattice_index];
        if ( (n_val > 0 && num_neighb < zlat) || (n_val == 0 && num_neighb > 0 ) )
        {
            if (n_vector_index == 0)
            {
                add_active_site(n_lattice_index);
            }
            else
            {
                add_pre_computed_rates(n_vector_index);
            }
        }
        else if (n_vector_index != 0)
        {
            remove_active_site(n_vector_index);
        }
    }
    if (!new_val)
    {
        // If an isolated atom was created, remove it
        for (int l = 0; l < zlat; l++)
        {
            int n_lattice_index = active_sites[vector_index].neighbors[l];
            if (sites[n_lattice_index] && !neighbs[n_lattice_index])
            {
                execute_move(vector_indices[n_lattice_index], zlat);
                remove_active_site(vector_indices[n_lattice_index]);
            }
        }
    }
}

void nanocrystal::kmc_step()
{
    // Create the partial sums list of the rates
    long double total_rate = 0.0;
    long double new_rate;
    std::vector<long double> process_rates(3);
    std::vector<int> count(3, 0);
    for (int i = 0; i < 3; i++)
    {
        for (int q = 0; q < poss_rates[i]; q++)
        {
            new_rate = total_rate + rate_counts[i][q] * pre_computed_rates[i][q];
            if (new_rate > total_rate)
            {
                total_rate = new_rate;
                partial_sums[i][count[i]] = total_rate;
                sums_order[i][count[i]] = q;
                count[i]++;
            }
        }
        process_rates[i] = total_rate;
    }
    // Choose the random process to occur in proportion to its rate
    long double cap = contdis(contgen) * total_rate;
    for (int i = 0; i < 3; i++)
    {
        if (cap < process_rates[i])
        {
            auto it = std::upper_bound(partial_sums[i].begin(), partial_sums[i].begin() + count[i], cap);
            int process_index = int(it - partial_sums[i].begin());
            int process = sums_order[i][process_index];
            
            // Choose the lattice position where that process occurs with uniform probability for all equivalent lattice sites
            int trial_index;
            int vector_index = 0;
            int process_label = 0;
            while (vector_index == 0)
            {
                trial_index = int(1 + contdis(contgen) * ((rate_counts[i][process] + num_open_slots_active_rate_indices[i][process]) - 1) + 0.5);
                vector_index = active_rate_indices[i][process][trial_index].first;
                process_label = active_rate_indices[i][process][trial_index].second;
            }
            // Change the lattice in accordance with the process
            execute_move(vector_index, process_label);
            time -= log( contdis(contgen) ) / total_rate;
            count_steps++;
            return;
        }
    }
}

void nanocrystal::kmc_run(unsigned long long num_steps, int config_rec, int data_rec, int seed, bool no_config)
{
    prepare_run(no_config, seed);
    contgen.seed(seed);
    contdis.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
    
    int config_steps = num_steps / config_rec;
    int data_steps = num_steps / data_rec;
    if (!no_config)
    {
        print_config();
    }

    count_steps = 0;
    while (count_steps < num_steps && num_atoms > 0)
    {
        kmc_step();
        
        // Uncomment below for debugging
        /*
        if (!brute_force_check())
        {
            std::cerr << "Fatal error\n";
        }
        */
        
        if ( count_steps % data_steps == 0 )
        {
            data_stream << count_steps << " " << time << " " << num_atoms << " " << energy << "\n";
            std::cout << count_steps << " " << time << " " << num_atoms << " " << energy << "\n";
        }
        if ( count_steps % config_steps == 0 && !no_config)
        {
            print_config();
        }
    }

    data_stream << count_steps << " " << time << " " << num_atoms << " " << energy << "\n";
    data_stream.close();
    
    if (!no_config)
    {
        config_stream.close();
    }

    dump_full_config();
    full_config_stream.close();
}

void nanocrystal::print_config()
{
    int val, num_neighbors;
    int lattice_index;
    int i, j, k;
    double x, y, z;

    for (int r = 1; r <= active_sites.size(); r++)
    {
        lattice_index = active_sites[r].lattice_index;
        
        if (lattice_index < N_tot)
        {
            val = sites[lattice_index];
            if (val == 1)
            {
                i = lattice_index / max_size_sq;
                j = (lattice_index % max_size_sq) / max_size;
                k = lattice_index % max_size;
                
                x = (1/sqrt_2) * (xtl_1_xyz[0] * i + xtl_2_xyz[0] * j + xtl_3_xyz[0] * k - (max_size / 2.0) * (xtl_1_xyz[0] + xtl_2_xyz[0] + xtl_3_xyz[0]));
                y = (1/sqrt_2) * (xtl_1_xyz[1] * i + xtl_2_xyz[1] * j + xtl_3_xyz[1] * k - (max_size / 2.0) * (xtl_1_xyz[1] + xtl_2_xyz[1] + xtl_3_xyz[1]));
                z = (1/sqrt_2) * (xtl_1_xyz[2] * i + xtl_2_xyz[2] * j + xtl_3_xyz[2] * k - (max_size / 2.0) * (xtl_1_xyz[2] + xtl_2_xyz[2] + xtl_3_xyz[2]));
            
                num_neighbors = neighbs[lattice_index];
                config_stream << val << " " << x << " " << y << " " << z << " " << num_neighbors << "\n";
            }
        }
    }
    config_stream << num_atoms << " " << count_steps << " " << time << " " << energy << "\n";
}

void nanocrystal::dump_full_config()
{
    auto size_int = sizeof(int);
    if (full_config_stream.is_open())
    {
        full_config_stream.write(reinterpret_cast<const char *>(&num_atoms), size_int);
        full_config_stream.write(reinterpret_cast<const char *>(&count_steps), sizeof(long long unsigned int));

        char temp_val;
        int char_size = sizeof(char);
        for (int i = 0; i < N_tot; i++)
        {
            temp_val = sites[i];
            if (temp_val)
            {
                full_config_stream.write(reinterpret_cast<const char *>(&temp_val), char_size);
                full_config_stream.write(reinterpret_cast<const char *>(&i), size_int);
            }
        }
        full_config_stream.close();
    }
}

bool nanocrystal::brute_force_check()
{
    // Brute force check that all rates are computed correctly (used only for debugging)
    int vind, lind, num_neighb, val;
    long double k_tot_copy = 0.0;
    
    std::array<std::vector<int>, 3> rate_count_copy;
    for (int i = 0; i < 3; i++)
    {
        rate_count_copy[i].resize(poss_rates[i]);
    }

    lind = 0;
    for (int i = 0; i < max_size; i++)
    {
        for (int j = 0; j < max_size; j++)
        {
            for (int k = 0; k < max_size; k++)
            {
                val = sites[lind];
                num_neighb = neighbs[lind];
        
                int cnum_neighb = 0;
                std::vector<int> cneighbors = get_nearest_neighbors(lind);
                for (int l = 0; l < zlat; l++)
                {
                    cnum_neighb += sites[cneighbors[l]];
                }
                if (cnum_neighb != num_neighb)
                {
                    std::cerr << "incorrect number of neighbors.\n";
                    return false;
                }
            
                if ( (val > 0 && num_neighb < zlat) || (val == 0 && num_neighb > 0) )
                {
                    vind = vector_indices[lind];
                    if (!vind)
                    {
                        std::cerr << "Site " << lind << " should be active and isn't\n";
                        return false;
                    }

                    if (!val)
                    {
                        if ( num_neighb > 0 && num_neighb < zlat )
                        {
                            int rate_num = num_neighb;
                            k_tot_copy += pre_computed_rates[0][rate_num];
                            rate_count_copy[0][rate_num] += 1;
                            
                            if (abs(exp((epsilon * num_neighb + mu) / (2*T)) - pre_computed_rates[0][rate_num]) > 1e-10)
                            {
                                std::cerr << "rate is incorrect!\n";
                            }
                            if (active_sites[vind].poss_rate_indices[0] == 0)
                            {
                                std::cerr << "This rate should be active for site " << lind << "\n";
                                return false;
                            }
                        }
                        else if (active_sites[vind].poss_rate_indices[0])
                        {
                            std::cerr << "This rate shouldn't be active for site " << lind << "\n";
                        }
                    }

                    else
                    {
                        if ( num_neighb > 0 && num_neighb < zlat )
                        {
                            int rate_num = num_neighb;
                            k_tot_copy += pre_computed_rates[1][rate_num];
                            rate_count_copy[1][rate_num] += 1;
                            if (abs(exp(-(epsilon * num_neighb + mu) / (2*T)) - pre_computed_rates[1][rate_num]) > 1e-10)
                            {
                                std::cerr << "rate is incorrect!\n";
                            }
                            if (active_sites[vind].poss_rate_indices[zlat] == 0)
                            {
                                std::cerr << "This rate should be active for site " << lind << "\n";
                                return false;
                            }
                        }
                        
                        else if (active_sites[vind].poss_rate_indices[zlat] > 0)
                        {
                            std::cerr << "This rate shouldn't be active for site " << lind << "\n";
                            return false;
                        }

                        for (int l = 0; l < zlat; l++)
                        {
                            int n_lat = active_sites[vind].neighbors[l];
                            if (!sites[n_lat])
                            {
                                int nnum_neighb = neighbs[n_lat];
                                int rate_num = nnum_neighb - num_neighb + 8;
                                k_tot_copy += pre_computed_rates[2][rate_num];
                                rate_count_copy[2][rate_num] += 1;
                                if (abs(k_D * exp((epsilon * (nnum_neighb - 1 - num_neighb) ) / (2*T)) - pre_computed_rates[2][rate_num]) > 1e-10)
                                {
                                    std::cerr << "rate is incorrect!\n";
                                }

                                if (active_sites[vind].poss_rate_indices[l] == 0)
                                {
                                    std::cerr << "This rate should be active for site " << lind << "\n";
                                    return false;
                                }
                            }
                            
                            else if (active_sites[vind].poss_rate_indices[l])
                            {
                                std::cerr << "This rate should not be active for site " << lind << "\n";
                                return false;
                            }
                        }
                    }
                }
                lind += 1;
            }
        }
    }
    double new_rate, curr_rate;
    curr_rate = 0;
    for (int i = 0; i < 3; i++)
    {
        for (int q = 0; q < poss_rates[i]; q++)
        {
            if (rate_counts[i][q] != rate_count_copy[i][q])
            {
                std::cerr << "Rates recorded incorrectly\n";
                return false;
            }
            new_rate = curr_rate + rate_counts[i][q] * pre_computed_rates[i][q];
            
            if (new_rate > curr_rate)
            {
                curr_rate = new_rate;
            }
        }
    }
    
    if (abs(curr_rate - k_tot_copy)/curr_rate > 1e-10)
    {
        std::cerr << abs(curr_rate - k_tot_copy)/curr_rate << " The rates have been calculated incorrectly\n";
        return false;
    }
    return true;
}

#endif /* BCSDbep_h */
