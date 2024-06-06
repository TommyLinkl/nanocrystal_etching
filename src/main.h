///
//  BCSDbep.h
//  Created by Sam Oaks-Leaf. 
//  A header describing a nanocrystal class designed for reshaping in the BEP framework for surface diffusion.

#ifndef BCSDirbligs_h
#define BCSDirbligs_h

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

class blocked_crystal
{
    class atom
    {
        public:
            unsigned int lattice_index, vector_index, surface_index;
            std::array<unsigned int, 13> poss_rate_indices;
            atom(unsigned int index_l, unsigned int index_v);
            atom();
    };
    private:
        // Three vectors set geometry of FCC lattice
        std::vector<double> xtl_1_xyz = {0, 1, 1}; // Distance is measured in units of sqrt(2) * radius of metal atom or 1/sqrt(2) * lattice constant
        std::vector<double> xtl_2_xyz = {1, 0, 1};
        std::vector<double> xtl_3_xyz = {1, 1, 0};
        
        double pi = 3.1415926535897;
        double sqrt_2 = sqrt(2.0);
        int poss_rates = 1734;
        
        // A matrix gives neighbors relative to a single site
        std::vector<std::vector<int> > order = { {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}, 
        
        {1, -1, 0}, {-1, 1, 0}, {0, 1, -1}, {0, -1, 1}, {-1, 0, 1}, {1, 0, -1} };

        std::vector<std::vector<int> > order2 = { {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}, 
        
        {1, -1, 0}, {-1, 1, 0}, {0, 1, -1}, {0, -1, 1}, {-1, 0, 1}, {1, 0, -1},
        
        {2, 0, 0}, {-2, 0, 0}, {0, 2, 0}, {0, -2, 0}, {0, 0, 2}, {0, 0, -2},
        
        {1, 1, 0}, {-1, -1, 0}, {1, 0, 1}, {-1, 0, -1}, {0, 1, 1}, {0, -1, -1},

        {2, -2, 0}, {-2, 2, 0}, {2, 0, -2}, {-2, 0, 2}, {0, 2, -2}, {0, -2, 2},

        {2, -1, 0}, {-2, 1, 0}, {2, 0, -1}, {-2, 0, 1}, {0, 2, -1}, {0, -2, 1}, {-1, 2, 0}, {1, -2, 0}, {0, -1, 2}, {0, 1, -2}, {-1, 0, 2}, {1, 0, -2}, 
        
        {1, -1, 1}, {-1, 1, -1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, -1}, {-1, 1, 1},

        {2, -1, -1}, {-2, 1, 1}, {-1, 2, -1}, {1, -2, 1}, {-1, -1, 2}, {1, 1, -2} };
        
        const double eps_11_conversion = 0.225; // Units of eV
        const double kb_to_ev = 8.6173e-5;
        double tol = 1e-9;
        //int resizing_factor = 10;
        int resizing_factor = 2.0;
        double increase_factor = 2.0;
        int lig_step = 50; 
        double max_time = INFINITY;
       
    public:

        // Parameters that define interactions and rates in model
        
        double T, mu, k_1, k_2, k_1D, k_2D, eps_11, eps_12, eps_22; // Energy measured in units of eps_11, time in units of rate_1 (usually unless k_1 set to zero)
        bool operational, count_h;
        // We track particle number, energy, total rate, and time elapsed throughout any simulation
        long double k_tot, time;
        double energy;
        unsigned int N_1, N_2;
        long long unsigned int count_steps = 0;
        double initial_size;

        //track the "width" of the nanoparticle in the 100 and 111 directions
        std::array<std::vector<int>, 3> h100_counts; //counts the number of atoms in each {100} layer
        std::array<std::vector<int>, 4> h111_counts;
        std::array<std::array<int, 2>, 3> h100; // mins and max layers with at least 1 atom present
        std::array<std::array<int, 2>, 4> h111;
        //int min_poss_111; // The minimum possible {111} plane. Will be negative, while the min possible {100} plane will always be 0.

        bool must_update_max_size = false;
        bool stationary = false;
        
        std::array<std::array<long double, 3>, 3> smn;
        std::array<long double, 3> com;
        long double asphericity, asphericity_i;
        long double shape_anisotropy, shape_anisotropy_i;

        unsigned int max_size;
        unsigned int max_est;
        unsigned int max_size_sq;
        unsigned int N_tot;
        
        int Ns_2, Ns, num_open_slots_surface_sites;
        std::vector<unsigned int> surface_sites;
        std::vector<unsigned int> copy_surface_sites;
        std::vector<unsigned int> open_slots_surface_sites;
        
        // Simulation Tools
        std::string filebase;
        std::string direc;
        int trial_number;

        std::mt19937 contgen;
        std::uniform_real_distribution<double> contdis;
        std::mt19937 liggen;
        std::uniform_int_distribution<int> ligdis;

        std::ofstream config_stream;
        std::ofstream inertia_stream;
        std::ofstream data_stream;
        std::ofstream full_config_stream;
        
        std::vector<double> partial_sums;
        std::vector<int> sums_order;
        std::vector<long double> all_rates; // Pre-computed Rates of evaporation and deposition of each species
        std::vector<unsigned int> rate_counts;

        std::array<std::vector<unsigned int>, 1734> active_rate_indices;
        std::array<std::vector<char>, 1734> active_rate_labels;
        std::array<std::vector<unsigned int>, 1734> open_slots_active_rate_indices;
        std::vector<unsigned int> num_open_slots_active_rate_indices;
        std::vector<unsigned int> max_size_active_rate_indices;
        
        std::vector<atom> active_sites;
        std::vector<unsigned int> open_slots_active_sites;
        unsigned int num_open_slots_active_sites;
        unsigned int max_size_active_sites;
        unsigned int size_active_sites;

        std::vector<char> sites;
        std::vector<char> neighb_1;
        std::vector<char> neighb_2;
        std::vector<unsigned int> vector_indices;

        double rho, rho_star;
        bool fixed_density;
        
        //Some performance checking tools
        int pv, nv, pn1, pn2;
    
        blocked_crystal(std::string, std::string, std::string, int);

        blocked_crystal(std::string, std::string, std::string, int, int);

        blocked_crystal(std::function<bool(int, int, double)> f, double r, double temp, double chem_pot_1, double e_scale_12, double e_scale_22, double quick_2, double squick_1, double squick_2, int r_max, std::string name, std::string dir, int trial_n);

        blocked_crystal(std::string);

        blocked_crystal();

        void add_active_site(unsigned int);

        void remove_active_site(unsigned int);

        void add_surface_site(unsigned int);

        void remove_surface_site(unsigned int, unsigned int);
        
        bool execute_move(unsigned int, int, int);

        void configure_surface(bool);

        void init_ligs(double, int);

        void add_ligs();
        
        void add_site_to_new_surface(unsigned int, int);

        void update_neighbors(unsigned int index, int old_val, int new_val);
        
        void kmc_step();
        
        void kmc_run_time(double total_time, int config_rec, int data_rec, unsigned int seed, unsigned int, bool);

        void kmc_run_mass_down(unsigned int N_1_end, int N_1_end_high, int config_rec, int data_rec, unsigned int seed, unsigned int, bool);

        void kmc_run_mass_up(unsigned int N_1_end, int N_1_end_low, int config_rec, int data_rec, unsigned int seed, unsigned int,  bool);

        void kmc_run_steps(unsigned int num_steps, int config_rec, int data_rec, unsigned int seed, unsigned int, bool);

        bool compare_rates(int i, int j);

        void print_surface_vec();

        void pre_compute_rates(double);

        void write_active_sites();

        void read_active_sites(std::ifstream& in_str);

        unsigned int add_rate(int, int, int);

        void remove_rate(int, int);

        void remove_all_rates(unsigned int);

        void add_all_rates(unsigned int);

        void add_to_gyration_tensor(int index);

        void remove_from_gyration_tensor(int index);

        void integrate_order_params(double);

        void print_gyration_tensor();

        void print_config();

        void print_everything();

        //void print_full_config();

        void prepare_run(bool, int);
        
        void print_xyz(int max_particles_viewed, std::ofstream& view_stream);

        void print_minimal_config(std::string);

        void add_to_hcounts(int);

        void remove_from_hcounts(int);

        void initialize_hcounts();

        //std::array<int, 3> get_ijk(int);

        //std::array<double, 3> get_xyz(int);

        double get_bond_energy(int lattice_index);

        void truncate_corners(int, int);

        bool brute_force_check();

        ~blocked_crystal(){};
        
};

bool blocked_crystal::brute_force_check()
{
    // Brute force check that rates are computed correctly
    int vind; 
    int ind = 0.0;
    int num_1, num_2, val;
    double k_tot_copy = 0.0;
    std::vector<int> rate_count_copy(poss_rates);

    int n_i, n_j, n_k, n_val, n_lat;
    for (int i = 0; i < max_size; i++)
    {
        for (int j = 0; j < max_size; j++)
        {
            for (int k = 0; k < max_size; k++)
            {
                val = sites[ind];
                num_1 = neighb_1[ind];
                num_2 = neighb_2[ind];
        
                int cnum_1 = 0, cnum_2 = 0;

                for (int l = 0; l < 12; l++)
                {
                    n_i = nnmod(i + order[l][0], max_size);
                    n_j = nnmod(j + order[l][1], max_size);
                    n_k = nnmod(k + order[l][2], max_size);
                    
                    
                    n_lat = n_i * max_size_sq + n_j * max_size + n_k;
                    n_val = sites[n_lat];
                    if (n_val == 1)
                    {
                        cnum_1++;
                    }
                    else if (n_val == 2)
                    {
                        cnum_2++;
                    }
                }
            
                if (cnum_1 != num_1 || cnum_2 != num_2)
                {
                    std::cerr << "incorrect number of neighbors.\n";
                    return false;
                }
            
                if ( (val > 0 && (num_1 + num_2) < 12) || (val == 0 && (num_1 + num_2) > 0 ) )
                {
                    vind = vector_indices[ind];
                    if ( vind == 0)
                    {
                        std::cerr << "Site " << ind << " should be active and isn't\n";
                        return false;
                    }

                    else if (vind >= max_size_active_sites)
                    {
                        std::cerr << "Vector index for site " << ind << " exceeds available indices\n";
                        return false;
                    }
                    if (val == 0)
                    {
                        if ( (num_1 + num_2) > 0 && (num_1 + num_2) < 12 )
                        {
                            //k_tot_copy += k_1 * exp((num_1 * eps_11 + num_2 * eps_12 + mu) / (2*T));
                            k_tot_copy += all_rates[6 * (17 * num_1 + num_2)];
                            rate_count_copy[6 * (17 * num_1 + num_2)] += 1;
                            rate_count_copy[6 * (17 * num_1 + num_2) + 1] += 1;
                            if (abs(k_1 * exp((eps_11 * num_1 + eps_12 * num_2 + mu) / (2*T)) - all_rates[6 * (17 * num_1 + num_2)]) > 10e-10)
                            {
                                std::cerr << "rate is incorrect!\n";
                            }
                            if (active_sites[vind].poss_rate_indices[0] == 0)
                            {
                                std::cerr << "This rate should be active for site " << ind << "\n";
                                return false;
                            }
                            k_tot_copy += k_2 * exp((num_1 * eps_12 + num_2 * eps_22 + mu) / (2*T));
                            if (active_sites[vind].poss_rate_indices[1] == 0)
                            {
                                std::cerr << "This rate should be active for site " << ind << "\n";
                                return false;
                            }
                        }
                        else if (active_sites[vind].poss_rate_indices[1] > 0 || active_sites[vind].poss_rate_indices[0] > 0)
                        {
                            std::cerr << "This rate shouldn't be active for site " << ind << "\n";
                        }
                    }

                    else 
                    {
                        if (val == 1)
                        {
                            if ( (num_1 + num_2) > 0 && (num_1 + num_2) < 12 )
                            {
                                //k_tot_copy += k_1 * exp(-(num_1 * eps_11 + num_2 * eps_12 + mu) / (2*T));
                                k_tot_copy += all_rates[6 * (17 * num_1 + num_2) + 2];
                                rate_count_copy[6 * (17 * num_1 + num_2) + 2] += 1;
                                if (abs(k_1 * exp(-(eps_11 * num_1 + eps_12 * num_2 + mu) / (2*T)) - all_rates[6 * (17 * num_1 + num_2) + 2]) > 10e-10)
                                {
                                    std::cerr << "rate is incorrect!\n";
                                }
                                if (active_sites[vind].poss_rate_indices[12] == 0)
                                {
                                    std::cerr << "This rate should be active for site " << ind << "\n";
                                    return false;
                                }   
                            }
                            
                            else if (active_sites[vind].poss_rate_indices[12] > 0)
                            {
                                std::cerr << "This rate shouldn't be active for site " << ind << "\n";
                                return false;
                            }

                            for (int l = 0; l < 12; l++)
                            {
                                n_i = nnmod(i + order[l][0], max_size);
                                n_j = nnmod(j + order[l][1], max_size);
                                n_k = nnmod(k + order[l][2], max_size);
                            
                                n_lat = n_i * max_size_sq + n_j * max_size + n_k;
                                int nnum_1 = neighb_1[n_lat];
                                int nnum_2 = neighb_2[n_lat];
                                if (sites[n_lat] == 0)
                                {
                                    //k_tot_copy += k_1D * exp((eps_11 * (nnum_1 - 1 - num_1) + eps_12 * (nnum_2 - num_2)) / (2*T));
                                    k_tot_copy += all_rates[6 * (17 * (nnum_1 - 1 - num_1 + 8) + (nnum_2 - num_2 + 8)) + 4];
                                    rate_count_copy[6 * (17 * (nnum_1 - 1 - num_1 + 8) + (nnum_2 - num_2 + 8)) + 4] += 1;
                                    if (abs(k_1D * exp((eps_11 * (nnum_1 - 1 - num_1) + eps_12 * (nnum_2 - num_2)) / (2*T)) - all_rates[6 * (17 * (nnum_1 - 1 - num_1 + 8) + (nnum_2 - num_2 + 8)) + 4]) > 10e-10)
                                    {
                                        std::cerr << "rate is incorrect!\n";
                                    }

                                    if (active_sites[vind].poss_rate_indices[l] == 0)
                                    {
                                        std::cerr << "This rate should be active for site " << ind << "\n";
                                        return false;
                                    }
                                }
                                
                                else if (active_sites[vind].poss_rate_indices[l] != 0)
                                {
                                    std::cerr << "This rate should not be active for site " << ind << "\n";
                                    return false;
                                }
                            }
                        }

                        else if (val == 2)
                        {
                            if ( (num_1 + num_2) > 0 && (num_1 + num_2) < 12 )
                            {
                                k_tot_copy += k_2 * exp(-(num_1 * eps_11 + num_2 * eps_12 + mu) / (2*T));
                                if (active_sites[vind].poss_rate_indices[12] == 0)
                                {
                                    std::cerr << "This rate should be active for site " << ind << "\n";
                                    return false;
                                }   
                            }
                            for (int l = 0; l < 12; l++)
                            {
                                n_i = nnmod(i + order[l][0], max_size);
                                n_j = nnmod(j + order[l][1], max_size);
                                n_k = nnmod(k + order[l][2], max_size);
                                
                                n_lat = n_i * max_size_sq + n_j * max_size + n_k;
                                int nnum_1 = neighb_1[n_lat];
                                int nnum_2 = neighb_2[n_lat];

                                if (sites[n_lat] == 0)
                                {
                                    k_tot_copy += k_2D * exp((eps_12 * (nnum_1 - num_1) + eps_22 * (nnum_2 - 1 - num_2)) / (2*T));
                                    
                                    if (active_sites[vind].poss_rate_indices[l] == 0)
                                    {
                                        std::cerr << "This rate should be active for site " << ind << "\n";
                                        return false;
                                    }
                                }
                                else if (active_sites[vind].poss_rate_indices[l] != 0)
                                {
                                    std::cerr << "This rate should not be active for site " << ind << "\n";
                                    return false;
                                }
                            }
                        }

                        else
                        {
                            std::cerr << "Site " << ind << " is not a valid type\n";
                            return false;
                        }
                    }
                }
                
                ind+=1;
            }
        }
    }
    double new_rate, curr_rate;
    curr_rate = 0;
    for (int q = 0; q < 1734; q++)
    {
        if (rate_counts[q] != rate_count_copy[q])
        {
            std::cerr << "Rates recorded incorrectly\n";
            return false;
        }
        new_rate = curr_rate + rate_counts[q] * all_rates[q];
        
        if (new_rate > curr_rate)
        {
            curr_rate = new_rate;
        }
    }
    
    if (abs(curr_rate - k_tot_copy)/curr_rate > pow(10,-12))
    {
        std::cerr << abs(curr_rate - k_tot_copy)/curr_rate << " The rates have been calculated incorrectly\n";
	    return false;
    }
    return true;
}

blocked_crystal::blocked_crystal(std::function<bool(int, int, double)> f, double r, double temp, double chem_pot_1, double e_scale_12, double e_scale_22, double quick_2, double squick_1, double squick_2, int r_max, std::string name, std::string dir, int trial_n)
{
    trial_number = trial_n;
    filebase = name;
    direc = dir;
    fixed_density = false;
    T = temp * kb_to_ev / eps_11_conversion; // converts T from kelvin to units of epsilon_11
    mu = chem_pot_1;
    k_1 = 1.0;
    k_2 = quick_2;
    eps_11 = 1.0; // All energies will be expressed in units of epsilon_11 in the code, connection to physical units (eV) hidden in "eps_11_conversion"
    eps_12 = e_scale_12;
    eps_22 = e_scale_22;
    initial_size = r;
    max_size = r_max;
    k_1D = squick_1;
    k_2D = squick_2;

    max_size_sq = max_size * max_size;
    N_tot = max_size * max_size_sq;
    //min_poss_111 = int(max_size_sq / sqrt(3) + 0.1);

    time = 0.0;
    count_steps = 0;
    
    if (initial_size >= max_size)
    {
        operational = false;
        throw std::invalid_argument("Initial size too large.");
    }
    
    // pre_comp_rates array will be organized as pre_comp_rates[i][j][k] = Boltzmann Factor for 0 --> i with
    // j neighbors of type 1 and k neighbors of type 2
    all_rates.resize(poss_rates);
    pre_compute_rates(1.0);
    shape_anisotropy_i = 0.0;
    asphericity_i = 0.0;
    double center = max_size / 2.0;
    
    N_1 = 0;
    N_2 = 0;

    // Determining Initial Configuration

    std::array<double, 3> trial_site;
    // Resize hcounts
    for (int i = 0; i < 3; i++)
    {
        h100_counts[i].resize(max_size * 2);
        h100[i][0] = 2 * max_size;
    }
    for (int i = 0; i < 4; i++)
    {
        h111_counts[i].resize(3 * max_size);
        h111[i][0] = 3 * max_size;
    }
    
    sites.reserve(N_tot);
    neighb_1.reserve(N_tot);
    neighb_2.reserve(N_tot);

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
                    N_1 += 1;
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
    configure_surface(true);
    operational = true;
    count_h = false;
}

blocked_crystal::blocked_crystal(std::string infilename) // For loading just a minimal config
{

    std::ifstream load_file(infilename, std::ios::in | std::ios::binary);
    //std::cout << "opening file...\n";
    if (load_file.is_open())
    {
        
        int double_size = sizeof(double);
        int ldouble_size = sizeof(long double);
        load_file.read(reinterpret_cast<char *>(&initial_size), double_size);
        load_file.read(reinterpret_cast<char *>(&T), double_size);
        load_file.read(reinterpret_cast<char *>(&mu), double_size);
        load_file.read(reinterpret_cast<char *>(&eps_11), double_size);
        load_file.read(reinterpret_cast<char *>(&eps_12), double_size);
        load_file.read(reinterpret_cast<char *>(&eps_22), double_size);
        load_file.read(reinterpret_cast<char *>(&k_1), double_size);
        load_file.read(reinterpret_cast<char *>(&k_2), double_size);
        load_file.read(reinterpret_cast<char *>(&k_1D), double_size);
        load_file.read(reinterpret_cast<char *>(&k_2D), double_size);
        load_file.read(reinterpret_cast<char *>(&time), ldouble_size);
        load_file.read(reinterpret_cast<char *>(&asphericity_i), ldouble_size);
        load_file.read(reinterpret_cast<char *>(&shape_anisotropy_i),ldouble_size);
        
        int int_size = sizeof(int);
        int uint_size = sizeof(unsigned int);
        load_file.read(reinterpret_cast<char *>(&max_size), uint_size);
        load_file.read(reinterpret_cast<char *>(&N_1), int_size);
        load_file.read(reinterpret_cast<char *>(&N_2), int_size);
        load_file.read(reinterpret_cast<char *>(&count_steps), sizeof(long long unsigned int));
        
        max_size_sq = max_size * max_size;
        N_tot = max_size * max_size_sq;
        
        sites.resize(N_tot);
        neighb_1.reserve(N_tot);
        neighb_2.reserve(N_tot);

        all_rates.resize(poss_rates);
        pre_compute_rates(1.0);

        char temp_val; 
        unsigned int temp_index;
        int char_size = sizeof(char);
        for (int i = 0; i < (N_1 + N_2); i++)
        {
            load_file.read(reinterpret_cast<char *>(&temp_val), char_size);
            load_file.read(reinterpret_cast<char *>(&temp_index), uint_size);
            sites[temp_index] = temp_val;
        }
        load_file.close();
    }
    else
    {
        throw std::invalid_argument("Can't load dump file\n");
    }
    configure_surface(true);
    operational = true;
}

void blocked_crystal::print_minimal_config(std::string filename)
{
    std::ofstream minimal_stream(filename, std::ios::out | std::ios::binary);
    if (minimal_stream.is_open())
    {
        int double_size = sizeof(double);
        int ldouble_size = sizeof(long double);
        minimal_stream.write(reinterpret_cast<const char *>(&initial_size), double_size);
        minimal_stream.write(reinterpret_cast<const char *>(&T), double_size);
        minimal_stream.write(reinterpret_cast<const char *>(&mu), double_size);
        minimal_stream.write(reinterpret_cast<const char *>(&eps_11), double_size);
        minimal_stream.write(reinterpret_cast<const char *>(&eps_12), double_size);
        minimal_stream.write(reinterpret_cast<const char *>(&eps_22), double_size);
        minimal_stream.write(reinterpret_cast<const char *>(&k_1), double_size);
        minimal_stream.write(reinterpret_cast<const char *>(&k_2), double_size);
        minimal_stream.write(reinterpret_cast<const char *>(&k_1D), double_size);
        minimal_stream.write(reinterpret_cast<const char *>(&k_2D), double_size);
        minimal_stream.write(reinterpret_cast<const char *>(&time), ldouble_size);
        minimal_stream.write(reinterpret_cast<const char *>(&asphericity_i), ldouble_size);
        minimal_stream.write(reinterpret_cast<const char *>(&shape_anisotropy_i),ldouble_size);

        int int_size = sizeof(int);
        int uint_size = sizeof(unsigned int);
        minimal_stream.write(reinterpret_cast<const char *>(&max_size), uint_size);
        minimal_stream.write(reinterpret_cast<const char *>(&N_1), int_size);
        minimal_stream.write(reinterpret_cast<const char *>(&N_2), int_size);
        minimal_stream.write(reinterpret_cast<const char *>(&count_steps), sizeof(long long unsigned int));   

        char temp_val; 
        int char_size = sizeof(char);
        for (unsigned int i = 0; i < N_tot; i++)
        {
            temp_val = sites[i];
            if (temp_val > 0)
            {
                minimal_stream.write(reinterpret_cast<const char *>(&temp_val), char_size);
                minimal_stream.write(reinterpret_cast<const char *>(&i), uint_size);
            }
        }
        minimal_stream.close();
    }
}

blocked_crystal::blocked_crystal(std::string filename, std::string label, std::string dir, int t_n)
{
    trial_number = t_n;
    std::ifstream load_file(filename, std::ios::in | std::ios::binary);
    //std::cout << "opening file...\n";
    if (load_file.is_open())
    {
        int string_size = sizeof(std::string);
        
        int double_size = sizeof(double);
        int ldouble_size = sizeof(long double);
        load_file.read(reinterpret_cast<char *>(&initial_size), double_size);
        load_file.read(reinterpret_cast<char *>(&T), double_size);
        load_file.read(reinterpret_cast<char *>(&mu), double_size);
        load_file.read(reinterpret_cast<char *>(&eps_11), double_size);
        load_file.read(reinterpret_cast<char *>(&eps_12), double_size);
        load_file.read(reinterpret_cast<char *>(&eps_22), double_size);
        load_file.read(reinterpret_cast<char *>(&k_1), double_size);
        load_file.read(reinterpret_cast<char *>(&k_2), double_size);
        load_file.read(reinterpret_cast<char *>(&k_1D), double_size);
        load_file.read(reinterpret_cast<char *>(&k_2D), double_size);
        
        int int_size = sizeof(int);
        int uint_size = sizeof(unsigned int);
        load_file.read(reinterpret_cast<char *>(&max_size), uint_size);
        
        filebase = label;
        direc = dir;

        max_size_sq = max_size * max_size;
        N_tot = max_size * max_size_sq;
        
        load_file.read(reinterpret_cast<char *>(&N_1), int_size);
        load_file.read(reinterpret_cast<char *>(&N_2), int_size);
        load_file.read(reinterpret_cast<char *>(&count_steps), sizeof(long long unsigned int));

        load_file.read(reinterpret_cast<char *>(&time), ldouble_size);
        load_file.read(reinterpret_cast<char *>(&asphericity), ldouble_size);
        load_file.read(reinterpret_cast<char *>(&asphericity_i), ldouble_size);
        load_file.read(reinterpret_cast<char *>(&shape_anisotropy), ldouble_size);
        load_file.read(reinterpret_cast<char *>(&shape_anisotropy_i),ldouble_size);
        load_file.read(reinterpret_cast<char *>(&k_tot),ldouble_size);
        
        //std::cout << "starting sites...\n";
        read_in_vec<char>(sites, load_file);
        //std::cout << "done with sites...\n";
        read_in_vec<char>(neighb_1, load_file);
        read_in_vec<char>(neighb_2, load_file);
        read_in_vec<unsigned int>(vector_indices, load_file);

        read_in_vec<double>(partial_sums, load_file);
        read_in_vec<int>(sums_order, load_file);
        read_in_vec<long double>(all_rates, load_file);
        read_in_vec<unsigned int>(rate_counts, load_file);

        load_file.read(reinterpret_cast<char *>(&max_size_active_sites), uint_size);
        //load_file.read(reinterpret_cast<char *>(&size_active_sites), uint_size);
        load_file.read(reinterpret_cast<char *>(&max_est), uint_size);

        unsigned int uitemp_val;
        max_size_active_rate_indices.reserve(poss_rates);
        num_open_slots_active_rate_indices.reserve(poss_rates);

        for (int i = 0; i < poss_rates; i++)
        {
            load_file.read(reinterpret_cast<char *>(&uitemp_val), uint_size);
            max_size_active_rate_indices.push_back(uitemp_val);
            read_in_vec<unsigned int>(active_rate_indices[i], load_file);
            read_in_vec<char>(active_rate_labels[i], load_file);
            
            load_file.read(reinterpret_cast<char *>(&uitemp_val), uint_size);
            num_open_slots_active_rate_indices.push_back(uitemp_val);
            read_in_vec<unsigned int>(open_slots_active_rate_indices[i], load_file);
        }
        read_active_sites(load_file);
        copy_surface_sites.resize(max_size_active_sites);
        read_in_vec<unsigned int>(surface_sites, load_file);
        load_file.read(reinterpret_cast<char *>(&Ns), int_size);
        load_file.read(reinterpret_cast<char *>(&Ns_2), int_size);
        load_file.read(reinterpret_cast<char *>(&num_open_slots_surface_sites), int_size);
        read_in_vec<unsigned int>(open_slots_surface_sites, load_file);
        load_file.read(reinterpret_cast<char *>(&num_open_slots_active_sites), uint_size);
        read_in_vec<unsigned int>(open_slots_active_sites, load_file);

        long double ldtemp_val;
        for (int i = 0; i < 3; i++)
        {
            load_file.read(reinterpret_cast<char *>(&ldtemp_val), ldouble_size);
            com[i] = ldtemp_val;

            for (int j = 0; j < 3; j++)
            {
                load_file.read(reinterpret_cast<char *>(&ldtemp_val), ldouble_size);
                smn[i][j] = ldtemp_val; 
            }
        }

        load_file.read(reinterpret_cast<char *>(&rho), double_size);
        load_file.read(reinterpret_cast<char *>(&rho_star), double_size);
        load_file.read(reinterpret_cast<char *>(&fixed_density), sizeof(bool));
        load_file.close();
        //std::cout << "loaded\n";
    }

    else
    {
        throw std::invalid_argument("Can't load dump file");
    } 
}

/*
blocked_crystal::blocked_crystal(std::string filename, std::string label, std::string dir, int t_n, int f_n)
{
    trial_number = t_n;
    std::ifstream load_file(filename, std::ios::in);
    int count_trajs = 1;
    int t_N_1, t_N_2;
    int temp_index, val, num_1, num_2, num_3, tail;

    if (load_file.is_open())
    {
        std::string trash;

        while (count_trajs < f_n)
        {
            for (int i = 0; i < 20; i++)
            {
                load_file >> trash;
            }
    
            for (int i = 0; i < 5; i++)
            {
                load_file >> trash;
            }

            load_file >> t_N_1 >> t_N_2 >> trash >> trash >> trash;
            
            for (int i = 0; i < (t_N_1 + t_N_2); i++)
            {
                load_file >> trash >> trash >> trash >> trash >> trash >> trash;
            }
            
            count_trajs += 1;
        }

        for (int i = 0; i < 10; i++)
        {
            load_file >> trash;
        }
    
        load_file >> T >> mu >> trash >> eps_12 >> eps_22 >> k_1 >> k_2 >> k_1D >> k_2D >> max_size;
        //std::cout << T << " " << mu << " " << eps_11 << " " << eps_12 << " " << eps_22 << " " << k_1 << " " << k_2 << " " << k_1D << " " << k_2D << " " << max_size << "\n"; 
        
        k_1D *= exp(mu/T);
        k_2D *= exp(mu/T);
        eps_11 = 1.0;

        pre_compute_rates();
        
        filebase = label;
        direc = dir;

        max_size_sq = max_size * max_size;
        N_tot = max_size * max_size_sq;

        for (int i = 0; i < 5; i++)
        {
            load_file >> trash;
        }

        sites.resize(N_tot);
        neighb_1.resize(N_tot);
        neighb_2.resize(N_tot);
        
        load_file >> N_1 >> N_2 >> time >> asphericity_i >> shape_anisotropy_i;
        
        for (int i = 0; i < (N_1 + N_2); i++)
        {
            load_file >> val >> temp_index >> num_1 >> num_2 >> num_3 >> tail;
            if (val == 1)
            {
                sites[temp_index] = 1;
            }

            else if (val == 2)
            {
                sites[temp_index] = 2;
            }

            else
            {
                sites[temp_index] = 3;
            }
        }

        load_file.close();
        operational = true;
        configure_surface(true);
    }

    else
    {
        operational = false;
        throw std::invalid_argument("Can't load config file.");
    }  
}
*/

blocked_crystal::blocked_crystal()
{
    time = 0;
    return;
}

blocked_crystal::atom::atom(unsigned int index_l, unsigned int index_v)
{
    lattice_index = index_l;
    vector_index = index_v;
    surface_index = 0;
    
    for (int i = 0; i < 13; i++)
    {
        poss_rate_indices[i] = 0;
    }
}

blocked_crystal::atom::atom()
{
    lattice_index = 4294967294;
    vector_index = 0;
    surface_index = 0;
    
    for (int i = 0; i < 13; i++)
    {
        poss_rate_indices[i] = 0;
    }
}

void blocked_crystal::integrate_order_params(double dt)
{
    long double a = -(smn[0][0] + smn[1][1] + smn[2][2]);
    long double b = ( smn[0][0] * smn[1][1] - smn[0][1] * smn[0][1] ) + ( smn[0][0] * smn[2][2] - smn[0][2] * smn[0][2] )  + ( smn[1][1] * smn[2][2] - smn[1][2] * smn[1][2] );
    long double c = -( smn[0][0] * smn[1][1] * smn[2][2] + 2 * smn[0][1] * smn[0][2] * smn[1][2] - smn[0][0] * smn[1][2] * smn[1][2] - smn[1][1] * smn[0][2] * smn[0][2] - smn[2][2] * smn[0][1] * smn[0][1]);

    long double p = -(a * a / 3.0 - b);
    long double q = (2 * a * a * a / 27.0 - a * b / 3.0 + c );
    long double l_1, l_2, l_3, l_x, l_y, l_z;

    if (abs(q) < tol)
    {
        l_x = abs(-a / 3.0);
        
        if (p <= 0)
        {
            l_y = abs(sqrt(-p) - a / 3.0);
            l_z = abs(-sqrt(-p) - a / 3.0);
        }

        else if (abs(p) < tol)
        {
            l_y = l_x;
            l_z = l_x;
        }

        else
        {
            std::cerr << "Error calculating asphericity " << filebase << " trial " << trial_number << "\n";
            asphericity_i += dt * asphericity;
            shape_anisotropy_i += dt * shape_anisotropy;
            return;
        }
    }

    else if (abs(p) < tol)
    {
        std::cerr << "Error calculating asphericity " << filebase << " trial " << trial_number << "\n";
        asphericity_i += dt * asphericity;
        shape_anisotropy_i += dt * shape_anisotropy;
        return;
    }
    
    
    else
    {
        long double big_A, big_B, phi;
        
        long double inner = -q * q / 4 - p * p * p / 27.0;
        
        big_A = -q / 2.0;
        
        if (abs(inner) < tol)
        {
            big_B = 0;
            phi = 0;
        }
        
        else if ( inner > 0 )
        {
            big_B = sqrt(-q * q / 4 - p * p * p / 27.0);
            phi = atan(big_B / big_A);
        }
        
        else
        {
            std::cerr << "Error calculating asphericity " << filebase << " trial " << trial_number << "\n";
            asphericity_i += dt * asphericity;
            shape_anisotropy_i += dt * shape_anisotropy;
            return;
        }
        
        long double mag = sqrt( big_A * big_A + big_B * big_B );
        
        long double mag_3 = pow(mag, 1.0 / 3.0);
        
        l_x = abs(mag_3 * 2 * cos( phi / 3.0 ) - a / 3.0);
        l_y = abs(mag_3 * 2 * cos( ( phi + 2 * pi ) / 3.0 ) - a / 3.0);
        l_z = abs(mag_3 * 2 * cos( ( phi + 4 * pi ) / 3.0 )  - a / 3.0);
        
    }

    if ( l_x > l_y )
    {
        if ( l_x > l_z )
        {
            l_3 = l_x;

            if ( l_y > l_z )
            {
                l_2 = l_y;
                l_1 = l_z;
            }
            else
            {
                l_2 = l_z;
                l_1 = l_y;
            }
        }
        else
        {
            l_3 = l_z;
            l_2 = l_x;
            l_1 = l_y;
        }
    }
    
    else if ( l_y > l_z )
    {
        l_3 = l_y;
        
        if ( l_z > l_x )
        {
            l_2 = l_z;
            l_1 = l_x;
        }
        else
        {
            l_2 = l_x;
            l_1 = l_z;
        }
    }
    
    else
    {
        l_3 = l_z;
        l_2 = l_y;
        l_1 = l_x;
    }
    asphericity = ( l_3 * l_3 - 0.5 * ( l_2 * l_2 + l_1 * l_1 ) );
    asphericity_i += dt * asphericity;
    shape_anisotropy = ( 1.5 * ( pow( l_1, 4) + pow(l_2, 4) + pow(l_3, 4) ) / pow( ( l_1 * l_1 + l_2 * l_2 + l_3 * l_3 ), 2) - 0.5 );
    shape_anisotropy_i += dt * shape_anisotropy;
    
    return;
}

void blocked_crystal::add_active_site(unsigned int site_index)
{
    unsigned int new_atom_index;

    if (num_open_slots_active_sites > 0)
    {
        new_atom_index = open_slots_active_sites[(num_open_slots_active_sites - 1)];
        open_slots_active_sites[(num_open_slots_active_sites - 1)] = 0;
        num_open_slots_active_sites--;
    }

    else
    {
        size_active_sites++;

        if (size_active_sites == max_size_active_sites)
        {
            max_size_active_sites *= increase_factor;
            active_sites.resize(max_size_active_sites);
            surface_sites.resize(max_size_active_sites);
            copy_surface_sites.resize(max_size_active_sites);
            open_slots_active_sites.resize(max_size_active_sites);
            open_slots_surface_sites.resize(max_size_active_sites);
        }

        new_atom_index = size_active_sites;
    }
    
    active_sites[new_atom_index] = atom(site_index, new_atom_index);
    
    vector_indices[site_index] = new_atom_index;
    int num_1 = neighb_1[site_index];
    
    if ( num_1 >= 2 && num_1 <= 5 )//num_1 <= 7 )
    {
        add_surface_site(new_atom_index);
    }

    add_all_rates(new_atom_index);
}

void blocked_crystal::add_surface_site(unsigned int vector_index)
{
    int site_index = active_sites[vector_index].lattice_index;
    
    unsigned int new_atom_index;
    if (num_open_slots_surface_sites > 0)
    {
        new_atom_index = open_slots_surface_sites[(num_open_slots_surface_sites - 1)];
        open_slots_surface_sites[(num_open_slots_surface_sites - 1)] = 0;
        num_open_slots_surface_sites--;
    }

    else
    {
        Ns++;  
        new_atom_index = Ns;
    }

    surface_sites[new_atom_index] = vector_index;
    active_sites[vector_index].surface_index = new_atom_index;
    
    if (sites[site_index] == 2)
    {
        Ns_2++;
    }
}
void blocked_crystal::remove_active_site(unsigned int vector_index)
{
    int lat_ind = active_sites[vector_index].lattice_index;
    int num_1 = neighb_1[lat_ind];
    int surf_ind = active_sites[vector_index].surface_index;
    
    if (surf_ind > 0)
    {
        remove_surface_site(surf_ind, vector_index);
    }

    num_open_slots_active_sites++;
    open_slots_active_sites[(num_open_slots_active_sites - 1)] = vector_index;
    vector_indices[active_sites[vector_index].lattice_index] = 0;

    active_sites[vector_index].lattice_index = N_tot;
    active_sites[vector_index].vector_index = 0;  
}


void blocked_crystal::remove_surface_site(unsigned int surf_ind, unsigned int vector_index)
{
    num_open_slots_surface_sites++;
    open_slots_surface_sites[(num_open_slots_surface_sites - 1)] = surf_ind;
    surface_sites[surf_ind] = 0;
    active_sites[vector_index].surface_index = 0;
    int lattice_index = active_sites[vector_index].lattice_index;
    
    if (sites[lattice_index] == 2)
    {
        Ns_2--;
    }
}

bool blocked_crystal::execute_move(unsigned int vector_index, int rate_num, int rate_label)
{
    unsigned int lattice_index = active_sites[vector_index].lattice_index;
   
    int old_val = sites[lattice_index];
    int num_1 = neighb_1[lattice_index];
    int num_2 = neighb_2[lattice_index];
  
    if (old_val == 0)
    {
        if (rate_label == 0)
        {
            nv = 1;
            remove_all_rates(vector_index);

            if (count_h)
            {
                add_to_hcounts(lattice_index);
            }

            add_to_gyration_tensor(lattice_index);
            energy -= (neighb_1[lattice_index] * eps_11 + neighb_2[lattice_index] * eps_12);

            sites[lattice_index] = 1;
            
            update_neighbors(lattice_index, 0, 1);
            add_all_rates(vector_index);
        }
        
        else if (rate_label == 1)
        {
            nv = 2;
            remove_all_rates(vector_index);
            
            N_2++;
            sites[lattice_index] = 2;
            energy -= (neighb_1[lattice_index] * eps_12 + neighb_2[lattice_index] * eps_22);
            if (active_sites[vector_index].surface_index > 0)
            {
                Ns_2++;
            }

            
            update_neighbors(lattice_index, 0, 2);
            add_all_rates(vector_index);
        }
        
        else
        {
            std::cerr << int(rate_num) << " " << int(rate_label) << "\n";
            std::cerr << "Illegal Move attempted in " << filebase << " trial " << trial_number << "\n";
            time = max_time;
            return false;
        }
    }

    else if (old_val == 1)
    {  
        if (rate_label == 12)
        {
            nv = 0;
            remove_all_rates(vector_index);

            energy += neighb_1[lattice_index] * eps_11 + neighb_2[lattice_index] * eps_12;
            if (count_h)
            {
                remove_from_hcounts(lattice_index);
            }
            remove_from_gyration_tensor(lattice_index);
            
            sites[lattice_index] = 0;
            
            update_neighbors(lattice_index, 1, 0);
            add_all_rates(vector_index);
        }

        else if (rate_label < 12)
        {
            nv = 1;
            unsigned int i, j, k, q, n_i, n_j, n_k;
            unsigned int n_lattice_index, n_vector_index;

            i = lattice_index / max_size_sq;
            q = lattice_index % max_size_sq;
            j = q / max_size;
            k = q % max_size;
            
            n_i = nnmod(i + order[rate_label][0], max_size);
            n_j = nnmod(j + order[rate_label][1], max_size);
            n_k = nnmod(k + order[rate_label][2], max_size);

            /*
            if (n_i < int(0.05 * max_size) + 1 || n_i > int(0.95 * max_size) || n_j < int(0.05 * max_size) + 1 || n_j > int(0.95 * max_size) || n_k < int(0.05 * max_size) + 1 || n_k > int(0.95 * max_size))
            {
                std::cerr << "Increasing max size in trial " << trial_number << "\n";
                must_update_max_size = true;
                return false;
            }
            */

            n_lattice_index = n_i * max_size_sq + n_j * max_size + n_k;
            energy += neighb_1[lattice_index] * eps_11 + neighb_2[lattice_index] * eps_12;
            remove_all_rates(vector_index);
        
            if (count_h)
            {
                remove_from_hcounts(lattice_index);
            }

            remove_from_gyration_tensor(lattice_index);
            
            sites[lattice_index] = 0;

            update_neighbors(lattice_index, 1, 0);

            //add_all_rates(vector_index);
            /*
            if (!brute_force_check())
            {
                time = max_time;
                std::cerr << "fatal error\n";
                return false;
            }
            */
         
            n_vector_index = vector_indices[n_lattice_index];

            if (n_vector_index > 0)
            {
                remove_all_rates(n_vector_index);
            }
            
            if (count_h)
            {
                add_to_hcounts(n_lattice_index);
            }

            add_to_gyration_tensor(n_lattice_index);
            
            sites[n_lattice_index] = 1;

            update_neighbors(n_lattice_index, 0, 1);
            if (n_vector_index == 0)
            {
                add_active_site(n_lattice_index);
            }
            else 
            {
                add_all_rates(n_vector_index);
            }
            
            /*
            if (!brute_force_check())
            {
                time = max_time;
                std::cerr << "fatal error\n";
                return false;
            }
            */
            
            energy -= neighb_1[n_lattice_index] * eps_11 + neighb_2[n_lattice_index] * eps_12;
           
        }
        else
        {
            std::cerr << int(rate_num) << " " << int(rate_label) << "\n";
            std::cerr << "Illegal Move attempted in " << filebase << " trial " << trial_number << "\n";
            time = max_time;
            return false;
        }
    }

    else if (old_val == 2)
    {
        if (rate_label == 12)
        {
            nv = 0;
            remove_all_rates(vector_index);
            
            sites[lattice_index] = 0;
            
            N_2--;
            energy += neighb_1[lattice_index] * eps_12 + neighb_2[lattice_index] * eps_22;
            
            if (active_sites[vector_index].surface_index > 0)
            {
                Ns_2--;
            }

            update_neighbors(lattice_index, 2, 0);
            add_all_rates(vector_index);
        }

        else if (rate_label < 12)
        {
            nv = 2;
            int i, j, k, q, n_i, n_j, n_k;
            unsigned int n_lattice_index, n_vector_index;

            i = lattice_index / max_size_sq;
            q = lattice_index % max_size_sq;
            j = q / max_size;
            k = q % max_size;

            n_i = nnmod(i + order[rate_label][0], max_size);
            n_j = nnmod(j + order[rate_label][1], max_size);
            n_k = nnmod(k + order[rate_label][2], max_size);

            /*
            if (n_i < int(0.05 * max_size) + 1 || n_i > int(0.95 * max_size) || n_j < int(0.05 * max_size) + 1 || n_j > int(0.95 * max_size) || n_k < int(0.05 * max_size) + 1 || n_k > int(0.95 * max_size))
            {
                std::cerr << "Increasing max size in trial " << trial_number << "\n";
                must_update_max_size = true;
                return false;
            }
            */
            
            n_lattice_index = n_i * max_size_sq + n_j * max_size + n_k;
            energy += neighb_1[lattice_index] * eps_12 + neighb_2[lattice_index] * eps_22;    
            remove_all_rates(vector_index);
            
            sites[lattice_index] = 0;
            
            if (active_sites[vector_index].surface_index > 0)
            {
                Ns_2--;
            }

            update_neighbors(lattice_index, 2, 0);
            
            //add_all_rates(vector_index);

            n_vector_index = vector_indices[n_lattice_index];
            
            if (n_vector_index > 0)
            {
                remove_all_rates(n_vector_index);
            }

            sites[n_lattice_index] = 2;
            
            if (active_sites[n_vector_index].surface_index > 0)
            {
                Ns_2++;
            }
    
            update_neighbors(n_lattice_index, 0, 2);   
            
            if (n_vector_index == 0)
            {
                add_active_site(n_lattice_index);
            }
            else 
            {
                add_all_rates(n_vector_index);
            }
            
            energy -= neighb_1[n_lattice_index] * eps_12 + neighb_2[n_lattice_index] * eps_22; 
        }

        else
        {
            std::cerr << int(rate_num) << " " << int(rate_label) << "\n";
            std::cerr << "Illegal Move attempted in " << filebase << " trial " << trial_number << "\n";
            time = max_time;
            return false;
        }
    }
    return true;
}

void blocked_crystal::add_to_gyration_tensor(int index)
{
    int i, j, k, q;
    std::array<double, 3> r_np1;
    double N = N_1 * 1.0;
    
    i = index / max_size_sq;
    q = index % max_size_sq;
    j = q / max_size;
    k = q % max_size;
    
    for (int s = 0; s < 3; s++)
    {
        r_np1[s] = (1/sqrt_2) * (xtl_1_xyz[s] * i + xtl_2_xyz[s] * j + xtl_3_xyz[s] * k - (max_size / 2) * (xtl_1_xyz[s] + xtl_2_xyz[s] + xtl_3_xyz[s]));
    }

    for (int m = 0; m < 3; m++)
    {
    for (int n = m; n < 3; n++)
    {
        smn[m][n] = ( N / (N + 1) ) * smn[m][n] + ( r_np1[m] - com[m] ) * ( r_np1[n] - com[n] ) / ( N + 1 ) * ( 1 - ( 1 / ( N + 1 ) ) );
        smn[n][m] = smn[m][n];
    }
    }

    for (int s = 0; s < 3; s++)
    {
        com[s] = ( N / ( N + 1 ) ) * com[s] + ( 1 / ( N + 1 ) ) * r_np1[s];
    }

    N_1 = N_1 + 1;
}

void blocked_crystal::remove_from_gyration_tensor(int index)
{
    int i, j, k, q;
    std::array<double, 3> r_n;
    double N = N_1 * 1.0;
    
    i = index / max_size_sq;
    q = index % max_size_sq;
    j = q / max_size;
    k = q % max_size;
    
    for (int s = 0; s < 3; s++)
    {
        r_n[s] = (1/sqrt_2) * (xtl_1_xyz[s] * i + xtl_2_xyz[s] * j + xtl_3_xyz[s] * k - (max_size / 2) * (xtl_1_xyz[s] + xtl_2_xyz[s] + xtl_3_xyz[s]));
    }

    for (int m = 0; m < 3; m++)
    {
        for (int n = m; n < 3; n++)
        {
            smn[m][n] = ( N / (N - 1) ) * smn[m][n] - ( r_n[m] - com[m] ) * ( r_n[n] - com[n] ) / ( N - 1 ) * ( 1 + ( 1 / ( N - 1 ) ) );
            smn[n][m] = smn[m][n];
        }
    }

    for (int s = 0; s < 3; s++)
    {
        com[s] = ( N / ( N - 1 ) ) * com[s] - ( 1 / ( N - 1 ) ) * r_n[s];
    }

    N_1 = N_1 - 1;
}
/*
void blocked_crystal::print_gyration_tensor()
{
    inertia_stream << N_1 << " " << N_2 << " " << time << "\n";
    
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            inertia_stream << smn[i][j] << " ";
        }
        inertia_stream << "\n";
    }
}
*/
void blocked_crystal::prepare_run(bool no_config, int seed_1)
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
        config_stream <<  "T mu eps_11 eps_12 eps_22 k_1 k_2 k_1D k_2D max_size\n";
        config_stream << T << " " << mu << " " << eps_11 << " " << eps_12 << " " << eps_22 << " " << k_1 << " " << k_2 << " " << k_1D << " " << k_2D  << " " << max_size << "\n";
        config_stream << "val index neighb_1 neighb_2 " << seed_1 << "\n";
    }

    data_stream <<  "T mu eps_11 eps_12 eps_22 k_1 k_2 k_1D k_2D max_size\n";
    data_stream << T << " " << mu << " " << eps_11 << " " << eps_12 << " " << eps_22 << " " << k_1 << " " << k_2 << " " << k_1D << " " << k_2D  << " " << max_size << "\n";
    data_stream <<  "N_1 N_2 time shape_anisotropy shape_anisotropy_i " << seed_1 << "\n";
    
    int double_size = sizeof(double);
    full_config_stream.write(reinterpret_cast<const char*>(&initial_size), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&T), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&mu), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&eps_11), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&eps_12), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&eps_22), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&k_1), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&k_2), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&k_1D), double_size);
    full_config_stream.write(reinterpret_cast<const char*>(&k_2D), double_size);

    int uint_size = sizeof(unsigned int);
    full_config_stream.write(reinterpret_cast<const char*>(&max_size), uint_size);
}

void blocked_crystal::remove_rate(int rate_num, int vector_index)
{
    /*
    if (rate_num >= poss_rates || rate_num < 0)
    {
        std::cout << "Attempting to remove non-existent process\n";
    }
    if (rate_counts[rate_num] == 0)
    {
        std::cout << "Attempting to remove non-existent process\n";
    }
    */
    
    rate_counts[rate_num]--;
    active_rate_indices[rate_num][vector_index] = 0;
    active_rate_labels[rate_num][vector_index] = 13;
    
    num_open_slots_active_rate_indices[rate_num]++;
    
    /*
    if (vector_index <= 0)
    {
        std::cout << "Illegal vector index " << vector_index << " !\n";
        time = max_time;
        return;
    }
    */
    
    open_slots_active_rate_indices[rate_num][num_open_slots_active_rate_indices[rate_num] - 1] = vector_index;
    if (num_open_slots_active_rate_indices[rate_num] == (max_size_active_rate_indices[rate_num] / resizing_factor))    
    {
        unsigned int count_swaps = 0;
        unsigned int swap_index;
        char swap_label;
        unsigned int next_open_rate_index;

        for (int j = rate_counts[rate_num] + num_open_slots_active_rate_indices[rate_num]; j > 0; j--)
        {
            swap_index = active_rate_indices[rate_num][j];
            swap_label = active_rate_labels[rate_num][j];
            
            if (swap_index > 0)
            {
                next_open_rate_index = open_slots_active_rate_indices[rate_num][count_swaps];
                
                if (next_open_rate_index < j)
                {
                    active_rate_indices[rate_num][next_open_rate_index] = swap_index;
                    active_rate_labels[rate_num][next_open_rate_index] = swap_label;
                    
                    active_sites[swap_index].poss_rate_indices[swap_label] = next_open_rate_index;
                    
                    active_rate_indices[rate_num][j] = 0;
                    active_rate_labels[rate_num][j] = poss_rates;
                }
                
                else
                {
                    j++;
                }

                open_slots_active_rate_indices[rate_num][count_swaps] = 0;
                count_swaps++;
            }

            if (count_swaps == num_open_slots_active_rate_indices[rate_num])
            {
                break;
            }
        }

        num_open_slots_active_rate_indices[rate_num] = 0;
    }
}


unsigned int blocked_crystal::add_rate(int rate_num, int vector_index, int rate_label)
{
    int new_index;
    /*
    if (rate_num >= poss_rates || rate_num < 0)
    {
        std::cout << "Attempting to add non-existent process\n";
    }
    */

    int num_slots = num_open_slots_active_rate_indices[rate_num];
    rate_counts[rate_num]++;

    if (num_slots > 0)
    {
        new_index = open_slots_active_rate_indices[rate_num][ (num_slots - 1) ];
        open_slots_active_rate_indices[rate_num][ (num_slots-1) ] = 0;
        num_open_slots_active_rate_indices[rate_num]--;
    }

    else
    {
        if (rate_counts[rate_num] + num_slots == max_size_active_rate_indices[rate_num])
        {
            max_size_active_rate_indices[rate_num] *= 2;
            active_rate_indices[rate_num].resize(max_size_active_rate_indices[rate_num]);
            active_rate_labels[rate_num].resize(max_size_active_rate_indices[rate_num]);
            open_slots_active_rate_indices[rate_num].resize(max_size_active_rate_indices[rate_num]);
            
            int new_slot_size = max_size_active_rate_indices[rate_num] / resizing_factor;
            
            open_slots_active_rate_indices[rate_num].resize(max_size_active_rate_indices[rate_num] / resizing_factor);
        }

        new_index = rate_counts[rate_num] + num_slots;
    }
    
    active_rate_indices[rate_num][new_index] = vector_index;
    active_rate_labels[rate_num][new_index] = rate_label;
    
    return new_index;
}

void blocked_crystal::remove_all_rates(unsigned int vector_index)
{
    unsigned int lat_index = active_sites[vector_index].lattice_index;
    
    /*
    if ( vector_index == 0 || vector_index > size_active_sites)
    {
        std::cout << "Illegal vector index " << vector_index << " !\n";
    }
    */

    int num_1 = neighb_1[lat_index];
    int num_2 = neighb_2[lat_index];
    int val = sites[lat_index];
    unsigned int curr_ind;

    if (val == 0)
    {
        curr_ind = active_sites[vector_index].poss_rate_indices[0];
        if (curr_ind > 0)
        {
            remove_rate(6 * (17 * num_1 + num_2), curr_ind);
            active_sites[vector_index].poss_rate_indices[0] = 0;
        }
        
        curr_ind = active_sites[vector_index].poss_rate_indices[1];
        if (curr_ind > 0)
        {
            remove_rate(6 * (17 * num_1 + num_2) + 1, curr_ind);
            active_sites[vector_index].poss_rate_indices[1] = 0;
        }
    }

    else if (val == 1)
    {
        int i, j, k;
        i = lat_index / max_size_sq;   
        j = (lat_index / max_size) % max_size;
        k = lat_index % max_size;
        curr_ind = active_sites[vector_index].poss_rate_indices[12];
        
        if (curr_ind > 0)
        {
            remove_rate(6 * (17 * num_1 + num_2) + 2, curr_ind);
            active_sites[vector_index].poss_rate_indices[12] = 0;
        }

        int n_i, n_j, n_k, n_index;
        for (int l = 0; l < 12; l++)
        {
            curr_ind = active_sites[vector_index].poss_rate_indices[l];
            if (curr_ind > 0)
            {
                n_i = nnmod(i + order[l][0], max_size);
                n_j = nnmod(j + order[l][1], max_size);
                n_k = nnmod(k + order[l][2], max_size);
                n_index = n_i * max_size_sq + n_j * max_size + n_k;
                remove_rate(6 * (17 * (neighb_1[n_index] - 1 - num_1 + 8) + (neighb_2[n_index] - num_2 + 8)) + 4, curr_ind);
                active_sites[vector_index].poss_rate_indices[l] = 0;
            }
        }
    }

    else if (val == 2)
    {
        int i, j, k;
        i = lat_index / max_size_sq;   
        j = (lat_index / max_size) % max_size;
        k = lat_index % max_size;
        
        curr_ind = active_sites[vector_index].poss_rate_indices[12];
        if (curr_ind > 0)
        {
            remove_rate(6 * (17 * num_1 + num_2) + 3, curr_ind);
            active_sites[vector_index].poss_rate_indices[12] = 0;
        }
        
        int n_i, n_j, n_k, n_index;
        for (int l = 0; l < 12; l++)
        {
            curr_ind = active_sites[vector_index].poss_rate_indices[l];
            if (curr_ind > 0)
            {
                n_i = nnmod(i + order[l][0], max_size);
                n_j = nnmod(j + order[l][1], max_size);
                n_k = nnmod(k + order[l][2], max_size);
                n_index = n_i * max_size_sq + n_j * max_size + n_k;
                remove_rate(6 * (17 * (neighb_1[n_index] - num_1 + 8) + (neighb_2[n_index] - 1 - num_2 + 8)) + 5, curr_ind);
                active_sites[vector_index].poss_rate_indices[l] = 0;
            }
        }
    }

    else
    {
        std::cerr << "Error in removing rates!\n";
	    time = max_time;
	    return;
    }
}

void blocked_crystal::add_all_rates(unsigned int vector_index)
{   
    unsigned int lat_index = active_sites[vector_index].lattice_index;
    int val = sites[lat_index];
    int num_1 = neighb_1[lat_index];
    int num_2 = neighb_2[lat_index];

    if (val == 0)
    {
        if ( (num_1 + num_2) > 0 && (num_1 + num_2) < 12 )
        {
            active_sites[vector_index].poss_rate_indices[0] = add_rate(6 * (17 * num_1 + num_2), vector_index, 0);
            active_sites[vector_index].poss_rate_indices[1] = add_rate(6 * (17 * num_1 + num_2) + 1, vector_index, 1);
        }
    }

    else if (val == 1)
    {
        if ((num_1 + num_2) > 0 && (num_1 + num_2) < 12)
        {
            active_sites[vector_index].poss_rate_indices[12] = add_rate(6 * (17 * num_1 + num_2) + 2, vector_index, 12);
        }

        int i, j, k, q, n_i, n_j, n_k;
        unsigned int n_index;

        i = lat_index / max_size_sq;
        q = lat_index % max_size_sq;
        j = q / max_size;
        k = q % max_size;
        
        for (int l = 0; l < 12; l++)
        {
            n_i = nnmod(i + order[l][0], max_size);
            n_j = nnmod(j + order[l][1], max_size);
            n_k = nnmod(k + order[l][2], max_size);
            
            /*
            if (n_i < int(0.05 * max_size) + 1 || n_i > int(0.95 * max_size) || n_j < int(0.05 * max_size) + 1 || n_j > int(0.95 * max_size) || n_k < int(0.05 * max_size) + 1 || n_k > int(0.95 * max_size))
            {
                std::cerr << "Reached boundary in " << trial_number << "\n";
                must_update_max_size = true;
                return;
            }
            */

            n_index = n_i * max_size_sq + n_j * max_size + n_k;
            if (sites[n_index] == 0)
            {
                active_sites[vector_index].poss_rate_indices[l] = add_rate(6 * (17 * (neighb_1[n_index] - 1 - num_1 + 8) + (neighb_2[n_index] - num_2 + 8) ) + 4, vector_index, l);
            }

            else
            {
                active_sites[vector_index].poss_rate_indices[l] = 0;
            }
        }
    }
    
    else if (val == 2)
    {
        if ((num_1 + num_2) > 0 && (num_1 + num_2) < 12)
        {
            active_sites[vector_index].poss_rate_indices[12] = add_rate(6 * (17 * num_1 + num_2) + 3, vector_index, 12);
        }
          
        int i, j, k, q, n_i, n_j, n_k;
        unsigned int n_index;

        i = lat_index / max_size_sq;
        q = lat_index % max_size_sq;
        j = q / max_size;
        k = q % max_size;

        for (int l = 0; l < 12; l++)
        {
            n_i = nnmod(i + order[l][0], max_size);
            n_j = nnmod(j + order[l][1], max_size);
            n_k = nnmod(k + order[l][2], max_size);
            
            /*
            if (n_i < int(0.05 * max_size) + 1 || n_i > int(0.95 * max_size) || n_j < int(0.05 * max_size) + 1 || n_j > int(0.95 * max_size) || n_k < int(0.05 * max_size) + 1 || n_k > int(0.95 * max_size))
            {
                std::cerr << "size too small " << trial_number << "\n";
                must_update_max_size = true;
            }
            */

            n_index = n_i * max_size_sq + n_j * max_size + n_k;

            if (sites[n_index] == 0)
            {
                active_sites[vector_index].poss_rate_indices[l] = add_rate(6 * (17 * (neighb_1[n_index] - num_1 + 8) + (neighb_2[n_index] - 1 - num_2 + 8) ) + 5, vector_index, l);
            }
            
            else
            {
                active_sites[vector_index].poss_rate_indices[l] = 0;
            }
        }
    }
}

void blocked_crystal::pre_compute_rates(double time_conversion)
{
    // Pre computes all possible rates scaling them by time_conversion factor
    long double boltz;
    int j;
    int i, r, num_1, num_2;
    for (int q = 0; q < 1734; q++) 
    {
        i = q / 6;
        r = q % 6;
        num_1 = (i / 17); 
        num_2 = (i % 17);
        
        if (r == 0)
        {
            all_rates[q] = k_1 * exp((eps_11 * num_1 + eps_12 * num_2 + mu )/(2 * T)) * time_conversion; // Deposition rates for metal
        }
        else if (r == 1)
        {
            all_rates[q] = k_2 * exp((eps_12 * num_1 + eps_22 * num_2 + mu )/(2 * T)) * time_conversion; // Depostion rates for ligand 
        }
        else if (r == 2)
        {
            all_rates[q] = k_1 * exp(-(eps_11 * num_1 + eps_12 * num_2 + mu )/(2 * T)) * time_conversion; // Sublimation rates for metal 
        }
        else if (r == 3)
        {
            all_rates[q] = k_2 * exp(-(eps_12 * num_1 + eps_22 * num_2 + mu )/(2 * T)) * time_conversion; // Sublimation rates for ligand 
        }
        else if (r >= 4)
        {
            num_1 -= 8; // Now difference in number of nearest neighbors, final - initial. 
            num_2 -= 8;

            if (r == 4)
            {
                all_rates[q] = k_1D * exp((eps_11 * num_1 + eps_12 * num_2)/(2 * T)) * time_conversion;
            }
            else
            {
                all_rates[q] = k_2D * exp((eps_12 * num_1 + eps_22 * num_2)/(2 * T)) * time_conversion;  
            }
        }
    }
}

void blocked_crystal::configure_surface(bool first)
{
    // Determining initial rates
    partial_sums.resize(poss_rates);
    sums_order.resize(poss_rates);
    rate_counts.resize(poss_rates);
    num_open_slots_active_rate_indices.resize(poss_rates);

    //max_est = max_size_sq / 6 + resizing_factor;
    max_est = max_size_sq + resizing_factor;
    unsigned int new_slot_size = max_est / resizing_factor;
    max_size_active_rate_indices.resize(poss_rates, max_est);

    for (int i = 0; i < 1734; i++)
    {
        active_rate_indices[i].resize(max_est);
        active_rate_labels[i].resize(max_est);
        open_slots_active_rate_indices[i].resize(new_slot_size);
    }

    max_size_active_sites = max_size_sq * 6;
    active_sites.resize(max_size_active_sites);
    surface_sites.resize(max_size_active_sites);
    copy_surface_sites.resize(max_size_active_sites);
    open_slots_active_sites.resize(max_size_active_sites);
    open_slots_surface_sites.resize(max_size_active_sites);
    vector_indices.resize(N_tot);

    num_open_slots_active_sites = 0;
    num_open_slots_surface_sites = 0;
    size_active_sites = 0;

    int val;
    unsigned int index = 0;
    index = 0;

    N_1 = 0;
    N_2 = 0;
    Ns = 0;
    Ns_2 = 0;

    asphericity = 0;
    shape_anisotropy = 0;

    for (int i = 0; i < 3; i++)
    {
        com[i] = 0;

        for (int j = 0; j < 3; j++)
        {
            smn[i][j] = 0;
        }
    }

    for (int i = 0; i < max_size; i++ )
    {
        for (int j = 0; j < max_size; j++ )
        {
            for (int k = 0; k < max_size; k++ )
            {
                val = sites[index];
                add_site_to_new_surface(index, val);
                index += 1;
            }
        }
    }

    for (int i = 0; i < N_tot; i++)
    {
        if ((sites[i] > 0 && (neighb_1[i] + neighb_2[i]) < 12) || (sites[i] == 0 && (neighb_1[i] + neighb_2[i]) > 0 ))
        {
            add_active_site(i);
        }
    }
}

void blocked_crystal::init_ligs(double r_star, int seed)
{
    k_2 = 0.0;
    pre_compute_rates(1.0);

    fixed_density = true;
    rho_star = r_star;
    liggen.seed(seed);

    int true_Ns = Ns - num_open_slots_surface_sites;
    int lat_ind, v_ind, old_val;
    std::copy_if( surface_sites.begin() + 1, surface_sites.begin() + Ns + 1, copy_surface_sites.begin(), [](int x) { return x > 0; });
    std::shuffle( copy_surface_sites.begin(), copy_surface_sites.begin() + true_Ns, liggen);
    int counter = 0;
    int targ_Ns_2 = int(rho_star * (Ns - num_open_slots_surface_sites));
    
    while( Ns_2 < targ_Ns_2 )
    {
        v_ind = copy_surface_sites[counter];
        lat_ind = active_sites[v_ind].lattice_index;
        old_val = sites[lat_ind];
        /*
        if (old_val != 2 && active_sites[v_ind].surface_index > 0)
        {
            if (old_val == 1)
            {
                remove_from_gyration_tensor(lat_ind);
            }
    
            remove_all_rates(v_ind);
            
            sites[lat_ind] = 2;
            
            N_2++;
            Ns_2++;
            
            add_all_rates(v_ind);
            update_neighbors(lat_ind, old_val, 2);
        }
        */
        if (old_val == 0 && active_sites[v_ind].surface_index > 0)
        {
            remove_all_rates(v_ind);

            sites[lat_ind] = 2;

            N_2++;
            Ns_2++;

            add_all_rates(v_ind);
            update_neighbors(lat_ind, 0, 2);
        }

        counter++;
        targ_Ns_2 = int(rho_star * (Ns - num_open_slots_surface_sites));
        //std::cout << N_1 << " " << N_2 << " " << Ns << " " << Ns_2 << " " << targ_Ns_2 << "\n";
    }
    

    /*
    for (int i = 0; i < targ_Ns_2; i++)
    {
        v_ind = ligdis(liggen);
        lat_ind = active_sites[v_ind].lattice_index;
        
        if (lat_ind < N_tot)
        {
            old_val = sites[lat_ind];
            num_1 = neighb_1[lat_ind];
        }

        while (lat_ind >= N_tot || old_val == 2 || num_1 < 3 || num_1 > 5)
        {
            v_ind = ligdis(liggen);
            lat_ind = active_sites[v_ind].lattice_index;
            
            if (lat_ind < N_tot)
            {
                old_val = sites[lat_ind];
                num_1 = neighb_1[lat_ind];
            }
        }
        
        if (old_val == 1)
        {
             remove_from_gyration_tensor(lat_ind);
        }
    
        remove_all_rates(v_ind);
        
        sites[lat_ind] = 2;
        N_2++;
        Ns_2++;
        
        add_all_rates(v_ind);
        update_neighbors(lat_ind, old_val, 2);
    }
    */
}

void blocked_crystal::add_ligs()
{
    int true_Ns = Ns - num_open_slots_surface_sites;
    int lat_ind, v_ind, old_val;
    int targ_Ns_2 = int(rho_star * true_Ns); 
    //std::cout << "Moving " << Ns_2 << " to " << targ_Ns_2 << "\n";
    int counter = 0;
    if ( Ns_2 < targ_Ns_2)
    {
        std::copy_if( surface_sites.begin() + 1, surface_sites.begin() + Ns + 1, copy_surface_sites.begin(), [](int x) { return x > 0; });
        std::shuffle(copy_surface_sites.begin(), copy_surface_sites.begin() + true_Ns, liggen);
        
        while (Ns_2 < targ_Ns_2 && counter < true_Ns)
        {
            //std::cout << "adding ligs" << "\n";
            v_ind = copy_surface_sites[counter];
            lat_ind = active_sites[v_ind].lattice_index;
            old_val = sites[lat_ind];
            /*
            if (old_val != 2 && active_sites[v_ind].surface_index > 0)
            {
                if (old_val == 1)
                {
                    remove_from_gyration_tensor(lat_ind);
                }
        
                remove_all_rates(v_ind);
                
                sites[lat_ind] = 2;
                
                N_2++;
                Ns_2++;
                
                add_all_rates(v_ind);
                update_neighbors(lat_ind, old_val, 2);
            }
            */
            if (old_val == 0 && active_sites[v_ind].surface_index > 0)
            {
                remove_all_rates(v_ind);
                
                sites[lat_ind] = 2;
                
                N_2++;
                Ns_2++;
                
                add_all_rates(v_ind);
                update_neighbors(lat_ind, 0, 2);
            }
            counter++;
            targ_Ns_2 = int(rho_star * (Ns - num_open_slots_surface_sites));
            //std::cout << N_1 << " " << N_2 << " " << Ns << " " << Ns_2 << " " << targ_Ns_2 << "\n";
        }
    }

    else if ( Ns_2 > targ_Ns_2)
    {
        std::copy_if( surface_sites.begin() + 1, surface_sites.begin() + Ns + 1, copy_surface_sites.begin(), [](int x) { return x > 0; });
        std::shuffle( copy_surface_sites.begin(), copy_surface_sites.begin() + true_Ns, liggen);
        //std::uniform_int_distribution<int> zero_or_one(0,1);
        while ( Ns_2 > targ_Ns_2 && counter < true_Ns)
        {
            //std::cout << "adding ligs" << "\n";
            v_ind = copy_surface_sites[counter];
            lat_ind = active_sites[v_ind].lattice_index;
            old_val = sites[lat_ind];
            
            /*
            int new_val;
            if (old_val == 2 && active_sites[v_ind].surface_index > 0)
            {
                remove_all_rates(v_ind);
                
                N_2--;
                Ns_2--;
                
                new_val = zero_or_one(liggen);
                
                sites[lat_ind] = new_val;
                
                if (new_val == 1)
                {
                    add_to_gyration_tensor(lat_ind);
                }
                
                add_all_rates(v_ind);
                update_neighbors(lat_ind, old_val, new_val);
            }
            */
            if (old_val == 2 && active_sites[v_ind].surface_index > 0)
            {
                remove_all_rates(v_ind);
                
                N_2--;
                Ns_2--;
                
                sites[lat_ind] = 0;
                
                add_all_rates(v_ind);
                update_neighbors(lat_ind, old_val, 0);
            }
            counter++;
            targ_Ns_2 = int(rho_star * (Ns - num_open_slots_surface_sites));
            //std::cout << N_1 << " " << N_2 << " " << Ns << " " << Ns_2 << " " << targ_Ns_2 << "\n";
        }
    }
    //std::cout << "Now we have " << Ns_2 << " with " << N_1 << " and " << N_2 << " after " << counter << " iterations\n";
    //copy_surface_sites.clear();
    /*
    int lat_ind = N_tot;
    int v_ind, old_val = 2, num_1 = 0;

    ligdis.param(std::uniform_int_distribution<int>::param_type(1, size_active_sites));
    for (int i = 0; i < add_2; i++)
    {
        v_ind = ligdis(liggen);
        lat_ind = active_sites[v_ind].lattice_index;
        
        if (lat_ind < N_tot)
        {
            old_val = sites[lat_ind];
            num_1 = neighb_1[lat_ind];
        }

        while (lat_ind >= N_tot || old_val == 2 || num_1 < 3 || num_1 > 5)
        {
            v_ind = ligdis(liggen);
            lat_ind = active_sites[v_ind].lattice_index;
            
            if (lat_ind < N_tot)
            {
                old_val = sites[lat_ind];
                num_1 = neighb_1[lat_ind];
            }
        }

        if (old_val == 1)
        {
            remove_from_gyration_tensor(lat_ind);
        }

        remove_all_rates(v_ind);
        
        sites[lat_ind] = 2;
        N_2++;
        Ns_2++;
        
        add_all_rates(v_ind);
        update_neighbors(lat_ind, old_val, 2);   
    }
    */
}

void blocked_crystal::add_site_to_new_surface(unsigned int index, int val)
{
    int num_1 = 0, num_2 = 0;
    unsigned int i, j, k, q, n_i, n_j, n_k;

    i = index / max_size_sq;
    q = index % max_size_sq;
    j = q / max_size;
    k = q % max_size;
    
    unsigned int n_index;
    int n_val;

    for (int m = 0; m < 12; m++)
    {
        n_i = nnmod(i + order[m][0], max_size);
        n_j = nnmod(j + order[m][1], max_size);
        n_k = nnmod(k + order[m][2], max_size);
        
        n_index = n_i * max_size_sq + n_j * max_size + n_k;

        n_val = sites[n_index];
        
        if (n_val == 1)
        {
            num_1 += 1;
        }

        else if (n_val == 2)
        {
            num_2 += 1;
        }
    }
        
    neighb_1.push_back(num_1);
    neighb_2.push_back(num_2);

    if (val == 1)
    {
        energy -= 0.5 * (neighb_1[index] * eps_11 + neighb_2[index] * eps_12);
        add_to_gyration_tensor(index);
        //std::cout << energy << " " << N_1 << " " << energy/N_1 << "\n";
    }

    else if (val == 2)
    {
        energy -= 0.5 * (neighb_1[index] * eps_12 + neighb_2[index] * eps_22);
        N_2 += 1;
    }
}

void blocked_crystal::update_neighbors(unsigned int lattice_index, int old_val, int new_val)
{
    unsigned int n_lattice_index, n_vector_index;
    int num_1, num_2;
    unsigned int i, j, k, q, n_i, n_j, n_k;
    int n_val;
    
    i = lattice_index / max_size_sq;
    q = lattice_index % max_size_sq;
    j = q / max_size;
    k = q % max_size;

    for (int l = 0; l < 54; l++)
    {
        n_i = nnmod(i + order2[l][0], max_size);
        n_j = nnmod(j + order2[l][1], max_size);
        n_k = nnmod(k + order2[l][2], max_size);

        /*
        if (n_i < int(0.05 * max_size) + 1 || n_i > int(0.95 * max_size) || n_j < int(0.05 * max_size) + 1 || n_j > int(0.95 * max_size) || n_k < int(0.05 * max_size) + 1 || n_k > int(0.95 * max_size))
        {
            std::cerr << "Increasing max size in trial " << trial_number << "\n";
            must_update_max_size = true;
            return;
        }
        */

        n_lattice_index = n_i * max_size_sq + n_j * max_size + n_k;
        n_vector_index = vector_indices[n_lattice_index];
        
        if (n_vector_index > 0)
        {
            remove_all_rates(n_vector_index);
        }
    }

    for (int l = 0; l < zlat; l++)
    {
        n_i = nnmod(i + order[l][0], max_size);
        n_j = nnmod(j + order[l][1], max_size);
        n_k = nnmod(k + order[l][2], max_size);

        n_lattice_index = n_i * max_size_sq + n_j * max_size + n_k;
    
        if (new_val == 1)
        {
            neighb_1[n_lattice_index] += 1;
        }
        
        else if (new_val == 2)
        {
            neighb_2[n_lattice_index] += 1;
        }

        if (old_val == 1)
        {
            neighb_1[n_lattice_index] -= 1;
        }

        else if (old_val == 2)
        {
            neighb_2[n_lattice_index] -= 1;
        }   
    }

    for (int l = 0; l < 54; l++)
    {
        n_i = nnmod(i + order2[l][0], max_size);
        n_j = nnmod(j + order2[l][1], max_size);
        n_k = nnmod(k + order2[l][2], max_size);

        n_lattice_index = n_i * max_size_sq + n_j * max_size + n_k;
        n_vector_index = vector_indices[n_lattice_index];
        n_val = sites[n_lattice_index];
        num_1 = neighb_1[n_lattice_index];
        num_2 = neighb_2[n_lattice_index];
        
        if ( (n_val > 0 && (num_1 + num_2) < 12) || (n_val == 0 && (num_1 + num_2) > 0 ) )
        {
            if (n_vector_index == 0)
            {
                add_active_site(n_lattice_index);
            }

            else
            {
                add_all_rates(n_vector_index);
                int surf_ind = active_sites[n_vector_index].surface_index;
                if (num_1 >= 2 && num_1 <= 5) 
                {
                    if (surf_ind == 0)
                    {
                        add_surface_site(n_vector_index);
                    }
                }
                else if (surf_ind > 0)
                {
                    remove_surface_site(surf_ind, n_vector_index);
                }
            }
        }
        
        else if (n_vector_index != 0)
        {
            remove_active_site(n_vector_index);   
        }
    }
}

void blocked_crystal::kmc_step()
{   
    // Create the partial sums list of the rates
    long double curr_rate = 0.0;
    long double new_rate;
    int count = 0;
    for (int q = 0; q < poss_rates; q++)
    {
        new_rate = curr_rate + rate_counts[q] * all_rates[q];
        if (new_rate > curr_rate)
        {
            curr_rate = new_rate;
            partial_sums[count] = curr_rate;
            sums_order[count] = q;
            count++;
        }
    }

    // Choose the random process to occur in proportion to its rate 
    k_tot = curr_rate;
    long double cap = contdis(contgen) * k_tot;
    auto it = std::upper_bound(partial_sums.begin(), partial_sums.begin() + count, cap);
    int process_index = int(it - partial_sums.begin());
    int process = sums_order[process_index];
    
    // Choose the lattice position where that process occurs with uniform probability for all equivalent lattice sites
    int vector_index = 0;
    int debug_count = 0;
    unsigned int trial_index;
    int process_label = 0;
    while (vector_index == 0)
    {
        trial_index = int(1 + contdis(contgen) * ((rate_counts[process] + num_open_slots_active_rate_indices[process]) - 1)  + 0.5);
        vector_index = active_rate_indices[process][trial_index];
        process_label = active_rate_labels[process][trial_index];
        debug_count++;
        if (debug_count == max_size_active_rate_indices[process] * 25)
        {
            std::cerr << "We got stuck in  " << filebase << " trial " << trial_number << " trying to do process " << r << " with count "  << rate_counts[r] << " and open slots " << num_open_slots_active_rate_indices[r] <<  " At time rate cap and sum = " << time <<  " " << k_tot << " " << cap << "\n";
            time = max_time;
            return;
        }
    }
    
    // Change the lattice in accordance with the process
    execute_move(vector_index, process, process_label);
    time -= log( contdis(contgen) ) / k_tot;
    count_steps++;
    return;
}

void blocked_crystal::kmc_run_time(double total_time, int config_rec, int data_rec, unsigned int seed, unsigned int seed_2, bool no_config)
{
    prepare_run(no_config, seed);
    
    contgen.seed(seed);
    liggen.seed(seed_2);
    contdis.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
    
    if (!no_config)
    {
        print_config();
        //print_everything();
    }
    
    double config_time = total_time / config_rec;
    double data_time = total_time  / data_rec;

    double start_time = time;
    double config_start_time = time;
    double data_start_time = time;

    double data_time_elapsed = 0.0;
    double config_time_elapsed = 0.0;
    data_stream << N_1 << " " << N_2 << " "  << count_steps << " " << time << " " << energy << " " << rho << " " << asphericity << " "  << asphericity_i << " " << shape_anisotropy << " " << shape_anisotropy_i << " ";
    for (int i = 0; i < 3; i++)
    {
        data_stream << h100[i][0] << " " << h100[i][1] << " ";
    }
    for (int i = 0; i < 4; i++)
    {
        data_stream << h111[i][0] << " " << h111[i][1] << " ";
    }
    data_stream << "\n";
    count_steps = 0;
    while ( (time - start_time) < total_time )
    {
        kmc_step();
        brute_force_check();

        if (must_update_max_size)
        {
            std::string msg = "Stopping trial " + std::to_string(trial_number) + " due to reached boundary\n";
            throw std::overflow_error(msg);
        }
        if (stationary)
        {
            std::string msg = "Stopping trial " + std::to_string(trial_number) + " due to effectively zero total rate\n";
            throw std::underflow_error(msg);   
        }

        data_time_elapsed = time - data_start_time;
        config_time_elapsed = time - config_start_time;

        if ( data_time_elapsed >= data_time && (time - start_time) < total_time )
        {
            data_stream << N_1 << " " << N_2 << " "  << count_steps << " " << time << " " << energy << " " << rho << " " << asphericity << " "  << asphericity_i << " " << shape_anisotropy << " " << shape_anisotropy_i << " ";
            for (int i = 0; i < 3; i++)
            {
                data_stream << h100[i][0] << " " << h100[i][1] << " ";
            }
            for (int i = 0; i < 4; i++)
            {
                data_stream << h111[i][0] << " " << h111[i][1] << " ";
            }
            data_stream << "\n";
            //std::cout << N_1 << " " << N_2 << " " <<  count_steps << " " << time << " " << asphericity << " "  << asphericity_i << " " << shape_anisotropy << " " << shape_anisotropy_i << "\n";
            data_start_time = time;
        }
    
        if ( config_time_elapsed >= config_time && (time - start_time) < total_time && !no_config)
        {
            print_config();
            config_start_time = time;
        }
	//count_steps++;
    }

    //print_full_config();
    data_stream << N_1 << " " << N_2 << " "  << count_steps << " " << time << " " << rho << " " << asphericity << " "  << asphericity_i << " " << shape_anisotropy << " " << shape_anisotropy_i << "\n";
    
    full_config_stream.close();
    if (!no_config)
    {
        print_config();
        config_stream.close();
    }
    data_stream.close();
}

void blocked_crystal::kmc_run_mass_down(unsigned int N_1_end, int N_1_end_high, int config_rec, int data_rec, unsigned int seed, unsigned int seed_2, bool no_config)
{
    prepare_run(no_config, seed);
    contgen.seed(seed);
    liggen.seed(seed_2);
    contdis.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
    
    int config_mass = (N_1 - N_1_end) / config_rec;
    int data_mass = (N_1 - N_1_end) / data_rec;

    int config_start_mass = N_1;
    int data_start_mass = N_1;

    int data_mass_elapsed = 0;
    int config_mass_elapsed = 0;
    
    if (!no_config)
    {
        print_config();
    }

    count_steps = 0;

    while ( N_1 > N_1_end && time < max_time && count_steps < N_1_end_high )
    {
        kmc_step();

        if (must_update_max_size)
        {
            std::cerr << "Stopping trial " << trial_number << " due to reached boundary\n";
            break;
        }
        if (stationary)
        {
            std::cerr << "Stopping trial " << trial_number << " due to reached stationary state\n";
            break;   
        }

        data_mass_elapsed = data_start_mass - N_1;
        config_mass_elapsed = config_start_mass - N_1;
    
        
        if ( data_mass_elapsed >= data_mass && time < max_time )
        {
            data_stream << N_1 << " " << N_2 << " "  << count_steps << " " << time << " " << rho << " " << asphericity << " "  << asphericity_i << " " << shape_anisotropy << " " << shape_anisotropy_i << "\n";
            //std::cout << N_1 << " " << N_2 << " " << count_steps << " " << time << " " << asphericity << " "  << asphericity_i << " " << shape_anisotropy << " " << shape_anisotropy_i << "\n";
            data_start_mass = N_1;
        }
    
        if ( config_mass_elapsed >= config_mass && time < max_time && !no_config)
        {
            print_config();
            config_start_mass = N_1;
        }
        //count_steps++;
    }

    data_stream << N_1 << " " << N_2 << " "  << count_steps << " " << time << " " << rho << " " << asphericity << " "  << asphericity_i << " " << shape_anisotropy << " " << shape_anisotropy_i << "\n";
    //print_full_config();

    full_config_stream.close();
    if (!no_config)
    {
        config_stream.close();
    }
    data_stream.close();

}

void blocked_crystal::kmc_run_mass_up(unsigned int N_1_end, int N_1_end_low, int config_rec, int data_rec, unsigned int seed, unsigned int seed_2, bool no_config)
{
    prepare_run(no_config, seed);
    contgen.seed(seed);
    liggen.seed(seed_2);
    //std::cout << lig_step <<"\n";
    contdis.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
    
    int config_mass = (N_1_end - N_1) / config_rec;
    int data_mass = (N_1_end - N_1) / data_rec;

    int config_start_mass = N_1;
    int data_start_mass = N_1;

    int data_mass_elapsed = 0;
    int config_mass_elapsed = 0;

    if (!no_config)
    {
        print_config();
    }
    count_steps = 0;
    while ( N_1 < N_1_end && time < max_time && count_steps < N_1_end_low)
    {
        kmc_step();

        if (must_update_max_size)
        {
            std::cerr << "Stopping trial " << trial_number << " due to reached boundary\n";
            break;
        }
        if (stationary)
        {
            std::cerr << "Stopping trial " << trial_number << " due to reached stationary state\n";
            break;   
        }

        data_mass_elapsed = N_1 - data_start_mass;
        config_mass_elapsed = N_1 - config_start_mass;
        
        if ( data_mass_elapsed >= data_mass && time < max_time )
        {
            data_stream << N_1 << " " << N_2 << " "  << count_steps << " " << time << " " << rho << " " << asphericity << " "  << asphericity_i << " " << shape_anisotropy << " " << shape_anisotropy_i << "\n";
            //std::cout << N_1 << " " << N_2 << " " << count_steps << " " << time << " " << asphericity << " "  << asphericity_i << " " << shape_anisotropy << " " << shape_anisotropy_i << "\n";
            data_start_mass = N_1;
        }
    
    
        if ( config_mass_elapsed >= config_mass && time < max_time && !no_config)
        {
            print_config();
            config_start_mass = N_1;
        }
        //count_steps++;
    }

    data_stream << N_1 << " " << N_2 << " "  << count_steps << " " << time << " " << rho << " " << asphericity << " "  << asphericity_i << " " << shape_anisotropy << " " << shape_anisotropy_i << "\n";
    
    //print_full_config();
    
    full_config_stream.close();
    if (!no_config)
    {
        config_stream.close();
    }
    data_stream.close();
}

void blocked_crystal::kmc_run_steps(unsigned int num_steps, int config_rec, int data_rec, unsigned int seed, bool no_config)
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
    while (count_steps < num_steps)
    {
        kmc_step();
        if ( count_steps % data_steps == 0 )
        {
            data_stream << count_steps << " " << time << " " << num_particles << " " << energy << "\n";
        }
        if ( count_steps % config_steps == 0 && !no_config)
        {
            print_config();
        }
    }

    data_stream << count_steps << " " << time << " " << num_particles << " " << energy << "\n";
    data_stream.close();
    
    if (!no_config)
    {
        config_stream.close();
    }

    print_minimal_config();
    full_config_stream.close();
}

void blocked_crystal::print_config()
{
    int val, num_neighbors;
    unsigned int lattice_index;
    int i, j, k;
    double x, y, z;

    for (unsigned int r = 1; r <= size_active_sites; r++)
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
            
                num_neighbors = neighb[lattice_index];
                config_stream << val << " " << x << " " << y << " " << z << " " << num_neighbors << "\n";
            }
        }
    }
    config_stream << N_1 << " " << N_2 << " " << count_steps << " " << time << " " << energy << "\n";
}
#endif /* BCSDbep_h */
