//  helper.h
//  Auxiliary function for ease of operating the nanocrystal class
//  Created by Sam Oaks-Leaf.

#ifndef helper_h
#define helper_h
#include<vector>
#include<cmath>
#include<array>

template <typename T>
void read_in_vec(std::vector<T>& vec, std::ifstream& in_stream)
{
    // Reads in vector from in_stream and stores it in vec.
    int int_size = sizeof(int);
    int length;
    in_stream.read(reinterpret_cast<char *>(&length), int_size);
    vec.reserve(length);
    T temp_value;
    int multiplicity;
    int T_size = sizeof(T);
    if (in_stream.is_open())
    {
        int count = 0;
        while (count < length)
        {
            in_stream.read(reinterpret_cast<char *>(&temp_value), T_size);
            in_stream.read(reinterpret_cast<char *>(&multiplicity), int_size);
            vec.insert(vec.end(), multiplicity, temp_value);
            count += multiplicity;
        }  
    }
}

template <typename T>
void write_out_vec(std::vector<T>& vec, std::ofstream& out_stream)
{
    // Writes vector vec to file through out_stream
    int int_size = sizeof(int);
    int length = vec.size();
    out_stream.write(reinterpret_cast<const char *>(&length), int_size);
    
    int T_size = sizeof(T);
    
    if (out_stream.is_open() && length > 0)
    {
        T temp_old = vec[0];
        out_stream.write(reinterpret_cast<const char *>(&temp_old), T_size);
        
        int count_multiplicity = 0;
        
        for ( const auto temp_new : vec )
        {
            
            if ( temp_new == temp_old )
            {
                count_multiplicity++;
            }
            
            else
            {
                out_stream.write(reinterpret_cast<const char *>(&count_multiplicity), int_size);
                count_multiplicity = 1;
                
                out_stream.write(reinterpret_cast<const char *>(&temp_new), T_size);
                temp_old = temp_new;
            }
        }
        out_stream.write(reinterpret_cast<const char *>(&count_multiplicity), int_size);
    }
}

inline int opp_pos(int pos)
{
    //Returns index of lattice vector pointing opposite to vector at index pos in the array "order"
    return pos + 1 - 2 * (pos % 2);
}

inline int nnmod( int i, int n )
{
    // Returns i mod n such that result is always non negative
    return (i % n + n) % n;
}

std::array<int, 3> get_ijk(int lattice_index, int max_size)
{
    // Gets crystal lattice vector from lattice index
    int i = lattice_index / (max_size * max_size);
    int j = (lattice_index / max_size) % max_size;
    int k = lattice_index % max_size;
    std::array<int, 3> r = {i, j, k};
    return r;
}

std::array<double, 3> get_xyz(int lattice_index, int max_size)
{
    // Gets cartesian lattice vector from lattice index
    std::array<int, 3> ijk = get_ijk(lattice_index, max_size);
    double x = (ijk[1] + ijk[2]) / 2.0;
    double y = (ijk[0] + ijk[2]) / 2.0;
    double z = (ijk[0] + ijk[1]) / 2.0;
    std::array<double, 3> r = {x, y, z};
    return r;
}

// A few functions used for initializing in certain shapes
bool sphere(int lattice_index, int max_size, double r)
{
    std::array<double, 3> vec = get_xyz(lattice_index, max_size);
    double x = vec[0] - max_size/2.0;
    double y = vec[1] - max_size/2.0;
    double z = vec[2] - max_size/2.0;

    if (x * x + y * y + z * z <= r * r)
    {
        return true;
    }
    return false;
}

bool octo(int lattice_index, int max_size, double r)
{
    std::array<double, 3> vec = get_xyz(lattice_index, max_size);
    double x = vec[0] - max_size/2.0;
    double y = vec[1] - max_size/2.0;
    double z = vec[2] - max_size/2.0;

    if (abs(x) + abs(y) + abs(z) <= r)
    {
        return true;
    }

    return false;
}

bool cube(int lattice_index, int max_size, double r)
{
    std::array<double, 3> vec = get_xyz(lattice_index, max_size);
    if (abs(vec[0] - max_size/2.0) < r && abs(vec[1] - max_size/2.0) < r && abs(vec[2] - max_size/2.0) < r)
    {
        return true;
    }

    return false;
}

bool cubocto(int lattice_index, int max_size, double r)
{
    return (cube(lattice_index, max_size, r) * octo(lattice_index, max_size, r));
}

#endif /* helper_h */
