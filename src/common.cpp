#include <string>
#include <unordered_map>
#include <memory>
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <sstream>
#include <random>
#include <chrono>
#include <algorithm>
#include <common.h>

using namespace std;

/**
 * Common utility functions used throughout the Monte Carlo code
 * This file contains random number generation and interpolation functions
 */

// Initialize high-quality random number generator with time-based seed

// Use Mersenne Twister algorithm with 64-bit state size for high quality randomness
std::mt19937_64 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());
std::uniform_real_distribution<double> unif(0.0, 1.0); // Uniform distribution between 0 and 1

/**
 * Generates a random number between 0 and 1
 * 
 * Uses the Mersenne Twister algorithm (mt19937_64) for high-quality
 * pseudo-random numbers needed in Monte Carlo simulations
 * 
 * @return Double precision random number in range [0,1)
 */
double randomnum() {
    return unif(rng); 
}

/**
 * Performs linear interpolation on tabulated data
 * 
 * Used for cross-section interpolation at arbitrary neutron energies.
 * Table is structured as a 2D vector where table[0] contains energies
 * and table[1] contains corresponding cross-section values.
 * 
 * @param table 2D vector with [0]=energies, [1]=values
 * @param en Neutron energy at which to interpolate
 * @return Interpolated cross-section value
 */
double interpolate(vector<vector<double>>& table, double en) {
    auto& energies = table[0];
    auto& values = table[1]; 

    // Return zero for energies outside the table range
    if (en < energies.front() || en > energies.back()) {
        return 0.0;
    }

    // Find position in the energy grid using binary search
    auto it = lower_bound(energies.begin(), energies.end(), en);
    int idx = distance(energies.begin(), it);

    // Handle boundary cases
    if (it == energies.begin()) return values.front();
    if (it == energies.end()) return values.back(); 

    // Perform linear interpolation between adjacent points
    int i = idx - 1;
    double x0 = energies[i], x1 = energies[i + 1];
    double y0 = values[i], y1 = values[i + 1];
    return y0 + (y1 - y0) * (en - x0) / (x1 - x0);
}