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
#include <xsread.h>
#include <geometry.h>
#include <interaction.h>
#include <deltaTracking.h>
#include <fileparsing.h>

/**
 * KhalDrogo - Monte Carlo Neutron Transport Code
 * 
 * Main program file implementing criticality source iteration algorithm
 * for k-effective calculations in complex reactor geometries.
 * 
 * The program uses delta tracking (Woodcock tracking) for efficient
 * neutron transport and implements a power iteration method to
 * converge the fission source distribution.
 */

using namespace std;

int main(int argc, char* argv[]) {
    // Data structures for reactor geometry and materials
    unordered_map<string, Surface> surfaces;    // Geometric surfaces (planes, cylinders, etc.)
    unordered_map<string, Cell> cells;          // Cells defined by surfaces
    unordered_map<string, Lattice> lattices;    // Regular arrangements of universes
    unordered_map<string, Universe> universes;  // Collections of cells
    vector<MaterialData> materials;             // Material compositions and properties
    
    // Simulation parameters
    int nn;          // Number of neutrons per generation
    Neutron source;  // Initial source parameters
    
    // Parse command line arguments for input file
    string filename = "input.txt";  // Default input file
    if (argc > 1) {
        filename = argv[1];  // Use file provided in command line
    }
    
    // Parse input file and load geometry and material data
    fileparsing(filename, surfaces, cells, lattices, universes, materials, nn, source);
    cout << "Input file loaded successfully" << endl;
    // Neutron population vectors for current and next generations
    vector<Neutron> oldGen; // Current generation neutrons
    vector<Neutron> newGen; // Next generation neutrons (from fission)

    // Initialize first neutron generation with isotropic directions
    double pi = 2*acos(0.0);
    double thetaN, phiN;
    for (int ii = 0; ii < nn; ii++){
        // Sample random direction isotropically
        thetaN = pi*randomnum();      // Polar angle [0,π]
        phiN = 2*pi*randomnum();      // Azimuthal angle [0,2π]
        
        // Calculate direction cosines
        source.u = sin(thetaN) * cos(phiN);  // x-direction
        source.v = sin(thetaN) * sin(phiN);  // y-direction
        source.w = cos(thetaN);              // z-direction
        
        // Add neutron to initial population
        oldGen.push_back(source);
    }

    // Calculate atomic densities for all materials
    for (MaterialData& matData : materials) {
        getAtomicDensity(matData.nuclides, matData.density, matData.enr);
    }
    cout << "Material properties calculated" << endl;
    
    // Criticality calculation parameters
    int maxCycles = 60;      // Total number of neutron generations to simulate
    int skipCycles = 10;     // Number of initial cycles to skip (source convergence)
    double k = 1;            // Multiplication factor for current cycle
    double keff = 0;         // Effective multiplication factor
    int active_cycles = 0;   // Number of active cycles used for k-effective averaging
    double keff_sum = 0.0;   // Sum of k-effective values for averaging
    vector<double> keff_values;  // Optional storage for statistical analysis
    double nu = 2.43;        // Initial estimate of neutrons per fission

    // Main Monte Carlo criticality power iteration loop
    for (int cycle = 0; cycle < maxCycles; cycle++) {
        // Clear next generation neutron population before starting new cycle        
        newGen.clear();
        
        // Track all neutrons and record interaction statistics
        int nfiss = 0;    // Number of fissions in this cycle
        int nscatt = 0;   // Number of scattering events
        int nabs = 0;     // Number of absorption events
        
        // Display progress information
        cout << "Progress " << (cycle+1)*100.0/maxCycles << "%" << endl;
        
        // Perform delta tracking for all neutrons in current generation
        // This will track each neutron until it is absorbed or leaves the system
        // Fission events generate neutrons for the next generation
        deltaTracking(materials, surfaces, cells, lattices, universes, 
                     oldGen, newGen, nfiss, nscatt, nabs, nu);

        // Calculate multiplication factor based on population ratio
        k = static_cast<double>(newGen.size()) / static_cast<double>(oldGen.size());
        
        // Alternative k-effective calculation based on fission reactions
        keff = nfiss*2.43/oldGen.size();
        
        // Adjust nu value to stabilize population (variance reduction technique)
        nu = nu/k;

        // Skip initial cycles to allow source distribution to converge
        if (cycle >= skipCycles) {
            // Accumulate k-effective values for active cycles
            keff_sum += keff;
            active_cycles++;
            keff_values.push_back(keff);  // Store for statistical analysis
        }

        // Safety check - terminate if all neutrons are lost
        if (newGen.empty()) {
            cout << "No neutrons left after cycle " << cycle << endl;
            break; // Exit if no neutrons are lost (subcritical system)
        }
        
        // Update current generation with neutrons from the new generation
        swap(oldGen, newGen);

    }
    // Calculate final k-effective average from active cycles
    double avg_keff = keff_sum / active_cycles;
    
    // Calculate standard deviation of k-effective (if enough active cycles)
    double std_dev = 0.0;
    if (active_cycles > 1) {
        double sum_squared_diff = 0.0;
        for (double k_val : keff_values) {
            sum_squared_diff += pow(k_val - avg_keff, 2);
        }
        std_dev = sqrt(sum_squared_diff / (active_cycles - 1));
    }
    
    // Display final results
    cout << "\nFinal results:" << endl;
    cout << "Average k-effective: " << avg_keff << endl;
    if (active_cycles > 1) {
        cout << "Standard deviation: " << std_dev << endl;
    }
    cout << "Active cycles: " << active_cycles << endl;

    return 0;
}