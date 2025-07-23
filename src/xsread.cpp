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

using namespace std;

/**
 * Reads neutron cross-section data from a file for a specific nuclide
 * 
 * @param name The name of the nuclide (e.g., "U235", "Pu239")
 * @return NuclideData structure containing all cross-section information
 */
NuclideData readXSfile(string name){
    NuclideData data;
    string filename = "XSdata/" + name + ".dat";
    ifstream file(filename);
    string line, dummy;
    int lineid = 0;
    int nlines, MT;
    data.name = name;
    
    // Read the cross-section data file line by line
    while (getline(file, line)){
        istringstream ss(line);
        
        // First line contains atomic mass number (A) and atomic weight (AW)
        if (lineid == 0){
            ss >> dummy >> dummy >> data.A >> data.AW;
            lineid++;
            continue;
        }
        
        // Second line contains number of header lines to skip
        if (lineid == 1){
            ss >> nlines;
            for (int i = 0; i < nlines; i++){
                getline(file, line);
            }
            lineid++;
            continue;
        }
        
        // Read different reaction types by MT number
        // MT: identifier for reaction type in ENDF format
        if (ss >> MT >> dummy >> nlines) {
            // MT=2: Elastic scattering cross section
            if (MT == 2){
                vector<vector<double>> scatt(2); // [0] = energies, [1] = cross-sections
                for (int i = 0; i < nlines; i++){
                    getline(file, line);
                    istringstream data_ss(line);
                    double energy, value;
                    data_ss >> energy >> value;
                    scatt[0].push_back(energy);
                    scatt[1].push_back(value);
                }
                data.scattering = scatt;
            // MT=18: Fission cross section
            } else if (MT == 18){ 
                vector<vector<double>> fiss(2);
                for (int i = 0; i < nlines; i++){
                    getline(file, line);
                    istringstream data_ss(line);
                    double energy, value;
                    data_ss >> energy >> value;
                    fiss[0].push_back(energy);
                    fiss[1].push_back(value);
                }
                data.fission = fiss;
            // MT=102: Radiative capture cross section
            } else if (MT == 102){
                vector<vector<double>> capt(2);
                for (int i = 0; i < nlines; i++){
                    getline(file, line);
                    istringstream data_ss(line);
                    double energy, value;
                    data_ss >> energy >> value;
                    capt[0].push_back(energy);
                    capt[1].push_back(value);
                }
                data.capture = capt;
            // Skip other reaction types
            } else {
                for (int i = 0; i < nlines; i++){
                    getline(file, line);
                }
            }
        }
    }
    file.close();
    return data;
}

/**
 * Calculates microscopic and macroscopic cross sections for a given neutron energy
 * 
 * @param nuclides Vector of nuclides in the material
 * @param en Neutron energy
 * @param matData Material data to store the calculated cross sections
 */
void getCS(vector<NuclideData>& nuclides, double en, MaterialData& matData) {
    matData.matXS = 0.0; // Initialize total macroscopic cross-section
    matData.matXSvec.clear(); // Clear the vector to avoid accumulation from previous calls
    
    // Loop through each nuclide in the material
    for(auto& nuclide: nuclides){
        // Calculate microscopic cross sections through interpolation
        nuclide.csscatt = interpolate(nuclide.scattering, en);
        
        // Handle fission (not all nuclides are fissionable)
        if (nuclide.fission.empty()){
            nuclide.csfiss = 0.0; // Default value if no fission data is available
        } else {
            nuclide.csfiss = interpolate(nuclide.fission, en);
        }
        
        nuclide.cscapt = interpolate(nuclide.capture, en);
        
        // Calculate total macroscopic cross section for this nuclide
        nuclide.xs_macr_tot = (nuclide.csscatt + nuclide.csfiss + nuclide.cscapt) * nuclide.atomic_density;
        
        // Add to material total cross section
        matData.matXS += nuclide.xs_macr_tot;
        matData.matXSvec.push_back(nuclide.xs_macr_tot);
    }
}

/**
 * Calculates atomic densities for all nuclides in a material
 * 
 * @param nuclides Vector of nuclides in the material
 * @param density Material density in g/cm³
 * @param enr Enrichment factor for U-235 (used only for uranium fuel)
 */
void getAtomicDensity(vector<NuclideData>& nuclides, double density, double enr) {
    const double NA = 6.0221409e23; // Avogadro's number
    
    // Calculate total molecular mass of the material
    double total_molecular_mass = 0.0;
    for (auto& nuclide : nuclides) {
        total_molecular_mass += nuclide.AW * nuclide.nnuclides;
    }
    
    // Calculate atomic density for each nuclide (atoms/barn-cm)
    for (auto& nuclide : nuclides) {
        // Special handling for uranium isotopes to account for enrichment
        if (nuclide.name == "U235"){
            // Apply enrichment factor to U-235
            nuclide.atomic_density = (density * NA * nuclide.nnuclides * enr / total_molecular_mass * 1e-24);
        } else if (nuclide.name == "U238") {
            // Apply depletion factor (1-enrichment) to U-238
            nuclide.atomic_density = (density * NA * nuclide.nnuclides * (1 - enr) / total_molecular_mass * 1e-24);
        } else {
            // Standard formula for other nuclides
            nuclide.atomic_density = (density * NA * nuclide.nnuclides) / total_molecular_mass * 1e-24;
        }
    }
}