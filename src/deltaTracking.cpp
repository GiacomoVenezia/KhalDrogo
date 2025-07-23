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

using namespace std;

/**
 * Implements the Delta Tracking algorithm for neutron transport
 *
 * Delta tracking (also called Woodcock tracking) is an efficient Monte Carlo
 * method that avoids explicit boundary crossing calculations by using
 * virtual collisions and rejection sampling based on the majorant cross-section.
 *
 * @param materials Vector of all material data
 * @param surfaces Map of all surfaces in the geometry
 * @param cells Map of all cells in the geometry
 * @param lattices Map of all lattices in the geometry
 * @param universes Map of all universes in the geometry
 * @param oldGen Vector of neutrons in current generation
 * @param newGen Vector to store fission neutrons for next generation
 * @param nfiss Counter for fission reactions
 * @param nscatt Counter for scattering reactions
 * @param nabs Counter for absorption reactions
 * @param nu Average number of neutrons produced per fission
 */
void deltaTracking(vector<MaterialData>& materials, unordered_map<string, Surface>& surfaces, unordered_map<string, Cell> cells,
    unordered_map<string, Lattice> lattices, unordered_map<string, Universe> universes, vector<Neutron>& oldGen, vector<Neutron>& newGen, int& nfiss , int& nscatt, int& nabs, double& nu) {
        // Physical constants
        double NA = 6.0221409e23;      // Avogadro's number [atoms/mol]
        double boltz = 8.6167e-11;     // Boltzmann constant [MeV/K]
        double Tmedium = 300*boltz;    // Medium temperature [MeV]
        double Tfiss = 1.2895;         // Average fission neutron energy [MeV]
        
        // Algorithm variables
        double l;                      // Path length
        int targ, interaction;
        double XS, M_target;           // Cross-section and target mass
        vector<double> xsvec;          // Reaction cross-sections
        
        // Process each neutron in the current generation
        for (Neutron& particle: oldGen){
            bool alive = true;
            // Continue tracking until neutron is absorbed, causes fission, or leaks
            while (alive){

                // Step 1: Calculate majorant cross-section (maximum across all materials)
                vector<double> XSmax;
                for (MaterialData& matData : materials) {
                    getCS(matData.nuclides, particle.En, matData);
                    XSmax.push_back(matData.matXS);
                }
                
                string material;
                double XSmaj = *max_element(XSmax.begin(), XSmax.end());

                // Step 2: Sample distance to next collision using majorant
                l = -log(randomnum())/XSmaj;
                
                // Step 3: Move neutron to collision site
                particle.x += l * particle.u;
                particle.y += l * particle.v;
                particle.z += l * particle.w;
                
                // Step 4: Find material at collision site
                material = universesearch(surfaces, cells, lattices, universes, {particle.x, particle.y, particle.z});
                
                // Check if neutron has leaked from the system
                if (material == "void") {
                    break;  // Exit tracking loop for this neutron
                }
                
                // Find properties of the material at collision site
                for (MaterialData& matData : materials){
                    if (matData.name == material){
                        XS = matData.matXS;  // Actual cross-section of the material
                        
                        // Sample target nuclide within material
                        targ = target(matData.matXSvec);
                        
                        // Get cross-sections for different reaction types
                        xsvec = {matData.nuclides[targ].csscatt, matData.nuclides[targ].csfiss, matData.nuclides[targ].cscapt};
                        
                        // Calculate target mass for scattering kinematics
                        M_target = (matData.nuclides[targ].AW*1e-3)/NA;
                        break;
                    }
                }
                
                // Step 5: Apply Delta Tracking rejection technique
                // Accept collision with probability XS/XSmaj (real vs. majorant)
                if (randomnum() < XS/XSmaj){
                    // Sample interaction type (0=scattering, 1=fission, 2=capture)
                    interaction = interact(xsvec);

                    if (interaction == 0) {
                        // Elastic scattering: update neutron direction and energy
                        scattering(particle, Tmedium, M_target);
                        nscatt++;
                    } else if (interaction == 1) {
                        // Fission: create new neutrons for next generation
                        fission(newGen, particle, nu, Tfiss);
                        nfiss++;
                        alive = false;  // Current neutron is absorbed
                    } else if (interaction == 2) {
                        // Radiative capture: neutron is absorbed
                        nabs++;
                        alive = false;
                    }
                }
                // If rejected (randomnum() >= XS/XSmaj), consider as virtual collision
                // and continue with the same neutron (delta tracking algorithm)
            }
        }
}
