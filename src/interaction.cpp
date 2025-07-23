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
using namespace std;

/**
 * Neutron interaction physics module
 * 
 * Contains functions for handling various neutron interactions:
 * - Scattering (elastic)
 * - Fission
 * - Radiative capture
 * 
 * Also includes routines for sampling energy distributions
 * and selecting interaction targets.
 */

/**
 * Determines the type of neutron interaction
 * based on reaction cross-section probabilities
 */
int interact(const vector<double>& csvec) {
    // Extract cross-sections for different reaction types
    double csscatt = csvec[0]; // Scattering cross-section
    double csfiss = csvec[1];  // Fission cross-section
    double cscapt = csvec[2];  // Capture cross-section

    // Calculate total cross-section
    double cstot = csscatt + csfiss + cscapt;

    // Sample random number to determine interaction type
    double eps = randomnum();

    // Return interaction type based on relative probabilities
    if (eps * cstot < csscatt) {
        return 0; // Scattering interaction (elastic)
    } else if (eps * cstot < csscatt + cscapt) {
        return 2; // Capture interaction (radiative)
    } else {
        return 1; // Fission interaction
    }
}

/**
 * Samples which target nuclide in a material undergoes interaction
 * 
 * @param cstot Vector of macroscopic cross-sections for each nuclide
 * @return Index of the selected nuclide
 */
int target(vector<double> cstot){
    // Generate random number
    double eps = randomnum();
    
    // Create cumulative sum vector for cross-sections
    vector<double> sums(cstot.size()+1);
    double sum = 0;
    sums[0] = 0;
    int target;
    
    // Calculate cumulative cross-sections
    for (int i = 0; i < cstot.size(); i++){
        sum = cstot[i] + sum;
        sums[i+1] = sum;
    }
    
    // Scale random number to total cross-section
    double epsxtot = sums[sums.size() - 1]*eps;
    
    // Find which nuclide the random number corresponds to
    for (int i = 0; i < sums.size()-1; i++){
        if (sums[i] < epsxtot && sums[i+1] >=epsxtot){
            target = i;
        } 
    }
    return target;
} 

/**
 * Samples neutron energy from Maxwell-Boltzmann distribution
 * 
 * Used for both thermal scattering and fission neutron energy sampling
 * 
 * @param Temp Temperature parameter in energy units (MeV)
 * @return Sampled energy in MeV
 */
double sampleEn(double Temp){
    // Marsaglia method for sampling from Maxwell-Boltzmann distribution
    double R = 2;
    double eps1;
    double eps2;
    
    // Generate random point in unit circle
    while (R > 1){
        eps1 = randomnum();
        eps2 = randomnum();
        R = (eps1*eps1)+(eps2*eps2);
    }
    
    // Additional random numbers for energy sampling
    double eps3 = randomnum();
    double eps4 = randomnum();
    
    // Calculate energy from Maxwell-Boltzmann distribution
    double newEn = -Temp*((eps1*eps1)*log(eps3)/R + log(eps4));
    return newEn;
}

/**
 * Handles fission reaction by creating secondary neutrons
 * 
 * @param newGen Vector to store the produced fission neutrons
 * @param particle Original neutron causing fission
 * @param nu Average number of neutrons produced per fission
 * @param Tfiss Temperature parameter for fission neutron spectrum (MeV)
 */
void fission(vector<Neutron>& newGen, Neutron particle ,double nu, double Tfiss){
    double pi = 2*acos(0.0);
    Neutron n;
    
    // Determine number of secondary neutrons using integer part of nu
    int secn = floor(nu);
    
    // Probabilistically add one more neutron to account for fractional part
    double eps = randomnum();
    if (eps < nu-secn){
        secn = secn + 1;
    }
    
    // Create each secondary neutron
    for (int i = 0; i < secn; i++){
        // Sample energy from fission spectrum (Maxwell-Boltzmann)
        double newEn = sampleEn(Tfiss);
        
        // Sample direction isotropically
        double thetaN = pi*randomnum();
        double phiN = 2*pi*randomnum();
        
        // Set neutron properties
        n.En = newEn;
        n.u = sin(thetaN)*cos(phiN);  // Direction cosine x
        n.v = sin(thetaN)*sin(phiN);  // Direction cosine y
        n.w = cos(thetaN);            // Direction cosine z
        n.x = particle.x;             // Position is same as parent
        n.y = particle.y;
        n.z = particle.z;
        
        // Add to next generation
        newGen.push_back(n);
    }
}

/**
 * Handles elastic scattering kinematics in the center of mass frame
 * 
 * Models both free gas thermal scattering for low energies and
 * standard elastic scattering for higher energies
 * 
 * @param particle Neutron undergoing scattering (updated in-place)
 * @param Tmedium Medium temperature in energy units (MeV)
 * @param Mtarg Mass of target nucleus (kg)
 */
void scattering(Neutron& particle, double Tmedium, double Mtarg) {
    double pi = 2*acos(0.0);
    double nmass = 1.6749286e-27; // Neutron mass [kg]
    
    // Convert neutron parameters to velocity components
    double sqrt2_nmassEn = sqrt(2*particle.En/ nmass);
    double vLx = sqrt2_nmassEn*particle.u;  // Lab velocity x-component
    double vLy = sqrt2_nmassEn*particle.v;  // Lab velocity y-component
    double vLz = sqrt2_nmassEn*particle.w;  // Lab velocity z-component
    
    // Free gas model threshold energy
    double Efg = 200e-6; // Thermal energy threshold [MeV]
    double Et = 0, VLx = 0, VLy = 0, VLz = 0;

    // For thermal neutrons, use free gas model with thermal motion of target
    if (particle.En < Efg){
        // Sample target nucleus velocity direction isotropically
        double thetaT = pi*randomnum();
        double phiT = 2*pi*randomnum();
        
        // Sample target energy from Maxwell-Boltzmann distribution
        Et = sampleEn(Tmedium);
        
        // Calculate target velocity components
        double sqrt2_Mtarg = sqrt(2 / Mtarg);
        VLx = sqrt2_Mtarg * sqrt(Et) * sin(thetaT) * cos(phiT);
        VLy = sqrt2_Mtarg * sqrt(Et) * sin(thetaT) * sin(phiT);
        VLz = sqrt2_Mtarg * sqrt(Et) * cos(thetaT);
    }
    
    // Calculate mass ratio (target/neutron)
    double A = Mtarg/nmass;

    // Calculate center-of-mass velocity components
    double VCMx = (vLx + A*VLx)/(1+A);
    double VCMy = (vLy + A*VLy)/(1+A);
    double VCMz = (vLz + A*VLz)/(1+A);

    // Convert to center-of-mass frame
    double vcx = vLx - VCMx;  // CM velocity x-component
    double vcy = vLy - VCMy;  // CM velocity y-component
    double vcz = vLz - VCMz;  // CM velocity z-component

    // Sample scattering angle in CM frame (isotropic for s-wave scattering)
    double mu = 1-2*randomnum();         // Cosine of polar angle
    double phi = 2*pi*randomnum();       // Azimuthal angle
    double sin_theta = sqrt(1 - mu*mu);  // Sine of polar angle

    // Calculate neutron speed in CM frame (conserved in elastic scattering)
    double vc_mag_sq = vcx * vcx + vcy * vcy + vcz * vcz;
    double vc_mag = sqrt(vc_mag_sq);

    // Calculate new velocity components in CM frame
    double vcnewx = vc_mag * sin_theta * cos(phi);
    double vcnewy = vc_mag * sin_theta * sin(phi);
    double vcnewz = vc_mag * mu;

    // Transform back to laboratory frame
    vLx = vcnewx + VCMx;
    vLy = vcnewy + VCMy;
    vLz = vcnewz + VCMz;

    // Calculate new neutron energy
    double vLmag_sq = vLx * vLx + vLy * vLy + vLz * vLz;
    particle.En = 0.5 * vLmag_sq * nmass;  // E = 1/2 m v²
    
    // Update neutron direction cosines
    double vLmag = sqrt(vLmag_sq);
    particle.u = vLx / vLmag;  // Direction cosine x
    particle.v = vLy / vLmag;  // Direction cosine y
    particle.w = vLz / vLmag;  // Direction cosine z
}
