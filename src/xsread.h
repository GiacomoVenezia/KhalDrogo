#ifndef XSREAD_H
#define XSREAD_H

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

struct NuclideData {
    string name;
    int nnuclides; // Number of nuclides
    double A;
    double AW;
    double csscatt;
    double csfiss;
    double cscapt;
    double xs_macr_tot;
    double atomic_density;
    vector<vector<double>> scattering;
    vector<vector<double>> fission;
    vector<vector<double>> capture;
};

struct MaterialData{
    string name;
    double density; // g/cm³
    double matXS; // Macroscopic cross-section
    double enr = 0.0; // Enrichment
    vector<double> matXSvec = {}; // Vector of macroscopic cross-sections for different nuclides
    vector<NuclideData> nuclides; // List of nuclides in the material
};

NuclideData readXSfile(string name);
void getCS(vector<NuclideData>& nuclides, double en, MaterialData& matData);
void getAtomicDensity(vector<NuclideData>& nuclides, double density, double enr);
#endif