#ifndef DELTA_TRACKING_H
#define DELTA_TRACKING_H

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

void deltaTracking(vector<MaterialData>& materials, unordered_map<string, Surface>& surfaces, unordered_map<string, Cell> cells,
    unordered_map<string, Lattice> lattices, unordered_map<string, Universe> universes, vector<Neutron>& oldGen, vector<Neutron>& newGen, int& nfiss , int& nscatt, int& nabs, double& nu);
#endif