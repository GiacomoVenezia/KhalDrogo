#ifndef FILEPARSING_H
#define FILEPARSING_H

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

void fileparsing(string filename, 
                 unordered_map<string, Surface>& surfaces,
                 unordered_map<string, Cell>& cells,
                 unordered_map<string, Lattice>& lattices,
                 unordered_map<string, Universe>& universes,
                 vector<MaterialData>& materials, int& nn, Neutron& source);

#endif