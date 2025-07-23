#ifndef GEOMETRY_H
#define GEOMETRY_H

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

struct Surface{
    // Contains the parameters to define a plane, sphere, cylinder, truncated cylinder, cuboid, or infinite prism
    double ax;
    double r;
    double cx;
    double cy;
    double cz;
    double plane1;
    double plane2;
    double plane3;
    double plane4;
    double plane5;
    double plane6;
    string type;
};

struct Cell{
    vector<string> def;
    string type, childname, material;
};
struct Lattice{
    double pitch, xbord1, xbord2;
    vector<vector<string>> unimat;
};

struct Universe{
    vector<string> cells;
};

bool testSurf(unordered_map<string, Surface> surfaces, vector<double> point, string key);

bool cellsearch(unordered_map<string, Surface> surfacess,
    unordered_map<string, Cell> cells,
    vector<double> point, string cellname);

string universesearch(unordered_map<string, Surface> surfaces,
    unordered_map<string, Cell> cells,
    unordered_map<string, Lattice> lattices,
    unordered_map<string, Universe> universes, vector<double> point);

#endif