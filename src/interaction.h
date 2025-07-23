#ifndef INTERACTION_H
#define INTERACTION_H

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

using namespace std;

struct Neutron {
    double En; // energy
    double u, v, w; // direction
    double x, y, z; // position
};

int interact(const vector<double>& csvec);
int target(vector<double> cstot);
double sampleEn(double Temp);
void fission(vector<Neutron>& newGen, Neutron particle ,double nu, double Tfiss);
void scattering(Neutron& particle, double Tmedium, double Mtarg);

#endif