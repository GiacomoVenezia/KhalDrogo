#ifndef COMMON_H
#define COMMON_H

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
double randomnum();
double interpolate(vector<vector<double>>& table, double en);
#endif