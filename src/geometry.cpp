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
#include <geometry.h>
#include <common.h>

using namespace std;

/**
 * Tests if a point is inside a geometric surface
 *
 * @param surfaces Map of all defined surfaces in the problem
 * @param point 3D coordinates of the point to test [x,y,z]
 * @param key Identifier of the surface to test against
 * @return true if point is inside the surface, false otherwise
 */
bool testSurf(unordered_map<string, Surface> surfaces, vector<double> point, string key){
    double c1, c2, c3, r;
    int ax;
    
    // Check point against different surface types
    if (surfaces[key].type == "PLANE"){
        // Plane defined by x_i = constant, where i is the axis (0=x, 1=y, 2=z)
        ax = surfaces[key].ax;
        return(point[ax] <= surfaces[key].plane1);
    } else if (surfaces[key].type == "SPHERE"){
        // Sphere defined by (x-x0)^2 + (y-y0)^2 + (z-z0)^2 < r^2
        c1 = surfaces[key].cx;
        c2 = surfaces[key].cy;
        c3 = surfaces[key].cz;
        return(surfaces[key].r > sqrt(pow(point[0]-c1,2) + pow(point[1]-c2,2) + pow(point[2]-c3,2)));
    } else if (surfaces[key].type == "CYLINDER"){
        // Infinite cylinder aligned with one axis
        ax = surfaces[key].ax;
        c1 = surfaces[key].cx;
        c2 = surfaces[key].cy;
        return(surfaces[key].r >= sqrt(pow(c1-point[ax-2],2)+pow(c2-point[ax-1],2)));
    } else if (surfaces[key].type == "CYL_TR"){
        // Truncated cylinder (bounded by planes perpendicular to its axis)
        ax = surfaces[key].ax;
        c1 = surfaces[key].cx;
        c2 = surfaces[key].cy;
        r = surfaces[key].r;
        return(r >= sqrt(pow(c1-point[ax-2],2)+pow(c2-point[ax-1],2)) && point[ax-1] <= surfaces[key].plane1 && point[ax-1] >= surfaces[key].plane2);
    } else if (surfaces[key].type == "CUBOID"){
        // Rectangular prism defined by min/max points
        return(point[0] >= surfaces[key].plane1 && point[0] <= surfaces[key].plane4 &&
               point[1] >= surfaces[key].plane2 && point[1] <= surfaces[key].plane5 &&
               point[2] >= surfaces[key].plane3 && point[2] <= surfaces[key].plane6);
    } else if (surfaces[key].type == "INF_PRISM"){
        // Infinite prism defined by 4 planes
        ax = surfaces[key].ax;
        return(point[ax-2] >= surfaces[key].plane1 && point[ax-2] <= surfaces[key].plane3 &&
               point[ax-1] >= surfaces[key].plane2 && point[ax-1] <= surfaces[key].plane4);
    }
    
    // Default return if surface type not recognized
    return false;
}

/**
 * Checks if a point is inside a cell defined by surface intersections
 *
 * @param surfaces Map of all defined surfaces in the problem
 * @param cells Map of all defined cells in the problem
 * @param point 3D coordinates of the point to test [x,y,z]
 * @param cellname Identifier of the cell to test
 * @return true if point is inside the cell, false otherwise
 */
bool cellsearch(unordered_map<string, Surface> surfaces,
    unordered_map<string, Cell> cells,
    vector<double> point, string cellname){

    // Get the cell definition (list of in/out operators and surfaces)
    vector<string> intersections = cells[cellname].def;
    int n = 0;
    
    // Check each surface condition (in/out) in the cell definition
    while (n < intersections.size()/2){
        string newK = intersections[n*2+1];  // Surface ID
        bool c = testSurf(surfaces, point, newK);  // Test if point is inside surface
        
        // Invert result if the operator is "out"
        if (intersections[n*2] == "out"){
            c = !c;
        }
        
        // If any condition fails, the point is not in the cell
        if (c == 0){
            break;
        }
        n++;
    }
    
    // If all conditions passed, the point is inside the cell
    if (n == intersections.size()/2){
        return true;
    }
    
    return false;
}

/**
 * Finds the material at a given point through recursive universe traversal
 * 
 * Traverses the geometry hierarchy (universes, cells, lattices) to determine
 * the material at a specific point in 3D space
 *
 * @param surfaces Map of all defined surfaces
 * @param cells Map of all defined cells
 * @param lattices Map of all defined lattices
 * @param universes Map of all defined universes
 * @param point 3D coordinates of the point to find material for [x,y,z]
 * @return Name of the material at the given point, or "void" if outside domain
 */
string universesearch(unordered_map<string, Surface> surfaces,
    unordered_map<string, Cell> cells,
    unordered_map<string, Lattice> lattices,
    unordered_map<string, Universe> universes, vector<double> point){
    string material;
    string univ = "U1";  // Start from top-level universe
    bool flag = true, flag2 = true;
    
    // Continue until material is found or point is determined to be outside domain
    while (flag){
        vector<string> cellsvec = universes[univ].cells;
        // Check each cell in current universe
        for (int ii=0; ii < cellsvec.size(); ii++){
            flag2 = cellsearch(surfaces, cells, point, cellsvec[ii]);
            
            // If point is in cell and cell contains material, we found the answer
            if (flag2 && cells[cellsvec[ii]].type == "MATERIAL"){
                flag = false;
                material = cells[cellsvec[ii]].material;
                break;
            // If point is in cell and cell contains another universe, continue search in that universe
            } else if (flag2 && cells[cellsvec[ii]].type == "UNIV"){
                univ = cells[cellsvec[ii]].childname;
                break;
            // If point is in cell and cell contains a lattice, find the specific universe in the lattice
            } else if (flag2 && cells[cellsvec[ii]].type == "LATTICE"){
                // Get lattice parameters
                double pitch = lattices[cells[cellsvec[ii]].childname].pitch;
                
                // Transform point to lattice local coordinates
                point = {point[0]-lattices[cells[cellsvec[ii]].childname].xbord1, 
                         point[1]-lattices[cells[cellsvec[ii]].childname].xbord1, 
                         point[2]};
                
                // Find lattice indices
                double i = floor(point[0]/pitch), j = floor(point[1]/pitch);
                
                // Transform point to cell local coordinates
                point = {point[0]-i*pitch, point[1]-j*pitch, point[2]};
                
                // Get universe from lattice matrix
                univ = lattices[cells[cellsvec[ii]].childname].unimat[i][j];
                break;
            // If point is not in any cell, it's outside the domain
            } else if (!flag2 && ii == cellsvec.size()-1){
                flag = false;
                material = "void";  // Return "void" for points outside geometry
                break;
            }
        }
    }
    return material;
}