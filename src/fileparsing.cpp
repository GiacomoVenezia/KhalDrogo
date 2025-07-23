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
#include <xsread.h>
#include <interaction.h>
#include <deltaTracking.h>
#include <fileparsing.h>

using namespace std;

/**
 * Parses an input file and populates the reactor geometry and material data structures
 *
 * @param filename Path to input file containing reactor description
 * @param surfaces Map to store surface definitions
 * @param cells Map to store cell definitions
 * @param lattices Map to store lattice definitions
 * @param universes Map to store universe definitions
 * @param materials Vector to store material data
 * @param nn Number of neutrons per generation
 * @param source Initial neutron source parameters
 */
void fileparsing(string filename, 
                 unordered_map<string, Surface>& surfaces,
                 unordered_map<string, Cell>& cells,
                 unordered_map<string, Lattice>& lattices,
                 unordered_map<string, Universe>& universes,
                 vector<MaterialData>& materials, int& nn, Neutron& source
                ) {
    // Open input file and read line by line
    ifstream file(filename);
    string line;
    while (getline(file, line)){
        istringstream ss(line);
        string surf, name, type;
        ss >> surf;

        // Process surface definitions (PLANE, SPHERE, CYLINDER, etc.)
        if (surf == "SURF"){
            ss >> name >> type;
            if (type == "PLANE"){
                double ax, pos;
                ss >> ax >> pos;
                Surface plane;
                plane.ax = ax;
                plane.plane1 = pos;
                plane.type = type;
                surfaces[name] = plane;
            } else if (type == "SPHERE"){
                double x, y, z, rad;
                ss >> x >> y >> z >> rad;
                Surface sphere;
                sphere.cx = x;
                sphere.cy = y;
                sphere.cz = z;
                sphere.r = rad;
                sphere.type = type;
                surfaces[name] = sphere;
            } else if (type == "CYLINDER"){
                double ax, c1, c2, rad;
                ss >> ax >> c1 >> c2 >> rad;
                Surface cylinder;
                cylinder.ax = ax;
                cylinder.cx = c1;
                cylinder.cy = c2;
                cylinder.r = rad;
                cylinder.type = type;
                surfaces[name] = cylinder;
            } else if (type == "CYL_TR"){
                double ax, c1, c2, rad, pos1, pos2;
                ss >> ax >> c1 >> c2 >> rad >> pos1 >> pos2;
                Surface tr_cyl;
                tr_cyl.ax = ax;
                tr_cyl.cx = c1;
                tr_cyl.cy = c2;
                tr_cyl.r = rad;
                tr_cyl.plane1 = pos1;
                tr_cyl.plane2 = pos2;
                tr_cyl.type = type;
                surfaces[name] = tr_cyl;
            } else if (type == "CUBOID"){
                double x1, y1, z1, x2, y2, z2;
                ss >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
                Surface cuboid;
                cuboid.plane1 = x1;
                cuboid.plane2 = y1;
                cuboid.plane3 = z1;
                cuboid.plane4 = x2;
                cuboid.plane5 = y2;
                cuboid.plane6 = z2;
                cuboid.type = type;
                surfaces[name] = cuboid;
            } else if (type == "INF_PRISM"){
                double ax, pos1, pos2, pos3, pos4;
                ss >> ax >> pos1 >> pos2 >> pos3 >> pos4;
                Surface inf_prism;
                inf_prism.ax = ax;
                inf_prism.plane1 = pos1;
                inf_prism.plane2 = pos2;
                inf_prism.plane3 = pos3;
                inf_prism.plane4 = pos4;
                inf_prism.type = type;
                surfaces[name] = inf_prism;
            }
        // Process cell definitions (regions bounded by surfaces)
        } else if (surf == "CELL"){
            ss >> name;
            Cell cell;
            string word;
            vector<string> words;
            // Store cell definition parameters (in/out surface relationships)
            while (ss >> word){
                words.push_back(word);
            }
            cell.def = words;
            getline(file, line);
            istringstream ss(line);
            ss >> type;
            cell.type = type;
            // Cell filled with material
            if (type == "MATERIAL"){
                string mat;
                ss >> mat;
                cell.material = mat;
                cells[name] = cell;
            // Cell filled with lattice structure
            } else if (type == "LATTICE"){
                double pitch;
                string univ, childname, disposition;
                vector<string> univs;
                ss >> childname >> pitch >> disposition;
                cell.childname = childname;
                cells[name] = cell;
                // Create alternating pattern of universes (checkerboard)
                if (disposition == "alternate"){
                    while (ss >> univ){
                        univs.push_back(univ);
                    }
                    getline(file, line);
                    istringstream ss(line);
                    double x1, x2;
                    ss >> x1 >> x2;
                    double width = x2-x1;
                    vector<vector<string>> univmat(width/(pitch), vector<string>(width/(pitch)));
                    // Fill matrix with alternating pattern of universes
                    for (int ii = 0; ii<width/(2*pitch); ii++){
                        for (int jj=0; jj<(x2-x1)/pitch; jj++){
                            univmat[ii*2][jj*2] = univs[0];
                            univmat[ii*2][jj*2+1] = univs[1];
                            univmat[ii*2+1][jj*2+1] = univs[1];
                            univmat[ii*2+1][jj*2] = univs[0];
                        }
                    }
                    Lattice lattice;
                    lattice.pitch = pitch;
                    lattice.unimat = univmat;
                    lattice.xbord1 = x1;
                    lattice.xbord2 = x2;
                    lattices[childname] = lattice;
                // Create uniform lattice with single universe type
                } else if (disposition == "mono"){
                    ss >> univ;
                    getline(file, line);
                    istringstream ss(line);
                    double x1, x2;
                    ss >> x1 >> x2;
                    double width = x2-x1;
                    // Fill matrix with the same universe in all positions
                    vector<vector<string>> univmat(width/(pitch), vector<string>(width/(pitch), univ));
                    Lattice lattice;
                    lattice.pitch = pitch;
                    lattice.unimat = univmat;
                    lattice.xbord1 = x1;
                    lattice.xbord2 = x2;
                    lattices[childname] = lattice;
                }
            // Cell filled with another universe
            } else if (type == "UNIV"){
                string childname, cellname;
                ss >> childname;
                cell.childname = childname;
                cells[name] = cell;
                vector<string> cellnames;
            while (ss >> cellname){
                cellnames.push_back(cellname);
            }
            Universe universe;
            universe.cells = cellnames;
            universes[childname] = universe;
            }
        // Process top-level universe definitions
        } else if (surf == "UNIV"){
            ss >> name;
            string cellname;
            vector<string> cellnames;
            // Collect all cells in this universe
            while (ss >> cellname){
                cellnames.push_back(cellname);
            }
            Universe universe;
            universe.cells = cellnames;
            universes[name] = universe;
        // Process standard material definitions
        } else if (surf == "DEFMAT"){
            string material;
            double density;
            MaterialData matData;
            ss >> material >> density;
            matData.name = material;
            matData.density = density;
            getline(file, line);
            istringstream ss(line);
            string isotope;
            int quantity;
            // Process each isotope in the material
            while (ss  >> quantity >> isotope) {
                NuclideData data = readXSfile(isotope);
                data.nnuclides = quantity;
                matData.nuclides.push_back(data);
            }
            materials.push_back(matData);
        // Process fuel material definitions with enrichment
        } else if (surf == "DEFFUEL"){
            string material;
            double density, enr;
            MaterialData matData;
            ss >> material >> density >> enr;
            matData.name = material;
            matData.density = density;
            matData.enr = enr; // Store enrichment for fuel
            getline(file, line);
            istringstream ss(line);
            string isotope;
            int quantity;
            // Process each isotope in the fuel material
            while (ss  >> quantity >> isotope) {
                NuclideData data = readXSfile(isotope);
                data.nnuclides = quantity;
                matData.nuclides.push_back(data);
                // Apply enrichment adjustments to uranium isotopes
                if (isotope == "U235") {
                    data.AW = data.AW*enr; // Adjust U235 for enrichment
                }
                if (isotope == "U238") {
                    data.AW = data.AW*(1 - enr); // Adjust U238 for depletion
                }
            }
            materials.push_back(matData);
        // Process simulation parameters
        } else if (surf == "GEN"){
            // Number of neutrons, initial energy and source position
            ss >> nn >> source.En >> source.x >> source.y >> source.z;
        }

    }
    file.close();
}