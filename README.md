# KhalDrogo - Monte Carlo Neutron Transport Code

## Description
KhalDrogo is a Monte Carlo neutron transport simulation code designed for reactor criticality analysis. The first version of this code was developed as part of nuclear engineering studies for the Monte Carlo course held by Pofessor J. Leppanen at Aalto University. The code implements a robust delta tracking algorithm for efficient neutron transport in complex geometries.

## Features
- Delta Tracking Algorithm: Efficient neutron transport through heterogeneous materials
- Constructive Solid Geometry (CSG): Hierarchical universe-based geometry with cells, surfaces, and lattices
- Cross-Section Processing**: Reads and processes neutron cross-section data from standard format files
- Physics Interactions: Models scattering, fission, and absorption reactions
- Criticality Source Simulation: Monte Carlo power iteration for fission source convergence
- Criticality Calculations: K-effective estimation with statistical error analysis

## Compilation
To compile the code:
```bash
make
```

## Usage
To run a simulation:
```bash
./drogo <input_file>
```

Example:
```bash
./drogo surfaceHWR.txt
```

## Input File Format
The input file defines:
- Universes: Hierarchical arrangement of cells
- Surfaces: Geometric boundaries (planes, cylinders, spheres, cuboids, infinite prisms and truncated cylinders)
- Cells: Regions defined by surface intersections
- Materials: Isotopic compositions and densities
- Lattices: Regular arrangements of universes
- Fuel: Isotopic composition, density and enrichment

## surfaceHWR.txt contains an example of input file simulating a 3D heavy water reactor

### Input File Keywords

#### Universe Definition
```
UNIV [universe_id] [cell_id1] [cell_id2] ...
```
Defines a universe with a unique identifier and the cells it contains.

#### Surface Definition
```
SURF [surface_id] [type] [parameters]
```
Supported surface types:
- `PLANE`: Plane - parameters: `ax pos` (ax + by = pos)
  - Example: `SURF S6 PLANE 0 0` (defines plane x = 0)
  - Example: `SURF S8 PLANE 1 0` (defines plane y = 0)

- `SPHERE`: Sphere - parameters: `x y z radius`
  - Example: `SURF S10 SPHERE 0 0 0 5.0` (sphere centered at origin with radius 5.0)

- `CYLINDER`: Cylinder - parameters: `ax c1 c2 radius`
  - where `ax` is the axis (0:x, 1:y, 2:z), c1 and c2 are center coordinates
  - Example: `SURF S4 CYLINDER 2 3 3 0.412` (cylinder along z-axis, centered at x=3,y=3 with radius 0.412)

- `CYL_TR`: Truncated cylinder - parameters: `ax c1 c2 radius pos1 pos2`
  - where pos1 and pos2 define the bounds along the cylinder axis
  - Example: `SURF S11 CYL_TR 2 0 0 5.0 0 10.0` (z-axis cylinder truncated between z=0 and z=10)

- `CUBOID`: Rectangular prism - parameters: `x1 y1 z1 x2 y2 z2`
  - Example: `SURF S1 CUBOID 5 5 5 205 205 205` (cuboid from (5,5,5) to (205,205,205))

- `INF_PRISM`: Infinite prism - parameters: `ax pos1 pos2 pos3 pos4`
  - Defines a prism with 4 bounding planes extending infinitely along one axis
  - Example: SURF S12 INF_PRISM 2 5 5 100 100 (infinite prism along axis z, boundary planes at x = 5, y = 5, x = 100, y = 100)

#### Cell Definition
```
CELL [cell_id] [operator] [surface_id1] [operator] [surface_id2] ...
[MATERIAL|LATTICE|UNIV] [material_name|lattice_id|universe_id]
```
Supported operators:
- `in`: Cell is inside the surface
- `out`: Cell is outside the surface

Cell types:
- Material cell: `MATERIAL [material_name]`
- Lattice cell: `LATTICE [lattice_id] [pitch] [type] [universe_id]`
- Universe cell: `UNIV [universe_id]`

#### Lattice Definition (under CELL with LATTICE type) 
```
LATTICE [lattice_id] [pitch] [type] [universe_id]
[xmin] [xmax]
```
- `pitch`: Lattice spacing
- `type`: `mono` for uniform lattice, `alternate` for alternating pattern (checkerboard)
- `universe_id`: Universe to fill the lattice with
- `xmin`, `xmax`: Lattice boundaries (lattice can be only cuboid shaped so the defined boundaries will act for x, y, z)

#### Material Definition
```
DEFMAT [material_name] [density]
[num_atoms1] [isotope1] [num_atoms2] [isotope2] ...
```
Defines a material with density and isotopic composition.

#### Fuel Definition
```
DEFFUEL [material_name] [density] [enrichment]
[num_atoms1] [isotope1] [num_atoms2] [isotope2] ...
```
Defines a fuel material with enrichment.

#### Simulation Parameters
```
GEN [num_particles] [initial_keff] [source_x] [source_y] [source_z]
```
- `num_particles`: Number of neutrons per generation
- `initial_keff`: Initial k-effective guess
- `source_x`, `source_y`, `source_z`: Initial source position coordinates

## Output
- K-effective with statistical error

## Under developement
- Parallelization of the delta tracking algorithm
- Collision flux exstimators to obtain reaction rates ecc.