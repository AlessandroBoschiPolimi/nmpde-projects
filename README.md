# TITLE

## Build Commands

**Linux**

To use dealii included in the AMSC constainer use

	cmake -S . -B build

otherwise, either
- export the following environment variable with the path to your dealii installation

		export DEAL_II_DIR=/path/to/dealii

- specify it just for the cmake command

		cmake -S . -B build -DDEAL_II_DIR=/path/to/dealii


Then build using

	cd build
	make

The resulting executable is build/PDE-06

**Windows**

ChatGPT keeps gaslighting me by saying that is possible to compile dealii on windows, upon further queries it turns out that all the methods he provided are no longer valid (or have never been). For this reason the following instruction don't make any sense, but nevertheless i still believe in a day (within the deadline of this project) when i'll be able to use an actual IDE, so i won't erase this.

Install dealii and then export an environment variable with the installation path

	set DEAL_II_DIR=C:\path\to\dealii

Then, open the folder with Visual Studio, it will automatically detect the CMakeLists.txt, then F5 to build and run, or ctrl+B to build.

Alternatively, to generate a Visual Studio solution under the build folder

	cmake -S . -B build -G "Visual Studio 17 2022"

As above, it's possible to specify the dealii path directly in the cmake command.

## Mesh Generation

To create the mesh run the following command outside the container

	gmsh scripts/mesh_ellipse.geo -3 -o mesh/mesh.msh

It requires gmsh installed

	sudo apt install gmsh

## Launch Command

To run the command run:
`./PDE-06 <work_file> <forcing_term: 0|1>`

The work file contains the execution parameters, which have to be provided in the following format.

	-----
	<NeoHooke | Guccione>                         // material type
	<file <filename> | cube | rod>                // mesh and optional input filename
	out <filename>                                // output filename
	it <integer>                                  // linear solver max iterations   
	N <csv, of, integers>                         // ids of Neumann boundaries
	<Neumann function> [parameters]               // see later
	D <boundary_id> <Dirichlet function>          // ids of Dirichlet boundaries, see later
	[D <boundary_id> <Dirichlet function> ...]
	<material parameters>

Note: comments are now supported by our parser, but it ignores only lines starting with '#'.

It's possible to submit multiple jobs in the same execution, each must start with a line with exactly 5 '-'.

The possible values for `Neumann function` are
- `bowl_pull_out` representing a force pulling in the direction normal to surface of an ellipsoid, and as parameters expects a double representing the scaling of the force.
- `todo`

The possible values for `Dirichlet function` are
- `zero`: homogeneous Dirichlet condition

Each Dirichlet boundary must be on a different line, since they can have different `Dirichlet function`.

The `material parameters` for NeoHooke are
- `C <double>`: controls resistance to isochoric (shape-changing) deformation
- `l <double>`: represents volumetric (compressibility) response

The `material parameters` for Guccione are
- `c <double>`: todo: meaning
- `b <9 space separated integers>`: 
- `a <double>`: 
- `anfun <int>,<int>,<int> <int>,<int>,<int> <int>,<int>,<int>`: 

There is very little validation on the input, pls be gentle.

Below an example using NeoHooke law, on the mesh provided in file ../mesh/mesh.msh

	------
	NeoHooke
	file ../mesh/mesh.msh
	out ../out/acorn
	it 10000
	N 1
	bowl_pull_out -0.5
	D 2 zero
	D 3 zero
	C 1
	l 2

	-----
	Guccione
	cube
	out ../out/guccione
	it 10000
	N 1
	cube_pull 0.5
	D 2 zero
	D 3 zero
	c 1
	b 1 2 3 4 5 6 7 8 9
	a 0.4
	anfun 1,0,0 0,1,0 0,0,1
