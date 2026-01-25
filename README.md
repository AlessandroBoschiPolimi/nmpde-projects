# Heart Mechanics Simplified Simulation

## Build Commands

### **Linux**

To use dealii included in the AMSC constainer use

	cmake -S . --preset linux-debug
	cmake -S . --preset linux-release

otherwise, either
- export the following environment variable with the path to your dealii installation

		export DEAL_II_DIR=/path/to/dealii

- specify it just for the cmake command

		cmake -S . --preset <preset> -DDEAL_II_DIR=/path/to/dealii


Then build using

	cd build/<debug | release>
	make -j

## Mesh Generation

To create the mesh run the following command outside the container

	cd mesh
	source generate.sh

It requires gmsh installed

	sudo apt install gmsh

## Launch Command

The resulting executable is `build/<debug | release>/PDE-06`, to run it

	mpirun -n <# processors> ./PDE-06 <work_file>

### Configuration File

The work file contains the execution parameters, which have to be provided in the following format.

	-----
	<NeoHooke | Guccione>                         // material type
	<file <filename> | cube [int] | rod>          // mesh and parameters, see later
	out <filename>                                // output filename
	it <integer>                                  // linear solver max iterations
	[new_damn <double>]                           // newton damping
	[f <Forcing term>]							  // see later
	N <csv, of, integers>                         // ids of Neumann boundaries
	<Neumann function> [parameters]               // see later
	[D <boundary_id> <Dirichlet function> ...]    // ids of Dirichlet boundaries, see later
	<material parameters>

Note: comments are now supported by our parser, but it ignores only lines starting with '#'.

It's possible to submit multiple jobs in the same execution, each must start with a line with exactly 5 '-'.

The optional parameters for the meshes are
- `file`: mandatory `filename`
- `cube`: optional `int`, representing the refinement of the mesh; the default is 1, corresponding to 6 subdivisions along each dimension; the parameter multiplies the subdivisions per dimension, so don't go past 3 if you want to see the results within a lifetime.

The possible values for `Forcing function` are
- `todo`

The possible values for `Neumann function` are
- `bowl_pull_out` representing a force pulling in the direction normal to surface of an ellipsoid, and as parameters expects a double representing the scaling of the force.
- `todo`

The possible values for `Dirichlet function` are
- `zero`: homogeneous Dirichlet condition

Each Dirichlet boundary must be on a different line, since they can have different `Dirichlet function`.

#### NeoHooke Settings

The `material parameters` for NeoHooke are
- `C <double>`: controls resistance to isochoric (shape-changing) deformation
- `l <double>`: represents volumetric (compressibility) response

#### Guccione Settings

The `material parameters` for Guccione are
- `c <double>`: todo: meaning
- `b <9 space separated integers>`: 
- `anfun <int>,<int>,<int> <int>,<int>,<int> <int>,<int>,<int>`: 

There is very little validation on the input, pls be gentle.

#### Options

After the `-----` division line some options can be added:

- `@skip`: to skip the job;

#### Examples
 
Below two examples

	# pressing a NeoHookean cube on one side, keeping the opposite fixed, 2 * default mesh refinement
	-----
	NeoHooke
	cube 2
	out ../out/example1
	it 2000
	N 4
	cube_push 0.4
	D 5 zero
	C 1
	l 2

	# example for Guccione
	-----
	Guccione
	rod
	out ../out/example2
	it 10000
	N 4
	rod_pull 0.004
	D 0 zero
	c 2
	b 8 4 4 4 2 2 4 2 2
	anfun 1,0,0 0,1,0 0,0,1
