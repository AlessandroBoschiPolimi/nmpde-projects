# <img src="./logo.png" alt="Description of image" style="width:60px;height:60px;"/> Heart Mechanics Simplified Simulation

<!--toc:start-->
- [Build Commands](#build-commands)
  - [Linux](#linux)
- [Mesh Generation](#mesh-generation)
- [Launch Command](#launch-command)
  - [Configuration File](#configuration-file)
    - [Forcing function](#forcing-function)
    - [Neumann function](#neumann-function)
    - [Dirichlet function](#dirichlet-function)
    - [NeoHooke Settings](#neohooke-settings)
    - [Guccione Settings](#guccione-settings)
    - [Options](#options)
    - [Examples](#examples)
- [Container](#container)
<!--toc:end-->

## Build Commands

### Linux

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

    [options]
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

> Comments are now supported by our parser, but it ignores only lines starting with '#'.
> All comments must be put before the division line and the options

It's possible to submit multiple jobs in the same execution, each must start with a line with exactly 5 '-'.

The optional parameters for the meshes are
- `file`: mandatory `filename`
- `cube`: optional `int`, representing the refinement of the mesh; the default is 1, corresponding to 6 subdivisions along each dimension; the parameter multiplies the subdivisions per dimension, so don't go past 3 if you want to see the results within a lifetime.

#### Forcing function

- `null` or an empty line: F = 0;
- rod
    - `bend_rod`: bends the rod (to use with the rod);
- cube
    - `cube_squeeze`: squeezes the cube;
    - `cube_tear`: applies a force from inside out, tearing the cube;
    - `cube_squeeze_z_lot`:
    - `cube_sqeeze_z-_little`:
- bowl

#### Neumann function

- bowl
    - `bowl_pull_out`: represents a force pulling in the direction normal to surface of an ellipsoid,
    and as parameters expects a double representing the scaling of the force;
    - `bowl_push_in`: applies a force with modulus 'tau' in the normal direction w.r.t. surface of an ellipsoid with semi-axes with the ratio 7:7:17. OUTGOING. 
    This is the force used in paper 5 problem 2, inside the acorn;
- cube
    - `cube_pull`: applies a force in the direction normal to the surface of the face of the cube,
    on the whole surface. OUTGOING;
    - `cube_push`: applies a force in the direction normal to the surface of the face of the cube,
    on the whole surface. INGOING;
- rod
    - `rod_pull`

#### Dirichlet function

- `zero`: homogeneous Dirichlet condition;
- `sin`: sinusoidal condition;

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

Before the `-----` division line some options can be added:

- `@skip`: to skip the job;
- `@passed`: to mark the job as passed. A special line will print when the job is run;
- `@name(<name>)`: name of the test. Spaces are not allowed;
- `@time`: ;

#### Examples
 
Below two examples

	# pressing a NeoHookean cube on one side, keeping the opposite fixed, 2 * default mesh refinement
    @name(Pressing_Cube_Neo)
    @passed
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
    @skip @passed
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

## Container

In order to compile this project the deal.II library is required.
As [suggested on deal.II](https://github.com/dealii/dealii.git) it is possible 
to use docker images in which dealii is already configured.

