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