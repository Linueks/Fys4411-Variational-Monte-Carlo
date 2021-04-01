# Project1 FYS4411 Computational Physics, Quantum Mechanical Systems Spring 2021 
Below is the readme from the forked repository which contained the code skeleton this project is based on. Through the project significant changes have been made to the existing code.

## Usage
No command line interface / reading variables from file has been implemented as of yet, so when changing values a recompile is needed (stupid I know), this would be something prioritized in future work on this code. Now there just three boolean flags in the main.cpp file which toggles SimpleHarmonicOscillator, InteractingHarmonicOscillator and whether or not to use numeric differentiation for the kinetic energy calculation. There is also a flag for toggling Importance sampling.



## Simple Variational Monte Carlo solve for FYS4411

Example class structure for the first VMC project of FYS4411 (spring 2021). You may, if you wish, fork this repository and make it the basis of your project. If you choose to do this, you will have to implement a lot of the crucial functions yourself. The relevant functions you need to implement are spread throughout the project, but they are all commented with a note saying what each function should do.

Please note that this is only a start, and when you have implemented all of these functions you will only have completed the first exercise. However, once this is done, you will have a very good basis for further work, and adding functionality will be easier with a good class structure.

If you want to write your own code from scratch, you are of course welcome to do so, and feel free to use this code as inspiration for your own class structure.

- If you choose to use this code as a basis for your work, the first thing you should do is fork it, pull it down to your computer, and make sure it compiles and runs. See the next section on how to compile and run the project. After this you should spend at least 10 minutes looking at the structure and familiarizing yourself with how the classes interact with eachother. 
- A good way to do this may be to simply start at the top of the main.cpp file, and go through all the calls to the System class functions. Consider also the base classes WaveFunction, Hamiltonian, and InitialState and see which functions are virtual (which functions NEED to be implemented by any sub-class).
- You can skip over the output function in the Sampler class and the entire Random class.


### Compiling and running the project
There are now several options you can use for compiling the project. If you use QT Creator, you can import this project into the IDE and point it to the `.pro`-file. If not, you can use CMake to create a Makefile for you which you can then run. You can install CMake through one of the Linux package managers, e.g., `apt install cmake`, `pacman -S cmake`, etc. For Mac you can install using `brew install cmake`. Other ways of installing are shown here: [https://cmake.org/install/](https://cmake.org/install/).

#### Compiling the project using CMake
In a Linux/Mac terminal this can be done by the following commands
```bash
# Create build-directory
mkdir build

# Move into the build-directory
cd build

# Run CMake to create a Makefile
cmake ../

# Make the Makefile using two threads
make -j2

# Move the executable to the top-directory
mv vmc ..
```
Or, simply run the script `compile_project` via
```bash
./compile_project
```
and the same set of commands are done for you. Now the project can be run by executing
```bash
./vmc
```
in the top-directory.

##### Cleaning the directory
Run `make clean` in the top-directory to remove the executable `vmc` and the `build`-directory.

##### Windows
Compilation of the project using Windows is still an open question to me, but please include a pull-request if you've got an example. CMake should be OS-independent, but `make` does not work on Windows.