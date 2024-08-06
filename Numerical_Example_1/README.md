# Numerical_Example_1

## Overview

This folder contains the files and scripts necessary for running the numerical simulation of the first numerical example using Kalman Particle Filtering (PF) and Finite State Projection (FFSP) methods. The following details explain the purpose and usage of each file and subdirectory within this folder.

## Folder Structure

### Numerical_Example_1/
- `.DS_Store`: A system file created by macOS, can be ignored.
- **Code_for_the_files/**
    - `Continuous_Evolution_1.m`: MATLAB script for continuous evolution processes.
    - `Evolution_Matrix_FFSP_2.m`: MATLAB script for evolution matrix computation.
    - `Evolution_Probability_FFSP_2.m`: MATLAB script for evolution probability computation.
    - `Exponential.m`: MATLAB script for exponential processes.
    - `FFSP_2.m`: MATLAB script for FFSP method.
    - `Gillespie_General_MAK.m`: MATLAB script for Gillespie algorithm simulation.
    - `Hidden_State.m`: MATLAB script for hidden state computation.
    - `Offsprings.m`: MATLAB script for offsprings in reactions.
    - `bimolecular.m`: MATLAB script for bimolecular reactions.
    - `indeces_FFSP.m`: MATLAB script for FFSP indices.
    - `jump.m`: MATLAB script for jump processes.
    - `next_reaction.m`: MATLAB script for next reaction method.
    - `particle_filter_1.m`: MATLAB script for particle filtering.
    - `propensity_1.m`: MATLAB script for propensity functions.
    - `propensity_1_bimolecular.m`: MATLAB script for propensity functions for bimolecular reactions.
    - `propensity_G.m`: MATLAB script for propensity functions.
    - `propensity_bimolecular.m`: MATLAB script for propensity functions for bimolecular reactions.
    - `sampling.m`: MATLAB script for sampling.
- `Filters_results.jpg`: Image file showing the results of the filters.
- `Kalman_PF_FFSP.m`: MATLAB script to run Kalman Particle Filtering and FFSP methods.
- `generate_directory_structure.m`: MATLAB script to generate the directory structure.
- `main.m`: Main MATLAB script to run the generate_directory_structure function.
- `variables_for_plotting.mat`: MATLAB file containing variables for plotting and to be opened with the Kalman_PF_FFSP.m file.
- `README.md`: This file.

## File Descriptions

### Code_for_the_files/
This directory contains all the MATLAB scripts necessary for simulating the first numerical example. Each script is responsible for different aspects of the simulation, from defining the reaction network to running particle filters and FFSP methods.

### Filters_results.jpg
An image file showing the results of the filters applied in the simulation.

### Kalman_PF_FFSP.m
A MATLAB script to run the Kalman Particle Filtering and Finite State Projection methods for the first numerical example.

### generate_directory_structure.m
A MATLAB script to generate the directory structure of the project.

### main.m
The main MATLAB script generate the project structure through generate_directory_structure.

### variables_for_plotting.mat
A MATLAB file containing variables for plotting the results of the simulations with Kalman_PF_FFSP.m.

###create_readme
script for creating the README.md file.

## Contributing

Please follow the standard Git workflow for contributing to this repository:
1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Commit your changes (`git commit -m 'Add new feature'`).
4. Push to the branch (`git push origin feature-branch`).
5. Create a new Pull Request.

---

Feel free to explore the repository, and if you have any questions or suggestions, please send an email to elena.dambrosio@bsse.ethz.ch. Happy coding!