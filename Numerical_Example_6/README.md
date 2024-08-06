# Numerical_Example_6

## Overview

This folder contains the files and scripts necessary for running the numerical simulation of the sixth numerical example. The following details explain the purpose and usage of each file and subdirectory within this folder.

## Folder Structure

### Numerical_Example_6/
- `.DS_Store`: A system file created by macOS, can be ignored.
- **Code_for_the_files/**
    - `.DS_Store`: A system file created by macOS, can be ignored.
    - `Evolution_Matrix_FFSP_2.m`: MATLAB script for evolution matrix computation for the FFSP algorithm.
    - `Evolution_Probability_FFSP_2.m`: MATLAB script for evolution probability computation in the ODE evolution.
    - `Exponential.m`: MATLAB script for simulating the exponential distribution.
    - `FFSP_2.m`: MATLAB script for FFSP method.
    - `Gillespie_General.m`: MATLAB script for general Gillespie algorithm simulation in this hybrid setting.
    - `Hidden_State.m`: MATLAB script for hidden state space computation.
    - `bimolecular.m`: MATLAB script for identifying bimolecular networks.
    - `indeces_FFSP.m`: MATLAB script for computing indices for defining the FFSP evolution matrix.
    - `next_reaction.m`: MATLAB script for next reaction method.
    - `propensity_1.m`: MATLAB script for propensity functions.
    - `propensity_1_bimolecular.m`: MATLAB script for propensity functions for bimolecular reactions.
    - `propensity_G.m`: MATLAB script for general propensity functions.
    - `propensity_G_FSP.m`: MATLAB script for propensity functions specific to FSP.
    - `propensity_bimolecular.m`: MATLAB script for propensity functions for bimolecular reactions.
    - `propensity_bimolecular_FSP.m`: MATLAB script for propensity functions for bimolecular reactions specific to FSP.
- **Figures/**
    - `Filtering_11_cell_tight.tif`: TIFF image of filtering results for the 11th cell.
    - `Filtering_1_cell_tight.tif`: TIFF image of filtering results for the 1st cell.
    - `Filtering_5_cell_tight.tif`: TIFF image of filtering results for the 5th cell.
- `README.md`: This file.
- `create_readme.m`: MATLAB script to generate the README.md file.
- `data45.mat`: MATLAB data file with the trajectories of the observed cells during the experiments.
- `generate_directory_structure.m`: MATLAB script to generate the directory structure.
- `main.m`: Main MATLAB script to create the project directory.
- `mean_cell_intensity_open_loop.m`: MATLAB script for calculating mean cell intensity in open-loop conditions.
- `path_of_cells_1.mat`: MATLAB data file with the paths of cells from FFSP run generated with the script 'script_for_estimations.m'.
- `project_structure.md`: Markdown file explaining the project structure.
- `script_for_estimations.m`: MATLAB script for running the filtering estimations of the hybrid experimental setup.
- `script_for_plotting.asv`: Autosave version of the plotting script.
- `script_for_plotting.m`: MATLAB script for generating plots of the filtering estimations .

## File Descriptions

### Code_for_the_files/
This directory contains all the MATLAB scripts necessary for simulating the sixth numerical example. Each script is responsible for different aspects of the simulation, from defining the reaction network to running particle filters and FFSP methods.

### Figures/
This directory contains TIFF images of filtering results for different numbers of cells.

### data45.mat
MATLAB data file with the trajectories of the observed cells during the experiments.

### generate_directory_structure.m
A MATLAB script to generate the directory structure of the project.

### main.m
The main MATLAB script to create the project directory.

### mean_cell_intensity_open_loop.m
A MATLAB script for calculating the mean cell intensity in open-loop conditions.

### path_of_cells_1.mat
MATLAB data file with the paths of cells from the FFSP run.

### project_structure.md
Markdown file explaining the structure of the project and its components.

## Contributing

Please follow the standard Git workflow for contributing to this repository:
1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Commit your changes (`git commit -m 'Add new feature'`).
4. Push to the branch (`git push origin feature-branch`).
5. Create a new Pull Request.

---

Feel free to explore the repository, and if you have any questions or suggestions, please open an issue. Happy coding!