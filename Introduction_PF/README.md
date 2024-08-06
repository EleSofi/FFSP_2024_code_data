# Introduction_PF

## Overview

This folder contains the files and scripts necessary for running the numerical simulation of the Introduction Particle Filtering (PF) example. The following details explain the purpose and usage of each file and subdirectory within this folder.

## Folder Structure

### Introduction_PF/
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
- **Figures/**
    - `cartoon_microscope_without_background.png`: Image file for visualization.
    - `conditional_N.jpg`: Image file showing conditional N.
    - `mRNA_Observation_Trajectory.jpg`: Image file of mRNA observation trajectory.
    - `mRNA_filters_results.ai`: Adobe Illustrator file of mRNA filters results.
    - `mRNA_filters_results.png`: Image file of mRNA filters results.
    - `mRNA_protein.ai`: Adobe Illustrator file of mRNA and protein results.
    - `mRNA_protein.png`: Image file of mRNA and protein results.
- `Linear_Circuit_PF_FFSP.asv`: Version control file for MATLAB script.
- `Linear_Circuit_PF_FFSP.m`: MATLAB script for linear circuit PF and FFSP simulation.
- `generate_directory_structure.m`: MATLAB script to generate the directory structure.
- `main.asv`: Version control file for main MATLAB script.
- `main.m`: MATLAB script to create the Readme file through generate_directory_structure.m.
- `project_structure.md`: Markdown file explaining the project structure.

## File Descriptions

### Code_for_the_files/
This directory contains all the MATLAB scripts necessary for simulating the Introduction Particle Filtering example. Each script is responsible for different aspects of the simulation, from defining the reaction network to running particle filters and FFSP methods.

### Figures/
This directory contains various plots and images for visualizing the parameters and results when defining the Introduction Particle Filtering Network:
- `cartoon_microscope_without_background.png`: Visual representation of a microscope.
- `conditional_N.jpg`: Visualization of conditional N.
- `mRNA_Observation_Trajectory.jpg`: mRNA observation trajectory plot.
- `mRNA_filters_results.ai`: Vector graphic of mRNA filters results.
- `mRNA_filters_results.png`: Image of mRNA filters results.
- `mRNA_protein.ai`: Vector graphic of mRNA and protein results.
- `mRNA_protein.png`: Image of mRNA and protein results.

### Linear_Circuit_PF_FFSP.m
A MATLAB script for simulating the linear circuit using particle filtering and FFSP methods.

### generate_directory_structure.m
A MATLAB script to generate the directory structure of the project.

### main.m
MATLAB script to create the Readme file through generate_directory_structure.m.

### particles_different_N.eps
An EPS file showing particles for different N, used for visualization purposes.

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

Feel free to explore the repository, and if you have any questions or suggestions, please write an email to elena.dambrosio@bsse.ethz.ch. Happy coding!