# VF-Unity-Parallelized
VF-Unity-Parallelized is a streamlined version of VirtualFlow, integrating the features of both VFVS (VirtualFlow for Virtual Screening) and VFLP (VirtualFlow for Ligand Preparation) into a cohesive workflow. Designed to operate seamlessly on SLURM systems, this workflow allows users to easily incorporate any docking software of their choice, ensuring maximum flexibility.

The core workflow involves users supplying a SMILES text file, the receptor of interest, and the docking parameters to facilitate large-scale docking simulations. Out-of-the-box, VF-Unity-Parallelized is configured to support docking with QuickVina 2.0 and Smina.

Execution within SLURM environments is highly optimized, with computations distributed in parallel across multiple CPUs and nodes. This design ensures efficient linear scaling relative to the number of molecules provided.


## Prerequisites
Please clone the repository using: 
```
git clone https://github.com/VirtualFlow/VFUparr.git
```
Please ensure that the following packages are installed: 
- [RDKit version 2021.09.5](https://www.rdkit.org/docs/Install.html)
- [Open Babel 3.1.0](https://openbabel.org/docs/dev/Installation/install.html)
- [Python 3.7.13](https://www.python.org/downloads/)
- [SELFIES](https://github.com/aspuru-guzik-group/selfies)


## File Navigator
* `DATA`: This directory is where users can place the receptor file and the corresponding executables for running the docking process.
* `OUTPUTS`: This directory is designated for storing the results of the docking simulations.
* `all.ctrl`: Contains all user-specifiable parameters required for the screening process, including the docking parameterization.
* `dataset_calc.py`: A Python script for running the docking on specified ligands.
* `submit.sh`:  A Slurm submission script for submitting an array of jobs for processing.

## Quick Start Guide

To get started with the docking simulations, follow the steps outlined below. These steps ensure that your configuration is correctly set up for your specific docking scenario:

1. **Configure Receptor Location:**
   - Open `all.ctrl` and specify the exact location of your receptor in the designated section.

2. **Set Docking Parameters:**
   - Within `all.ctrl`, enter the appropriate `CENTER-X/Y/Z` and `SIZE-X/Y/Z` coordinates to define your docking area.

3. **Specify SMILES List Path:**
   - In `all.ctrl`, input the path to your SMILES list file. This file is crucial for defining the molecular inputs for the simulation.
   - Ensure the file adheres to the format: each line contains a SMILES representation followed by a comma and the molecule ID (e.g., `C[C@@H](N)C(=O)O, Molecule1`). Lines must be separated by newline characters to distinguish between different molecular entries.

4. **Slurm Cluster Account:**
   - In `submit.sh`, replace `TODO` in `#SBATCH --account=TODO` with your actual Slurm cluster account name to ensure proper job submission.

5. **Job Submission Configuration:**
   - Adjust the number of jobs to submit for your docking calculation in `submit.sh` by modifying `#SBATCH --array=1-999` accordingly. Ensure this number matches the `MAX_NUM_JOBS` parameter set in `all.ctrl`.

6. **Executable Permissions:**
   - Make sure the docking executables have the correct executable permissions by running `chmod 777 ./DATA/qvina`.

7. **Submit Your Job:**
   - Finally, submit your job to the Slurm cluster with the command: `sbatch submit.sh`.

By following these steps, you'll be properly set up to conduct your docking simulations. Ensure all paths and parameters are double-checked for accuracy before submitting your job.

### Analyzing Your Jobs

Upon completion of docking calculations, the results will be systematically saved in the working directory of your repository. Look for files named following the pattern `OUTPUT_*_*.txt`, where each represents a different output from your simulations. These text files are comprehensive, containing vital information for each molecule processed:

- **SMILES String:** The unique identifier for the chemical structure of the molecule.
- **Docking Score:** A numerical value indicating the predicted affinity between the receptor and the ligand. A score of 10,000 indicates a failed docking calculation. 
- **Molecule ID:** A specific identifier assigned to the molecule for easy reference.
- **Input Ligand Location:** The path to the file containing the input ligand used in the docking simulation.
- **Docking Pose File Location:** The path to the file showing the preferred orientation (pose) of the ligand when bound to the protein receptor.

This organized output allows for efficient analysis and interpretation of your docking simulations, enabling a deeper understanding of the interaction between molecules and their potential efficacy.


### Contributing
If you are interested in contributing to VirtualFlow, whether it is to report a bug or to extend VirtualFlow with your own code, please see the file [CONTRIBUTING.md](CONTRIBUTING.md) and the file [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).

### License
The project ist distributed under the GNU GPL v2.0. Please see the file [LICENSE](LICENSE) for more details. 


### Citation
Gorgulla, Christoph, et al. "VirtualFlow 2.0-The Next Generation Drug Discovery Platform Enabling Adaptive Screens of 69 Billion Molecules." bioRxiv (2023): 2023-04.
