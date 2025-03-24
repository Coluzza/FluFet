# FluFet

FluFet is a simulation project designed to prepare simulations for LAMMPS to study the effects of immobilization, absorption, and flow on various surfaces with polymer grafting functionalization. The project supports different degrees of hydrophobicity or hydrophilicity and includes polymer brushes with ends of different charges (positive, negative, neutral) that can also act as donors or acceptors of hydrogen bonds.

## Features

- **Surface Types**:
  - Spherical
  - Cylindrical
  - Planar

- **Polymer Grafting Functionalization**:
  - Supports various degrees of hydrophobicity or hydrophilicity
  - Polymer brushes with ends of different charges (positive, negative, neutral)
  - Capability to act as hydrogen bond donors or acceptors

- **Simulation Preparation**:
  - Prepares input files for LAMMPS
  - Configurable parameters for different simulation scenarios
 
## Installation

To install and run FluFet, you need to have the following prerequisites:

- [LAMMPS](https://lammps.sandia.gov/) installed on your system
- A C compiler (e.g., GCC)

Clone the repository:

```bash
git clone https://github.com/Coluzza/FluFet.git
cd FluFet
```

## Usage

To prepare a simulation, you need to configure the parameters in the `GraftFlow_Brush_general.c` file. Here are the key options you can modify:

- **Surface Type**: Choose between spherical, cylindrical, or planar surfaces.
- **Degree of Hydrophobicity/Hydrophilicity**: Adjust the degree to match your simulation requirements.
- **Polymer Brush Charges**: Define the charges of the polymer brush ends (positive, negative, neutral).
- **Hydrogen Bonding**: Set the polymer brush ends to act as hydrogen bond donors or acceptors.

This script automates the preparation and submission of simulation jobs for LAMMPS on various HPC systems. It handles different simulation types and configurations, creating the necessary job scripts and submitting them to the appropriate job scheduler.

## Features

- **Automated Preparation**: Creates and configures job scripts for different HPC systems.
- **Parameter Handling**: Reads and processes simulation parameters from a `param.dat` file.
- **Job Submission**: Submits jobs to the HPC scheduler with dependencies, ensuring the correct order of execution.

## Usage

### Prerequisites

- Ensure you have LAMMPS installed and accessible on your HPC system.
- Ensure you have the necessary modules loaded (e.g., CUDA, OpenMPI, etc.).
- Prepare the input files (`pdb`, `param.dat`, `pka`, and `mask_linker.dat`).

### Command

To run the script, use the following command:

```bash
./GRAFTBRUSHFLOW_Run_HPC.sh <PDB_FILE> <PARAM_FILE> <PKA_FILE> <HPC_SYSTEM> <MASK_LINKER>
```

### Arguments

- **PDB_FILE**: Path to the PDB file.
- **PARAM_FILE**: Path to the parameter file (`param.dat`).
- **PKA_FILE**: Path to the pKa file.
- **HPC_SYSTEM**: HPC system to use (`HYPERION`, `KAROLINA`, `NOTS`, or `COGITATORE`).
- **MASK_LINKER**: Path to the mask linker file.

### Example

```bash
./GRAFTBRUSHFLOW_Run_HPC.sh protein.pdb param.dat protein.pka HYPERION mask_linker.dat
```

## `param.dat` Structure

The `param.dat` file defines the parameters for the simulation. Below is the structure and the expected parameters:

```plaintext
# Simulation type and geometry
Simul_Type = <type>  # e.g., BRUSH, FLOW
Geometry = <geometry>  # e.g., CYLINDER, SPHERE, PLANE

# Simulation parameters
Seed = <seed>  # Random seed for the simulation
Polymer_Fraction = <fraction>  # Fraction of polymer grafting
Polymer_Length = <length>  # Length of the polymer chains
Final_Flow_Temp = <temperature>  # Final temperature for flow simulations
E_Scale = <scale>  # Energy scale factor
N_Proteins = <number>  # Number of proteins in the simulation
pH = <ph_value>  # pH value for the simulation
Pressure = <pressure>  # Pressure for the simulation
```

### Example `param.dat`

```plaintext
Simul_Type = BRUSH
Geometry = CYLINDER

Seed = 12345
Polymer_Fraction = 0.2
Polymer_Length = 50
Final_Flow_Temp = 300
E_Scale = 1.0
N_Proteins = 100
pH = 7.0
Pressure = 1.0
```
