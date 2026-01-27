
# SSE_H2SQ
J0, J1, J2, Implemented as of now

# IT_QMC

This repository contains the setup and workflow for running Quantum Monte Carlo (QMC) simulations for the IT_QMC project. The simulations can be executed either directly from the terminal or through a Jupyter notebook interface.

---

## Project Setup

To initialize the project and generate parameter-specific folders:

1. Change directory to `build`:
   `cd build`

2. Make the launch script executable:
   `chmod +x launch_dir.sh`

3. Run the script:
   `./launch_dir.sh`

This process creates a `files/` directory.  
Each subfolder inside `files/` corresponds to a unique parameter set and contains separate object files and executables for that configuration.

---

## Option 1: Terminal-Based Workflow

Use this method to manually compile and run simulations for a chosen parameter set.

1. Navigate to a specific parameter folder, for example:  
   `cd files/betaVp_100/fangleVp_0/L4/J2_0.00/J3_0.00/M1`

2. Clean previous builds:  
   `make clean`

3. Compile the code:  
   `make`

4. Run the executable:  
   `./main.exe`

This option is useful when you want direct control over compilation and execution, especially for debugging or testing individual cases.

---

## Option 2: Jupyter Notebook Workflow

This method runs the same simulations using a Jupyter notebook, which is more convenient for automation and batch runs.

1. Open the notebook file:  
   `H2SQ.ipynb`

2. Inside the notebook, define simulation parameters as a Python dictionary.

3. Execute the cells to call the `qmc_runner`, which internally performs the same steps as the terminal workflow.

This approach is recommended for running multiple parameter sets efficiently and for combining simulations with data analysis.

---

## Summary

Two equivalent ways are available to run the simulations:

- **Terminal workflow** – manual build and execution for a single parameter set.  
- **Jupyter workflow** – automated and efficient execution using Python.

Both approaches generate identical results; the Jupyter method simply provides a higher-level interface.

---

## Notes

- The `files/` directory is auto-generated and contains all parameter-specific runs.
- Each parameter folder is self-contained and can be built and executed independently.
- Use the terminal workflow for debugging and development.
- Use the Jupyter workflow for large-scale or repeated experiments.

---

End of README.