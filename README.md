# mg: A Multigrid Solver for the Poisson Equation

This document outlines the structure, compilation, execution, and implementation details for the multigrid solver for the Poisson equation.

## Folder Structure

The project is organized as follows:

```
mg/
├── .vscode/                   // Any files pertaining to VS Code
│   └── settings.json
├── data/                      // Contains data after running the programs
│   ├── comparison.csv
│   ├── residuals.csv
│   └── summary.csv
├── figs/                      // Contains figures
│   ├── q1_perf-vs-lmax.png
│   ├── q1_residuals-lmax2.png
│   ├── q1_residuals-lmax3.png
│   ├── q1_residuals-lmax4.png
│   ├── q1_residuals-lmax5.png
│   ├── q1_residuals-lmax6.png
│   ├── q1_residuals-lmax7.png
│   └── q2_comparison.png
├── include/                   // Header files for...
│   ├── common.h               // common things between the first and second question,
│   └── multigrid.h            // and multigrid-specific functions
├── scripts/
│   └── plot.py
│── src/                       // All source files for...
│   ├── common/                // the multigrid-specific functions,
│   │   └── multigrid.c
│   ├── q1/                    // the main driver for the first question,
│   │   └── main.c
│   └── q2/                    // and the main driver for the second question
│       └── main.c
├── .clang-format
├── .editorconfig
├── Makefile
├── README.md
└── requirements.txt
```

A more detailed description of top-level files and folders follows:

- data/:
- figs/:
- include/:
- scripts/:
- src/:
- README.md: This report file.

## How to Compile and Run

1. **Prerequisites**:

   - A C/C++ compiler (e.g., GCC or G++).
   - Python packages as given in requirements.txt. Can be installed by running (within the root directory):

      ```bash
      pip install -r requirements.txt
      ```

2.  **Compile**:

    - Run the make command within the root directory. This will compile the C source files and create two executables - the first is q1 in bin/, and the second is q2 also in bin/, both of which are for answering the first and second question, respectively.

3.  **Run**:

    - Execute the compiled programs from the terminal using make run; this will run both programs and also carry out the plotting. This will generate the required plots in the figs/ directory.

## Code Explanation

We explain how our implementation follows the pseudocode given while also explaining the intricacies of the programs.

### smooth Function

### residual Function

### restriction Function

### prolongate Function

### gauss Function

### vcycle Function

## Figure Explanations

### Performance vs. lmax

### lmax Residual Histories

### Performance Difference for Different Multigrid Strategies
