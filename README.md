## Overview
This repository contains Stochastic Series Expansion (SSE) Quantum Monte Carlo (QMC) and Exact Diagonalization (ED) simulators targeting the long range XXZ model in the presence of external fields. 

## Dependencies

### Runtime Dependencies
- Python 3.10+  
- Julia 1.12+   
- CMake 3.24+  
- Clang 17+ or GCC 14+ (C++ Compiler)  

CMake and a C++ compiler are required for building the long range QMC backend. Julia is needed for our ED implementation

### Python Packages
This package also requires a number of Python libraries:  
- Numpy  
- Scipy  
- Matplotlib  

These are installed automatically by pip as part of the Heisenberg Ion package install

### C++ Dependencies
- spdlog  
- HDF5  

The spdlog dependency is fetched and compiled automatically by CMake while building the C++ source.
The HDF5 package needs to installed. This can be done as follows: 
```
brew install hdf5
```

## Getting Started
To use this package, first clone the Heisenberg Ion repository:  
```
git clone https://github.com/OpenQuantumDesign/oqd-heisenberg-ion.git
```

Then, install the package (after navigating to the local project root):  
```
pip install .
```

Execute an example long range QMC calculation with the following:
```
oqd-heisenberg-ion -i examples/input_files/qmc_long_range.txt
```

Note that the above commands expects the current working directory to be the project root i.e. the local directory where the src and example folders are located. The example input file being used in this calculation can be found at:  
```
examples/input_files
```

The above calculation creates a folder titled ```results``` in the project root. A different output folder can be provided as a command line override as follows: 
```
oqd-heisenberg-ion -i examples/input_files/qmc_long_range.txt -o root_folder [output folder path]
```

Details about configuring inputs for each of the simulators can be found in the [documentation](https://openquantumdesign.github.io/oqd-heisenberg-ion/). Also, the example notebook exhibits different approaches for specifying the input parameters. 

### Documentation 
Documentation is developed using MkDocs. To install documentation-related dependencies, use:

```
pip install -e ".[docs]"
```

The documentation server can be deployed locally using: 

```
cp -r examples/ docs/examples/
mkdocs serve
```