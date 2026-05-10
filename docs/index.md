---
hide:
    - navigation
---

# Docs

!!! Note
    Heisenberg Ion is still under active development so breaking changes are possible.
    See [Features in Development](feat_dev.md) for details.

## Introduction
Welcome to Heisenberg Ion, our open source quantum many-body physics simulator targeting lattice systems that can be natively realized in trapped ion architectures. We currently support Exact Diagonalization (ED) and Quantum Monte Carlo (QMC) Stochastic Series Expansion (SSE) engines for long range anisotropic Heisenberg models in the presence of external fields. 

While Heisenberg Ion is primarily a python package, it uses a C++ engine for large-scale QMC simulations. 
Our ED calculator uses a Julia backend. 
Both of these simulators are wrapped in a lightweight Python preprocessor responsible for driving the required engine.  

## Dependencies

### Runtime Dependencies
- Python 3.10+  
- Julia 1.12+  
- CMake 3.24+  
- Clang 17+ or GCC 14+ (C++ Compiler)  

CMake and a C++ compiler are required for building the long range QMC source code 

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

Then, install the package after navigating to the local project root:  
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

Details about configuring inputs for each of the simulators can be found in the [Input Specification](user_guide/input_specs.md) section. Also, see the [User Guide](user_guide/overview.md) and the example notebook in the repository for user-level documentation. A description of the SSE algorithm used in this package can be found in the [Algorithms](algorithms/sse.md) section.

### Documentation 
Documentation is developed using MkDocs. To install documentation-related dependencies, use:

```
pip install -e ".[docs]
```

The documentation server can be deployed locally using: 

```
cp -r examples/ docs/examples/
mkdocs serve
```