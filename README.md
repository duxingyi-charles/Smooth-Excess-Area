# Optimizing Global Injectivity for Constrained Parameterization

![](figure/teaser.jpeg)

[Xingyi Du](https://duxingyi-charles.github.io/), [Danny M. Kaufman](https://research.adobe.com/person/danny-kaufman/), [Qingnan Zhou](https://research.adobe.com/person/qingnan-zhou/), [Shahar Kovalsky](https://shaharkov.github.io/), [Yajie Yan](https://yajieyan.github.io/), [Noam Aigerman](https://research.adobe.com/person/noam-aigerman/), [Tao Ju](https://www.cse.wustl.edu/~taoju/)
<br/>*ACM Transaction on Graphics (Proceedings of SIGGRAPH 2021)*<br/>

[`Project Page`](https://duxingyi-charles.github.io/publication/optimizing-global-injectivity-for-constrained-parameterization/)
[`Dataset`](https://doi.org/10.5281/zenodo.5547887)

## Abstract

Injective parameterizations of triangulated meshes are critical across applications but remain challenging to compute. Existing algorithms to find injectivity either require initialization from an injective starting state, which is currently only possible without positional constraints, or else can only prevent triangle inversion, which is insufficient to ensure injectivity. Here we present, to our knowledge, the first algorithm for recovering a globally injective parameterization from an arbitrary non-injective initial mesh subject to stationary constraints. These initial meshes can be inverted, wound about interior vertices and/or overlapping. Our algorithm in turn enables globally injective mapping for meshes with arbitrary positional constraints. Our key contribution is a new energy, called smooth excess area (**SEA**), that measures non-injectivity in a map. This energy is well-defined across both injective and non-injective maps and is smooth almost everywhere, making it readily minimizable using standard gradient-based solvers starting from a non-injective initial state. Importantly, we show that maps minimizing SEA are guaranteed to be locally injective and almost globally injective, in the sense that the overlapping area can be made arbitrarily small. Analyzing SEA’s behavior over a new benchmark set designed to test injective mapping, we find that optimizing SEA successfully recovers globally injective maps for 85% of the benchmark and obtains locally injective maps for 90%. In contrast, state-of-the-art methods for removing triangle inversion obtain locally injective maps for less than 6% of the benchmark, and achieve global injectivity (largely by chance as prior methods are not designed to recover it) on less than 4%.

## SEA-QN

![](figure/overview.png)

Here we release SEA-QN, a program that computes injective mapping under positional constraints by minimize SEA (Smooth Excess Area) energy using quasi-Newton method.

The above figure illustrates **what SEA does**. It takes as input a source mesh, a group of positional constraints and a non-injective initial map, and outputs an injective map satisfying the positional constraints.

This program has been tested on macOS xx.15.5 (Apple Clang xx.0.3) and Windows 10 (visual studio 2019 and 2022).

## Build

### Mac

We use [NLopt](https://nlopt.readthedocs.io/en/latest/) (version 2.6.1)'s L-BFGS quasi-Newton implementation.

The easiest way to build on Mac is to run the script, which installs NLopt using [homebrew](https://brew.sh/) and compiles the program.

    ./build_mac.sh

The program `SEA_QN` will be generated in the `build` subdirectory.

### Windows

To build the program, you can use CMake to generate a visual studio project from CMakeLists.txt.

## How to use

The executable `SEA_QN` asks for 3 arguments: path to an input data file, path to a solver options file, and path to the file to store the result.

    ./SEA_QN [input_file] [solver_options_file] [result_file]

An example is provided in the `example` subdirectory. Test it by:

    ./SEA_QN example/input example/solver_options example/my_result

The result will be written to `example/my_result`.

In the 3 arguments, `input_file` is mandatory, while the rest two are optional. If `solver_options_file` is not specified, `SEA_QN` will look for a file named `solver_options` in the same directory as the binary. If that file is not found, the program will fall back to default options. If `result_file` is not given, results will be written to a file named `result` in the directory of the binary.


## File format

### input_file

_Input file_ contains vertices and faces(triangles/tetrahedrons) information about the source mesh and initial embedding, as well as the indices of constrained vertices (called handles, usually are just boundary vertices). Vertices are indexed from 0.


    [num_sourceVert] [dimension_sourceVert]
    ... (num_sourceVert * dimension_sourceVert) Matrix ...
    [num_initVert]   [dimension_initVert]
    ... (num_initVert * dimension_initVert) Matrix ...
    [num_simplex]    [simplex_size]
    ... (num_simplex * simplex_size) Matrix ...
    [num_handles]
    ... (num_handles * 1) Matrix ...
 
 See `example/input` for a concrete example.
 
:bell:  **Important**: Since TLC aims at constrained embedding problem, the user should at least provide the indices of boundary vertices as handles in the `input_file`, or provide them in a `handleFile` as described below. 
To make this easier, we provide a script to generate a `handleFile` containing boundary vertex indices for a given input mesh. See below for usage.

 
:tada: **It's possible to use your own mesh formats.** We provide two python scripts in directory `IO` to convert common mesh formats to our `input_file` format.
 
 To use the two scripts, make sure to install [meshio](https://github.com/nschloe/meshio) with
 
     pip install meshio
 
 To convert triangle meshes to our input format, run
 
    ./convert_input_2D.py [inputObjFile] [handleFile] [outFile]
 
 Currently, we only support OBJ file with initial mesh as uv coordinates. Check out our [dataset](https://github.com/duxingyi-charles/Locally-Injective-Mappings-Benchmark) for some concrete OBJ and handle files.
 The generated `outFile` will have the format of our `input_file`.
 
 For your convenience, we also provide a script in directory `IO` to generate a `handleFile` containing all the boundary vertex indices for a given input mesh. The script works for both triangle/tetrahedron mesh.
   
     ./extract_boundary_vert.py [inputMeshFile] [outputHandleFile] 
 
 To convert tetrahedron rest(source) and initial meshes to our input format, run
 
    ./convert_input_3D.py [restFile] [initFile] [handleFile] [outFile]
 
 All tet-mesh formats supported by `meshio` should be handled by this script. We have tested the VTK format. For more examples in VTK format, please check out our [dataset](https://github.com/duxingyi-charles/Locally-Injective-Mappings-Benchmark).
    
 
### solver_options_file

_Solver options file_ contains parameters for TLC energy, options for NLopt solver, and a list of intermediate status to record during optimization.


    form
    [harmonic OR Tutte]
    alphaRatio
    [val]
    alpha
    [val]
    ftol_abs
    [val]
    ftol_rel
    [val]
    xtol_abs
    [val]
    xtol_rel
    [val]
    algorithm
    [LBFGS]
    maxeval
    [val]
    stopCode
    [none OR all_good]
    record
    vert    [0 OR 1]
    energy  [0 OR 1]
    minArea [0 OR 1]

The following table explains each option in details.
We **recommend** using the default values (especially "form", "alphaRatio" and "alpha") as they are most successful in our experiments. 
 
See `example\solver_options` for a concrete example.

|                | possible values  | default value | explanation                                                                                                                    |
|----------------|------------------|---------------|--------------------------------------------------------------------------------------------------------------------------------|
| form           | harmonic, Tutte  | Tutte         | two forms of TLC energy (see paper for details)                                                                                |
| alphaRatio     | [0, inf)         | 1e-6          | Specify the ratio of content (area or volume) between rest mesh and target domain. Default value 1e-6 is recommended.          |
| alpha          | (-inf, inf)      | -1            | If negative, alpha will be computed from alphaRatio. If non-negative, alpha will overwrite the value computed from alphaRatio. |
| ftol_abs       | (-inf, inf)      | 1e-8          | Absolute energy change stop threshold. Negative value means disabled.                                                          |
| ftol_rel       | (-inf, inf)      | 1e-8          | Relative energy change stop threshold. Negative value means disabled.                                                          |
| xtol_abs       | (-inf, inf)      | 1e-8          | Absolute variable change stop threshold. Negative value means disabled.                                                        |
| xtol_rel       | (-inf, inf)      | 1e-8          | Relative variable change stop threshold. Negative value means disabled.                                                        |
| algorithm      | LBFGS            | LBFGS         | Quasi-Newton method.                                                                                                           |
| maxeval        | positive integer | 10000         | max number of iterations stop threshold.                                                                                        |
| stopCode       | none, all_good   | all_good      | Custom stop criteria. "all_good": optimization will stop when there are no inverted elements.                                   |
| record:vert    | 0, 1             | 0             | 1: record target mesh vertices at each iteration.                                                                              |
| record:energy  | 0, 1             | 0             | 1: record TLC energy at each iteration.                                                                                        |
| record:minArea | 0, 1             | 0             | 1: record smallest simplex signed content (area or volume) at each iteration.                                                  |
      
   

### result_file

_Result file_ stores the vertices of result mesh, and also intermediate records as specified in solver options file.

    
    name dims
    data
    ...

See `example\result` for a concrete example.

We provide a script to convert a `result_file` to a mesh file in directory `IO`.

Usage

    ./get_result_mesh.py [inputFile] [resultFile] [outFile]

For example, 

    ./get_result_mesh.py example/input example/result result.vtk


## Dataset

**Download**:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5547887.svg)](https://doi.org/10.5281/zenodo.5547887)

We release the large benchmark dataset of triangle meshes used to evaluate our method and compare with existing methods. The dataset includes _1791_ triangular mesh examples. Each example includes a source mesh, up to 20 positional constraints, and a non-injective initial mesh generated by ARAP method. The dataset comes with both inputs and results of our method.
