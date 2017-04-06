### Package MESH
**M**ultilayer **E**lectromagnetic **S**olver for **H**eat transfer.

A program for computing electromagnetic far-field and near-field heat transfer for periodic, layered
structures, developed by [Kaifeng Chen](http://web.stanford.edu/~kfchen/) (<kfrancischen@gmail.com>) of the
[Fan group](http://web.stanford.edu/group/fan/) in the Stanford Electrical Engineering Department.

The program is built upon C++ and wrapped with Lua (>= 5.2), with OpenMP and MPI support. It is enabled with heat flux calculation in both far and near field for planar, grating and pattern geometries. The source code can be downloaded at [Github](https://github.com/kfrancischen/MESH). This document will cover the basic ideas behind MESH, complete descriptions of the Lua API and C++ API, and a few concrete examples created either to illustrate the simple usage or to reproduce some of results from existing literatures. The documents are organized as follows:

#### Overview
* [Equations & Features](features.md): detailed equations computed and features offered by the package.
* [Installation](installation.md): the installation steps on different platforms and clusters.

#### Lua API
* [Base Class](LuaAPI/baseClass.md): the base class the code is built upon.
* [SimulationPlanar](LuaAPI/planar.md): the inherited class for planar geometries.
* [SimulationGrating](LuaAPI/grating.md): the inherited class for 1D grating geometries.
* [SimulationPattern](LuaAPI/pattern.md): the inherited class for 2D pattern geometries.

#### C++ API
* [C++ classes and functions](C++API/classAndFunction.md): the C++ interface for all the functions.

#### Examples
* [Tutorial Example](Examples/tutorial.md): a simple tutorial about setting up the inputs and running the input script.
* [Single Plane Far-field](Examples/planeFarField.md): example about thermal radiation to the far field.
* [Two Planes Near-field](Examples/planeNearField.md): example about near-field heat transfer between two planes.
* [Anisotropic Material Near-field](Examples/anisotropic.md): example about simulations involving anisotropic material.
* [Iterate Over Gaps Near-field](Examples/iterate.md): example about iterating over different gaps for near-field heat transfer between two plates.
* [Two Gratings Near-field](Examples/gratingNearField.md): example the calculation between two grating structures in the near-field regime.
* [Two Rectangle Patterns Near-field](Examples/rectangleNearField.md): example about the near-field heat transfer between two rectangle patterned structures.
* [Mixed Patterns Near-field](Examples/mixedNearField.md): example about the near-field heat transfer between two structures with different kinds of patterns (rectangle and circle).
* [MPI Example](Examples/MPI.md): example about how to use Lua wrapped MPI interface for fast parallelization on distributed clusters.

#### Other information
* [Developer](develop.md): information for developers on how to contribute to MESH.
* [About](about.md): information about the author, license, copyrights and contact addresses.
