### Package MESH
**M**ultilayer **E**lectromagnetic **S**olver for **H**eat transfer.

A program for computing electromagnetic far-field and near-field heat transfer periodic, layered
structures, developed by Kaifeng Chen (<kfrancischen@gmail.com>) of the
[Fan group](http://web.stanford.edu/group/fan/) in the Stanford Electrical Engineering Department.

The program is built upon C++ and wrapped with Lua (>= 5.2), with OpenMP and MPI support. It is enabled with heat flux calculation in both far and near field, along with implementations about structure visualization with [POVRay](http://www.povray.org/). This documents will cover the basic idea over MESH, a complete description of the Lua API and C++ API, and a few concrete examples either created to illustrate the simple usage or reproducing some results from existing literatures. The documents are organized as follows:

#### Overview
* [Features](features.md): detailed features and equations offered by the package.
* [Installation](installation.md): the installation steps on different platforms and clusters.

#### Lua API
* [Base Class](LuaAPI/baseClass.md): the base class the code is built upon.
* [SimulationPlanar](LuaAPI/planar.md): the inherited class computing planar geometries
* [SimulationGrating](LuaAPI/grating.md): the inherited class computing 1D grating geometries.
* [SimulationPattern](LuaAPI/pattern.md): the inherited class computing 2D pattern geometries

#### C++ API
* [C++ classes and functions](C++API/classAndFunction.md): the C++ interface for all the functions.

#### Examples
* [Tutorial Example](Examples/tutorial.md): a simple tutorial about how to setup the inputs and run the input script
* [Single Plane Far-field](Examples/planeFarField.md): example about how to compute thermal radiation to the far field.
* [Two Planes Near-field](Examples/planeNearField.md): example about how to compute near-field heat transfer between two planes.
* [Anisotropic Material Near-field](Examples/anisotropic.md): example about how to use anisotropic material
* [Iterate Over Gaps Near-field](Examples/iterate.md): example about how to iterate over different gaps for near-field heat transfer between two plates.
* [Two Gratings Near-field](Examples/gratingNearField.md): example about how to setup the calculation between two grating structures in the near-field regime.
* [Two Rectangle Patterns Near-field](Examples/rectangleNearField.md): example about the near-field heat transfer between two rectangle patterned structures
* [Mixed Patterns Near-field](Examples/mixedNearField.md): example about the near-field heat transfer between two structures with different kinds of patterns (rectangle and circle).
* [MPI Example](Examples/MPI.md): example about how to use Lua wrapped MPI interface for fast parallelization on distributed clusters.

#### Other information
* [Developer](develop.md): information for developers on how to contribute to MESH
* [About](about.md): Information about the author, license, copyrights and contact addresses.
