#### Prerequisites
MESH comes with two versions: one has OpenMP support if the system has OpenMP libraries; one has MPI support. Both of them requires the following Prerequisites:

* Lapack and blas (or blas mutants, such as openblas, atlas, or mkl). For MacOS this is not necessary.
* Lua version >= 5.2 (version 5.3.x is preferred).

In order to have the MPI version of MESH installed, one needs to install MPI besides the above two libraries.

MESH can be downloaded by
```bash
git clone git@github.com:kfrancischen/MESH.git
unzip MESH_master.zip -d MESH
cd MESH
```

#### Installation on Linux

If the above libraries are installed on the system default directory, The vanilla version of MESH can be simply installed by
```bash
make
```

and the MPI version can be installed by
```bash
make meshMPI
```

Executables can be found in directory `build/`. The executables are called `mesh` and `meshMPI`.

For a customized installation, please change the paths for compilers in `Makefile.Linux`.

#### Installation on MacOS
Default clang compiler does not support OpenMP. One can either follow the same steps exactly the same as the installation on Linux machines (then no OpenMP is not supported), or installing gcc with OpenMP first by
```bash
brew install gcc --without-multilib
```

and then change the compiles from cc and c++ in `Makefile.Darwin` to corresponding GNU compilers, and add `-fopenmp` in the `CFLAGS` and `CXXFLAGS`.

#### Installation on Windows
//TODO, haven't tried yet

#### Location on Clusters
The mesh has been built on `hera`, `comet` and `stampede`. The executables are in the directory
```bash
/home/kfchen/MESH/build/
```

One can add this to one's own path by adding
```bash
export PATH="$PATH:/home/kfchen/MESH/build/"
```
in `.bashrc`.
