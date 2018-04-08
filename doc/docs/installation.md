#### Prerequisites
MESH comes with two versions: one has OpenMP support if the system has OpenMP libraries; one has MPI support. Both of them require the following prerequisites:

* [Lapack](http://www.netlib.org/lapack/) and [blas](http://www.netlib.org/blas/) (or blas mutants, such as [openblas](http://www.openblas.net/), [atlas](http://math-atlas.sourceforge.net/), or [mkl](https://software.intel.com/en-us/intel-mkl)). For MacOS this is not necessary.
* [Lua](https://www.lua.org/) (optional) version >= 5.2 (version 5.3.x is preferred). Required for Lua versions of MESH.
* [Python](https://www.python.org/) (optional) version 2 and 3 are both fine.

In order to have the MPI version of MESH installed, one needs to install MPI besides the above two libraries.

MESH can be downloaded by
```bash
git clone git@github.com:kfrancischen/MESH.git
cd MESH
```
Below are the instructions on install the Lua version:

### Installation on Linux

#### installation of CPP openmp version (no Lua required)

The CPP-only version can be installed by
```bash
make -f Makefile.without_lua
```
This command line will generate static MESH libraries called `libmesh.a` in `build/` folder, and can be used to compile along with C code.

#### installation of Lua vanilla version (no MPI required)

If the above libraries are installed on the system default directory, The vanilla version of MESH can be simply installed by
```bash
make mesh
```
Installing this version does not require MPI packages. The executable can be found in directory `build/`. The executable is called `mesh`.


#### installation of Lua MPI version (Lua required)

The MPI version can be installed by
```bash
make meshMPI
```

The executables can be found in directory `build/`. The executable is called `meshMPI`.

For a customized installation, please change the paths for compilers in `Makefile.Linux`. For OpenMP support, one can add (replace $4$ with the maximum number of cores supported by the computer)
```bash
export OMP_NUM_THREADS=4
```
to the `.bashrc` file.

#### installation of Python version (Python required)
To Install Python version, please modify `gensetup.py.sh`, and then
```bash
make meshPy
```

### Installation on MacOS
Default clang compiler does not support OpenMP. One can either follow the same steps exactly the same as the installation on Linux machines (then no OpenMP is not supported), or installing gcc with OpenMP first by
```bash
brew install gcc --without-multilib
```

and then change the compiles from cc and c++ in `Makefile.Darwin` to corresponding GNU compilers, and add `-fopenmp` in the `CFLAGS` and `CXXFLAGS`. The installation options are exactly the same as the installation on Linux.

### Installation on Windows
The easiest way to install mesh on windows is to download [Ubuntu on Windows](https://msdn.microsoft.com/en-us/commandline/wsl/about) and follow the instructions in the Linux installation part. With [Ubuntu on Windows](https://msdn.microsoft.com/en-us/commandline/wsl/about), one can do `sudo apt get` for the required packages.

### Location on Clusters
The mesh has been built on `hera`, `comet` and `stampede`. The executables are in the directory
```bash
/home/kfchen/MESH/build/
```

One can add this to one's own path by adding
```bash
export PATH="$PATH:/home/kfchen/MESH/build/"
```
in `.bashrc` (on stampede it is `.profile` for `bash`).

In addition, on `hera`, please add the following line to `.bashrc`
```bash
export LD_LIBRARY_PATH="LD_LIBRARY_PATH:/home/kfchen/mkl/mkl/lib/intel64"
```
On `comet`, please add the following to your job submission file
```bash
module purge
module load gnu
module load gnutools
module load mkl
module load openmpi_ib
```
and on `stampede`:
```bash
module purge
module load gcc
module load mvapich2
module load mkl
```
For a job that use more than $24$ cores for `comet` and $16$ cores for `stampede`, MPI version should be used. On `hera`, the OpenMP version is recommended.

(Update for `stampede2`, 09/27/2017):      
For `stampede2` please change `module load mvapich2` to `module load impi`

### Installation on clusters (not recommended)
If one wants to install MESH on his/her own directory, one `hera` please use
```bash
make -f Makefile.hera
```
On `comet` please install Lua at the same directory as MESH and type
```bash
module purge
module load gnu
module load gnutools
module load mkl
module load openmpi_ib
make -f Makefile.comet
```

On `stampede` please install Lua at the same directory as MESH and type
```bash
module purge
module load gcc
module load mvapich2
module load mkl
make -f Makefile.stampede
```

On clusters, both OpenMP version and MPI version will be generated.


Currently the Python version hasn't been tested on `comet` and `stampede`. For Python MPI support, one needs [mpi4py](http://mpi4py.readthedocs.io/en/stable/)