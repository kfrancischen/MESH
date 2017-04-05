The C++ interface have the same functions as the Lua interface except for the following changes.

* An instance of a class is initiated to smart pointers. For example to instantiate a `SimulationPlanar` object, the right way to do is
```cpp
Ptr<SimulationPlanar> s = SimulationPlanar::instanceNew();
```
This also applies to the instantiation of `SimulationGrating` and `SimulationPattern` objects. The advantage of using smart pointers is that there is no need to do manual gabage collection using `delete`.

* Function names: the function names starts with lower case whereas Lua interface starts with a capital letter.

* extra inputs in function `setMaterial`. In C++, this function is called as    
```cpp
void setMaterial(const std::string name, double** epsilon, const std::string type)
```      
  where type is one of "scalar", "diagonal" and "tensor". The 2D array `epsilon` then depends on the type. For scalar, the epsilon is of dimension $\text{# omega}\times 2$. For diagonal, the epsilon is of dimension $\text{# omega}\times 6$, and for tensor, the epsilon is of dimension $\text{# omega}\times 10$.

* changes in obtaining physical constants. In C++ there is no need to initiate or call any function to retrieve the constants. A constant, for example $q$, can be directly obtained by using constant.q.

* Extra Makefile is needed. An example for the Makefile for a `main.cpp` file is         
```bash
CFLAGS=-std=c++11 -O3 -ffast-math -march=native -fopenmp  
MESHPATH=../../          
INCLUDES=-I$(MESHPATH)/src        
ARMAINCLUDE=-I$(MESHPATH)/src/arma -DARMA_DONT_USE_WRAPPER -DARMA_NO_DEBUG         
LIBS=-L$(MESHPATH)/build -lmesh -lopenblas -llapack        
CXX=g++       
all:     
 $(CXX) $(CFLAGS) $(INCLUDES) ${ARMAINCLUDE} main.cpp -o main $(LIBS)
```
  The `MESHPATH` should point to where the `src` folder is placed.