The SimulationGrating class can be initiated in python script by
```python
from MESH import SimulationGrating
s = SimulationGrating()
```

```python
SetNumOfG(nG)
```
* Arguments:
    1. nG: [int], the number of total Fourier components in both directions. Note this nG might not be the true nG used in the simulation.

* Output: None


```python
GetNumOfG()
```
* Arguments: None

* Output: number of G. If function `InitSimulation` has been called, then this function returns the true nG used in the simulation, otherwise return the user input nG.


!!! note
    The following functions are added and specific to `SimulationGrating` object.

```python
SetLattice(p1)
```
* Arguments:
    1. p1: [double], the periodicity in $x$ direction in SI unit

* Output: None

```python
SetLayerPatternGrating(layer name, material name, center, width)
```
* Arguments:
    1. layer name: [string], the layer that this grating will be embedded. Such layer should already exist in the simulation, otherwise an error message will be printed out.
    2. material name: [string],  the material used as the grating. Such material should already exist in the simulation, otherwise an error message will be printed out.
    3. center: [double], the center of the grating in SI unit.
    4. width: [double], the width of the grating in SI unit.

* Output: None
