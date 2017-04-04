The SimulationGrating class can be initiated in Lua script by
```lua
s = SimulationGrating.new()
```

All of the function provided in [base class](baseClass.md) can be used except for the following changes.

!!! warning
    The following function needs to be called specifically for `SimulationGrating` object.

```lua
Setperiodicity(p1)
```
The above function can only be called without the second argument `p2`.

!!! note
    The following functions are added and specific to `SimulationGrating` object.


```lua
SetLayerPatternGrating(layer name, material name, center, width)
```
* Arguments:
    1. layer name: [string], the layer that this grating will be embedded. Such layer should already exist in the simulation, otherwise an error message will be printed out.
    2. material name: [string],  the material used as the grating. Such material should already exist in the simulation, otherwise an error message will be printed out.
    3. center: [double], the center of the grating in SI unit.
    4. width: [double], the width of the grating in SI unit.

* Output: None

```lua
OptUseAdaptive()
```
* Arguments: None

* Note: this function will use the [spatial adaptive resolution method](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.125404) to compute the Fourier transform of the dielectric.