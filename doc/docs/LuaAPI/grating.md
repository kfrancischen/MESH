The SimulationGrating class can be initiated in Lua script by
```lua
s = SimulationGrating.new()
```

!!! note
    The following functions are added and specific to `SimulationGrating` object.

```lua
SetLattice(p1)
```
* Arguments:
    1. p1: [double], the periodicity in $x$ direction in SI unit

* Output: None

```lua
SetLayerPatternGrating(layer name, material name, center, width)
```
* Arguments:
    1. layer name: [string], the layer that this grating will be embedded. Such layer should already exist in the simulation, otherwise an error message will be printed out.
    2. material name: [string],  the material used as the grating. Such material should already exist in the simulation, otherwise an error message will be printed out.
    3. center: [double], the center of the grating in SI unit.
    4. width: [double], the width of the grating in SI unit.

* Output: None
