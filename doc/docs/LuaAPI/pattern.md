The SimulationPattern class can be initiated in Lua script by
```lua
s = SimulationPattern.new()
```

All of the function provided in [base class](baseClass.md) can be used.


!!! note
    The following functions are added and specific to `SimulationPattern` object.

```lua
SetLayerPatternRectangle(layer name, material name, {centerx, centery}, {widthx, widthy})
```
* Arguments:
    1. layer name: [string], the layer that this rectangle pattern will be embedded. Such layer should already exist in the simulation, otherwise an error message will be printed out.
    2. material name: [string],  the material used as the rectangle pattern. Such material should already exist in the simulation, otherwise an error message will be printed out.
    3. {centerx, centery}: [double table], the centers of the rectangle pattern in $x$ and $y$ direction, respectively, in SI unit.
    4. {widthx, widthy}: [double table], the widths of the rectangle pattern in $x$ and $y$ direction, respectively, in SI unit.

* Output: None

```lua
SetLayerPatternCircle(layer name, material name, {centerx, centery}, radius)
```
* Arguments:
    1. layer name: [string], the layer that this circle pattern will be embedded. Such layer should already exist in the simulation, otherwise an error message will be printed out.
    2. material name: [string],  the material used as the circle pattern. Such material should already exist in the simulation, otherwise an error message will be printed out.
    3. {centerx, centery}: [double table], the centers of the circle pattern in $x$ and $y$ direction, respectively, in SI unit.
    4. radius: [double], the radius of the circle pattern in SI unit.

* Output: None
