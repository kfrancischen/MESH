The SimulationPattern class can be initiated in Lua script by
```lua
s = SimulationPattern.new()
```

All of the function provided in [base class](baseClass.md) can be used.


!!! note
    The following functions are added and specific to `SimulationPattern` object.

```lua
SetLattice(xLen, yLen, angle)
```
* Arguments:
    1. xLen: [double], the length of the periodicity in $x$ direction in SI unit.
    2. yLen: [double], the length of the periodicity in the other direction in SI unit.
    3. angle: the angle between the two lattice axises, in degree.

* Output: None

* Note: MESH handles the case when one of the real space lattice is aligned with $x$ axis (with length `xLen`), and the other axis should be in the upper plane. The angle between the two axises are specified by `angle`. See the following figure:
![coordinate system](coordinate.png)

```lua
GetReciprocalLattice()
```
* Arguments: None

* Outputs: {{$x_1$, $y_1$}, {$x_2$, $y_2$}}. Two tables of the reciprocal lattices in SI unit.

```lua
SetLayerPatternRectangle(layer name, material name, {centerx, centery}, angle, {widthx, widthy})
```
* Arguments:
    1. layer name: [string], the layer that this rectangle pattern will be embedded. Such layer should already exist in the simulation, otherwise an error message will be printed out.
    2. material name: [string],  the material used as the rectangle pattern. Such material should already exist in the simulation, otherwise an error message will be printed out.
    3. angle: the rotated angle with respect to the positive $x$ direction in a counterclockwise manner, in degree.
    4. {centerx, centery}: [double table], the centers of the rectangle pattern in $x$ and $y$ direction, respectively, in SI unit.
    5. {widthx, widthy}: [double table], the widths of the rectangle pattern in $x$ and $y$ direction, respectively, in SI unit.

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

```lua
SetLayerPatternEllipse(layer name, material name, {centerx, centery}, angle, {a, b})
```
* Arguments:
    1. layer name: [string], the layer that this ellipse pattern will be embedded. Such layer should already exist in the simulation, otherwise an error message will be printed out.
    2. material name: [string],  the material used as the ellipse pattern. Such material should already exist in the simulation, otherwise an error message will be printed out.
    3. angle: the rotated angle with respect to the positive $x$ direction in a counterclockwise manner, in degree.
    4. {centerx, centery}: [double table], the centers of the ellipse pattern in $x$ and $y$ direction, respectively, in SI unit.
    5. {a, b}: [double], the half widths the ellipse pattern in SI unit.

* Output: None

* Note: the ellipse with angle$=0$ is written as:
    $$\frac{(x-x_c)^2}{a^2}+\frac{(y-y_c)^2}{b^2}=1$$

```lua
SetLayerPatternPolygon(layer name, material name, {centerx, centery}, angle, { {v1_x, v2_x}, ..., {vn_x, vn_y} })
```
* Arguments:
    1. layer name: [string], the layer that this polygon pattern will be embedded. Such layer should already exist in the simulation, otherwise an error message will be printed out.
    2. material name: [string],  the material used as the polygon pattern. Such material should already exist in the simulation, otherwise an error message will be printed out.
    3. {centerx, centery}: [double table], the centers of the polygon pattern in $x$ and $y$ direction, respectively, in SI unit.
    4. angle: the rotated angle with respect to the positive $x$ direction in a counterclockwise manner, in degree.
    5. { {v1_x, v2_x}, ..., {vn_x, vn_y} }: [nested double table], coordinates of the vertices relative to the center in SI unit.

* Output: None
