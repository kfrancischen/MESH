MESH is written in a inheritance manner, so most of the functions in the base class can be directly accessed by subclasses. Usage of MESH involves writing a Lua script to call into various parts of MESH. Here we describe all of the MESH base class functions that can be called within the Lua environment.

For detailed function calls of a given geometry, please read the pages [SimulationPlanar](planar.md), [SimulationGrating](grating.md) and [SimulationPattern](pattern.md) for the geometries you are simulating.

!!! note
    The base class is just a wrapper for most of the functions, but it cannot be initiated in a Lua script. The only instances that can be initiated are the classes corresponding to different dimensions.

```lua
AddMaterial(material name, input file)
```
* Arguments:
    1. material name: [string], the name of the material added to the simulation. Such name is unique and if there already exists a material with the same name, an error message will be printed out.
    2. input file: [string], a file that contains the dielectric properties of the corresponding material. For scalar dielectric, the input file should be formatted as     
    ```    
    omega eps_r eps_i      
    ```    
    For diagonal dielectric, the format is    
    ```    
    omega eps_xx_r eps_xx_i eps_yy_r eps_yy_i eps_zz_r  eps_zz_i      
    ```    
    For tensor dielectric, the format is    
    ```    
    omega eps_xx_r eps_xx_i eps_xy_r eps_xy_i eps_yx_r eps_yx_i eps_yy_r eps_yy_i eps_zz_r  eps_zz_i      
    ```     
    The omega needs to be aligned for all the materials in the simulation.

* Output: None

```lua
SetMaterial(material name, new epsilon)
```
* Arguments:
    1. material name: [string], the name of the material whose epsilon will be changed. This material should already exist in the simulation, otherwise an error will be printed out.
    2. new epsilon: [nested table], the length equals the number of omega, and per row is the epsilon values with the same format as the input of `AddMaterial` function. i.e. for scalar case the length will be $2$, for diagonal case the length is $6$ and for tensor case the length is $10$.

* Output: None

```lua
AddLayer(layer name, thickness, material name)
```
* Arguments:
    1. layer name: [string], the name of the layer. Similarly, the names for layers are unique, and if such name already exists in the simulation, an error message will be printed out.
    2. thickness: [double], the thickness of the new layer in SI unit.
    3. material name: [string], the material that is used as the background of the layer. This material should already exist in the simulation, otherwise an error message will be printed out.

* Output: None

* Note: this new added layer will be placed on top of all the previous layers.

```lua
SetLayer(layer name, thickness, material name)
```
* Arguments:
    1. layer name: [string], the layer whose thickness and background will be changed. Such layer needs to already exist in the simulation, otherwise an error message will be printed out.
    2. thickness: [double], the new thickness of the layer.
    3. material name: [string], the material for the background of the layer. If such material does not exist, an error message will be printed out.

* Output: None

```lua
SetLayerThickness(layer name, thickness)
```
* Arguments:
    1. layer name: [string],  the layer whose thickness will be changed. Such layer needs to already exist in the simulation, otherwise an error message will be printed out.
    2. thickness: [double], the new thickness of the layer.

* Output: None  

```lua
AddLayerCopy(new layer name, original layer name)
```
* Arguments:
    1. new layer name: [string], the new layer that is copied from the original layer.  Such name cannot already exist in the simulation, otherwise an error message will be printed out.
    2. original layer name: [string], the original layer from whom everything is copied. If this layer does not exist in the simulation, an error will be printed out.

* Output: None

* Note: this function only copies the structure information, for example any pattern of the original layer, but does not copy any thermal information. For example, even the original layer is set as a source, the copied layer is still not a source. In addition, this new added layer will be placed on top of all the previous layers.

```lua
DeleteLayer(layer name)
```
* Arguments:
    1. layer name: [string], the name of the layer that will be deleted. Such layer should already be in the system, otherwise an error will be printed out.

* Output: None

```lua
SetPeriodicity(p1, p2)
```
* Arguments:
    1. p1: [double], the periodicity in $x$ direction in SI unit.
    2. p2: [double, optional for grating geometry], the periodicity in $y$ direction in SI unit.

* Output: None

```lua
SetGx(nGx)
```
* Arguments:
    1. nGx: [int], the number of positive Fourier components in $x$ direction. The total number of G is thus 2nGx+1.

* Output: None

```lua
SetGy(nGy)
```
* Arguments:
    1. nGy: [int], the number of positive Fourier components in $y$ direction. The total number of G is thus 2nGy+1.

* Output: None


```lua
SetSourceLayer(layer name)
```
* Arguments:
    1. layer name: [string], the name of the layer that is designated as the source layer. Such layer should already exist in the system, otherwise an error will be printed out.

* Output: None

* Note: a system can have more than $1$ source layers

```lua
SetProbeLayer(layer name)
```
* Arguments:
    1. layer name: [string], the name of the layer that is designated as the probe layer of the flux. Such layer should already exist in the system, otherwise an error will be printed out.

* Output: None

* Note: a system can have only one probe layer. Setting another layer as the probe layer will overwrite the previous one. In addition, the source layer should be above all the source layers in the real geometry.

```lua
SetThread(nthread)
```
* Arguments:
    1. nthreads: [int], number of threads used in OpenMP.

* Output: None
* Note: this function only works in an OpenMP setup.

```lua
SetKxIntegral(points, end)
```
* Arguments:
    1. points: [int], number of points in the integration
    2. end: [double, optional for grating and pattern geometries], the end of the integral over $k_x$. This end should be a normalized number with respect to $\omega/c$.

* Output: None

* Note: this function is essential doing
$$\int_{-\text{end}\cdot \omega/c}^{\text{end}\cdot \omega/c}dk_x$$ where the integral is evaluated as a summation of `points` points. In the case when `end` is not given, the upper bound of the integral will be $(G_x/2)/(\omega/c)$.

```lua
SetKyIntegral(points, end)
```
* Arguments:
    1. points: [int], number of points in the integration
    2. end: [double, optional for pattern geometries], the end of the integral over $k_y$. This end should be a normalized number with respect to $\omega/c$. In the case when `end` is not given, the upper bound of the integral will be $(G_y/2)/(\omega/c)$.

* Output: None

* Note: this function is essential doing
$$\int_{-\text{end}\cdot \omega/c}^{\text{end}\cdot \omega/c}dk_y$$ where the integral is evaluated as a summation of `points` points.

```lua
SetKxIntegralSym(points, end)
```
* Arguments:
    1. points: [int], number of points in the integration
    2. end: [double, optional for grating and pattern geometries], the end of the integral over $k_x$. This end should be a normalized number with respect to $\omega/c$.

* Output: None

* Note: this function is essential doing
$$2\times \int_{0}^{\text{end}\cdot \omega/c}dk_x$$ where the integral is evaluated as a summation of `points` points. In the case when `end` is not given, the upper bound of the integral will be $(G_x/2)/(\omega/c)$.

```lua
SetKyIntegralSym(points, end)
```
* Arguments:
    1. points: [int], number of points in the integration
    2. end: [double, optional for pattern geometries], the end of the integral over $k_y$. This end should be a normalized number with respect to $\omega/c$.

* Output: None

* Note: this function is essential doing
$$2\times \int_{0}^{\text{end}\cdot \omega/c}dk_y$$ where the integral is evaluated as a summation of `points` points. In the case when `end` is not given, the upper bound of the integral will be $(G_y/2)/(\omega/c)$.

```lua
BuildRCWA()
```
* Arguments: None

* Note: this function builds up the matrices for the dielectric.

```lua
IntegrateKxKy()
```
* Arguments: None

* Output: None

* Note: this function integrates over $k_x$ and $k_y$ based on the integral properties set by the user. So the function can only be called after the $k_x$ and $k_y$ integrals are configured, and RCWA matrices are built.

```lua
IntegrateKxKyMPI(rank, size)
```
* Arguments:
    1. rank: [int], the rank of the thread
    2. size: [int], the total size of the MPI run.

* Output: None

*  Note: this function can only be called during MPI. For an example of a funtion call, please refer to [MPI example](../Examples/MPI.md).

```lua
GetNumOfOmega()
```
* Arguments: None

* Output: [int], the number of total omega points computed in the simulation.

```lua
GetPhi()
```
* Arguments: None

* Output: [table of double], the $\Phi(\omega)$ values obtained from the simulation.

* Note: can only be called after $k_x$ and $k_y$ are integrated.

```lua
GetOmega()
```
* Arguments: None

* Output: [table of double], the omega values computed in the simulation.

```lua
GetPhiAtKxKy(omega index, kx, ky)
```
* Arguments:
    1. omega index: [int], the index of the omega value where $\Phi(\omega[\text{index}], k_x, k_y)$ is evaluated.
    2. kx: [double], the $k_x$ value where $\Phi(\omega[\text{index}], k_x, k_y)$ is evaluated. It is a normalized by $\omega[\text{index}]/c$.
    3. ky: [double], the $k_y$ value where $\Phi(\omega[\text{index}], k_x, k_y)$ is evaluated. It is a normalized by $\omega[\text{index}]/c$.

* Output: [double], the value of $\Phi(\omega, k_x, k_y)$.

```lua
OutputSysInfo()
```
* Arguments: None

* Output: the function prints out a system description to screen.

```lua
OutputStructurePOVRay(file name)
```
* Arguments:
    file name: [string], the output file that contains the POVRay object that describes the system. The file needs to have extension ".pov".
* Output: None


Also, MESH provides some options for printing intermediate information and methods for Fourier transform of the dielectric.


```lua
OptUseNaiveRule()
```
* Arguments: None

* Note: this function tells the RCWA to use the simplest closed form Fourier transform for the dielectric.

```lua
OptUseInverseRule()
```
* Arguments: None

* Note: this function tells the RCWA to use the inverse rule of the Fourier transform for the dielectric.

```lua
OptPrintIntermediate()
```
* Arguments: None

* Note: this function prints intermediate $\Phi(\omega, k_x, k_y)$ when function `IntegrateKxKy()` or `IntegrateKxKyMPI(rank, size)` is called.
