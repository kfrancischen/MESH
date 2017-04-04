The SimulationPlanar class can be initiated in Lua script by
```lua
s = SimulationPlanar.new()
```

All of the function provided in [base class](baseClass.md) can be used except for the following changes.


!!! failure
    The following function are unavailable and cannot be called for a `SimulationPlanar` object.

```lua
SetGx()
SetGy()
Setperiodicity(p1, p2)
```

!!! note
    The following functions are added and specific to `SimulationPlanar` object.

```lua
SetKParallel(end)
```
* Arguments:
    1. end: [double], the end of the $k_{\parallel}$ integral. It is a normalized number with respect to $\omega/c$.

* Output: None

* Note: this function is essentially doing
$$ \int_{0}^{\text{end}\cdot\omega/c}dk_{\parallel}$$
where the integral is evaluated either using Gauss-legendre methdod or Gauss-kronrod adaptive method.

```lua
GetPhiAtKParallel(omega index, k parallel value)
```
* Arguments:
    1. omega index: [int], the index of omega where $\Phi(\omega(\text[index]), k_{\parallel})$ is evaluated.
    2. k parallel: [double], the $k_{\parallel}$ value where $\Phi(\omega(\text[index]), k_{\parallel})$ is evaluated. It should be a normalized value with respect to $\omega(\text[index])/c$.

* Output: None

```lua
IntegrateKPrarallel()
```
* Arguments: None

* Note: before using this function, make sure the flux integral can be reduced to a $k_{\parallel}$ integral. In principle, the cases when all the materials possess only scalar dielectric or diagonal forms
$$ \begin{pmatrix}
\epsilon_{1} & 0 & 0\\
0 & \epsilon_{1} & 0\\
0 & 0 & \epsilon_{2}
\end{pmatrix}$$ can use this function.

```lua
OptUseQuadgl(degree)
```
* Arguments:
    1. degree: [int, optional], the degree of legendre polynomial in evaluating the integral. If degree is not given, it is set to be $1024$.

* Output: None

```lua
OptUseQuadgk()
```
* Arguments: None

* Note: the simulation will use this function if the integral option is not specified.