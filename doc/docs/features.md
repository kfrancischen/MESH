### Physical Quantities Computed By MESH

By classifying the geometries into three different categories: planar geometry where the basic components are plates, grating geometry which contains at least one layer has a grating along $x$ direction, and patterns where at least one layer has either rectangle or circular patterns. The details about these three types are discussed below.

In general, the total heat transfer between bodies with temperatures $T_1$ and $T_2$ is written as

$$P=\int d\omega [\Theta(T_1, \omega)-\Theta(T_2,\omega)]\Phi(\omega)$$

where the quantity $\Phi(\omega)$ characterize the strength of heat transfer between the two bodies involved.

#### Heat transfer in planar geometries

For planar geometries (implemented as [SimulationPlanar](LuaAPI/planar.md) object),  the quantity $\Phi(\omega)$ computed by integrating the $k_x$ and $k_y$ vector over the whole $k$ space, i.e.

$$\Phi(\omega)=\int_{-\infty}^{\infty}dk_x\int_{-\infty}^{\infty}dk_y \Phi(k_x,k_y, \omega)$$

And in the case of an isotropic material in $x$ and $y$ direction, the above integral can be reduced to

$$\Phi(\omega)=\int_{0}^{\infty}dk_{\parallel} \Phi(k_{\parallel}, \omega)$$

MESH directly provides function to either compute $\Phi(\omega, k_x, k_y)$, $\Phi(\omega, k_{\parallel})$ and $\Phi(\omega)$ directly.

#### Heat transfer in grating geometries

In the case of a grating geometry (implemented as [SimulationGrating](LuaAPI/grating.md) object), the heat transfer rate is written as

$$\Phi(\omega)=\int_{-G/2}^{G/2}dk_x\int_{-\infty}^{\infty}dk_y \Phi(k_x,k_y, \omega)$$

Here again, MESH gives access to a few things

* settings of the integral: for isotropic geometries $\int_{-G/2}^{G/2}$ and  $\int_{-\infty}^{\infty}$ can be reduced to twice of the integral over the positive axis.
* directly computation of $\Phi(k_x,k_y, \omega)$ and $\Phi(\omega)$, and allows to print out $\Phi(k_x,k_y, \omega)$ in the process of obtaining $\Phi(\omega)$.

#### Heat transfer in pattern geometries

In the case of a pattern geometry (implemented as [SimulationPattern](LuaAPI/pattern.md) object), the heat transfer rate is written as

$$\Phi(\omega)=\int_{-G_x/2}^{G_x/2}dk_x\int_{-G_y/2}^{-G_y/2}dk_y \Phi(k_x,k_y, \omega)$$

Again similar to the grating geometries, integral settings and printing of intermediate $\Phi(k_x,k_y, \omega)$ are supported. In addition, multiple different kinds of patterns can exist in one layer. However, the code now only supports a rectangle lattice and non-interleaving patterns. In MESH, the only two patterns supported are rectangle and circle. However, extending the code to non-trivial lattice and adding supports to more patterns are in principle doable.

#### Supports over both scalar dielectric and tensor dielectric

MESH also supports material whose dielectric can be scalar, diagonal, or a tensor. These three different types are encapsulated into a Union type in MESH. Currently for tensor types, MESH can only deal with the case when epsilon in the $z$ direction can be decomposed from $x$ and $y$ directions, i.e.

$$ \overleftrightarrow{\epsilon}=\begin{pmatrix}
\epsilon_{xx} & \epsilon_{xy} & 0\\
\epsilon_{yx} & \epsilon_{yy} & 0\\
0 & 0 & \epsilon_{zz}
\end{pmatrix}$$

#### A Comprehensive Lua Wrapper for Users

Lua wrapper for C code is widely used because of easy implementation and high readability. For example, [S4](https://web.stanford.edu/group/fan/S4/) which is also an implementation of RCWA for periodic geometries, also utilizes Lua as its front. Here MESH not only gives the user basic functionalities of computing physical quantities related to heat transfer, but also allows users to directly use MPI in a lua script, so that users have better control over the simulation at run time. To show the advantages of directly revealing MPI interface to the users, a concrete example of lua MPI interface is explained in the [example section](Examples/MPI.md). MESH also provides a vanilla version that can be built without MPI. In that case, OpenMP is used if the system has OpenMP libraries built in.