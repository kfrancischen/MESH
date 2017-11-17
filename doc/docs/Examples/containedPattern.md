The following example creates two contained circle patterns in the layer of `SiBottom`, resulting in a cylindrical shell made of vacuum.

```lua
s = SimulationPattern.new();
s:SetLattice(1e-7, 1e-7, 90);
s:SetNumOfG(250);
--s:OptSetLatticeTruncation("Parallelogramic")
s:AddMaterial("Si", "Si.txt");
s:AddMaterial("Vacuum", "Vacuum.txt");
s:AddLayer("SiBottom", 0, "Si");
s:SetLayerPatternCircle("SiBottom", "Vacuum", {50e-9, 50e-9}, 45e-9)
s:SetLayerPatternCircle("SiBottom", "Si", {50e-9, 50e-9}, 30e-9)
s:AddLayer("VacGap", 1e-7, "Vacuum");
s:AddLayer("SiTop", 0, "Si");
s:SetSourceLayer("SiBottom");
s:SetProbeLayer("VacGap");

s:OptPrintIntermediate();
s:SetKxIntegralSym(100);
s:SetKyIntegralSym(100);
s:InitSimulation();
s:OutputSysInfo();
s:IntegrateKxKy();
```
The output from the function `GetLayerPatternRealization` results in the following figure
![Fourier transform](fourier.png)
The resulting $\Phi(\omega)$ is
![composite](composite.png)

The corresponding Python version is
```python
from MESH import SimulationPattern

s = SimulationPattern()
s.SetLattice(1e-6, 1e-6, 90)
s.SetNumOfG(440)

s.AddMaterial("Si", "Si.txt")
s.AddMaterial("Vacuum", "Vacuum.txt")

s.AddLayer("SiBottom", 0, "Si")
s.AddLayer("VacGap", 1e-6, "Vacuum")
s.AddLayer("SiTop", 0, "Si")
s.SetLayerPatternRectangle("SiTop", "Vacuum", (2.5e-7, 2.5e-7), 0, (5e-7, 5e-7))
s.SetLayerPatternCircle("SiTop", "Vacuum", (7.5e-7, 7.5e-7), 2.5e-7)


s.SetSourceLayer("SiBottom")
s.SetProbeLayer("VacGap")
s.OutputSysInfo()

s.OptPrintIntermediate()
s.SetKxIntegralSym(20, 60)
s.SetKyIntegralSym(20, 60)
s.InitSimulation()
s.IntegrateKxKy()
```