This file is similar to the [Two Planes Near-field](planeNearField.md) example, except for a loop over the vacuum gap separations.

```lua
s = SimulationPlanar.new()
s:AddMaterial("GaAs", "GaAs.txt");
s:AddMaterial("Vacuum", "Vacuum.txt");
s:AddMaterial("PEC", "PEC.txt");

s:AddLayer("PECBottom", 0, "PEC");
s:AddLayer("GaAsBottom", 1e-6, "GaAs");
s:AddLayer("VacGap", 1e-8, "Vacuum");
s:AddLayerCopy("GaAsTop", "GaAsBottom");
s:AddLayerCopy("PECTop", "PECBottom");


s:SetSourceLayer("GaAsBottom");
s:SetProbeLayer("VacGap");
s:OptUseQuadgk();
s:SetThread(4);
s:SetKParallelIntegral(500);

for i = 10, 100, 10 do
  s:SetLayerThickness("VacGap", i * 1e-9);
  s:InitSimulation();
  s:IntegrateKParallel();
  phi = s:GetPhi();
  omega = s:GetOmega();
  for j = 1,s:GetNumOfOmega(), 1 do
    print(string.format("%e", omega[j]).."\t"..string.format("%e", phi[j]));
  end
end
```

One can use the same way described here in one's own simulation for a scan of thickness of layers.

The corresponding Python version is
```python
from MESH import SimulationPlanar

s = SimulationPlanar()
s.AddMaterial("GaAs", "GaAs.txt")
s.AddMaterial("Vacuum", "Vacuum.txt")
s.AddMaterial("PEC", "PEC.txt")

s.AddLayer("PECBottom", 0, "PEC")
s.AddLayer("GaAsBottom", 1e-6, "GaAs")
s.AddLayer("VacGap", 1e-8, "Vacuum")
s.AddLayerCopy("GaAsTop", "GaAsBottom")
s.AddLayerCopy("PECTop", "PECBottom")


s.SetSourceLayer("GaAsBottom")
s.SetProbeLayer("VacGap")
s.OptUseQuadgk()
s.SetThread(4)
s.SetKParallelIntegral(500)

for i in range(10, 110, 10):
  s.SetLayerThickness("VacGap", i * 1e-9)
  s.InitSimulation()
  s.IntegrateKParallel()



```