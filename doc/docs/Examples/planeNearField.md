This is the complete code for the [tutorial example](tutorial.md)
```lua
s = SimulationPlanar.new()
s:AddMaterial("GaAs", "GaAs.txt")
s:AddMaterial("PEC", "PEC.txt")
s:AddMaterial("Vacuum", "Vacuum.txt")

s:AddLayer("PECBottom", 0, "PEC");
s:AddLayer("GaAsBottom", 1e-6, "GaAs");
s:AddLayer("VacGap", 1e-8, "Vacuum");
s:AddLayerCopy("GaAsTop", "GaAsBottom");
s:AddLayerCopy("PECTop", "PECBottom");

s:SetSourceLayer("GaAsBottom");
s:SetProbeLayer("VacGap");
s:OptUseQuadgk();
s:SetKParallelIntegral(10);
s:SetThread(4);

s:InitSimulation();
s:IntegrateKParallel();

phi = s:GetPhi();
omega = s:GetOmega();
for i = 1,s:GetNumOfOmega(), 1 do
  print(string.format("%e", omega[i]).."\t"..string.format("%e", phi[i]));
end
```
This code can be run by:
```bash
mesh main.lua > output.txt
```
and the results will be contained in the `output.txt` file.

The python version of this is
```python
from MESH import SimulationPlanar


s = SimulationPlanar()
s.AddMaterial(material_name = 'GaAs', file_name = 'GaAs.txt')
s.AddMaterial(material_name = 'Vacuum', file_name = 'Vacuum.txt')
s.AddMaterial(material_name = 'PEC', file_name = 'PEC.txt')

s.AddLayer(layer_name = 'PECBottom', thickness = 0, material_name = 'PEC')
s.AddLayer(layer_name = 'GaAsBottom', thickness = 1e-6, material_name = 'GaAs')
s.AddLayer(layer_name = 'VacGap', thickness = 1e-8, material_name = 'Vacuum')
s.AddLayerCopy(layer_name = 'GaAsTop', copy_layer_name = 'GaAsBottom')
s.AddLayerCopy(layer_name = 'PECTop', copy_layer_name = 'PECBottom')

s.SetSourceLayer(layer_name = 'GaAsBottom')
s.SetProbeLayer(layer_name = 'VacGap')
s.OptUseQuadgk()
s.SetKParallelIntegral(integral_end = 10)
s.SetThread(num_thread = 4)

s.InitSimulation()
s.IntegrateKParallel()

phi = s.GetPhi()
omega = s.GetOmega()
for i in range(s.GetNumOfOmega()):
  print omega[i], '\t', phi[i]
```