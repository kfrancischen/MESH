This is a simple example calculating far-field radiation from a $1~\mu m$ thick GaAs slab to vacuum. Here one can set the upper bound of $k_{\parallel}$ integration to $1$ because it is far field.
```lua
s = SimulationPlanar.new()
s:SetThread(4);
s:AddMaterial("GaAs", "GaAs.txt");
s:AddMaterial("Vacuum", "Vacuum.txt");
s:AddMaterial("PEC", "PEC.txt");
s:AddLayer("PECLayer", 0, "PEC");
s:AddLayer("GaAsLayer", 1e-6, "GaAs");
s:AddLayer("VacuumLayer", 0, "Vacuum");

s:SetSourceLayer("GaAsLayer");
s:SetProbeLayer("VacuumLayer");

s:OptUseQuadgk();
s:SetKParallelIntegral(1);
s:InitSimulation();
s:IntegrateKParallel();
phi = s:GetPhi();
omega = s:GetOmega();
for i = 1,s:GetNumOfOmega(), 1 do
  print(string.format("%e", omega[i]).."\t"..string.format("%e", phi[i]));
end
```

The Python version of this example is
```python
from MESH import SimulationPlanar

s = SimulationPlanar()
s.AddMaterial(material_name = 'GaAs', file_name = 'GaAs.txt')
s.AddMaterial(material_name = 'Vacuum', file_name = 'Vacuum.txt')
s.AddMaterial(material_name = 'PEC', file_name = 'PEC.txt')

s.AddLayer(layer_name = 'PECLayer', thickness = 0, material_name = 'PEC')
s.AddLayer(layer_name = 'GaAsLayer', thickness = 1e-6, material_name = 'GaAs')
s.AddLayer(layer_name = 'VacuumLayer', thickness = 0, material_name = 'Vacuum')

s.SetSourceLayer(layer_name = 'GaAsLayer')
s.SetProbeLayer(layer_name = 'VacuumLayer')

s.OptUseQuadgk()
s.InitSimulation()
s.SetKParallelIntegral(1);
s.IntegrateKParallel()
phi = s.GetPhi()
omega = s.GetOmega()
for i in range(s.GetNumOfOmega()):
  print omega[i],'\t', phi[i]
```