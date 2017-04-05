This is a simple example calculating far-field radiation from a $1~\um m$ thick GaAs slab to vacuum. Here one can set the upper bound of $k_{\parallel}$ integration to $1$ because it is far field.
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
s:BuildRCWA();
s:SetKParallelIntegral(1);
s:IntegrateKParallel();
phi = s:GetPhi();
omega = s:GetOmega();
for i = 1,s:GetNumOfOmega(), 1 do
  print(string.format("%e", omega[i]).."\t"..string.format("%e", phi[i]));
end
```