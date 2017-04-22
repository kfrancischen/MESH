-- this file is a lua realization of the corresponding cpp file
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
for i = 1, s:GetNumOfOmega() do
  print(phi[i])
end
