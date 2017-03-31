-- this is a lua realization of the corresponding cpp file
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
s:SetKParallelIntegral(500)

for i = 10, 100, 10 do
  s:SetLayerThickness("VacGap", i * 1e-9);
  s:BuildRCWA();
  s:IntegrateKParallel();
end


