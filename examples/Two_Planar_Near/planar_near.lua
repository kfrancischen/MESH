-- this is a lua realization of the corresponding cpp file
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

s:BuildRCWA();
s:IntegrateKParallel();
s:SaveToFile("test_lua.txt");

numOfOmega = s:GetNumOfOmega();
phi = s:GetPhi();
omega = s:GetOmega();
for i = 1,s:GetNumOfOmega(), 1 do
  print(string.format("%e", omega[i]).."\t"..string.format("%e", phi[i]));
end