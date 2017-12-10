-- this is a lua realization of the corresponding cpp file
s = SimulationGrating.new()
s:SetLattice(1e-6);
s:AddMaterial("Au", "fullGold.txt");
s:AddMaterial("Vacuum", "fullVacuum.txt");

s:AddLayer("BottomAir", 0, "Vacuum");
s:AddLayer("GoldSubstrateBottom", 0.5e-6, "Au");
s:AddLayer("GoldGratingBottom", 5e-6,"Au");
s:SetLayerPatternGrating("GoldGratingBottom", "Vacuum", 0.5e-6, 0.2e-6);

s:AddLayer("VacGap", 1e-6, "Vacuum");
s:AddLayerCopy("GoldGratingTop", "GoldGratingBottom");
s:AddLayerCopy("GoldSubstrateTop", "GoldSubstrateBottom");
s:AddLayerCopy("TopAir", "BottomAir");

s:SetSourceLayer("GoldSubstrateBottom");
s:SetSourceLayer("GoldGratingBottom");
s:SetProbeLayer("VacGap");
s:OutputSysInfo();

s:OptPrintIntermediate("test");
s:SetThread(4);
s:SetNumOfG(101);
s:SetKxIntegralSym(500);
s:SetKyIntegralSym(200, 5);
s:InitSimulation();
s:IntegrateKxKy();
--print(string.format("%e", s:GetPhiAtKxKy(0, 0, 0)));
