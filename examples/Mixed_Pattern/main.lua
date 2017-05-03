
s = SimulationPattern.new();
s:SetLattice(1e-6, 1e-6, 90);
s:SetNumOfG(440);

s:AddMaterial("Si", "Si.txt");
s:AddMaterial("Vacuum", "Vacuum.txt");

s:AddLayer("SiBottom", 0, "Si");
s:AddLayer("VacGap", 1e-6, "Vacuum");
s:AddLayer("SiTop", 0, "Si");
s:SetLayerPatternRectangle("SiTop", "Vacuum", {2.5e-7, 2.5e-7}, {5e-7, 5e-7});
s:SetLayerPatternCircle("SiTop", "Vacuum", {7.5e-7, 7.5e-7}, 2.5e-7);


s:SetSourceLayer("SiBottom");
s:SetProbeLayer("VacGap");
s:OutputSysInfo();

s:OptPrintIntermediate();
s:SetKxIntegralSym(20, 60);
s:SetKyIntegralSym(20, 60);
s:InitSimulation();
s:IntegrateKxKy();

