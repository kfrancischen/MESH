-- this is a lua realization of the corresponding cpp file
constants = Constants();

function thetaDerivative(omega, T)
  local theta = constants.h_bar * omega / (math.exp(constants.h_bar * omega/constants.k_B/T) - 1);
  return math.pow(theta, 2) * math.exp(constants.h_bar * omega/constants.k_B/T) /constants.k_B / math.pow(T, 2);
end


s1 = SimulationPlanar.new();
s1:AddMaterial("Si", "Si.txt");
s1:AddMaterial("Vacuum", "Vacuum.txt");
s1:AddLayer("SiBottom", 0, "Si");
s1:AddLayer("VacGap", 20e-9, "Vacuum");
s1:AddLayerCopy("SiTop", "SiBottom");

s1:SetSourceLayer("SiBottom");
s1:SetProbeLayer("VacGap");

s1:SetKParallelIntegral(500);
s1:OptUseQuadgk();
s1:BuildRCWA();
s1:IntegrateKParallel();
phi = s1:GetPhi();
omega = s1:GetOmega();
for i = 1, s1:GetNumOfOmega() do
  local E = constants.h_bar * omega[i] / constants.q;
  local spectrum = thetaDerivative(omega[i], 300) * phi[i] * constants.q / constants.h_bar;
  print(spectrum/1e5);
end

f = 0.98;
width = math.sqrt(f * 50e-9 * 50e-9);
s2 = SimulationPattern.new();
s2:SetPeriodicity(50e-9, 50e-9);
s2:SetGx(10);
s2:SetGy(10);

s2:AddMaterial("Si", "Si.txt");
s2:AddMaterial("Vacuum", "Vacuum.txt");
s2:AddLayer("SiBottom", 0, "Si");
s2:SetLayerPatternRectangle("SiBottom", "Vacuum", {25e-9, 25e-9}, {width, width});
s2:AddLayer("VacGap", 20e-9, "Vacuum");
s2:AddLayerCopy("SiTop", "SiBottom");

s2:SetSourceLayer("SiBottom");
s2:SetProbeLayer("VacGap");
s2:GetSysInfo();

s2:OptPrintIntermediate();
s2:SetKxIntegralSym(20, 60);
s2:SetKyIntegralSym(20, 60);
s2:BuildRCWA();
s2:IntegrateKxKy();

