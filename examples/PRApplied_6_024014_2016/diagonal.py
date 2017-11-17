from MESH import SimulationPlanar
s = SimulationPlanar()
s.AddMaterial("hBN", "hBN.txt")
s.AddMaterial("MCT", "MCT.txt")
s.AddMaterial("PEC", "PEC.txt")
s.AddMaterial("Vacuum", "Vacuum.txt")

s.AddLayer("PECBottom", 0, "PEC");
s.AddLayer("hBNLayer", 5e-6, "hBN");
s.AddLayer("VacLayer", 1e-8, "Vacuum");
s.AddLayer("MCTLayer", 5e-6, "MCT");
s.AddLayerCopy("PECTop", "PECBottom");

s.SetSourceLayer("hBNLayer");
s.SetProbeLayer("VacLayer");

s.SetKParallelIntegral(100);
s.SetThread(4);
s.OptUseQuadgk();
s.InitSimulation();
s.IntegrateKParallel();
phi = s.GetPhi()
omega = s.GetOmega()
for i in range(s.GetNumOfOmega()):
  print omega[i], '\t', phi[i]