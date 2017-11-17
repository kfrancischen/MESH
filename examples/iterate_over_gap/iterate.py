from MESH import SimulationPlanar

s = SimulationPlanar()
s.AddMaterial("GaAs", "GaAs.txt")
s.AddMaterial("Vacuum", "Vacuum.txt")
s.AddMaterial("PEC", "PEC.txt")

s.AddLayer("PECBottom", 0, "PEC")
s.AddLayer("GaAsBottom", 1e-6, "GaAs")
s.AddLayer("VacGap", 1e-8, "Vacuum")
s.AddLayerCopy("GaAsTop", "GaAsBottom")
s.AddLayerCopy("PECTop", "PECBottom")


s.SetSourceLayer("GaAsBottom")
s.SetProbeLayer("VacGap")
s.OptUseQuadgk()
s.SetThread(4)
s.SetKParallelIntegral(500)

for i in range(10, 110, 10):
  s.SetLayerThickness("VacGap", i * 1e-9)
  s.InitSimulation()
  s.IntegrateKParallel()


