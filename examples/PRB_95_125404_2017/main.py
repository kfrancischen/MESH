from MESH import SimulationGrating

s = SimulationGrating()
s.SetLattice(0.5e-6)
s.SetNumOfG(101)
s.OptPrintIntermediate()

s.AddMaterial("Au", "fullGold.txt")
s.AddMaterial("Vacuum", "fullVacuum.txt")

s.AddLayer("GoldBottomSubstrate", 0, "Au")
s.AddLayer("GoldBottomGrating", 4.7e-6,"Au")
s.SetLayerPatternGrating("GoldBottomGrating", "Vacuum", 0.1e-6, 0.2e-6)
s.AddLayer("VacGap", 1e-6, "Vacuum")
s.AddLayerCopy("GoldLayerTopGrating", "GoldBottomGrating")
s.AddLayerCopy("GoldLayerTopSubstrate", "GoldBottomSubstrate")

s.SetSourceLayer("GoldBottomSubstrate")
s.SetSourceLayer("GoldBottomGrating")
s.SetProbeLayer("VacGap")

s.SetKxIntegralSym(500)
s.SetKyIntegralSym(200, 5)
s.OutputSysInfo()
s.SetThread(4)
s.InitSimulation()
s.IntegrateKxKy()
