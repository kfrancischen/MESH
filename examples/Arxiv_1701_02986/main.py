from MESH import SimulationPlanar, SimulationPattern, Constants
import math

consts = Constants()

def thetaDerivative(omega, T):
  theta = consts['h_bar'] * omega / (math.exp(consts['h_bar'] * omega / consts['k_B']/T) - 1)
  return theta ** 2 * math.exp(consts['h_bar'] * omega /consts['k_B'] / T) / consts['k_B'] / T ** 2

s1 = SimulationPlanar()
s1.AddMaterial("Si", "Si.txt")
s1.AddMaterial("Vacuum", "Vacuum.txt")
s1.AddLayer("SiBottom", 0, "Si")
s1.AddLayer("VacGap", 20e-9, "Vacuum")
s1.AddLayerCopy("SiTop", "SiBottom")

s1.SetSourceLayer("SiBottom")
s1.SetProbeLayer("VacGap")

s1.SetKParallelIntegral(500)
s1.OptUseQuadgk()
s1.InitSimulation()
s1.IntegrateKParallel()
phi = s1.GetPhi()
omega = s1.GetOmega()
for i in range(s1.GetNumOfOmega()):
  E = consts['h_bar'] * omega[i] / consts['q']
  spectrum = thetaDerivative(omega[i], 300) * phi[i] * consts['q'] / consts['h_bar'];
  print spectrum/1e5

f = 0.98
width = math.sqrt(f * 50e-9 * 50e-9)
s2 = SimulationPattern()
s2.SetLattice(50e-9, 50e-9, 90)
s2.SetNumOfG(440)

s2.AddMaterial("Si", "Si.txt")
s2.AddMaterial("Vacuum", "Vacuum.txt")
s2.AddLayer("SiBottom", 0, "Si")
s2.SetLayerPatternRectangle("SiBottom", "Vacuum", (25e-9, 25e-9), 0, (width, width))
s2.AddLayer("VacGap", 20e-9, "Vacuum")
s2.AddLayerCopy("SiTop", "SiBottom")

s2.SetSourceLayer("SiBottom")
s2.SetProbeLayer("VacGap")
s2.OutputSysInfo()

s2.OptPrintIntermediate()
s2.SetKxIntegralSym(20, 60)
s2.SetKyIntegralSym(20, 60)
s2.InitSimulation()
s2.IntegrateKxKy()

