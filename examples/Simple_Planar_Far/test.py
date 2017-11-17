from MESH import SimulationPlanar

s = SimulationPlanar()
s.AddMaterial(material_name = 'GaAs', file_name = 'GaAs.txt')
s.AddMaterial(material_name = 'Vacuum', file_name = 'Vacuum.txt')
s.AddMaterial(material_name = 'PEC', file_name = 'PEC.txt')

s.AddLayer(layer_name = 'PECLayer', thickness = 0, material_name = 'PEC')
s.AddLayer(layer_name = 'GaAsLayer', thickness = 1e-6, material_name = 'GaAs')
s.AddLayer(layer_name = 'VacuumLayer', thickness = 0, material_name = 'Vacuum')

s.SetSourceLayer(layer_name = 'GaAsLayer')
s.SetProbeLayer(layer_name = 'VacuumLayer')

s.OptUseQuadgk()
s.InitSimulation()
s.SetKParallelIntegral(1);
s.IntegrateKParallel()
phi = s.GetPhi()
omega = s.GetOmega()
for i in range(s.GetNumOfOmega()):
  print omega[i],'\t', phi[i]