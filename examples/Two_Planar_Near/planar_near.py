from MESH import SimulationPlanar


s = SimulationPlanar()
s.AddMaterial(material_name = 'GaAs', file_name = 'GaAs.txt')
s.AddMaterial(material_name = 'Vacuum', file_name = 'Vacuum.txt')
s.AddMaterial(material_name = 'PEC', file_name = 'PEC.txt')

s.AddLayer(layer_name = 'PECBottom', thickness = 0, material_name = 'PEC')
s.AddLayer(layer_name = 'GaAsBottom', thickness = 1e-6, material_name = 'GaAs')
s.AddLayer(layer_name = 'VacGap', thickness = 1e-8, material_name = 'Vacuum')
s.AddLayerCopy(layer_name = 'GaAsTop', copy_layer_name = 'GaAsBottom')
s.AddLayerCopy(layer_name = 'PECTop', copy_layer_name = 'PECBottom')

s.SetSourceLayer(layer_name = 'GaAsBottom')
s.SetProbeLayer(layer_name = 'VacGap')
s.OptUseQuadgk()
s.SetKParallelIntegral(integral_end = 10)
s.SetThread(num_thread = 4)

s.InitSimulation()
s.IntegrateKParallel()

phi = s.GetPhi()
omega = s.GetOmega()
for i in range(s.GetNumOfOmega()):
  print omega[i], '\t', phi[i]