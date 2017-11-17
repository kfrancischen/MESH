from MESH import SimulationGrating

s = SimulationGrating()
s.SetLattice(1e-6)
s.AddMaterial("Au", "fullGold.txt")
s.AddMaterial("Vacuum", "fullVacuum.txt")

s.AddLayer("BottomAir", 0, "Vacuum")
s.AddLayer("GoldSubstrateBottom", 0.5e-6, "Au")
s.AddLayer("GoldGratingBottom", 5e-6,"Au")
s.SetLayerPatternGrating("GoldGratingBottom", "Vacuum", 0.5e-6, 0.2e-6)

s.AddLayer("VacGap", 1e-6, "Vacuum")
s.AddLayerCopy("GoldGratingTop", "GoldGratingBottom")
s.AddLayerCopy("GoldSubstrateTop", "GoldSubstrateBottom")
s.AddLayerCopy("TopAir", "BottomAir")

s.SetSourceLayer("GoldSubstrateBottom")
s.SetSourceLayer("GoldGratingBottom")
s.SetProbeLayer("VacGap")

s.OptPrintIntermediate()
s.SetNumOfG(50)
s.SetKxIntegralSym(2)
s.SetKyIntegralSym(2, 5)
s.InitSimulation()


"""
-- this part is for MPI
-- in principle there is no need to change this part for your simulation
"""
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0:
  s.OutputSysInfo()

numOfOmega = s.GetNumOfOmega()

if rank == 0:
  s.IntegrateKxKyMPI(rank, size);
  phi_master = s.GetPhi()
  phi = list(phi_master)
  for i in range(1, size):
    phi_slave = comm.recv(source = i, tag = 11)
    for j in range(numOfOmega):
      phi[j] = phi[j] + phi_slave[j]

  for i in range(numOfOmega):
    print phi[i]
  
else:
  s.IntegrateKxKyMPI(rank, size)
  phi_slave = s.GetPhi()
  comm.send(phi_slave, dest = 0, tag = 11)
