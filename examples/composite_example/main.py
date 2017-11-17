from MESH import SimulationPattern

s = SimulationPattern()
s.SetLattice(1e-7, 1e-7, 90)
s.SetNumOfG(300)
s.AddMaterial("Si", "Si.txt")
s.AddMaterial("Vacuum", "Vacuum.txt")
s.AddLayer("SiBottom", 0, "Si")
s.SetLayerPatternCircle("SiBottom", "Vacuum", (50e-9, 50e-9), 45e-9)
s.SetLayerPatternCircle("SiBottom", "Si", (50e-9, 50e-9), 30e-9)
s.AddLayer("VacGap", 1e-7, "Vacuum")
s.AddLayer("SiTop", 0, "Si")
s.SetSourceLayer("SiBottom")
s.SetProbeLayer("VacGap")

s.OptPrintIntermediate()
s.SetKxIntegralSym(100)
s.SetKyIntegralSym(100)
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
