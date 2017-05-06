s = SimulationPattern.new();
s:SetLattice(1e-7, 1e-7, 90);
s:SetNumOfG(300);
--s:OptSetLatticeTruncation("Parallelogramic")
s:AddMaterial("Si", "Si.txt");
s:AddMaterial("Vacuum", "Vacuum.txt");
s:AddLayer("SiBottom", 0, "Si");
s:SetLayerPatternCircle("SiBottom", "Vacuum", {50e-9, 50e-9}, 45e-9)
s:SetLayerPatternCircle("SiBottom", "Si", {50e-9, 50e-9}, 30e-9)
--s:SetLayerPatternCircle("SiBottom", "Vacuum", {50e-9, 50e-9}, 15e-9)
s:AddLayer("VacGap", 1e-7, "Vacuum");
s:AddLayer("SiTop", 0, "Si");
--s:AddLayerCopy("SiTop", "SiBottom")
s:SetSourceLayer("SiBottom");
s:SetProbeLayer("VacGap");

s:OptPrintIntermediate();
s:SetKxIntegralSym(100);
s:SetKyIntegralSym(100);
s:InitSimulation();
--s:OutputSysInfo();

----------------------------------------------------------------
-- this part is for MPI
-- in principle there is no need to change this part for your simulation
----------------------------------------------------------------
-- start MPI
local sizeb = buffer.new_buffer(buffer.sizeof(buffer.int))
local rankb = buffer.new_buffer(buffer.sizeof(buffer.int))

MPI.Init()
MPI.Comm_rank(MPI.COMM_WORLD, rankb)
MPI.Comm_size(MPI.COMM_WORLD, sizeb)

local size = buffer.get_typed(sizeb, buffer.int, 0)
local rank = buffer.get_typed(rankb, buffer.int, 0)

if rank == 0 then
  s:OutputSysInfo();
end

status = MPI.Status()
numOfOmega = s:GetNumOfOmega()
omega = s:GetOmega();
-- rank 0 is the master node
if rank == 0 then
  s:IntegrateKxKyMPI(rank, size);
  phi_master = s:GetPhi();
  -- master collects values from slave
  for i = 1,size - 1 do
    local phi_local = buffer.new_buffer(numOfOmega * buffer.sizeof(buffer.double));
    MPI.Recv(phi_local, numOfOmega, MPI.DOUBLE, i, 0, MPI.COMM_WORLD, status);
    for j = 1, numOfOmega do
      phi_master[j] = phi_master[j] + buffer.get_typed(phi_local, buffer.double, j - 1);
    end
  end
  -- output all the phi values from the master
  for i = 1,numOfOmega do
    print(string.format("%e", omega[i]).."\t"..string.format("%e", phi_master[i]));
  end
-- rank 1-size are the slave nodes
else
  s:IntegrateKxKyMPI(rank, size);
  local phi_slave = s:GetPhi();
  -- slave nodes send phi values back to master
  local phi_ = buffer.new_buffer(numOfOmega * buffer.sizeof(buffer.double));
  for i = 1, numOfOmega do
    buffer.set_typed(phi_, buffer.double, i - 1, phi_slave[i]);
  end
  MPI.Send(phi_, numOfOmega, MPI.DOUBLE, 0, 0, MPI.COMM_WORLD);
end

MPI.Finalize();
