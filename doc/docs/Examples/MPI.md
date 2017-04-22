This is an MPI example, rewriting the same example in [Two Gratings Near-field](gratingNearField.md)
```lua
s = SimulationGrating.new();
s:SetPeriodicity(1e-6);
s:AddMaterial("Au", "fullGold.txt");
s:AddMaterial("Vacuum", "fullVacuum.txt");

s:AddLayer("BottomAir", 0, "Vacuum");
s:AddLayer("GoldSubstrateBottom", 0.5e-6, "Au");
s:AddLayer("GoldGratingBottom", 5e-6,"Au");
s:SetLayerPatternGrating("GoldGratingBottom", "Vacuum", 0.5e-6, 0.2e-6);

s:AddLayer("VacGap", 1e-6, "Vacuum");
s:AddLayerCopy("GoldGratingTop", "GoldGratingBottom");
s:AddLayerCopy("GoldSubstrateTop", "GoldSubstrateBottom");
s:AddLayerCopy("TopAir", "BottomAir");

s:SetSourceLayer("GoldSubstrateBottom");
s:SetSourceLayer("GoldGratingBottom");
s:SetProbeLayer("VacGap");

s:OptPrintIntermediate();
s:SetGx(50);
s:SetKxIntegralSym(500);
s:SetKyIntegralSym(200, 5);
s:InitSimulation();
```
The above are the same as the standard Lua script. The following code configures the MPI.
```lua
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
```
The above code uses master-slave manner, where rank $0$ is the master node that collects all the $\Phi(\omega)$ from the other slave nodes, and then sum them together. The code can be run by:
```bash
mpirun -np 40 meshMPI main.lua
```