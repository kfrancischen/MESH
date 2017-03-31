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
s:BuildRCWA();

-- start MPI
local sizeb = buffer.new_buffer(buffer.sizeof(buffer.int))
local rankb = buffer.new_buffer(buffer.sizeof(buffer.int))

MPI.Init()
MPI.Comm_rank(MPI.COMM_WORLD, rankb)
MPI.Comm_size(MPI.COMM_WORLD, sizeb)

local size = buffer.get_typed(sizeb, buffer.int, 0)
local rank = buffer.get_typed(rankb, buffer.int, 0)

if rank == 0 then
  s:GetSysInfo();
end

status = MPI.Status()
numOfOmega = s:GetNumOfOmega()

if rank == 0 then
  print(rank);
  s:IntegrateKxKyMPI(rank, size);
  phi_master = s:GetPhi();
  for i = 1,size - 1 do
    local phi_local = buffer.new_buffer(buffer.sizeof(numOfOmega * buffer.double));
    MPI.Recv(phi_local, numOfOmega, MPI.DOUBLE, 0, MPI.COMM_WORLD, status);
    for j = 1, numOfOmega do
      phi_master[j] = phi_master[j] + phi_local[j];
    end
  end
  for i = 1,numOfOmega do
    print(string.format("%e", phi_master[i]));
  end
else
  print(rank);
  s:IntegrateKxKyMPI(rank, size);
  local phi_slave = s:GetPhi();
  MPI.Send(phi_slave, numOfOmega, MPI.DOUBLE, 0, 0, MPI.COMM_WORLD);
end

MPI.Finalize();
