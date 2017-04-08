This is an example that has been used in the paper [Phys. Rev. Applied 6, 024014, 2016](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.6.024014), where `hBN` has a diagonal dielectric function as (top $5$ rows)
```tex
2.127000e+14	1.022020e+01	3.983602e-01	1.022020e+01	3.983602e-01	4.251609e+00	3.151862e-01
2.128640e+14	1.023858e+01	3.932224e-01	1.023858e+01	3.932224e-01	4.248610e+00	3.151862e-01
2.130281e+14	1.025696e+01	3.880846e-01	1.025696e+01	3.880846e-01	4.253492e+00	3.151862e-01
2.131921e+14	1.027534e+01	3.829468e-01	1.027534e+01	3.829468e-01	4.258538e+00	3.151862e-01
2.133561e+14	1.029372e+01	3.778090e-01	1.029372e+01	3.778090e-01	4.263583e+00	3.151862e-01
```
These values correspond to $\omega, \text{real}(\epsilon_{xx}), \text{imag}(\epsilon_{xx}), \text{real}(\epsilon_{yy}), \text{imag}(\epsilon_{yy}), \text{real}(\epsilon_{zz}), \text{imag}(\epsilon_{zz})$, respectively. And the following is a script computing the flux between MCT and hBN in the presence of a $10~\text{nm}$ vacuum gap.
```lua
s = SimulationPlanar.new()
s:AddMaterial("hBN", "hBN.txt")
s:AddMaterial("MCT", "MCT.txt")
s:AddMaterial("PEC", "PEC.txt")
s:AddMaterial("Vacuum", "Vacuum.txt")

s:AddLayer("PECBottom", 0, "PEC");
s:AddLayer("hBNLayer", 5e-6, "hBN");
s:AddLayer("VacLayer", 1e-8, "Vacuum");
s:AddLayer("MCTLayer", 5e-6, "MCT");
s:AddLayerCopy("PECTop", "PECBottom");

s:SetSourceLayer("hBNLayer");
s:SetProbeLayer("VacLayer");

s:SetKParallelIntegral(100);
s:SetThread(4);
s:OptUseQuadgk();
s:BuildRCWA();
s:IntegrateKParallel();
phi = s:GetPhi();
omega = s:GetOmega();
for i = 1,s:GetNumOfOmega(), 1 do
  print(string.format("%e", omega[i]).."\t"..string.format("%e", phi[i]));
end
```
The output of the Lua file is the same as the output from the original MATLAB code. The output from the MATLAB code is [Electron_5000_10_5000.mat](https://github.com/kfrancischen/MESH/tree/master/examples/PRApplied_6_024014_2016).