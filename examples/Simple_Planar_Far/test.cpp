#include "setup.h"

int main(){
  // initializing material
  Ptr<SimulationPlanar> s = SimulationPlanar::instanceNew();
  s->addMaterial("GaAs", "GaAs.txt");
  s->addMaterial("Vacuum", "Vacuum.txt");
  s->addMaterial("PEC", "PEC.txt");

  // initializing layer, from bottom to top
  s->addLayer("PECLayer", 0, "PEC");
  s->addLayer("GaAsLayer", 1e-6, "GaAs");
  s->addLayer("VacuumLayer", 0, "Vacuum");

  // set source and probe
  s->setSourceLayer("GaAsLayer");
  s->setProbeLayer("VacuumLayer");

  s->setKParallelIntegral(1.0);
  s->optUseQuadgk();
  s->setThread(4);
  s->buildRCWA();
  s->integrateKParallel();

  double* phi = s->getPhi();
  for(int i = 0; i < s->getNumOfOmega(); i++){
    std::cout << phi[i] << std::endl;
  }
  return 0;
}
