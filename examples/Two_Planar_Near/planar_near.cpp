#include "setup.h"

int main(){
  // initializing simulation
  Ptr<SimulationPlanar> s = SimulationPlanar::instanceNew();
  s->addMaterial("GaAs", "GaAs.txt");
  s->addMaterial("Vacuum", "Vacuum.txt");
  s->addMaterial("PEC", "PEC.txt");

  // add layer from bottom to the top
  s->addLayer("PECBottom", 0, "PEC");
  s->addLayer("GaAsBottom", 1e-6, "GaAs");
  s->addLayer("VacGap", 1e-8, "Vacuum");
  s->addLayerCopy("GaAsTop", "GaAsBottom");
  s->addLayerCopy("PECTop", "PECBottom");

  // set source and probe
  s->setSourceLayer("GaAsBottom");
  s->setProbeLayer("VacGap");

  // set simulation
  s->setKParallelIntegral(10);
  s->useQuadgk();
  s->build();
  s->setThread(4);
  s->runNaive();
  s->saveToFile("test.txt");
  double* Phi = s->getPhi();
  for(int i = 0; i < s->getNumOfOmega(); i++){
    std::cout << Phi[i] << std::endl;
  }
  return 0;
}
