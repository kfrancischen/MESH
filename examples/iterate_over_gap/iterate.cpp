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
  s->setKParallelIntegral(500);
  s->optUseQuadgk();
  s->setThread(4);
  // this function will interate the gap from 10-100 nm
  for(int i = 10; i < 100; i += 10){
    s->setLayerThickness("VacGap", i * 1e-9);
    s->buildRCWA();
    s->integrateKParallel();
    s->saveToFile("test" + std::to_string(i) + ".txt");
  }
  return 0;
}
