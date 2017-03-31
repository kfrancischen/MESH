#include "setup.h"

int main(){
  // initializing simulation
  Ptr<SimulationPlanar> s = SimulationPlanar::instanceNew();
  s->addMaterial("hBN", "hBN.txt");
  s->addMaterial("MCT", "MCT.txt");
  s->addMaterial("Vacuum", "Vacuum.txt");
  s->addMaterial("PEC", "PEC.txt");

  // add layer from bottom to the top
  s->addLayer("PECBottom", 0, "PEC");
  s->addLayer("hBNLayer", 5e-6, "hBN");
  s->addLayer("VacGap", 1e-8, "Vacuum");
  s->addLayer("MCTLayer", 5e-6, "MCT");
  s->addLayerCopy("PECTop", "PECBottom");

  // set source and probe
  s->setSourceLayer("hBNLayer");
  s->setProbeLayer("VacGap");

  // set simulation
  s->setKParallelIntegral(100);
  s->optUseQuadgk();
  s->setThread(4);
  s->buildRCWA();
  s->integrateKParallel();
  return 0;
}
