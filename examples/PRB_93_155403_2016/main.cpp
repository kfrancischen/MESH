#include "setup.h"

int main(){
  // initialize simulation
  Ptr<SimulationGrating> s = SimulationGrating::instanceNew();
  s->setPeriodicity(1e-6);
  s->setGx(50);
  // initialize materials
  s->addMaterial("Au", "fullGold.txt");
  s->addMaterial("Vacuum", "fullVacuum.txt");
  // initialize layers
  s->addLayer("BottomAir", 0, "Vacuum");
  s->addLayer("GoldSubstrateBottom", 0.5e-6, "Au");
  s->addLayer("GoldGratingBottom", 5e-6,"Au");
  s->setLayerPatternGrating("GoldGratingBottom", "Vacuum", 0.5e-6, 0.2e-6);

  s->addLayer("VacGap", 1e-6, "Vacuum");
  s->addLayerCopy("GoldGratingTop", "GoldGratingBottom");
  s->addLayerCopy("GoldSubstrateTop", "GoldSubstrateBottom");
  s->addLayerCopy("TopAir", "BottomAir");
  // set source and probe
  s->setSourceLayer("GoldGratingBottom");
  s->setSourceLayer("GoldSubstrateBottom");

  s->setProbeLayer("VacGap");

  // set output

  s->setKxIntegralSym(500);
  s->setKyIntegralSym(200, 5);
  //s->useInverseRule();
  s->getSysInfo();
  s->optPrintIntermediate();
  s->setThread(4);
  s->build();
  s->run();
  s->saveToFile("gold_to_vac.txt");
  // std::cout << s->getPhiAtKxKy(0, 0, 0) << std::endl;
  return 0;
}
