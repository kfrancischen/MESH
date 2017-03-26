#include "setup.h"

int main(){
  // initializing simulation
  Ptr<SimulationGrating> s = SimulationGrating::instanceNew();
  s->setPeriodicity(0.5e-6);
  s->setGx(50);
  // initialize materials
  s->addMaterial("Au", "fullGold.txt");
  s->addMaterial("Vacuum", "fullVacuum.txt");

  // initialize layers
  s->addLayer("GoldBottomSubstrate", 0, "Au");
  s->addLayer("GoldBottomGrating", 4.7e-6,"Au");
  s->setLayerPatternGrating("GoldBottomGrating", "Vacuum", 0.1e-6, 0.2e-6);
  s->addLayer("VacGap", 1e-6, "Vacuum");
  s->addLayerCopy("GoldLayerTopGrating", "GoldBottomGrating");
  s->addLayerCopy("GoldLayerTopSubstrate", "GoldBottomSubstrate");

  s->setSourceLayer("GoldBottomSubstrate");
  s->setSourceLayer("GoldBottomGrating");
  s->setProbeLayer("VacGap");

  // setoutput simulation
  s->setKxIntegralSym(500);
  s->setKyIntegralSym(200, 5);
  s->getSysInfo();
  s->setThread(4);
  s->build();
  s->run();
  s->saveToFile("gold_to_vac.txt");

  //std::cout << s->getPhiAtKxKy(0, 0.2, 0.2) << std::endl;
  return 0;
}
