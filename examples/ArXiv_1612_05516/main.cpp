#include "setup.h"

int main(){

  // initializing material
  Ptr<FileLoader> fileLoader = FileLoader::instanceNew();
  fileLoader->load("fullGold.txt");
  Ptr<Material> Gold = Material::instanceNew("Au", fileLoader->getOmegaList(), fileLoader->getEpsilonList(), fileLoader->getNumOfOmega());
  fileLoader->load("fullVacuum.txt");
  Ptr<Material> Vacuum = Material::instanceNew("Vacuum", fileLoader->getOmegaList(), fileLoader->getEpsilonList(), fileLoader->getNumOfOmega());

  // initializing layer
  Ptr<Layer> GoldLayerBottomSub = Layer::instanceNew("GoldLayerBottomSub", Gold, 0);
  Ptr<Layer> GoldLayerBottomGrating = Layer::instanceNew("GoldLayerBottomGrating", Gold, 4.7e-6);

  GoldLayerBottomSub->setIsSource();
  GoldLayerBottomGrating->addGratingPattern(Vacuum, 0.1e-6, 0.2e-6);
  GoldLayerBottomGrating->setIsSource();

  Ptr<Layer> vacGap = Layer::instanceNew("vacGap", Vacuum, 1e-6);

  Ptr<Layer> GoldLayerTopSub = Layer::instanceNew("GoldLayerTopSub", Gold, 0);
  Ptr<Layer> GoldLayerTopGrating = Layer::instanceNew("GoldLayerTopGrating", Gold, 4.7e-6);
  GoldLayerTopGrating->addGratingPattern(Vacuum, 0.1e-6, 0.2e-6);

  // initializing structure
  Ptr<Structure> structure = Structure::instanceNew();
  structure->setPeriodicity(0.5e-6);
  structure->addLayer(GoldLayerBottomSub);
  structure->addLayer(GoldLayerBottomGrating);
  structure->addLayer(vacGap);
  structure->addLayer(GoldLayerTopGrating);
  structure->addLayer(GoldLayerTopSub);

  // initializing simulation
  Ptr<SimulationGrating> s = SimulationGrating::instanceNew();
  s->addStructure(structure);
  s->setGx(50);
  s->setTargetLayerByLayer(vacGap);
  s->setOutputFile("gold_to_vac.txt");
  s->setKxIntegralSym(500);
  s->setKyIntegral(200, 5);
  s->build();
  s->run();

  //std::cout << s->getPhiAtKxKy(0, 0.2, 0.2) << std::endl;
  return 0;
}
