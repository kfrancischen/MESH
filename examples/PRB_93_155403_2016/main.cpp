#include "setup.h"

int main(){

  // initializing material
  int numOfOmega = 1;
  double omega[numOfOmega];
  dcomplex epsilon[numOfOmega];
  fileLoader("fullGold.txt", omega, epsilon, numOfOmega);
  Ptr<Material> Gold = Material::instanceNew("Au", omega, epsilon, numOfOmega);
  fileLoader("fullVacuum.txt", omega, epsilon, numOfOmega);
  Ptr<Material> Vacuum = Material::instanceNew("Vacuum", omega, epsilon, numOfOmega);
  // initializing layer
  Ptr<Layer> vacLayer = Layer::instanceNew("VacLayer", Vacuum, 0);
  Ptr<Layer> GoldLayerBottomSub = Layer::instanceNew("GoldLayerBottomSub", Gold, 0.5e-6);
  Ptr<Layer> GoldLayerBottomGrating = Layer::instanceNew("GoldLayerBottomGrating", Gold, 5e-6);
  //GoldLayerBottomSub->setIsSource();
  GoldLayerBottomGrating->addGratingPattern(Vacuum, 0.4e-6, 0.6e-6);
  GoldLayerBottomGrating->setIsSource();

  Ptr<Layer> vacGap = Layer::instanceNew("VacGap", Vacuum, 1e-6);

  Ptr<Layer> GoldLayerTopSub = Layer::instanceNew("GoldLayerTopSub", Gold, 0.5e-6);
  Ptr<Layer> GoldLayerTopGrating = Layer::instanceNew("GoldLayerTopGrating ", Gold, 5e-6);
  GoldLayerTopGrating->addGratingPattern(Vacuum, 0.4e-6, 0.6e-6);

  // initializing structure
  Ptr<Structure> structure = Structure::instanceNew();
  structure->setPeriodicity(1e-6);
  structure->addLayer(vacLayer);
  structure->addLayer(GoldLayerBottomSub);
  structure->addLayer(GoldLayerBottomGrating);
  structure->addLayer(vacGap);
  structure->addLayer(GoldLayerTopGrating);
  structure->addLayer(GoldLayerTopSub);
  structure->addLayer(vacLayer);

  // initializing simulation
  Ptr<SimulationGrating> s = SimulationGrating::instanceNew();
  s->addStructure(structure);
  s->setGx(50);
  s->setTargetLayerByLayer(vacGap);
  s->setOutputFile("gold_to_vac.txt");
  s->setKxIntegralSym(10);
  s->setKyIntegral(10, 1);
  s->build();
  s->run();

  //std::cout << s->getPhiAtKxKy(0, 0, 0) << std::endl;
  return 0;
}
