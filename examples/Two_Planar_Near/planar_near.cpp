#include "setup.h"

int main(){
  // initializing material
  Ptr<FileLoader> fileLoader = FileLoader::instanceNew();
  fileLoader->load("GaAs.txt");
  Ptr<Material> GaAs = Material::instanceNew("GaAs", fileLoader->getOmegaList(), fileLoader->getEpsilonList(), fileLoader->getNumOfOmega());
  fileLoader->load("Vacuum.txt");
  Ptr<Material> Vacuum = Material::instanceNew("Vacuum", fileLoader->getOmegaList(), fileLoader->getEpsilonList(), fileLoader->getNumOfOmega());
  fileLoader->load("PEC.txt");
  Ptr<Material> PEC = Material::instanceNew("PEC", fileLoader->getOmegaList(), fileLoader->getEpsilonList(), fileLoader->getNumOfOmega());

  // initializing layer
  Ptr<Layer> PECLayer = Layer::instanceNew("PECLayer", PEC, 0);
  Ptr<Layer> GaAsBottom = Layer::instanceNew("GaAsBottom", GaAs, 1e-6);
  Ptr<Layer> VacLayer = Layer::instanceNew("VacLayer", Vacuum, 1e-8);
  Ptr<Layer> GaAsTop = Layer::instanceNew("GaAsTop", GaAs, 1e-6);
  GaAsBottom->setIsSource();

  // initializing structure
  Ptr<Structure> structure = Structure::instanceNew();
  structure->addLayer(PECLayer);
  structure->addLayer(GaAsBottom);
  structure->addLayer(VacLayer);
  structure->addLayer(GaAsTop);
  structure->addLayer(PECLayer);

  // initializing simulation
  Ptr<SimulationPlanar> s = SimulationPlanar::instanceNew();
  s->addStructure(structure);
  s->setTargetLayerByLayer(VacLayer);
  s->setKxIntegral(10);
  s->setOutputFile("test.txt");
  s->build();
  s->run();

  return 0;
}
