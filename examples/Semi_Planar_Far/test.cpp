#include "setup.h"

int main(){
  // initializing material
  int numOfOmega = 199;
  double omega[numOfOmega];
  dcomplex epsilon[numOfOmega];
  fileLoader("GaAs.txt", omega, epsilon, numOfOmega);
  Ptr<Material> GaAs = Material::instanceNew("GaAs", omega, epsilon, numOfOmega);
  fileLoader("Vacuum.txt", omega, epsilon, numOfOmega);
  Ptr<Material> Vacuum = Material::instanceNew("Vacuum", omega, epsilon, numOfOmega);
  fileLoader("PEC.txt", omega, epsilon, numOfOmega);
  Ptr<Material> PEC = Material::instanceNew("PEC", omega, epsilon, numOfOmega);

  // initializing layer
  Ptr<Layer> PECLayer = Layer::instanceNew(PEC, 0);
  Ptr<Layer> GaAsLayer = Layer::instanceNew(GaAs, 1e-6);
  Ptr<Layer> VaccumLayer = Layer::instanceNew(Vacuum, 0);
  GaAsLayer->setIsSource();

  // initializing structure
  Ptr<Structure> structure = Structure::instanceNew();
  structure->addLayer(PECLayer);
  structure->addLayer(GaAsLayer);
  structure->addLayer(VaccumLayer);

  // initializing simulation
  Ptr<SimulationPlanar> s = SimulationPlanar::instanceNew();
  s->addStructure(structure);
  s->setTargetLayerByLayer(VaccumLayer);
  s->setKxIntegral(1.0);
  s->setOutputFile("test_output.txt");
  s->build();
  s->run();

  return 0;
}
