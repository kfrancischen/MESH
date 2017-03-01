#include "Rcwa.h"
#include "System.h"
#include "Mesh.h"

using namespace SYSTEM;
using namespace MESH;
int main(){

  // initializing material
  int numOfOmega = 1;
  double omega[numOfOmega];
  dcomplex epsilon[numOfOmega];
  fileLoader("GaAs.txt", omega, epsilon, numOfOmega);
  Material* GaAs = new Material("GaAs", omega, epsilon, numOfOmega);
  fileLoader("Vacuum.txt", omega, epsilon, numOfOmega);
  Material* Vacuum = new Material("Vacuum", omega, epsilon, numOfOmega);

  // initializing layer
  Layer* GaAsBottom = new Layer(GaAs, 0);
  GaAsBottom->addGratingPattern(Vacuum, 5e-7, 1e-6);
  Layer* VacLayer = new Layer(Vacuum, 0);
  GaAsBottom->setIsSource();

  // initializing structure
  Structure* structure = new Structure();
  structure->setPeriodicity(1e-6);
  structure->addLayer(GaAsBottom);
  structure->addLayer(VacLayer);

  // initializing simulation
  SimulationGrating* s = new SimulationGrating();
  s->addStructure(structure);
  s->setGx(2);
  s->setTargetLayerByLayer(VacLayer);
  s->setOutputFile("grating1.txt");
  s->setKxIntegralSym(10);
  s->setKyIntegral(10, 1);
  s->build();
  s->run();
  //std::cout << s->getPhiAtKxKy(0, 0, 0) << std::endl;

  return 0;
}
