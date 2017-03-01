#include "Rcwa.h"
#include "System.h"
#include "Mesh.h"

using namespace SYSTEM;
using namespace MESH;
int main(){
  // initializing material
  int numOfOmega = 199;
  double omega[numOfOmega];
  dcomplex epsilon[numOfOmega];
  fileLoader("GaAs.txt", omega, epsilon, numOfOmega);
  Material* GaAs = new Material("GaAs", omega, epsilon, numOfOmega);
  fileLoader("Vacuum.txt", omega, epsilon, numOfOmega);
  Material* Vacuum = new Material("Vacuum", omega, epsilon, numOfOmega);
  fileLoader("PEC.txt", omega, epsilon, numOfOmega);
  Material* PEC = new Material("PEC", omega, epsilon, numOfOmega);

  // initializing layer
  Layer* PECBottom = new Layer(PEC, 0);
  Layer* PECTop = new Layer(PEC, 0);
  Layer* GaAsBottom = new Layer(GaAs, 1e-6);
  Layer* VacLayer = new Layer(Vacuum, 1e-8);
  Layer* GaAsTop = new Layer(GaAs, 1e-6);
  GaAsBottom->setIsSource();

  // initializing structure
  Structure* structure = new Structure();
  structure->addLayer(PECBottom);
  structure->addLayer(GaAsBottom);
  structure->addLayer(VacLayer);
  structure->addLayer(GaAsTop);
  structure->addLayer(PECTop);

  // initializing simulation
  SimulationPlanar* s = new SimulationPlanar();
  s->addStructure(structure);
  s->setTargetLayerByLayer(VacLayer);
  s->setKxIntegral(10);
  s->setOutputFile("test.txt");
  s->build();
  s->run();

  return 0;
}
