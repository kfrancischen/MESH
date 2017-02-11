#define ARMA_DONT_USE_WRAPPER
#include "Rcwa.h"
#include "System.h"
#include "Mesh.h"

#include <armadillo>
using namespace SYSTEM;
using namespace MESH;
int main(){

  // initializing material
  int numOfOmega = 2;
  dcomplex epsilonGaAs[numOfOmega] = {dcomplex(1, -1), dcomplex(1, -1)};
  dcomplex epsilonVac[numOfOmega] = {dcomplex(1, 1e-10), dcomplex(1, 1e-10)};
  double omegaList[numOfOmega] = {1.0e15, 2.0e15};

  Material* GaAs = new Material("GaAs", epsilonGaAs, omegaList, numOfOmega);
  Material* Vacuum = new Material("Vacuum", epsilonVac, omegaList, numOfOmega);

  // initializing layer
  Layer* GaAsLayer = new Layer(GaAs, 0, ISSOURCE_);
  Layer* VacLayer = new Layer(Vacuum, 0);

  // initializing structure
  Structure* structure = new Structure();
  structure->addLayer(GaAsLayer);
  structure->addLayer(VacLayer);

  std::cout << structure->getNumOfLayer() << std::endl;

  // initializing simulation
  SimulationPlanar* simulation = new SimulationPlanar();
  simulation->addStructure(structure);
  simulation->setTargetLayerByLayer(VacLayer);
  return 0;
}
