#define ARMA_DONT_USE_WRAPPER
#include "Rcwa.h"
#include "System.h"
#include "Mesh.h"

#include <armadillo>
using namespace SYSTEM;
using namespace MESH;
int main(){
  RCWA::RCWAMatrix S;
  RCWA::initSMatrix(2, &S);
//  std::cout << S << std::endl << S(arma::span(0, 1), arma::span(0, 1));
  RCWA::RCWAMatrix T = S;
  //T(arma::span(0, 1), arma::span(0, 1)) = arma::ones<RCWA::RCWAMatrix>(2, 2);
  //std::cout << S << std::endl << S(arma::span(0, 1), arma::span(0, 1));

  RCWA::RCWAMatrices S_vec;
  S_vec.push_back(S);
  S_vec.push_back(S);

  S_vec[0] *= 2;
  std::cout << S_vec[1]<< std::endl << S_vec[0] << std::endl;

  std::string name = "abc";
  Material* material = new Material(name);
  Layer* layer = new Layer(material);
  Structure* structure = new Structure();
  structure->addLayer(layer);
  std::cout << structure->getNumOfLayer() << std::endl;

  SimulationPlanar* simulation = new SimulationPlanar();
  simulation->addStructure(structure);
  return 0;
}
