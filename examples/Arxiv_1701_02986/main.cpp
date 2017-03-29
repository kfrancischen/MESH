#include "setup.h"

static double thetaDerivative(double omega, double T){
  double theta = constant.h_bar * omega / (std::exp(constant.h_bar * omega/constant.k_B/T) - 1);
  return std::pow(theta, 2) * std::exp(constant.h_bar * omega/constant.k_B/T) /constant.k_B / std::pow(T, 2);
}

int main(){
  // first replicate the result in blue curve in Fig2
  Ptr<SimulationPlanar> s_1 = SimulationPlanar::instanceNew();
  s_1->addMaterial("Si", "Si.txt");
  s_1->addMaterial("Vacuum", "Vacuum.txt");
  s_1->addLayer("SiBottom", 0, "Si");
  s_1->addLayer("VacGap", 20e-9, "Vacuum");
  s_1->addLayerCopy("SiTop", "SiBottom");

  s_1->setSourceLayer("SiBottom");
  s_1->setProbeLayer("VacGap");

  s_1->setKParallelIntegral(500);
  s_1->optUseQuadgk();
  s_1->buildRCWA();
  s_1->integrateKParallel();
  double* phi = s_1->getPhi();
  double* omega = s_1->getOmega();
  for(int i = 0; i < s_1->getNumOfOmega(); i++){
    double E = constant.h_bar * omega[i] / constant.q;
    double spectrum = thetaDerivative(omega[i], 300) * phi[i] * constant.q / constant.h_bar;
    std::cout << E << "\t" << spectrum / 1e5 << std::endl;
  }

  // calculate the result in the case of f = 0.98 in Fig2
  double f = 0.98;
  double width = std::sqrt(f * 50e-9 * 50e-9);

  Ptr<SimulationPattern> s_2 = SimulationPattern::instanceNew();
  s_2->setPeriodicity(50e-9, 50e-9);
  s_2->setGx(10);
  s_2->setGy(10);

  s_2->addMaterial("Si", "Si.txt");
  s_2->addMaterial("Vacuum", "Vacuum.txt");
  s_2->addLayer("SiBottom", 0, "Si");
  s_2->setLayerPatternRectangle("SiBottom", "Vacuum", 25e-9, 25e-9, width, width);
  s_2->addLayer("VacGap", 20e-9, "Vacuum");
  s_2->addLayerCopy("SiTop", "SiBottom");

  s_2->setSourceLayer("SiBottom");
  s_2->setProbeLayer("VacGap");
  s_2->getSysInfo();

  s_2->optPrintIntermediate();
  s_2->setKxIntegralSym(20, 60);
  s_2->setKyIntegralSym(20, 60);
  s_2->buildRCWA();
  s_2->integrateKxKy();

  return 0;
}
