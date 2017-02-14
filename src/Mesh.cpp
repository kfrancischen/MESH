/* Copyright (C) 2016-2018, Stanford University
 * This file is part of MESH
 * Written by Kaifeng Chen (kfchen@stanford.edu)
 *
 * MESH is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * MESH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "Mesh.h"

namespace MESH{
  /*==============================================
  This function loads data from disk
  @args:
  fileName: the name of the input file
  omega: the input omega list
  epsilon: the input epsilon list
  size: the size of omega
  ==============================================*/
  void fileLoader(std::string fileName, double* omega, dcomplex* epsilon, int size){
    std::ifstream inputFile(fileName);
    std::string line;
    double real, imag;
    for(size_t i = 0; i < size; i++){
      inputFile >> omega[i] >> real >> imag;
      epsilon[i] = dcomplex(real, -imag);
    }
    inputFile.close();
  }

  /*==============================================
  This function saves data from disk
  @args:
  fileName: the name of the output file
  omega: the output omega list
  epsilon: the output epsilon list
  size: the size of omega
  ==============================================*/
  void fileSaver(std::string fileName, double* omega, double* Phi, int size){
    std::ofstream outputFile(fileName);
    for(size_t i = 0; i < size; i++){
      outputFile << omega[i] << "\t" << Phi[i] << std::endl;
    }
    outputFile.close();
  }

  /*==============================================
  This function wraps the data for quad_legendre
  @args:
  kx: the kx value (normalized)
  data: wrapper for all the arguments wrapped in wrapper
  ==============================================*/
  double wrapperFun(double kx, void* data){
    ArgWrapper wrapper = *(ArgWrapper*) data;
    return kx * poyntingFlux(
      wrapper.omega,
      &(wrapper.thicknessList),
      kx,
      0,
      &(wrapper.EMatrices),
      &(wrapper.grandImaginaryMatrices),
      &(wrapper.dielectricMatrixInverse),
      &(wrapper.Gx_mat),
      &(wrapper.Gy_mat),
      &(wrapper.sourceList),
      wrapper.targetLayer,
      1
    );
  }
  /*======================================================
  Implementaion of the parent simulation super class
  =======================================================*/
  Simulation::Simulation() : nGx_(0), nGy_(0), numOfCore_(1), numOfOmega_(0), structure_(nullptr),
  Phi_(nullptr), omegaList_(nullptr), kxStart_(0), kxEnd_(0), kyStart_(0), kyEnd_(0), numOfKx_(0), numOfKy_(0)
  {
    period_ = new double[2];
    period_[0] = 0;
    period_[1] = 0;
  }

  Simulation::~Simulation(){
    delete structure_;
    delete[] Phi_;
    delete[] omegaList_;
    delete[] period_;
  }

  /*==============================================
  This function adds structure to the simulation
  @args:
  structure: the structure of the simulation
  ==============================================*/
  void Simulation::addStructure(Structure* structure){
    structure_ = structure;
  }

  /*==============================================
  This function enables MPI
  @args:
  numOfCore: the number of cores in the simulation
  ==============================================*/
  void Simulation::enableMPI(int numOfCore){
    numOfCore_ = numOfCore;
  }

  /*==============================================
  This function sets number of positive Gx
  @args:
  Gx: number of positive Gx
  ==============================================*/
  void Simulation::setGx(int Gx){
    nGx_ = Gx;
  }

  /*==============================================
  This function sets number of positive Gy
  @args:
  Gy: number of positive Gy
  ==============================================*/
  void Simulation::setGy(int Gy){
    nGy_ = Gy;
  }

  /*==============================================
  This function sets the target layer by layer index
  @args:
  index: target layer index
  ==============================================*/
  void Simulation::setTargetLayerByIndex(int index){
    targetLayer_ = index;
  }

  /*==============================================
  This function sets the target layer by layer
  @args:
  index: target layer
  ==============================================*/
  void Simulation::setTargetLayerByLayer(Layer* layer){
    for(size_t i = 0; i < structure_->getNumOfLayer(); i++){
      if(structure_->getLayerByIndex(i) == layer){
        targetLayer_ = i;
        return;
      }
    }
  }

  /*==============================================
  This function cleans up the simulation
  ==============================================*/
  void Simulation::resetSimulation(){
    for(size_t i = 0; i < EMatricesVec_.size(); i++){
      EMatricesVec_[i].clear();
      grandImaginaryMatricesVec_.clear();
      dielectricMatrixInverseVec_[i].clear();
    }
    EMatricesVec_.clear();
    grandImaginaryMatricesVec_.clear();
    dielectricMatrixInverseVec_.clear();
  }

  /*==============================================
  This function gets the structure
  ==============================================*/
  Structure* Simulation::getStructure(){
    return structure_;
  }

  /*==============================================
  This function gets the omega array
  ==============================================*/
  double* Simulation::getOmegaList(){
    return omegaList_;
  }

  /*==============================================
  This function gets the \Phi
  ==============================================*/
  double* Simulation::getPhi(){
    return Phi_;
  }

  /*==============================================
  This function gets the periodicity
  ==============================================*/
  double* Simulation::getPeriodicity(){
    return period_;
  }

  /*==============================================
  This function gets the Phi at given kx and ky
  @args:
  omegaIndex: the index of omega
  kx: the kx value, normalized
  ky: the ky value, normalized
  @note
  used by grating and patterning
  ==============================================*/
  double Simulation::getPhiAtKxKy(int omegaIdx, double kx, double ky){
    int N = getN(nGx_, nGy_);
    return POW3(omegaList_[omegaIdx] / datum::c_0) / POW3(datum::pi) / 2.0 *
      poyntingFlux(omegaList_[omegaIdx] / datum::c_0, &thicknessListVec_, kx, ky, &(EMatricesVec_[omegaIdx]),
      &(grandImaginaryMatricesVec_[omegaIdx]), &(dielectricMatrixInverseVec_[omegaIdx]), &Gx_mat_, &Gy_mat_,
      &sourceList_, targetLayer_,N);
  }

  /*==============================================
  This function builds up the matrices
  ==============================================*/
  void Simulation::build(){
    period_ = structure_->getPeriodicity();
    Layer* firstLayer = structure_->getLayerByIndex(0);
    Material* backGround = firstLayer->getBackGround();
    numOfOmega_ = backGround->getNumOfOmega();
    omegaList_ = backGround->getOmegaList();
    Phi_ = new double[numOfOmega_];
    EMatricesVec_.resize(numOfOmega_);
    grandImaginaryMatricesVec_.resize(numOfOmega_);
    dielectricMatrixInverseVec_.resize(numOfOmega_);

    RCWAMatricesVec dielectricMatrixVec(numOfOmega_), dielectricImMatrixVec(numOfOmega_);
    int numOfLayer = structure_->getNumOfLayer();
    int N = getN(nGx_, nGy_);

    RCWAMatrix onePadding1N = eye<RCWAMatrix>(N, N);

    for(size_t i = 0; i < numOfLayer; i++){
      Layer* layer = structure_->getLayerByIndex(i);
      dcomplex* epsilon = (layer->getBackGround())->getEpsilon();

      switch (layer->getPattern()) {
        /*************************************
        /* if the pattern is a plane */
        case PLANAR_:{
          for(size_t j = 0; j < numOfOmega_; j++){
            dielectricMatrixVec[j].push_back(onePadding1N * epsilon[j]);
            dielectricMatrixInverseVec_[j].push_back(onePadding1N * dcomplex(1, 0) / epsilon[j]);
            dielectricImMatrixVec[j].push_back(onePadding1N * (epsilon[j]).imag());
          }
          getGMatrices(nGx_, nGy_, period_, &Gx_mat_, &Gy_mat_, NO_);
          break;
        }
        /*************************************
        /* if the pattern is a grating (1D) */
        /************************************/
        case GRATING_:{
          dcomplex IMAG_I = dcomplex(0, 1);
          dcomplex* epsilonBG = backGround->getEpsilon();
          RCWAMatrix G_row(1, N), G_col(N, 1);
          for(size_t i = 0; i < N; i++){
            G_row(0, i) = -i * 2 * datum::pi / period_[0];
            G_col(i, 0) = -G_row(0, i);
          }

          RCWAMatrix G_mat = toeplitz(G_col, G_row);
          int numOfMaterial = layer->getNumOfMaterial();

          RCWAVector centerVec(numOfMaterial), widthVec(numOfMaterial);
          int count = 0;

          for (size_t j = 0; j < numOfOmega_; j++) {
            RCWAMatrix dielectricMatrix(N, N, fill::zeros), dielectricImMatrix(N, N, fill::zeros);
            count = 0;

            for(const_PatternIter it = layer->getArg1Begin(); it != layer->getArg1End(); it++){
              centerVec(i) = (it->first + it->second) / 2;
              widthVec(i) = it->second - it->first;
              count++;
            }
            count = 0;
            for(const_MaterialIter it = layer->getVecBegin(); it != layer->getVecEnd(); it++){
              dcomplex epsilon = (*it)->getEpsilonAtIndex(j);
              dielectricMatrix += exp(IMAG_I * G_mat * centerVec(count)) * (epsilon - epsilonBG[j])
                * widthVec(count) % sinh(G_mat / (2*datum::pi) * widthVec(count));

              dielectricImMatrix += exp(IMAG_I * G_mat * centerVec(count)) * (epsilon - epsilonBG[j])
                * widthVec(count) % sinh(G_mat / (2*datum::pi) * widthVec(count));
              count++;
            }

            dielectricMatrix /= period_[0];
            dielectricImMatrix /= period_[0];
            dielectricMatrix += epsilonBG[j] * eye<RCWAMatrix>(N, N);
            dielectricImMatrix += epsilonBG[j].imag() * eye<RCWAMatrix>(N, N);

            dielectricMatrixVec[j].push_back(dielectricMatrix);
            dielectricImMatrixVec[j].push_back(dielectricImMatrix); 
            dielectricMatrixInverseVec_[j].push_back(dielectricMatrix.i());

          }
          break;
        }

        /*************************************
        /* if the pattern is a circle (2D) */
        /************************************/
        case CIRCLE_:{
          // TODO
          break;
        }

        /*************************************
        /* if the pattern is a rectangle (2D) */
        /************************************/
        case RECTANGLE_:{
          // TODO
          break;
        }

        default: break;
      }
    }

    for(size_t i = 0; i < numOfOmega_; i++){
      getEMatrices(
        &(EMatricesVec_[i]),
        &(dielectricMatrixVec[i]),
        numOfLayer,
        N
      );

      getGrandImaginaryMatrices(
        &(grandImaginaryMatricesVec_[i]),
        &(dielectricImMatrixVec[i]),
        numOfLayer,
        N
      );
    }

    thicknessListVec_ = zeros<RCWAVector>(numOfLayer);
    double* thicknessList = structure_->getThicknessList();
    sourceList_.resize(numOfLayer);
    for(size_t i = 0; i < numOfLayer; i++){
      thicknessListVec_(i) = thicknessList[i];
      sourceList_[i] = (structure_->getLayerByIndex(i))->checkIsSource();
    }
  }

  /*==============================================
  This function computes the flux
  ==============================================*/
  void Simulation::run(){
    //TODO
  }

  /*======================================================
  Implementaion of the class on planar simulation
  =======================================================*/
  void SimulationPlanar::setKxIntegral(double end){
    kxStart_ = 0;
    numOfKx_ = 0;
    kxEnd_ = end;
  }

  /*==============================================
  This function gets the flux at a given kx
  @args:
  omegaIndex: the index of omega
  kx: the kx value, normalized
  ==============================================*/
  double SimulationPlanar::getPhiAtKx(int omegaIdx, double kx){
    return POW3(omegaList_[omegaIdx] / datum::c_0) / POW2(datum::pi) * kx *
      poyntingFlux(omegaList_[omegaIdx] / datum::c_0, &thicknessListVec_, kx, 0, &(EMatricesVec_[omegaIdx]),
      &(grandImaginaryMatricesVec_[omegaIdx]), &(dielectricMatrixInverseVec_[omegaIdx]), &Gx_mat_, &Gy_mat_,
      &sourceList_, targetLayer_,1);
  }

  /*==============================================
  This function integrates kx using gauss_legendre method
  ==============================================*/
  void SimulationPlanar::run(){
    if(numOfCore_ == 1){
      ArgWrapper wrapper;
      wrapper.thicknessList = thicknessListVec_;
      wrapper.Gx_mat = Gx_mat_;
      wrapper.Gy_mat = Gy_mat_;
      wrapper.sourceList = sourceList_;
      wrapper.targetLayer = targetLayer_;

      for(size_t omegaIdx = 0; omegaIdx < numOfOmega_; omegaIdx++){
        wrapper.omega = omegaList_[omegaIdx] / datum::c_0;
        wrapper.EMatrices = EMatricesVec_[omegaIdx];
        wrapper.grandImaginaryMatrices = grandImaginaryMatricesVec_[omegaIdx];
        wrapper.dielectricMatrixInverse = dielectricMatrixInverseVec_[omegaIdx];
        Phi_[omegaIdx] = POW3(omegaList_[omegaIdx] / datum::c_0) / POW2(datum::pi) *
          gauss_legendre(DEGREE, wrapperFun, &wrapper, kxStart_, kxEnd_);
      }

    }
    else{
      // TODO: implement MPI
    }
  }

  /*======================================================
  Implementaion of the class on 1D grating simulation
  =======================================================*/
  void SimulationGrating::setKxIntegral(double points){
    numOfKx_ = points;
    kxStart_ = -datum::pi / period_[0];
    kxEnd_ = -kxStart_;
  }

  /*==============================================
  This function set the integral of kx when the system is symmetric in x direction
  @args:
  points: number of kx points
  ==============================================*/
  void SimulationGrating::setKxIntegralSym(double points){
    numOfKx_ = points;
    kxStart_ = 0;
    kxEnd_ = datum::pi / period_[0];
    prefactor_ *= 2;
  }

  /*==============================================
  This function set the integral of ky
  @args:
  points: number of ky points
  end: the upperbound of the integral
  ==============================================*/
  void SimulationGrating::setKyIntegral(double points, double end){
    kyStart_ = 0;
    numOfKy_ = points;
    kyEnd_ = end;
  }

  /*======================================================
  Implementaion of the class on 2D patterning simulation
  =======================================================*/
  void SimulationPattern::setKxIntegral(double points){
    kxStart_ = -datum::pi / period_[0];
    numOfKx_ = points;
    kxEnd_ = -kxStart_;
  }

  /*==============================================
  This function set the integral of kx when the system is symmetric in x direction
  @args:
  points: number of kx points
  ==============================================*/
  void SimulationPattern::setKxIntegralSym(double points){
    kxStart_ = 0;
    numOfKx_ = points;
    kxEnd_ = datum::pi / period_[0];
    prefactor_ *= 2;
  }

  /*==============================================
  This function set the integral of ky
  @args:
  points: number of ky points
  end: the upperbound of the integral
  ==============================================*/
  void SimulationPattern::setKyIntegral(double points){
    kyStart_ = -datum::pi / period_[1];
    numOfKy_ = points;
    kyEnd_ = -kyStart_;
  }

  /*==============================================
  This function set the integral of ky when the system is symmetric in y direction
  @args:
  points: number of ky points
  ==============================================*/
  void SimulationPattern::setKyIntegralSym(double points){
    kyStart_ = 0;
    numOfKy_ = points;
    kyEnd_ = datum::pi / period_[1];
    prefactor_ *= 2;
  }

}