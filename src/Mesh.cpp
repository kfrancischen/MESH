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
  void filerLoader(std::string fileName, double* omega, dcomplex* epsilon, int size){
    // TODO
  }

  void saveData(std::string fileName, double* omega, double* Phi, int size){
    // TODO
  }

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
    period_[0] = 0;
    period_[1] = 0;
  }

  Simulation::~Simulation(){
    delete Phi_;
    delete structure_;
    delete omegaList_;
  }

  void Simulation::addStructure(Structure* structure){
    structure_ = structure;
  }


  void Simulation::enableMPI(int numOfCore){
    numOfCore_ = numOfCore;
  }

  void Simulation::setPeriod(double p1, double p2){
    period_[0] = p1;
    period_[1] = p2;
  }

  void Simulation::setTargetLayerByIndex(int index){
    targetLayer_ = index;
  }

  void Simulation::setTargetLayerByLayer(Layer* layer){
    for(size_t i = 0; i < structure_->getNumOfLayer(); i++){
      if(structure_->getLayerByIndex(i) == layer){
        targetLayer_ = i;
        return;
      }
    }
  }

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

  Structure* Simulation::getStructure(){
    return structure_;
  }

  double* Simulation::getOmegaList(){
    return omegaList_;
  }

  double* Simulation::getPhi(){
    return Phi_;
  }

  double* Simulation::getPeriod(){
    return period_;
  }


  double Simulation::getPhiAtKxKy(int omegaIndex, double kx, double ky){
    // TODO
  }

  void Simulation::build(){
    Layer* firstLayer = structure_->getLayerByIndex(0);
    Material* backGround = firstLayer->getBackGround();
    numOfOmega_ = backGround->getNumOfOmega();
    omegaList_ = backGround->getOmegaList();
    Phi_ = new double[numOfOmega_];
    EMatricesVec_.reserve(numOfOmega_);
    grandImaginaryMatricesVec_.reserve(numOfOmega_);
    dielectricMatrixInverseVec_.reserve(numOfOmega_);


    RCWAMatricesVec dielectricMatrixVec(numOfOmega_), dielectricImMatrixVec(numOfOmega_);
    int numOfLayer = structure_->getNumOfLayer();
    int N = getN(nGx_, nGy_);

    RCWAMatrix onePadding1N = eye<RCWAMatrix>(N, N);

    for(size_t i = 0; i < numOfLayer; i++){
      Layer* layer = structure_->getLayerByIndex(i);
      dcomplex* epsilon = (layer->getBackGround())->getEpsilon();

      switch (layer->getPattern()) {

        case PLANAR_:{
          for(size_t j = 0; j < numOfOmega_; j++){
            dielectricMatrixVec[j].push_back(onePadding1N * epsilon[j]);
            dielectricMatrixInverseVec_[j].push_back(onePadding1N * dcomplex(1, 0) / epsilon[j]);
            dielectricImMatrixVec[j].push_back(onePadding1N * (epsilon[j]).imag());
          }
          getGMatrices(nGx_, nGy_, period_, &Gx_mat_, &Gy_mat_, NO_);
          break;
        }

        case CIRCLE_:{
          // TODO
          break;
        }

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
        numOfLayer
      );

    }

    thicknessListVec_ = zeros<RCWAVector>(numOfLayer);
    sourceList_.reserve(numOfLayer);
    double* thicknessList = structure_->getThicknessList();
    for(size_t i = 0; i < numOfLayer; i++){
      thicknessListVec_(i) = thicknessList[i];
      sourceList_[i] = (structure_->getLayerByIndex(i))->checkIsSource();
    }
  }

  void Simulation::run(){
    //TODO
  }

  /*======================================================
  Implementaion of the class on planar simulation
  =======================================================*/
  void SimulationPlanar::setKxIntegral(double points, double end){
    kxStart_ = 0;
    numOfKx_ = points;
    kxEnd_ = end;
  }

  double SimulationPlanar::getPhiAtKx(int omegaIndex, double kx){
    // TODO
  }

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

        Phi_[omegaIdx] = POW2(omegaList_[omegaIdx] / datum::c_0) / POW2(datum::pi) *
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

  void SimulationGrating::setKxIntegralSym(double points){
    numOfKx_ = points;
    kxStart_ = 0;
    kxEnd_ = datum::pi / period_[0];
    prefactor_ *= 2;
  }

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

  void SimulationPattern::setKxIntegralSym(double points){
    kxStart_ = 0;
    numOfKx_ = points;
    kxEnd_ = datum::pi / period_[0];
    prefactor_ *= 2;
  }

  void SimulationPattern::setKyIntegral(double points){
    kyStart_ = -datum::pi / period_[1];
    numOfKy_ = points;
    kyEnd_ = -kyStart_;
  }

  void SimulationPattern::setKyIntegralSym(double points){
    kyStart_ = 0;
    numOfKy_ = points;
    kyEnd_ = datum::pi / period_[1];
    prefactor_ *= 2;
  }

}