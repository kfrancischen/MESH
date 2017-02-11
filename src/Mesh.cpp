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

#define POW2(x) pow(x, 2)
#define _2PI datum::pi * 2
namespace MESH{
  void filerLoader(std::string fileName, double* omega, dcomplex* epsilon, int size){
    // TODO
  }

  void saveData(std::string fileName, double* omega, double* fluxSpectrum, int size){
    // TODO
  }

  double wrapperFun(double kx, void* data){
    ArgWrapper* wrapper = (ArgWrapper*) data;
    return kx * poyntingFlux(
      wrapper->omega,
      &wrapper->thicknessList,
      kx,
      0,
      &wrapper->EMatrices,
      &wrapper->grandImaginaryMatrices,
      &wrapper->dielectricMatrixInverse,
      &wrapper->Gx_mat,
      &wrapper->Gy_mat,
      &wrapper->sourceList,
      wrapper->targetLayer,
      0
    );
  }
  /*======================================================
  Implementaion of the parent simulation super class
  =======================================================*/
  Simulation::Simulation() : nGx_(0), nGy_(0), numOfCore_(1), numOfOmega_(0), structure_(nullptr),
  fluxSpectrum_(nullptr), omegaList_(nullptr), kxStart_(0), kxEnd_(0), kyStart_(0), kyEnd_(0), numOfKx_(0), numOfKy_(0)
  {
    period_[0] = 0;
    period_[1] = 0;
  }

  Simulation::~Simulation(){
    delete fluxSpectrum_;
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

  double* Simulation::getFluxSpectrum(){
    return fluxSpectrum_;
  }

  double* Simulation::getPeriod(){
    return period_;
  }

  /*======================================================
  Implementaion of the class on planar simulation
  =======================================================*/
  void SimulationPlanar::setKxIntegral(double points, double end){
    kxStart_ = 0;
    numOfKx_ = points;
    kxEnd_ = end;
  }

  void SimulationPlanar::build(){
    Layer* firstLayer = structure_->getLayerByIndex(0);
    Material* backGround = firstLayer->getBackGround();
    numOfOmega_ = backGround->getNumOfOmega();
    omegaList_ = backGround->getOmegaList();
    fluxSpectrum_ = new double[numOfOmega_];
    EMatricesVec_.reserve(numOfOmega_);
    grandImaginaryMatricesVec_.reserve(numOfOmega_);
    dielectricMatrixInverseVec_.reserve(numOfOmega_);

    RCWAMatricesVec dielectricMatrixVec(numOfOmega_), dielectricImMatrixVec(numOfOmega_);

    int numOfLayer = structure_->getNumOfLayer();
    for(size_t i = 0; i < numOfLayer; i++){
      Layer* layer = structure_->getLayerByIndex(i);
      dcomplex* epsilon = (layer->getBackGround())->getEpsilon();
      for(size_t j = 0; j < numOfOmega_; j++){
        double real = (epsilon[j]).real();
        double imag = (epsilon[j]).imag();
        dielectricMatrixVec[j].push_back(RCWAMatrix(real, imag));
        dielectricMatrixInverseVec_[j].push_back(RCWAMatrix(1.0 / real, 1.0 / imag));
        dielectricImMatrixVec[j].push_back(RCWAMatrix(0, imag));
      }
    }
    getGMatrices(nGx_, nGy_, period_, &Gx_mat_, &Gy_mat_, NO_);
    int N = getN(nGx_, nGy_);
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
  }

  void SimulationPlanar::run(){
    int numOfLayer = structure_->getNumOfLayer();

    RCWAVector thicknessListVec = zeros<RCWAVector>(numOfLayer);
    double* thicknessList = structure_->getThicknessList();
    SourceList sourceList(numOfLayer);

    for(size_t i = 0; i < numOfLayer; i++){
      thicknessListVec(i) = thicknessList[i];
      sourceList[i] = (structure_->getLayerByIndex(i))->checkSource();
    }

    if(numOfCore_ == 1){
      ArgWrapper* wrapper;
      wrapper->thicknessList = thicknessListVec;
      wrapper->Gx_mat = Gx_mat_;
      wrapper->Gy_mat = Gy_mat_;
      wrapper->sourceList = sourceList;
      wrapper->targetLayer = targetLayer_;

      for(size_t omegaIdx = 0; omegaIdx < numOfOmega_; omegaIdx++){
        wrapper->omega = omegaList_[omegaIdx] / datum::c_0;
        wrapper->EMatrices = EMatricesVec_[omegaIdx];
        wrapper->grandImaginaryMatrices = grandImaginaryMatricesVec_[omegaIdx];
        wrapper->dielectricMatrixInverse = dielectricMatrixInverseVec_[omegaIdx];
        fluxSpectrum_[omegaIdx] = POW2(omegaList_[omegaIdx] / datum::c_0) / POW2(datum::pi) *
          gauss_legendre(DEGREE, wrapperFun, wrapper, kxStart_, kxEnd_);
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
    kxStart_ = -_2PI / period_[0];
    kxEnd_ = -kxStart_;
  }

  void SimulationGrating::setKxIntegralSym(double points){
    numOfKx_ = points;
    kxStart_ = 0;
    kxEnd_ = _2PI / period_[0];
    prefactor_ *= 2;
  }

  void SimulationGrating::setKyIntegral(double points, double end){
    kyStart_ = 0;
    numOfKy_ = points;
    kyEnd_ = end;
  }

  void SimulationGrating::build(){
    // TODO
  }

  void SimulationGrating::run(){
    // TODO
  }

  /*======================================================
  Implementaion of the class on 2D patterning simulation
  =======================================================*/
  void SimulationPattern::setKxIntegral(double points){
    kxStart_ = -_2PI / period_[0];
    numOfKx_ = points;
    kxEnd_ = -kxStart_;
  }

  void SimulationPattern::setKxIntegralSym(double points){
    kxStart_ = 0;
    numOfKx_ = points;
    kxEnd_ = _2PI / period_[0];
    prefactor_ *= 2;
  }

  void SimulationPattern::setKyIntegral(double points){
    kyStart_ = -_2PI / period_[1];
    numOfKy_ = points;
    kyEnd_ = -kyStart_;
  }

  void SimulationPattern::setKyIntegralSym(double points){
    kyStart_ = 0;
    numOfKy_ = points;
    kyEnd_ = _2PI / period_[1];
    prefactor_ *= 2;
  }

  void SimulationPattern::build(){
    // TODO
  }

  void SimulationPattern::run(){
    // TODO
  }
}