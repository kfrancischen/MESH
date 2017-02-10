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
#include "gauss_legendre.h"

namespace MESH{
  void filerLoader(std::string fileName, double* omega, dcomplex* epsilon, int size){

  }

  void saveData(std::string fileName, double* omega, double* fluxSpectrum, int size){

  }
  /*======================================================
  Implementaion of the parent simulation super class
  =======================================================*/
  Simulation::Simulation() : nGx_(0), nGy_(0), numOfCore_(1), numOfPoints_(0), structure_(nullptr),
  fluxSpectrum_(nullptr), omegaList_(nullptr)
  {
  }

  Simulation::~Simulation(){
    delete fluxSpectrum_;
    delete structure_;
    delete omegaList_;
  }

  void Simulation::addStructure(Structure* structure){
    structure_ = structure;
  }

  void Simulation::setNumOfPoints(int numOfPoints){
    numOfPoints_ = numOfPoints;
    fluxSpectrum_ = new double[numOfPoints_];
    omegaList_ = new double[numOfPoints_];
  }

  void Simulation::setOmega(double* omegaList){
    for(size_t i = 0; i < numOfPoints_; i++){
      omegaList_[i] = omegaList[i];
    }
  }

  void Simulation::enableMPI(int numOfCore){
    numOfCore_ = numOfCore;
  }

  void Simulation::setPeriod(double p1, double p2){
    period_[0] = p1;
    period_[1] = p2;
  }

  void Simulation::resetSimulation(){
    for(size_t i = 0; i < EMatricesVec.size(); i++){
      EMatricesVec[i].clear();
      dielectricImMatrixVec[i].clear();
      dielectricMatrixVec[i].clear();
      dielectricMatrixInverseVec[i].clear();
    }
    EMatricesVec.clear();
    dielectricImMatrixVec.clear();
    dielectricMatrixVec.clear();
    dielectricMatrixInverseVec.clear();
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
  void SimulationPlanar::initMatrices(){
    // TODO
  }

  void SimulationPlanar::setKxIntegral(double start, double points, double end){
    // nothing to do
  }

  void SimulationPlanar::setKxIntegral(double start, double points){
    // nothing to do
  }

  void SimulationPlanar::setKyIntegral(double start, double points){
    // nothing to do
  }

  void SimulationPlanar::run(){
    // TODO
  }
  /*======================================================
  Implementaion of the class on 1D grating simulation
  =======================================================*/
  void SimulationGrating::initMatrices(){
    // TODO
  }

  void SimulationGrating::setKxIntegral(double start, double points, double end){
    // TODO
  }

  void SimulationGrating::setKxIntegral(double start, double points){
    // nothing to do
  }

  void SimulationGrating::setKyIntegral(double start, double points){
    // TODO
  }

  void SimulationGrating::run(){
    // TODO
  }

  /*======================================================
  Implementaion of the class on 2D patterning simulation
  =======================================================*/
  void SimulationPattern::initMatrices(){
    // TODO
  }

  void SimulationPattern::setKxIntegral(double start, double points, double end){
    // nothing to do
  }

  void SimulationPattern::setKxIntegral(double start, double points){
    // TODO
  }

  void SimulationPattern::setKyIntegral(double start, double points){
    // TODO
  }

  void SimulationPattern::run(){
    // TODO
  }
}