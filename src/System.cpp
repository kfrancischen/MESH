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
#include "System.h"

namespace SYSTEM{
  /*======================================================
  Implementaion of the Material class
  =======================================================*/
  Material::Material(std::string name, dcomplex* epsilonList, double* omegaList, int numOfOmega):
  name_(name), numOfOmega_(numOfOmega){
    epsilonList_ = new dcomplex[numOfOmega_];
    std::copy(epsilonList, epsilonList + numOfOmega_, epsilonList_);
    omegaList_ = new double[numOfOmega_];
    std::copy(omegaList, omegaList + numOfOmega_, omegaList_);
  }

  Material::Material(std::string name) :
  name_(name), epsilonList_(nullptr), omegaList_(nullptr), numOfOmega_(0){}


  Material::Material(const Material &material){
    numOfOmega_ = material.numOfOmega_;
    name_ = material.name_;
    epsilonList_ = new dcomplex[numOfOmega_];
    std::copy(material.epsilonList_, material.epsilonList_ + numOfOmega_, epsilonList_);
    omegaList_ = new double[numOfOmega_];
    std::copy(material.omegaList_, material.omegaList_ + numOfOmega_, omegaList_);
  }

  Material::~Material(){
    delete epsilonList_;
    delete omegaList_;
  }

  std::string Material::getName(){
    return name_;
  }

  dcomplex* Material::getEpsilon(){
    return epsilonList_;
  }

  double* Material::getOmegaList(){
    return omegaList_;
  }

  int Material::getNumOfOmega(){
    return numOfOmega_;
  }

  void Material::setName(const std::string name){
    name_ = name;
  }

  void Material::setOmega(const double* omegaList, int numOfOmega){
    numOfOmega_ = numOfOmega;
    omegaList_ = new double[numOfOmega_];
    std::copy(omegaList, omegaList + numOfOmega_, omegaList_);
  }

  void Material::setEpsilon(const dcomplex* epsilonList, int numOfOmega){
    numOfOmega_ = numOfOmega;
    epsilonList_ = new dcomplex[numOfOmega_];
    std::copy(epsilonList, epsilonList + numOfOmega_, epsilonList_);
  }

  /*======================================================
  Implementaion of the Layer class
  =======================================================*/
  Layer::Layer(Material* material, double thickness, SOURCE source) :
    thickness_(thickness), source_(source){
    backGround_ = material;
  }

  Layer::Layer(Material* material) :
    source_(ISNOTSOURCE_), thickness_(0){
    //backGround_ = new Material(*material);
    backGround_ = material;
  }

  Layer::Layer() :
    backGround_(nullptr), source_(ISNOTSOURCE_), thickness_(0){}

  Layer::~Layer(){
    for(MaterialIter it = materialVec_.begin(); it != materialVec_.end(); it++){
      delete(*it);
    }
  }

  Layer::Layer(const Layer& layer){
    backGround_ = new Material(*(layer.backGround_));
    pattern_ = layer.pattern_;
    source_ = layer.source_;
    thickness_ = layer.thickness_;
    for(const_MaterialIter it = layer.materialVec_.cbegin(); it != layer.materialVec_.cend(); it++){
      Material* material = new Material(*(*it));
      materialVec_.push_back(material);
    }

    for(const_PatternIter it = layer.args1_.cbegin(); it != layer.args1_.cend(); it++){
      args1_.push_back(*it);
    }

    for(const_PatternIter it = layer.args2_.cbegin(); it != layer.args2_.cend(); it++){
      args2_.push_back(*it);
    }

  }

  void Layer::setBackGround(Material *material){
    //backGround_ = new Material(*material);
    backGround_ = material;
  }

  void Layer::setThickness(double thickness){
    thickness_ = thickness;
  }

  void Layer::isSource(){
    source_ = ISSOURCE_;
  }

  void Layer::isNotSource(){
    source_ = ISNOTSOURCE_;
  }

  SOURCE Layer::checkSource(){
    return source_;
  }

  Material* Layer::getBackGround(){
    return backGround_;
  }

  Material* Layer::getMaterialByName(std::string name){
    for(const_MaterialIter it = this->getVecBegin(); it != this->getVecEnd(); it++){
      if(name.compare((*it)->getName()) == 0) return *it;
    }
    return nullptr;
  }

  int Layer::getNumOfMaterial(){
    return materialVec_.size();
  }

  double Layer::getThickness(){
    return thickness_;
  }

  const_MaterialIter Layer::getVecBegin(){
    return materialVec_.cbegin();
  }

  const_MaterialIter Layer::getVecEnd(){
    return materialVec_.cend();
  }

  const_PatternIter Layer::getArg1Begin(){
    return args1_.cbegin();
  }

  const_PatternIter Layer::getArg2Begin(){
    return args2_.cbegin();
  }

  const_PatternIter Layer::getArg1End(){
    return args1_.cend();
  }

  const_PatternIter Layer::getArg2End(){
    return args2_.cend();
  }

  void Layer::addPattern(Material * material, double args1[2], double args2[2], std::string pattern){
    pattern_ = pattern;
    materialVec_.push_back(material);
    args1_.push_back(std::make_pair(args1[0], args1[1]));
    args2_.push_back(std::make_pair(args2[0], args2[2]));
  }

  void Layer::addPattern(Material * material, double args1[2]){
    materialVec_.push_back(material);
    args1_.push_back(std::make_pair(args1[0], args1[1]));
  }

  /*======================================================
  Implementaion of the structure class
  =======================================================*/
  Structure::Structure() : nGx_(0), nGy_(0){}

  Structure::~Structure(){
    for(LayerIter it = layerMap_.begin(); it != layerMap_.end(); it++){
      delete(it->second);
    }
  }

  Structure::Structure(Structure& structure){
    nGx_ = structure.nGx_;
    nGy_ = structure.nGy_;
    for(const_LayerIter it = structure.layerMap_.cbegin(); it != structure.layerMap_.cend(); it++){
      layerMap_.insert(LayerMap::value_type(it->first, it->second));
    }
  }


  void Structure::addLayer(Layer* layer){
    int size = layerMap_.size();
    layerMap_.insert(LayerMap::value_type(size, layer));
  }

  Layer* Structure::getLayerByIndex(int index){
    return layerMap_.at(index);
  }

  int Structure::getNumOfLayer(){
    return layerMap_.size();
  }

  double* Structure::getThicknessList(){
    int numOfLayer = this->getNumOfLayer();
    double* thicknessList = new double[numOfLayer];
    int count = 0;
    for(const_LayerIter it = this->getMapBegin(); it != this->getMapEnd(); it++){
      thicknessList[count] = (it->second)->getThickness();
    }
    return thicknessList;
  }

  const_LayerIter Structure::getMapBegin(){
    return layerMap_.cbegin();
  }

  const_LayerIter Structure::getMapEnd(){
    return layerMap_.cend();
  }

  void Structure::setGx(int nGx){
    nGx_ = nGx;
  }

  void Structure::setGy(int nGy){
    nGy_ = nGy;
  }

  int Structure::getGx(){
    return nGx_;
  }

  int Structure::getGy(){
    return nGy_;
  }

}