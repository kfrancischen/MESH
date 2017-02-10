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
  Material* Layer::getBackGround(){
    return backGround_;
  }

  int Layer::getNumOfMaterial(){
    return materialVec_.size();
  }


  const_MaterialIter Layer::getVecBegin(){
    return materialVec_.cbegin();
  }

  const_MaterialIter Layer::getVecEnd(){
    return materialVec_.cend();
  }

  void Layer::addPattern(Material * material, double args1[2], double args2[2], std::string pattern){
    pattern_ = pattern;
    materialVec_.push_back(material);
    args1_[0] = args1[0];
    args1_[1] = args1[1];
    args2_[0] = args2[0];
    args2_[1] = args2[1];
  }


  Material::Material(std::string name, dcomplex* epsilonList, double* omegaList, size_t size):
  name_(name), size_(size){
    epsilonList_ = new dcomplex[size_];
    std::copy(epsilonList, epsilonList + size_, epsilonList_);
    omegaList_ = new double[size_];
    std::copy(omegaList, omegaList + size_, omegaList_);
  }

  Material::Material(std::string name) :
  name_(name), epsilonList_(nullptr), omegaList_(nullptr), size_(0){}


  Material::Material(const Material &material){
    size_ = material.size_;
    name_ = material.name_;
    epsilonList_ = new dcomplex[size_];
    std::copy(material.epsilonList_, material.epsilonList_ + size_, epsilonList_);
    omegaList_ = new double[size_];
    std::copy(material.omegaList_, material.omegaList_ + size_, omegaList_);
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

  double* Material::getOmega(){
    return omegaList_;
  }

  void Material::setName(const std::string name){
    name_ = name;
  }

  void Material::setOmega(const double* omega, size_t size){
    size_ = size;
    omegaList_ = new double[size_];
    std::copy(omega, omega + size_, omegaList_);
  }

  void Material::setEpsilon(const dcomplex* epsilon, size_t size){
    size_ = size;
    epsilonList_ = new dcomplex[size_];
    std::copy(epsilon, epsilon + size_, epsilonList_);
  }

  /*======================================================
  Implementaion of the Layer class
  =======================================================*/
  Layer::Layer(Material* material){
    backGround_ = new Material(*material);
  }

  Layer::Layer() : backGround_(nullptr){}

  Layer::~Layer(){
    for(MaterialIter it = materialVec_.begin(); it != materialVec_.end(); it++){
      delete(*it);
    }
  }

  Layer::Layer(const Layer& layer){
    backGround_ = new Material(*(layer.backGround_));
    pattern_ = layer.pattern_;
    args1_[0] = layer.args1_[0];
    args1_[1] = layer.args1_[1];
    args2_[0] = layer.args2_[0];
    args2_[1] = layer.args2_[1];
    for(const_MaterialIter it = layer.materialVec_.cbegin(); it != layer.materialVec_.cend(); it++){
      Material* material = new Material(*(*it));
      materialVec_.push_back(material);
    }
  }

  void Layer::setBackGround(Material *material){
    backGround_ = new Material(*material);
  }
  /*======================================================
  Implementaion of the structure class
  =======================================================*/
  Structure::Structure(){}

  Structure::~Structure(){
    for(LayerIter it = layerMap_.begin(); it != layerMap_.end(); it++){
      delete(it->second);
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


}