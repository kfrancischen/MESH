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
#include <iostream>
namespace SYSTEM{
  /*==============================================*/
  // Implementaion of the Material class
  /*==============================================*/
  Material::Material(
    const std::string name,
    const double* omegaList,
    const EPSILON& epsilonList,
    const int numOfOmega): NamedInterface(name)
    , numOfOmega_(numOfOmega){
    epsilonList_.epsilonVals = new EpsilonVal[numOfOmega_];
    epsilonList_.type_ = epsilonList.type_;
    if(epsilonList.type_ == SCALAR_){
      for(int i = 0; i < numOfOmega_; i++){
        epsilonList_.epsilonVals[i].scalar[0] = epsilonList.epsilonVals[i].scalar[0];
        epsilonList_.epsilonVals[i].scalar[1] = epsilonList.epsilonVals[i].scalar[1];
      }
      //std::copy(epsilonList.epsilonVals, epsilonList.epsilonVals + numOfOmega_, epsilonList_);
    }
    else if(epsilonList.type_ == DIAGONAL_){
      for(int i = 0; i < numOfOmega_; i++){
        for(int j = 0; j < 6; j++){
          epsilonList_.epsilonVals[i].diagonal[j] = epsilonList.epsilonVals[i].diagonal[j];
        }
      }
    }
    else{
      for(int i = 0; i < numOfOmega_; i++){
        for(int j = 0; j < 10; j++){
          epsilonList_.epsilonVals[i].tensor[j] = epsilonList.epsilonVals[i].tensor[j];
        }
      }
    }
    omegaList_ = new double[numOfOmega_];
    std::copy(omegaList, omegaList + numOfOmega_, omegaList_);
  }

  /*==============================================*/
  // This is a thin wrapper for the usage of smart pointer
  /*==============================================*/
  Ptr<Material> Material::instanceNew(
    const std::string name,
    const double* omegaList,
    const EPSILON& epsilonList,
    const int numOfOmega
  ){
    return new Material(name, omegaList, epsilonList, numOfOmega);
  }

  /*==============================================*/
  // destructor
  /*==============================================*/
  Material::~Material(){
    delete[] epsilonList_.epsilonVals;
    epsilonList_.epsilonVals = nullptr;
    delete[] omegaList_;
    omegaList_ = nullptr;
  }
  /*==============================================*/
  // function return the name of the material
  /*==============================================*/
  std::string Material::getName(){
    return this->name();
  }
  /*==============================================*/
  // function check whether the given material has a tensor dielectric
  /*==============================================*/
  EPSTYPE Material::getType(){
    return epsilonList_.type_;
  }

  /*==============================================*/
  // function return the epsilon at a specific index
  /*==============================================*/
  EpsilonVal Material::getEpsilonAtIndex(const int index){
    if(index >= numOfOmega_){
      std::cerr << std::to_string(index) + ": out of range!" << std::endl;
      throw UTILITY::RangeException(std::to_string(index) + ": out of range!");
    }
    return epsilonList_.epsilonVals[index];
  }

  /*==============================================*/
  // function return the omega list
  /*==============================================*/
  double* Material::getOmegaList(){
    return omegaList_;
  }
  /*==============================================*/
  // function return the number of omega points
  /*==============================================*/
  int Material::getNumOfOmega(){
    return numOfOmega_;
  }
  /*==============================================*/
  // function set the omega of the material
  // @args:
  // omegaList: the omegaList of the material
  // numOfOmega: the number of omega points
  /*==============================================*/
  void Material::setOmega(const double* &omegaList, const int numOfOmega){
    numOfOmega_ = numOfOmega;
    omegaList_ = new double[numOfOmega_];
    std::copy(omegaList, omegaList + numOfOmega_, omegaList_);
  }
  /*==============================================*/
  // function set the epsilon values of the material
  // @args:
  // epsilonList: the epsilonList of the material
  // numOfOmega: the number of omega points
  /*==============================================*/
  void Material::setEpsilon(const EPSILON& epsilonList, const int numOfOmega){
    numOfOmega_ = numOfOmega;
    epsilonList_.epsilonVals = new EpsilonVal[numOfOmega_];
    if(epsilonList.type_ == SCALAR_){
      for(int i = 0; i < numOfOmega_; i++){
        epsilonList_.epsilonVals[i].scalar[0] = epsilonList.epsilonVals[i].scalar[0];
        epsilonList_.epsilonVals[i].scalar[1] = epsilonList.epsilonVals[i].scalar[1];
      }
      //std::copy(epsilonList.epsilonVals, epsilonList.epsilonVals + numOfOmega_, epsilonList_);
    }
    else if(epsilonList.type_ == DIAGONAL_){
      for(int i = 0; i < numOfOmega_; i++){
        for(int j = 0; j < 6; j++){
          epsilonList_.epsilonVals[i].diagonal[j] = epsilonList.epsilonVals[i].diagonal[j];
        }
      }
    }
    else{
      for(int i = 0; i < numOfOmega_; i++){
        for(int j = 0; j < 10; j++){
          epsilonList_.epsilonVals[i].tensor[j] = epsilonList.epsilonVals[i].tensor[j];
        }
      }
    }
  }

  /*==============================================*/
  // Implementaion of the Layer class
  /*==============================================*/
  Layer::Layer(const string name, const Ptr<Material>& material, const double thickness) :
    NamedInterface(name), thickness_(thickness), source_(ISNOTSOURCE_){
    backGround_ = material;
  }

  /*==============================================*/
  // This is a thin wrapper for the usage of smart pointer
  /*==============================================*/
  Ptr<Layer> Layer::instanceNew(
    const string name,
    const Ptr<Material>& material,
    const double thickness
  ){
    return new Layer(name, material, thickness);
  }
  /*==============================================*/
  // Plain layer constructor
  /*==============================================*/
  Layer::Layer(const string name) : NamedInterface(name),
    thickness_(0), backGround_(nullptr), source_(ISNOTSOURCE_){}

  /*==============================================*/
  // This is a thin wrapper for the usage of smart pointer
  /*==============================================*/
  Ptr<Layer> Layer::instanceNew(const string name)
  {
    return new Layer(name);
  }

  /*==============================================*/
  // destructor
  /*==============================================*/
  Layer::~Layer(){
  }

  /*==============================================*/
  // function return a copy of the layer
  // @arg:
  // name: the name of the new layer
  /*==============================================*/

  Ptr<Layer> Layer::layerCopy(const string name){
    Ptr<Layer> newLayer = Layer::instanceNew(name);
    newLayer->setBackGround(this->getBackGround());
    newLayer->setThickness(this->getThickness());

    const_MaterialIter itMat = this->getMaterialsBegin();
    const_PatternIter itPattern = this->getPatternsBegin();
    for(int count = 0; (itMat + count) != this->getMaterialsEnd(); count++){
      Pattern pattern = *(itPattern + count);
      switch(pattern.type_){
        case GRATING_:{
          newLayer->addGratingPattern(*(itMat + count), pattern.arg1_.first, pattern.arg1_.second);
          break;
        }
        case RECTANGLE_:{
          const double arg1[2] = {pattern.arg1_.first, pattern.arg1_.second};
          const double arg2[2] = {pattern.arg2_.first, pattern.arg2_.second};
          newLayer->addRectanlgePattern(*(itMat + count), arg1, arg2);
          break;
        }
        case CIRCLE_:{
          const double arg[2] = {pattern.arg1_.first, pattern.arg2_.first};
          const double radius = pattern.arg1_.second;
          newLayer->addCirclePattern(*(itMat + count), arg, radius);
        }
        default: break;
      }
    }
    return newLayer;
  }
  /*==============================================*/
  // set background by material
  // @args:
  // material: the background material
  /*==============================================*/
  void Layer::setBackGround(const Ptr<Material>& material){
    //backGround_ = new Material(*material);
    backGround_ = material;
    if(material->getType() == TENSOR_) hasTensor_++;
  }
  /*==============================================*/
  // set thickness of the layer
  // @args:
  // thickness: the thickness of the layer
  /*==============================================*/
  void Layer::setThickness(const double thickness){
    thickness_ = thickness;
  }
  /*==============================================*/
  // set the layer to be the source
  /*==============================================*/
  void Layer::setIsSource(){
    source_ = ISSOURCE_;
  }
  /*==============================================*/
  // check whether the layer is a source
  /*==============================================*/
  bool Layer::checkIsSource(){
    if(source_ == ISSOURCE_){
      return true;
    }
    return false;
  }
  /*==============================================*/
  // set the layer contains a material with tensor dielectric
  /*==============================================*/
  void Layer::containTensor(bool val){
    if(val) hasTensor_++;
    else hasTensor_--;
  }
  /*==============================================*/
  // check whether the layer contains a material with tensor dielectric
  /*==============================================*/
  bool Layer::hasTensor(){
    return hasTensor_ != 0;
  }
  /*==============================================*/
  // check whether the layer contains a given material
  /*==============================================*/
  bool Layer::hasMaterial(const Ptr<Material>& material){
    if((backGround_->name()).compare(material->getName()) == 0) return true;
    for(const_MaterialIter it = this->getMaterialsBegin(); it != this->getMaterialsEnd(); it++){
      if(*it == material) return true;
    }
    return false;
  }
  /*==============================================*/
  // get the background material
  /*==============================================*/
  Ptr<Material> Layer::getBackGround(){
    if(backGround_ == nullptr){
      std::cerr << "Backgroud not set yet!" << std::endl;
      throw UTILITY::NullPointerException("Backgroud not set yet!");
    }
    return backGround_;
  }
  /*==============================================*/
  // get the material by the name
  /*==============================================*/
  Ptr<Material> Layer::getMaterialByName(const std::string name){
    for(const_MaterialIter it = this->getMaterialsBegin(); it != this->getMaterialsEnd(); it++){
      if(name.compare((*it)->getName()) == 0) return *it;
    }
    return NULL;
  }
  /*==============================================*/
  // get the number of materials other than background
  /*==============================================*/
  int Layer::getNumOfMaterial(){
    return materialVec_.size();
  }
  /*==============================================*/
  // get the thickness of the layer
  /*==============================================*/
  double Layer::getThickness(){
    return thickness_;
  }
  /*==============================================*/
  // function return the name of the material
  /*==============================================*/
  std::string Layer::getName(){
    return this->name();
  }
  /*==============================================*/
  // material iterator begin
  /*==============================================*/
  const_MaterialIter Layer::getMaterialsBegin(){
    return materialVec_.cbegin();
  }
  /*==============================================*/
  // material iterator end
  /*==============================================*/
  const_MaterialIter Layer::getMaterialsEnd(){
    return materialVec_.cend();
  }
  /*==============================================*/
  // pattern parameter iterator begin
  /*==============================================*/
  const_PatternIter Layer::getPatternsBegin(){
    return patternVec_.cbegin();
  }
  /*==============================================*/
  // pattern parameter iterator end
  /*==============================================*/
  const_PatternIter Layer::getPatternsEnd(){
    return patternVec_.cend();
  }
  /*==============================================*/
  // add a rectangle pattern
  // @args:
  // material: the material used for this part of the pattern
  // args1: the position of centers (x,y)
  // args2: the widths in x and y directions
  /*==============================================*/
  void Layer::addRectanlgePattern(
    const Ptr<Material>& material,
    const double args1[2],
    const double args2[2]
  ){

    materialVec_.push_back(material);
    if(material->getType() == TENSOR_) hasTensor_++;
    Pattern pattern;
    pattern.arg1_ = std::make_pair(args1[0], args1[1]);
    pattern.arg2_ = std::make_pair(args2[0], args2[1]);
    pattern.type_ = RECTANGLE_;
    patternVec_.push_back(pattern);
  }

  /*==============================================*/
  // add a circular pattern
  // @args:
  // material: the material used for this part of the pattern
  // args1: the position of centers (x,y)
  // radius: the radius of the circle
  /*==============================================*/
  void Layer::addCirclePattern(
    const Ptr<Material>& material,
    const double args[2],
    const double radius
  ){
    materialVec_.push_back(material);
    if(material->getType() == TENSOR_) hasTensor_++;
    Pattern pattern;
    pattern.arg1_ = std::make_pair(args[0], radius);
    pattern.arg2_ = std::make_pair(args[1], radius);
    pattern.type_ = CIRCLE_;
    patternVec_.push_back(pattern);
  }
  /*==============================================*/
  // add a grating pattern
  // @args:
  // material: the material used for this part of the pattern
  // start: the start position of the this pattern
  // end: the end position of this pattern
  /*==============================================*/
  void Layer::addGratingPattern(
    const Ptr<Material>& material,
    const double center,
    const double width
  ){
    materialVec_.push_back(material);
    if(material->getType() == TENSOR_) hasTensor_++;
    Pattern pattern;
    pattern.arg1_ = std::make_pair(center, width);
    pattern.arg2_ = std::make_pair(0, 0);
    pattern.type_ = GRATING_;
    patternVec_.push_back(pattern);
  }

  /*==============================================*/
  // Implementaion of the structure class
  /*==============================================*/
  Structure::Structure(){
    period_[0] = 0;
    period_[1] = 0;
  }
  /*==============================================*/
  // This is a thin wrapper for the usage of smart pointer
  /*==============================================*/
  Ptr<Structure> Structure::instanceNew(){
    return new Structure();
  }
  /*==============================================*/
  // destructor
  /*==============================================*/
  Structure::~Structure(){
    //delete[] period_;
    //period_ = nullptr;
  }
  /*==============================================*/
  // copy constructor
  /*==============================================*/
  Structure::Structure(const Structure& structure){
    for(const_LayerIter it = structure.layerMap_.cbegin(); it != structure.layerMap_.cend(); it++){
      layerMap_.insert(LayerMap::value_type(it->first, it->second));
    }
  }
  /*==============================================*/
  // function adding a material to the structure
  // @args:
  // material: the added material
  /*==============================================*/
  void Structure::addMaterial(const Ptr<Material>& material){
    materialMap_.insert(MaterialMap::value_type(material->getName(), material));
  }
  /*==============================================*/
  // function adding a layer to the structure
  // @args:
  // layer: the added layer
  /*==============================================*/
  void Structure::addLayer(const Ptr<Layer>& layer){
    int size = layerMap_.size();
    layerMap_.insert(LayerMap::value_type(size, layer));
  }
  /*==============================================*/
  // function setting the periodicity of a layer
  // @args:
  // p1: periodicity in x direction
  // p2: periodicity in y direction
  /*==============================================*/
  void Structure::setPeriodicity(const double p1, const double p2){
    period_[0] = p1;
    period_[1] = p2;
  }
  /*==============================================*/
  // function deleting a layer in the structure by its name
  // @args:
  // name: the name of the layer
  /*==============================================*/
  void Structure::deleteLayerByName(const string name){
    for(const_LayerIter it = this->getLayersBegin(); it != this->getLayersEnd(); it++){
      if(!name.compare((it->second)->getName())){
        this->deleteLayer(it);
        break;
      }
    }
    this->reorganizeLayers();
  }
  /*==============================================*/
  // function deleting a layer in the structure by its pointer
  // @args:
  // layer: the pointer of the layer
  /*==============================================*/
  void Structure::deleteLayerByLayer(const Ptr<Layer>& layer){
    this->deleteLayerByName(layer->getName());
  }
  /*==============================================*/
  // function getting a layer by its index
  // @args:
  // index: the index of the wanted layer
  /*==============================================*/
  Ptr<Layer> Structure::getLayerByIndex(const int index){
    if(index >= this->getNumOfLayer()){
      return NULL;
    }
    return layerMap_.at(index);
  }

  /*==============================================*/
  // function getting a layer by its name
  // @args:
  // name: the name of the wanted layer
  /*==============================================*/
  Ptr<Layer> Structure::getLayerByName(const std::string name){
    for(const_LayerIter it = this->getLayersBegin(); it != this->getLayersEnd(); it++){
      if(!name.compare((it->second)->getName())){
        return it->second;
      }
    }
    return NULL;
  }
  /*==============================================*/
  // function getting the number of layers
  /*==============================================*/
  int Structure::getNumOfLayer(){
    return layerMap_.size();
  }
  /*==============================================*/
  // function getting the thickness list
  /*==============================================*/
  void Structure::getThicknessList(double* &thicknessList){
    int count = 0;
    for(const_LayerIter it = this->getLayersBegin(); it != this->getLayersEnd(); it++){
      thicknessList[count] = (it->second)->getThickness();
      count++;
    }
    return;
  }

  /*==============================================*/
  // layer iterator begin
  /*==============================================*/
  const_LayerIter Structure::getLayersBegin(){
    return layerMap_.cbegin();
  }
  /*==============================================*/
  // layer iterator end
  /*==============================================*/
  const_LayerIter Structure::getLayersEnd(){
    return layerMap_.cend();
  }

  /*==============================================*/
  // erase one layer of the system
  // @args:
  // it: the iterator of the Layermap
  /*==============================================*/
  void Structure::deleteLayer(const_LayerIter it){
    layerMap_.erase(it);
  }
  /*==============================================*/
  // function rehashing the layermap
  /*==============================================*/
  void Structure::reorganizeLayers(){
    LayerMap newMap;
    for(const_LayerIter it = this->getLayersBegin(); it != this->getLayersEnd(); it++){
      newMap.insert(LayerMap::value_type(newMap.size(), it->second));
    }
    layerMap_.clear();
    layerMap_ = newMap;
  }
}