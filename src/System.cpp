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
  void Material::setOmega(const double* omegaList, const int numOfOmega){
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
          epsilonList_.epsilonVals[i].diagonal[i] = epsilonList.epsilonVals[i].diagonal[i];
        }
      }
    }
    else{
      for(int i = 0; i < numOfOmega_; i++){
        for(int j = 0; j < 10; j++){
          epsilonList_.epsilonVals[i].tensor[i] = epsilonList.epsilonVals[i].tensor[i];
        }
      }
    }
  }

  /*==============================================*/
  // Implementaion of the Layer class
  /*==============================================*/
  Layer::Layer(const string name, const Ptr<Material>& material, const double thickness) :
    NamedInterface(name), thickness_(thickness), pattern_(PLANAR_), source_(ISNOTSOURCE_){
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
    backGround_(nullptr), source_(ISNOTSOURCE_), thickness_(0), pattern_(PLANAR_){}

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

    const_MaterialIter itMat = this->getVecBegin();
    const_PatternIter itArg1 = this->getArg1Begin();
    const_PatternIter itArg2 = this->getArg2Begin();
    switch (this->getPattern()) {
      case GRATING_:{
        for(int count = 0; (itMat + count) != this->getVecEnd(); count++){
          const std::pair<double, double> arg1Pair = *(itArg1 + count);
          newLayer->addGratingPattern(*(itMat + count), arg1Pair.first, arg1Pair.second);
        }
        break;
      }
      case RECTANGLE_:{
        for(int count = 0; (itMat + count) != this->getVecEnd(); count++){
          const std::pair<double, double> arg1Pair = *(itArg1 + count);
          const std::pair<double, double> arg2Pair = *(itArg2 + count);
          const double arg1[2] = {arg1Pair.first, arg1Pair.second};
          const double arg2[2] = {arg2Pair.first, arg2Pair.second};
          newLayer->addRectanlgePattern(*(itMat + count), arg1, arg2);
        }
        break;
      }
      case CIRCLE_:{
        for(int count = 0; (itMat + count) != this->getVecEnd(); count++){
          const std::pair<double, double> arg1Pair = *(itArg1 + count);
          const std::pair<double, double> arg2Pair = *(itArg2 + count);
          const double arg[2] = {arg1Pair.first, arg2Pair.first};
          const double radius = arg1Pair.second;
          newLayer->addCirclePattern(*(itMat + count), arg, radius);
        }
        break;
      }
      default:break;
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
    hasTensor_ = val;
  }
  /*==============================================*/
  // check whether the layer contains a material with tensor dielectric
  /*==============================================*/
  bool Layer::hasTensor(){
    return hasTensor_;
  }
  /*==============================================*/
  // get the background material
  /*==============================================*/
  Ptr<Material> Layer::getBackGround(){
    if(backGround_ == nullptr){
      throw UTILITY::NullPointerException("Backgroud not set yet!");
    }
    return backGround_;
  }
  /*==============================================*/
  // get the material by the name
  /*==============================================*/
  Ptr<Material> Layer::getMaterialByName(const std::string name){
    for(const_MaterialIter it = this->getVecBegin(); it != this->getVecEnd(); it++){
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
  // set the pattern of the layer
  /*==============================================*/
  PATTEN Layer::getPattern(){
    return pattern_;
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
  const_MaterialIter Layer::getVecBegin(){
    return materialVec_.cbegin();
  }
  /*==============================================*/
  // material iterator end
  /*==============================================*/
  const_MaterialIter Layer::getVecEnd(){
    return materialVec_.cend();
  }
  /*==============================================*/
  // pattern parameter iterator begin
  /*==============================================*/
  const_PatternIter Layer::getArg1Begin(){
    return args1_.cbegin();
  }
  /*==============================================*/
  // pattern parameter iterator end
  /*==============================================*/
  const_PatternIter Layer::getArg2Begin(){
    return args2_.cbegin();
  }
  /*==============================================*/
  // pattern parameter iterator begin
  /*==============================================*/
  const_PatternIter Layer::getArg1End(){
    return args1_.cend();
  }
  /*==============================================*/
  // pattern parameter iterator end
  /*==============================================*/
  const_PatternIter Layer::getArg2End(){
    return args2_.cend();
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
    pattern_ = RECTANGLE_;
    materialVec_.push_back(material);
    args1_.push_back(std::make_pair(args1[0], args1[1]));
    args2_.push_back(std::make_pair(args2[0], args2[1]));

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
    pattern_ = CIRCLE_;
    materialVec_.push_back(material);
    args1_.push_back(std::make_pair(args[0], radius));
    args2_.push_back(std::make_pair(args[1], radius));
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
    pattern_ = GRATING_;
    materialVec_.push_back(material);
    args1_.push_back(std::make_pair(center, width));
  }

  /*==============================================*/
  // Implementaion of the structure class
  /*==============================================*/
  Structure::Structure(){
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
  // function adding a layer to the structure
  // @args:
  // layer: the added layer
  /*==============================================*/
  void Structure::addLayer(const Ptr<Layer>& layer){
    int size = layerMap_.size();
    layerMap_.insert(LayerMap::value_type(size, layer));
  }
  /*==============================================*/
  // function deleting a layer in the structure by its name
  // @args:
  // name: the name of the layer
  /*==============================================*/
  void Structure::deleteLayerByName(const string name){
    for(const_LayerIter it = this->getMapBegin(); it != this->getMapEnd(); it++){
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
    for(const_LayerIter it = this->getMapBegin(); it != this->getMapEnd(); it++){
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
  void Structure::getThicknessList(double* thicknessList){
    int count = 0;
    for(const_LayerIter it = this->getMapBegin(); it != this->getMapEnd(); it++){
      thicknessList[count] = (it->second)->getThickness();
      count++;
    }
    return;
  }
  /*==============================================*/
  // layer iterator begin
  /*==============================================*/
  const_LayerIter Structure::getMapBegin(){
    return layerMap_.cbegin();
  }
  /*==============================================*/
  // layer iterator end
  /*==============================================*/
  const_LayerIter Structure::getMapEnd(){
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
    for(const_LayerIter it = this->getMapBegin(); it != this->getMapEnd(); it++){
      newMap.insert(LayerMap::value_type(newMap.size(), it->second));
    }
    layerMap_.clear();
    layerMap_ = newMap;
  }
}