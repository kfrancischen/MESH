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
  Material::Material(
    const std::string name,
    const double* omegaList,
    const dcomplex* epsilonList,
    const int numOfOmega):
  name_(name), numOfOmega_(numOfOmega){
    epsilonList_ = new dcomplex[numOfOmega_];
    std::copy(epsilonList, epsilonList + numOfOmega_, epsilonList_);
    omegaList_ = new double[numOfOmega_];
    std::copy(omegaList, omegaList + numOfOmega_, omegaList_);
  }
  /*======================================================
  This is a thin wrapper for the usage of smart pointer
  =======================================================*/
  Ptr<Material> Material::instanceNew(
    const std::string name,
    const double* omegaList,
    const dcomplex* epsilonList,
    const int numOfOmega
  ){
    return new Material(name, omegaList, epsilonList, numOfOmega);
  }

  /*======================================================
  Material constructor by name
  =======================================================*/
  Material::Material(const std::string name) :
  name_(name), epsilonList_(nullptr), omegaList_(nullptr), numOfOmega_(0){}

  /*======================================================
  This is a thin wrapper for the usage of smart pointer
  =======================================================*/
  Ptr<Material> Material::instanceNew(
    const std::string name
  ){
    return new Material(name);
  }
  /*======================================================
  destructor
  =======================================================*/
  Material::~Material(){
    delete[] epsilonList_;
    delete[] omegaList_;
  }
  /*======================================================
  function return the name of the material
  =======================================================*/
  std::string Material::getName(){
    return name_;
  }
  /*======================================================
  function return the epsilon list of the material
  =======================================================*/
  dcomplex* Material::getEpsilon(){
    return epsilonList_;
  }

  /*======================================================
  function return the epsilon at a specific index
  =======================================================*/
  dcomplex Material::getEpsilonAtIndex(const int index){
    if(index >= numOfOmega_){
      throw UTILITY::RangeException(std::to_string(index) + ": out of range!");
    }
    return epsilonList_[index];
  }

  /*======================================================
  function return the omega list
  =======================================================*/
  double* Material::getOmegaList(){
    return omegaList_;
  }
  /*======================================================
  function return the number of omega points
  =======================================================*/
  int Material::getNumOfOmega(){
    return numOfOmega_;
  }
  /*======================================================
  function set the name of the material
  @args:
  name: the name of the material
  =======================================================*/
  void Material::setName(const std::string name){
    name_ = name;
  }
  /*======================================================
  function set the omega of the material
  @args:
  omegaList: the omegaList of the material
  numOfOmega: the number of omega points
  =======================================================*/
  void Material::setOmega(const double* omegaList, int numOfOmega){
    numOfOmega_ = numOfOmega;
    omegaList_ = new double[numOfOmega_];
    std::copy(omegaList, omegaList + numOfOmega_, omegaList_);
  }
  /*======================================================
  function set the epsilon values of the material
  @args:
  epsilonList: the epsilonList of the material
  numOfOmega: the number of omega points
  =======================================================*/
  void Material::setEpsilon(const dcomplex* epsilonList, int numOfOmega){
    numOfOmega_ = numOfOmega;
    epsilonList_ = new dcomplex[numOfOmega_];
    std::copy(epsilonList, epsilonList + numOfOmega_, epsilonList_);
  }

  /*======================================================
  Implementaion of the Layer class
  =======================================================*/
  Layer::Layer(const Ptr<Material> material, const double thickness) :
    thickness_(thickness), pattern_(PLANAR_), source_(ISNOTSOURCE_){
    backGround_ = material;
  }

  /*======================================================
  This is a thin wrapper for the usage of smart pointer
  =======================================================*/
  Ptr<Layer> Layer::instanceNew(
    const Ptr<Material> material,
    const double thickness
  ){
    return new Layer(material, thickness);
  }
  /*======================================================
  Layer constructor with back ground material
  =======================================================*/
  Layer::Layer(const Ptr<Material> material) :
    source_(ISNOTSOURCE_), thickness_(0), pattern_(PLANAR_){
    //backGround_ = new Material(*material);
    backGround_ = material;
  }
  /*======================================================
  This is a thin wrapper for the usage of smart pointer
  =======================================================*/
  Ptr<Layer> Layer::instanceNew(
    const Ptr<Material> material
  ){
    return new Layer(material);
  }
  /*======================================================
  Plain layer constructor
  =======================================================*/
  Layer::Layer() :
    backGround_(nullptr), source_(ISNOTSOURCE_), thickness_(0), pattern_(PLANAR_){}

  /*======================================================
  This is a thin wrapper for the usage of smart pointer
  =======================================================*/
  Ptr<Layer> Layer::instanceNew()
  {
    return new Layer();
  }

  /*======================================================
  destructor
  =======================================================*/
  Layer::~Layer(){
  }
  /*======================================================
  copy constructor
  =======================================================*/
  Layer::Layer(const Layer& layer){
  }
  /*======================================================
  set background by material
  @args:
  material: the background material
  =======================================================*/
  void Layer::setBackGround(const Ptr<Material> material){
    //backGround_ = new Material(*material);
    backGround_ = material;
  }
  /*======================================================
  set thickness of the layer
  @args:
  thickness: the thickness of the layer
  =======================================================*/
  void Layer::setThickness(const double thickness){
    thickness_ = thickness;
  }
  /*======================================================
  set the layer to be the source
  =======================================================*/
  void Layer::setIsSource(){
    source_ = ISSOURCE_;
  }
  /*======================================================
  set the layer not to be the source
  =======================================================*/
  void Layer::setIsNotSource(){
    source_ = ISNOTSOURCE_;
  }
  /*======================================================
  check whether the layer is a source
  =======================================================*/
  SOURCE Layer::checkIsSource(){
    return source_;
  }
  /*======================================================
  get the background material
  =======================================================*/
  Ptr<Material> Layer::getBackGround(){
    if(backGround_ == nullptr){
      throw UTILITY::NullPointerException("Backgroud not set yet!");
    }
    return backGround_;
  }
  /*======================================================
  get the material by the name
  =======================================================*/
  Ptr<Material> Layer::getMaterialByName(const std::string name){
    for(const_MaterialIter it = this->getVecBegin(); it != this->getVecEnd(); it++){
      if(name.compare((*it)->getName()) == 0) return *it;
    }
    return nullptr;
  }
  /*======================================================
  get the number of materials other than background
  =======================================================*/
  int Layer::getNumOfMaterial(){
    return materialVec_.size();
  }
  /*======================================================
  get the thickness of the layer
  =======================================================*/
  double Layer::getThickness(){
    return thickness_;
  }
  /*======================================================
  set the pattern of the layer
  =======================================================*/
  PATTEN Layer::getPattern(){
    return pattern_;
  }
  /*======================================================
  material iterator begin
  =======================================================*/
  const_MaterialIter Layer::getVecBegin(){
    return materialVec_.cbegin();
  }
  /*======================================================
  material iterator end
  =======================================================*/
  const_MaterialIter Layer::getVecEnd(){
    return materialVec_.cend();
  }
  /*======================================================
  pattern parameter iterator begin
  =======================================================*/
  const_PatternIter Layer::getArg1Begin(){
    return args1_.cbegin();
  }
  /*======================================================
  pattern parameter iterator end
  =======================================================*/
  const_PatternIter Layer::getArg2Begin(){
    return args2_.cbegin();
  }
  /*======================================================
  pattern parameter iterator begin
  =======================================================*/
  const_PatternIter Layer::getArg1End(){
    return args1_.cend();
  }
  /*======================================================
  pattern parameter iterator end
  =======================================================*/
  const_PatternIter Layer::getArg2End(){
    return args2_.cend();
  }
  /*======================================================
  add a rectangle pattern
  @args:
  material: the material used for this part of the pattern
  args1: the position of centers (x,y)
  args2: the widths in x and y directions
  =======================================================*/
  void Layer::addRectanlgePattern(
    const Ptr<Material> material,
    const double args1[2],
    const double args2[2]
  ){
    pattern_ = RECTANGLE_;
    materialVec_.push_back(material);
    args1_.push_back(std::make_pair(args1[0], args1[1]));
    args2_.push_back(std::make_pair(args2[0], args2[1]));
  }

  /*======================================================
  add a circular pattern
  @args:
  material: the material used for this part of the pattern
  args1: the position of centers (x,y)
  radius: the radius of the circle
  =======================================================*/
  void Layer::addCirclePattern(
    const Ptr<Material> material,
    const double args[2],
    const double radius
  ){
    pattern_ = CIRCLE_;
    materialVec_.push_back(material);
    args1_.push_back(std::make_pair(args[0], radius));
    args2_.push_back(std::make_pair(args[1], radius));

  }
  /*======================================================
  add a grating pattern
  @args:
  material: the material used for this part of the pattern
  start: the start position of the this pattern
  end: the end position of this pattern
  =======================================================*/
  void Layer::addGratingPattern(
    const Ptr<Material> material,
    const double start,
    const double end
  ){
    pattern_ = GRATING_;
    materialVec_.push_back(material);
    args1_.push_back(std::make_pair(start, end));
  }

  /*======================================================
  Implementaion of the structure class
  =======================================================*/
  Structure::Structure(){
    period_ = new double[2];
    period_[0] = 0;
    period_[1] = 0;
  }
  /*======================================================
  This is a thin wrapper for the usage of smart pointer
  =======================================================*/
  Ptr<Structure> Structure::instanceNew(){
    return new Structure();
  }
  /*======================================================
  destructor
  =======================================================*/
  Structure::~Structure(){
    delete[] period_;
    period_ = nullptr;
  }
  /*======================================================
  copy constructor
  =======================================================*/
  Structure::Structure(const Structure& structure){
    for(const_LayerIter it = structure.layerMap_.cbegin(); it != structure.layerMap_.cend(); it++){
      layerMap_.insert(LayerMap::value_type(it->first, it->second));
    }
  }
  /*======================================================
  set periodicity of the system
  @args:
  p1: periodicity in x direction
  p2: periodicity in y direction, for 1D it is 0
  =======================================================*/
  void Structure::setPeriodicity(const double p1, const double p2){
    period_[0] = p1;
    period_[1] = p2;
  }
  /*======================================================
  function adding a layer to the structure
  @args:
  layer: the added layer
  =======================================================*/
  void Structure::addLayer(const Ptr<Layer> layer){
    int size = layerMap_.size();
    layerMap_.insert(LayerMap::value_type(size, layer));
  }
  /*======================================================
  function getting a layer by its index
  @args:
  index: the index of the wanted layer
  =======================================================*/
  Ptr<Layer> Structure::getLayerByIndex(const int index){
    if(index >= this->getNumOfLayer()){
      throw UTILITY::RangeException(std::to_string(index) + ": out of range!");
    }
    return layerMap_.at(index);
  }
  /*======================================================
  function getting the number of layers
  =======================================================*/
  int Structure::getNumOfLayer(){
    return layerMap_.size();
  }
  /*======================================================
  function getting the thickness list
  =======================================================*/
  void Structure::getThicknessList(double* thicknessList){
    int count = 0;
    for(const_LayerIter it = this->getMapBegin(); it != this->getMapEnd(); it++){
      thicknessList[count] = (it->second)->getThickness();
      count++;
    }
    return;
  }
  /*======================================================
  layer iterator begin
  =======================================================*/
  const_LayerIter Structure::getMapBegin(){
    return layerMap_.cbegin();
  }
  /*======================================================
  layer iterator end
  =======================================================*/
  const_LayerIter Structure::getMapEnd(){
    return layerMap_.cend();
  }
  /*======================================================
  get the periodicity of the system
  =======================================================*/
  double* Structure::getPeriodicity(){
    return period_;
  }

}