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
  // function returning the POV for one pattern
  // @args:
  // type: the type of the pattern
  /*==============================================*/
  std::string Layer::getPOVRayForPattern(const Pattern pattern){
    double t = thickness_ * 1e6;
    std::string output = "";
    switch (pattern.type_) {
      case GRATING_:{
        output += std::string("box{\n\t<-1, -1, 0>, <1, 1, ") + std::to_string(t) + std::string(">, 1\n");
        output += std::string("\tscale +x*") + std::to_string(pattern.arg1_.first * 1e6 /2) + std::string("\n\tscale +x*0");
        output += std::string("\ttranslate +x*\n") + std::to_string(pattern.arg2_.first * 1e6) + std::string("\ttranslate +y*0\n");
        break;
      }
      case RECTANGLE_:{
        output += std::string("box{\n\t<-1, -1, 0>, <1, 1, ") + std::to_string(t) + std::string(">, 1\n");
        output += std::string("\tscale +x*") + std::to_string(pattern.arg1_.first * 1e6 /2) + std::string("\n\tscale +x*") + std::to_string(pattern.arg1_.second * 1e6/2);
        output += std::string("\ttranslate +x*\n") + std::to_string(pattern.arg2_.first * 1e6) + std::string("\ttranslate +y*\n") + std::to_string(pattern.arg2_.second * 1e6);
        break;
      }
      case CIRCLE_:{
        output += std::string("cylinder{\n\t<0,0,0>, <0,0,") + std::to_string(t) + std::string(">, \n");
        output += std::string("\ttranslate +x*\n") + std::to_string(pattern.arg1_.first * 1e6) + std::string("\ttranslate +y*\n") + std::to_string(pattern.arg2_.first * 1e6);
        break;
      }
      default: break;
    }
    output += "}\n";
    return output;
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
  void Structure::getThicknessList(double* thicknessList){
    int count = 0;
    for(const_LayerIter it = this->getLayersBegin(); it != this->getLayersEnd(); it++){
      thicknessList[count] = (it->second)->getThickness();
      count++;
    }
    return;
  }
  /*==============================================*/
  // function saving the structure to a POVRay file
  // @args:
  // outfile: the output file name, should end with .pov
  /*==============================================*/
  void Structure::getPOVRay(const std::string outfile){
    double charsize = std::max(period_[0], period_[1]) * 1e6;
    // get Voronoi defining points
  	double vorpts[8];
  	vorpts[0] = period_[0] * 1e6;
  	vorpts[1] = 0;
  	vorpts[2] = period_[1] * 1e6;
  	vorpts[3] = 0;
  	vorpts[4] = vorpts[0] + vorpts[2];
  	vorpts[5] = 0;
  	vorpts[6] = vorpts[0] - vorpts[2];
  	vorpts[7] = 0;
  	if(std::abs(vorpts[6]) < std::abs(vorpts[4])){
  		vorpts[4] = vorpts[6];
  	}
    // output preample
    std::string outputString = "// -w320 -h240\n"
      "\n"
      "#version 3.7;\n"
      "\n"
      "#include \"colors.inc\"\n"
      "#include \"textures.inc\"\n"
      "#include \"shapes.inc\"\n"
      "\n"
      "global_settings {max_trace_level 5 assumed_gamma 1.0}\n"
      "\n"
      "camera {\n"
      "\tlocation<";

    outputString += std::to_string(-3.0*charsize) + std::string(",")
      + std::to_string(6.0*charsize) + std::string(",")
      + std::to_string(-9.0*charsize);

    outputString += ">\n\tdirection <0, 0,  2.25>\n"
      "\tright x*1.33\n"
      "\tlook_at <0,0,0>\n"
      "}\n"
      "\n"
      "#declare Dist=80.0;\n"
      "light_source {< -25, 50, -50> color White\n"
      "\tfade_distance Dist fade_power 2\n"
      "}\n"
      "light_source {< 50, 10,  -4> color Gray30\n"
      "\tfade_distance Dist fade_power 2\n"
      "}\n"
      "light_source {< 0, 100,  0> color Gray30\n"
      "\tfade_distance Dist fade_power 2\n"
      "}\n"
      "\n"
      "sky_sphere {\n"
      "\tpigment {\n"
      "\t\tgradient y\n"
      "\t\tcolor_map {\n"
      "\t\t\t[0, 1  color White color White]\n"
      "\t\t}\n"
      "\t}\n"
      "}\n"
      "\n"
      "#declare Xaxis = union{\n"
      "\tcylinder{\n"
      "\t\t<0,0,0>,<0.8,0,0>,0.05\n"
      "\t}\n"
      "\tcone{\n"
      "\t\t<0.8,0,0>, 0.1, <1,0,0>, 0\n"
      "\t}\n"
      "\ttexture { pigment { color Red } }\n"
      "}\n"
      "#declare Yaxis = union{\n"
      "\tcylinder{\n"
      "\t\t<0,0,0>,<0,0.8,0>,0.05\n"
      "\t}\n"
      "\tcone{\n"
      "\t\t<0,0.8,0>, 0.1, <0,1,0>, 0\n"
      "\t}\n"
      "\ttexture { pigment { color Green } }\n"
      "}\n"
      "#declare Zaxis = union{\n"
      "\tcylinder{\n"
      "\t<0,0,0>,<0,0,0.8>,0.05\n"
      "\t}\n"
      "\tcone{\n"
      "\t\t<0,0,0.8>, 0.1, <0,0,1>, 0\n"
      "\t}\n"
      "\ttexture { pigment { color Blue } }\n"
      "}\n"
      "#declare Axes = union{\n"
      "\tobject { Xaxis }\n"
      "\tobject { Yaxis }\n"
      "\tobject { Zaxis }\n"
      "}\n";

    // needs to collect material
    for(MaterialMap::const_iterator it = materialMap_.cbegin(); it != materialMap_.cend(); it++){
      std::string name = it->first;
      if(!name.compare("air") || !name.compare("vacuum") || !name.compare("Air") || !name.compare("Vacuum")){
        outputString += std::string("#declare Material_") + name
          + std::string("= texture{ pigment{ color transmit 1.0 } }\n");
      }
      else{
        outputString += std::string("#declare Material_") + name
          + std::string("= texture{ pigment{ rgb <")
          + std::to_string((double)rand()/(double)RAND_MAX) + std::string(",")
          + std::to_string((double)rand()/(double)RAND_MAX) + std::string(",")
          + std::to_string((double)rand()/(double)RAND_MAX) + std::string(">}\n");
      }
    }
    // needs to collect layer info
    int layerCounter = 0;
    for(int i = 0; i < this->getNumOfLayer(); i++){
      Ptr<Layer> layer = this->getLayerByIndex(i);
      std::string name = layer->getName();
      outputString += std::string("#declare Layer_") + std::to_string(layerCounter++) + std::string(" = union{\n");
      outputString += std::string("\tdifference{\n") + std::string("\t\tintersection{\n");

      for(int j = 0; j < 3; j++){
        double dist = 0.5 * hypot(vorpts[2*i], vorpts[2*i+1]);
        outputString += std::string("\t\t\tplane{ <") + std::to_string(vorpts[2*i]) + std::string(",")
          + std::to_string(vorpts[2*i+1]) + std::string(",0>, ") + std::to_string(dist) + std::string(" }\n");
        outputString += std::string("\t\t\tplane{ <") + std::to_string(-vorpts[2*i]) + std::string(",")
          + std::to_string(-vorpts[2*i+1]) + std::string(",0>, ") + std::to_string(dist) + std::string(" }\n");
      }
      outputString += std::string("\t\t\tplane{ <0,0,-1>, 0 }\n") + std::string("\t\t\tplane{ <0,0,1>, ")
        + std::to_string(layer->getThickness()) + std::string("}\n\t\t}\n");

      outputString += std::string("// nshapes = ") + std::to_string(layer->getNumOfMaterial()) + std::string("%d\n");

      for(const_PatternIter it = layer->getPatternsBegin(); it != layer->getPatternsEnd(); it++){
        outputString += layer->getPOVRayForPattern((*it));
      }
      outputString += std::string("\t}\n");
    }
    // Output postamble
    outputString += std::string("#declare Layers = union {\n");
    layerCounter = 0;
    for(LayerMap::const_iterator it = layerMap_.cbegin(); it != layerMap_.cend(); it++){
      Ptr<Layer> layer = it->second;
      std::string name = layer->getName();
      std::string comment = (0 < layerCounter && layerCounter + 1 < this->getNumOfLayer()) ? "" : "//";
      outputString += std::string("\t") + comment + std::string("object{ Layer_") + std::to_string(layerCounter++) + std::string(" }\n");
    }

    outputString += "}\n"
      "\n"
      "Axes\n"
      "Layers\n";
    // output to file
    std::ofstream file(outfile);
    file << outputString;
    file.close();

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