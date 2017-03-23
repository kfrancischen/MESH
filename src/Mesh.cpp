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
  /*==============================================*/
  // Constructor of the FileLoader class
  /*==============================================*/
  FileLoader::FileLoader(const int numOfOmega){
    if(numOfOmega != 0){
      preSet_ = true;
      numOfOmega_ = numOfOmega;
    }
    else{
      numOfOmega_ = 0;
    }
    omegaList_ = nullptr;
    //epsilonList_ = nullptr;
  }
  /*==============================================*/
  // This is a thin wrapper for the usage of smart pointer
  /*==============================================*/
  Ptr<FileLoader> FileLoader::instanceNew(){
    return new FileLoader();
  }
  /*==============================================*/
  // This is a thin wrapper for the usage of smart pointer
  /*==============================================*/
  Ptr<FileLoader> FileLoader::instanceNew(const int numOfOmega){
    return new FileLoader(numOfOmega);
  }
  /*==============================================*/
  // Function checking a line contains more than 3 spaces
  /*==============================================*/
  EPSTYPE checkType(std::string line){
    int numOfSpace = 0;
    for(int i = 0; i < line.length(); i++){
      if(isspace(line.at(i))) numOfSpace++;
    }
    switch (numOfSpace) {
      case 2: return SCALAR_;
      case 6: return DIAGONAL_;
      case 10: return TENSOR_;
      default:{
        throw UTILITY::UnknownTypeException("Input type wrong!");
      }
    }
    return SCALAR_;
  }
  /*==============================================*/
  // Function reads a file from disk
  // @args
  // fileName: the name of the file
  /*==============================================*/
  void FileLoader::load(const std::string fileName){
    std::ifstream inputFile(fileName);
    if(!inputFile.good()){
      throw UTILITY::FileNotExistException(fileName + " not exists!");
    }
    std::string line;
    int count = 0;
    EPSTYPE type = SCALAR_;
    while(std::getline(inputFile, line)){
      count++;
      type = checkType(line);
    }
    inputFile.close();

    if(!preSet_){
      if(numOfOmega_ != 0 && numOfOmega_ != count){
        throw UTILITY::StorageException(fileName + " wrong length!");
      }
      numOfOmega_ = count;
    }
    else{
      if(numOfOmega_ > count){
        throw UTILITY::RangeException("not enought omega points!");
      }
    }
    if(omegaList_ == nullptr){
      omegaList_ = new double[numOfOmega_];
    }

    if(epsilonList_.epsilonVals == nullptr){
      epsilonList_.epsilonVals = new EpsilonVal[numOfOmega_];
      epsilonList_.type_ = type;
    }
    std::ifstream inputFile2(fileName);
    for(int i = 0; i < numOfOmega_; i++){
      inputFile2 >> omegaList_[i];
      if(epsilonList_.type_ == SCALAR_){
        double imag;
        inputFile2 >> epsilonList_.epsilonVals[i].scalar[0] >> imag;
        epsilonList_.epsilonVals[i].scalar[1] = -imag;
      }
      else if(epsilonList_.type_ == DIAGONAL_){
        double real1, imag1, real2, imag2, real3, imag3;
        inputFile2 >> real1 >> imag1 >> real2 >> imag2 >> real3 >> imag3;
        epsilonList_.epsilonVals[i].diagonal[0] = real1;
        epsilonList_.epsilonVals[i].diagonal[1] = -imag1;
        epsilonList_.epsilonVals[i].diagonal[2] = real2;
        epsilonList_.epsilonVals[i].diagonal[3] = -imag2;
        epsilonList_.epsilonVals[i].diagonal[4] = real3;
        epsilonList_.epsilonVals[i].diagonal[5] = -imag3;
      }
      else{
        double real1, imag1, real2, imag2, real3, imag3, real4, imag4, real5, imag5;
        inputFile2 >> real1 >> imag1 >> real2 >> imag2 >> real3 >> imag3 >> real4 >> imag4 >> real5 >> imag5;
        epsilonList_.epsilonVals[i].tensor[0] = real1;
        epsilonList_.epsilonVals[i].tensor[1] = -imag1;
        epsilonList_.epsilonVals[i].tensor[2] = real2;
        epsilonList_.epsilonVals[i].tensor[3] = -imag2;
        epsilonList_.epsilonVals[i].tensor[4] = real3;
        epsilonList_.epsilonVals[i].tensor[5] = -imag3;
        epsilonList_.epsilonVals[i].tensor[6] = real4;
        epsilonList_.epsilonVals[i].tensor[7] = -imag4;
        epsilonList_.epsilonVals[i].tensor[8] = real5;
        epsilonList_.epsilonVals[i].tensor[9] = -imag5;
      }
    }
    inputFile2.close();
  }

  /*==============================================*/
  // Function returning the omega values
  /*==============================================*/
  double* FileLoader::getOmegaList(){
    return omegaList_;
  }
  /*==============================================*/
  // Function returning the epsilon values
  /*==============================================*/
  EPSILON FileLoader::getEpsilonList(){
    return epsilonList_;
  }
  /*==============================================*/
  // Function returning number of omega points
  /*==============================================*/
  int FileLoader::getNumOfOmega(){
    return numOfOmega_;
  }
  /*==============================================*/
  // Class destructor
  /*==============================================*/
  FileLoader::~FileLoader(){
    delete [] omegaList_;
    omegaList_ = nullptr;
    delete [] epsilonList_.epsilonVals;
    epsilonList_.epsilonVals = nullptr;
  }
  /*==============================================*/
  // This function wraps the data for quad_gaussian_kronrod
  // @args:
  // kx: the kx value (normalized)
  // wrapper: wrapper for all the arguments wrapped in wrapper
  /*==============================================*/
  static double wrapperFunQuadgk(const double kx, ArgWrapper* data){
    ArgWrapper wrapper = *data;
    return kx * poyntingFlux(
      wrapper.omega,
      wrapper.thicknessList,
      kx,
      0,
      wrapper.EMatrices,
      wrapper.grandImaginaryMatrices,
      wrapper.eps_zz_Inv,
      wrapper.Gx_mat,
      wrapper.Gy_mat,
      wrapper.sourceList,
      wrapper.targetLayer,
      1
    );
  }
  /*==============================================*/
  // This function wraps the data for quad_legendre
  // @args:
  // kx: the kx value (normalized)
  // data: wrapper for all the arguments wrapped in wrapper
  /*==============================================*/
  static double wrapperFunQuadgl(const double kx, void* data){
    return wrapperFunQuadgk(kx, (ArgWrapper*) data);
  }
  /*======================================================*/
  // Implementaion of the parent simulation super class
  /*=======================================================*/
  Simulation::Simulation() : nGx_(0), nGy_(0), numOfOmega_(0), structure_(nullptr),
  Phi_(nullptr), omegaList_(nullptr), kxStart_(0), kxEnd_(0), kyStart_(0), kyEnd_(0), numOfKx_(0), numOfKy_(0)
  {
    output_ = "output.txt";
    targetLayer_ = -1;
    dim_ = NO_;
    structure_ = Structure::instanceNew();
    fileLoader_ = FileLoader::instanceNew();
  }
  /*==============================================*/
  // Class destructor
  /*==============================================*/
  Simulation::~Simulation(){
    delete[] Phi_;
    Phi_ = nullptr;
  }

  /*==============================================*/
  // This function saves data to disk
  /*==============================================*/
  void Simulation::saveToFile(){
    std::ofstream outputFile(output_);
    for(int i = 0; i < numOfOmega_; i++){
      outputFile << omegaList_[i] << "\t" << Phi_[i] << std::endl;
    }
    outputFile.close();
  }

  /*==============================================*/
  // This function set the output file name
  // @args:
  // name: the output file name
  /*==============================================*/
  void Simulation::setOutputFile(std::string name){
    output_ = name;
  }
  /*==============================================*/
  // This function set the periodicity
  // @args:
  // p1: the periodicity along x direction
  // p2: the periodicity along x direction
  /*==============================================*/
  void Simulation::setPeriodicity(const double p1, const double p2){
    period_[0] = p1;
    period_[1] = p2;
  }
  /*==============================================*/
  // This function adds material to the system
  // @args:
  // name: the name of the material
  // infile: the input file containing the properties of the material
  /*==============================================*/
  void Simulation::addMaterial(const std::string name, const std::string infile){
    if(materialInstanceMap_.find(name) != materialInstanceMap_.cend()){
      throw UTILITY::NameInUseException(name + ": Material already exist!");
      return;
    }

    fileLoader_->load(infile);
    Ptr<Material> material = Material::instanceNew(name,
     fileLoader_->getOmegaList(), 
     fileLoader_->getEpsilonList(),
     fileLoader_->getNumOfOmega()
    );
    materialInstanceMap_.insert(MaterialInstanceMap::value_type(name, material));
  }

  /*==============================================*/
  // This function reset the dielectric of a material
  // @args:
  // name: the name of the material
  // epsilon: the new epsilon values
  // type: the type of the epsilon, one of ['scalar', 'diagonal', 'tensor']
  // Note: the material should have already existed in the system
  // If type == 'scalar', epsilon should have size [numOfOmega][2]
  // If type == 'diagonal', epsilon should have size [numOfOmega][6]
  // If type == 'tensor', epsilon should have size [numOfOmega][10]
  /*==============================================*/
  void Simulation::setMaterial(const std::string name, const double** epsilon, const std::string type){
    if(materialInstanceMap_.find(name) == materialInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(name + ": Material does not exist!");
      return;
    }
    Ptr<Material> material = materialInstanceMap_.find(name)->second;
    int numOfOmega = material->getNumOfOmega();
    EPSILON newEpsilon;
    newEpsilon.epsilonVals = new EpsilonVal[numOfOmega];
    // if a scalar
    if(type == "scalar"){
      for(int i = 0; i < numOfOmega; i++){
        newEpsilon.epsilonVals[i].scalar[0] = epsilon[i][0];
        newEpsilon.epsilonVals[i].scalar[1] = epsilon[i][1];
      }
      newEpsilon.type_ = SCALAR_;
    }
    // if a diagonal
    else if(type == "diagonal"){
      for(int i = 0; i < numOfOmega; i++){
        for(int j = 0; j < 6; j++){
          newEpsilon.epsilonVals[i].scalar[j] = epsilon[i][j];
        }
      }
      newEpsilon.type_ = DIAGONAL_;
    }
    // if a tensor
    else if(type == "tensor"){
      for(int i = 0; i < numOfOmega; i++){
        for(int j = 0; j < 10; j++){
          newEpsilon.epsilonVals[i].scalar[j] = epsilon[i][j];
        }
      }
      newEpsilon.type_ = TENSOR_;
    }
    // error
    else{
      throw UTILITY::AttributeNotSupportedException("Please choose 'type' from 'scalar', 'diagonal' or 'tensor'!");
    }

    material->setEpsilon(newEpsilon, numOfOmega);
    delete [] newEpsilon.epsilonVals;
    newEpsilon.epsilonVals = nullptr;
  }
  /*==============================================*/
  // This function adds a new layer to the system
  // @args:
  // name: the name of the layer
  // thick: the thickness of the layer
  // materialName: the name of the material
  /*==============================================*/
  void Simulation::addLayer(const std::string name, const double thick, const std::string materialName){
    if(materialInstanceMap_.find(materialName) == materialInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(materialName + ": Material does not exist!");
      return;
    }
    if(layerInstanceMap_.find(name) != layerInstanceMap_.cend()){
      throw UTILITY::NameInUseException(name + ": Layer already exists!");
      return;
    }
    Ptr<Material> material = materialInstanceMap_.find(materialName)->second;
    Ptr<Layer> layer = Layer::instanceNew(name, material, thick);
    structure_->addLayer(layer);
    layerInstanceMap_.insert(LayerInstanceMap::value_type(name, layer));
  }
  /*==============================================*/
  // This function change the background material of a layer and its thickness
  // @args:
  // name: the name of the layer
  // thick: the new thickness
  // materialName: the new background
  /*==============================================*/
  void Simulation::setLayer(const std::string name, const double thick, const std::string materialName){
    if(materialInstanceMap_.find(materialName) == materialInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(materialName + ": Material does not exist!");
      return;
    }
    if(layerInstanceMap_.find(name) == layerInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(name + ": Layer does not exist!");
      return;
    }
    Ptr<Material> material = materialInstanceMap_.find(materialName)->second;
    Ptr<Layer> layer = layerInstanceMap_.find(name)->second;
    layer->setBackGround(material);
    layer->setThickness(thick);
  }
  /*==============================================*/
  // This function change the thickness
  // @args:
  // name: the name of the layer
  // thick: the new thickness
  /*==============================================*/
  void Simulation::setLayerThickness(const std::string name, const double thick){
    if(layerInstanceMap_.find(name) == layerInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(name + ": Layer does not exist!");
      return;
    }
    Ptr<Layer> layer = layerInstanceMap_.find(name)->second;
    layer->setThickness(thick);
  }
  /*==============================================*/
  // This function make a copy of an existing layer
  // @args:
  // name: the name of the copied layer
  // originalName: the name of the original layer
  /*==============================================*/
  void Simulation::addLayerCopy(const std::string name, const std::string originalName){
    if(layerInstanceMap_.find(originalName) == layerInstanceMap_.cend()){
      throw UTILITY::NameInUseException(originalName + ": Layer does not exist!");
      return;
    }
    if(layerInstanceMap_.find(name) != layerInstanceMap_.cend()){
      throw UTILITY::NameInUseException(name + ": cannot add a layer that already exist!");
      return;
    }   
    Ptr<Layer> originalLayer = layerInstanceMap_.find(originalName)->second;
    Ptr<Layer> newLayer = originalLayer->layerCopy(name);
    layerInstanceMap_.insert(LayerInstanceMap::value_type(name, newLayer));
    structure_->addLayer(newLayer);
  }
  /*==============================================*/
  // This function deletes an existing layer
  // @args:
  // name: the name of the layer
  /*==============================================*/
  void Simulation::deleteLayer(const std::string name){
    if(layerInstanceMap_.find(name) == layerInstanceMap_.cend()){
      throw UTILITY::NameInUseException(name + ": Layer does not exist!");
      return;
    }
    structure_->deleteLayerByName(name);
  }

  /*==============================================*/
  // This function add grating to a layer
  // @args:
  // layerName: the name of the layer
  // materialName: the name of the material
  // center: the center of the grating
  // width: the width of the grating
  /*==============================================*/
  void Simulation::setLayerPatternGrating(
    const std::string layerName,
    const std::string materialName,
    const double center,
    const double width
  ){
    if(materialInstanceMap_.find(materialName) == materialInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(materialName + ": Material does not exist!");
      return;
    }
    if(layerInstanceMap_.find(layerName) == layerInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(layerName + ": Layer does not exist!");
      return;
    }
    Ptr<Material> material = materialInstanceMap_.find(materialName)->second;
    Ptr<Layer> layer = layerInstanceMap_.find(layerName)->second;
    layer->addGratingPattern(material, center, width);
  }
  /*==============================================*/
  // This function add rectangle pattern to a layer
  // @args:
  // layerName: the name of the layer
  // materialName: the name of the material
  // centerx: the center of the rectangle in x direction
  // centery: the center of the rectangle in y direction
  // widthx: the width of the rectangle in x direction
  // widthy: the width of the rectangle in y direction
  /*==============================================*/  
  void Simulation::setLayerPatternRectangle(
    const std::string layerName,
    const std::string materialName,
    const double centerx,
    const double centery,
    const double widthx,
    const double widthy
  ){
    if(materialInstanceMap_.find(materialName) == materialInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(materialName + ": Material does not exist!");
      return;
    }
    if(layerInstanceMap_.find(layerName) == layerInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(layerName + ": Layer does not exist!");
      return;
    }
    Ptr<Material> material = materialInstanceMap_.find(materialName)->second;
    Ptr<Layer> layer = layerInstanceMap_.find(layerName)->second;
    double arg1[2] = {centerx, centery};
    double arg2[2] = {widthx, widthy};
    layer->addRectanlgePattern(material, arg1, arg2);
  }
  /*==============================================*/
  // This function add circle pattern to a layer
  // @args:
  // layerName: the name of the layer
  // materialName: the name of the material
  // centerx: the center of the circle in x direction
  // centery: the center of the circle in y direction
  // radius: the radius of the circle
  /*==============================================*/ 
  void Simulation::setLayerPatternCircle(
    const std::string layerName,
    const std::string materialName,
    const double centerx,
    const double centery,
    const double radius
  ){
    if(materialInstanceMap_.find(materialName) == materialInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(materialName + ": Material does not exist!");
      return;
    }
    if(layerInstanceMap_.find(layerName) == layerInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(layerName + ": Layer does not exist!");
      return;
    }
    Ptr<Material> material = materialInstanceMap_.find(materialName)->second;
    Ptr<Layer> layer = layerInstanceMap_.find(layerName)->second;
    double arg1[2] = {centerx, centery};
    layer->addCirclePattern(material, arg1, radius);
  }
  /*==============================================*/
  // This function sets a layer as the source
  // @args:
  // name: the name of the source layer
  /*==============================================*/
  void Simulation::setSourceLayer(const std::string name){
    if(layerInstanceMap_.find(name) == layerInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(name + ": Layer does not exist!");
      return;
    }
    Ptr<Layer> layer = layerInstanceMap_.find(name)->second;
    layer->setIsSource();
  }
  /*==============================================*/
  // This function sets the probe layer
  // @args:
  // name: the name of the probe layer
  /*==============================================*/
  void Simulation::setProbeLayer(const std::string name){
    if(layerInstanceMap_.find(name) == layerInstanceMap_.cend()){
      throw UTILITY::IllegalNameException(name + ": Layer does not exist!");
      return;
    }
    Ptr<Layer> layer = layerInstanceMap_.find(name)->second;
    this->setTargetLayerByLayer(layer);
  }
  /*==============================================*/
  // This function sets number of positive Gx
  // @args:
  // Gx: number of positive Gx
  /*==============================================*/
  void Simulation::setGx(const int Gx){
    nGx_ = Gx;
  }

  /*==============================================*/
  // This function sets number of positive Gy
  // @args:
  // Gy: number of positive Gy
  /*==============================================*/
  void Simulation::setGy(const int Gy){
    nGy_ = Gy;
  }

  /*==============================================*/
  // This function sets the target layer by layer
  // @args:
  // index: target layer
  /*==============================================*/
  void Simulation::setTargetLayerByLayer(const Ptr<Layer>& layer){
    for(int i = 0; i < structure_->getNumOfLayer(); i++){
      if(structure_->getLayerByIndex(i) == layer){
        targetLayer_ = i;
        return;
      }
    }
  }

  /*==============================================*/
  // This function cleans up the simulation
  /*==============================================*/
  void Simulation::resetSimulation(){
    for(size_t i = 0; i < EMatricesVec_.size(); i++){
      EMatricesVec_[i].clear();
      grandImaginaryMatrixVec_.clear();
      eps_zz_Inv_MatrixVec_[i].clear();
    }
    EMatricesVec_.clear();
    grandImaginaryMatrixVec_.clear();
    eps_zz_Inv_MatrixVec_.clear();
  }

  /*==============================================*/
  // This function rebuilds the simulation
  /*==============================================*/
  void Simulation::rebuild(){
    this->resetSimulation();
    this->build();
  }
  /*==============================================*/
  // This function gets the structure
  /*==============================================*/
  Ptr<Structure> Simulation::getStructure(){
    return structure_;
  }
  
  /*==============================================*/
  // This function gets the Phi at given kx and ky
  // @args:
  // omegaIndex: the index of omega
  // kx: the kx value, normalized
  // ky: the ky value, normalized
  // @note
  // used by grating and patterning
  // N: the number of total G
  /*==============================================*/
  double Simulation::getPhiAtKxKy(const int omegaIdx, const double kx, const double ky){
    if(omegaIdx >= numOfOmega_){
      throw UTILITY::RangeException(std::to_string(omegaIdx) + ": out of range!");
    }
    int N = getN(nGx_, nGy_);
    return POW3(omegaList_[omegaIdx] / datum::c_0) / POW3(datum::pi) / 2.0 *
      poyntingFlux(omegaList_[omegaIdx] / datum::c_0, thicknessListVec_, kx, ky, EMatricesVec_[omegaIdx],
      grandImaginaryMatrixVec_[omegaIdx], eps_zz_Inv_MatrixVec_[omegaIdx], Gx_mat_, Gy_mat_,
      sourceList_, targetLayer_,N);
  }


  /*==============================================*/
  // This function builds up the matrices
  /*==============================================*/
  void Simulation::build(){
    // check each layer whether a tensor exist
    for(const_LayerIter it = structure_->getMapBegin(); it != structure_->getMapEnd(); it++){
      Ptr<Layer> layer = it->second;
      layer->containTensor(false);
      if(layer->getBackGround()->getType() == TENSOR_){
        layer->containTensor(true);
        break;
      }
      for(const_MaterialIter m_it = layer->getVecBegin(); m_it != layer->getVecEnd(); it++){
        if((*m_it)->getType() == TENSOR_){
          layer->containTensor(true);
          break;
        }
      }
    }
    // essential, get the shared Gx_mat_ and Gy_mat_
    getGMatrices(nGx_, nGy_, period_, Gx_mat_, Gy_mat_, dim_);
    // get constants
    Ptr<Layer> firstLayer = structure_->getLayerByIndex(0);
    Ptr<Material> backGround = firstLayer->getBackGround();
    numOfOmega_ = backGround->getNumOfOmega();
    omegaList_ = backGround->getOmegaList();

    Phi_ = new double[numOfOmega_];
    EMatricesVec_.resize(numOfOmega_);
    grandImaginaryMatrixVec_.resize(numOfOmega_);
    eps_zz_Inv_MatrixVec_.resize(numOfOmega_);

    RCWAMatricesVec eps_xx_MatrixVec(numOfOmega_), eps_xy_MatrixVec(numOfOmega_), eps_yx_MatrixVec(numOfOmega_), eps_yy_MatrixVec(numOfOmega_);
    RCWAMatricesVec im_eps_xx_MatrixVec(numOfOmega_), im_eps_xy_MatrixVec(numOfOmega_), im_eps_yx_MatrixVec(numOfOmega_), im_eps_yy_MatrixVec(numOfOmega_), im_eps_zz_MatrixVec(numOfOmega_);
    int numOfLayer = structure_->getNumOfLayer();
    int N = getN(nGx_, nGy_);

    for(int i = 0; i < numOfLayer; i++){
      Ptr<Layer> layer = structure_->getLayerByIndex(i);
      switch (layer->getPattern()) {
        /*************************************/
        // if the pattern is a plane
        /*************************************/
        case PLANAR_:{

          FMM::transformPlanar(
            eps_xx_MatrixVec,
            eps_xy_MatrixVec,
            eps_yx_MatrixVec,
            eps_yy_MatrixVec,
            eps_zz_Inv_MatrixVec_,
            im_eps_xx_MatrixVec,
            im_eps_xy_MatrixVec,
            im_eps_yx_MatrixVec,
            im_eps_yy_MatrixVec,
            im_eps_zz_MatrixVec,
            layer,
            N
          );

          break;
        }
        /*************************************/
        // if the pattern is a grating (1D)
        /************************************/
        case GRATING_:{
          if(options_.FMMRule == SPATIALADAPTIVE_){
            // to be filled
          }

          else{
            FMM::transformGratingNaive(
              eps_xx_MatrixVec,
              eps_xy_MatrixVec,
              eps_yx_MatrixVec,
              eps_yy_MatrixVec,
              eps_zz_Inv_MatrixVec_,
              im_eps_xx_MatrixVec,
              im_eps_xy_MatrixVec,
              im_eps_yx_MatrixVec,
              im_eps_yy_MatrixVec,
              im_eps_zz_MatrixVec,
              layer,
              N,
              period_[0],
              options_.FMMRule == INVERSERULE_
            );
          }
          break;
        }

        /*************************************/
        // if the pattern is a rectangle (2D)
        /************************************/
        case RECTANGLE_:{
          FMM::transformRectangle(
            eps_xx_MatrixVec,
            eps_xy_MatrixVec,
            eps_yx_MatrixVec,
            eps_yy_MatrixVec,
            eps_zz_Inv_MatrixVec_,
            im_eps_xx_MatrixVec,
            im_eps_xy_MatrixVec,
            im_eps_yx_MatrixVec,
            im_eps_yy_MatrixVec,
            im_eps_zz_MatrixVec,
            layer,
            nGx_,
            nGy_,
            period_,
            options_.FMMRule == INVERSERULE_
          );
          break;
        }

      /*************************************/
      // if the pattern is a circle (2D)
      /************************************/
        case CIRCLE_:{
          FMM::transformCircleTensor(
            eps_xx_MatrixVec,
            eps_xy_MatrixVec,
            eps_yx_MatrixVec,
            eps_yy_MatrixVec,
            eps_zz_Inv_MatrixVec_,
            im_eps_xx_MatrixVec,
            im_eps_xy_MatrixVec,
            im_eps_yx_MatrixVec,
            im_eps_yy_MatrixVec,
            im_eps_zz_MatrixVec,
            layer,
            N,
            period_
          );
          break;
        }
        default: break;
      }
    }


    for(int i = 0; i < numOfOmega_; i++){
      getEMatrices(
        EMatricesVec_[i],
        eps_xx_MatrixVec[i],
        eps_xy_MatrixVec[i],
        eps_yx_MatrixVec[i],
        eps_yy_MatrixVec[i],
        numOfLayer,
        N
      );

      getGrandImaginaryMatrices(
        grandImaginaryMatrixVec_[i],
        im_eps_xx_MatrixVec[i],
        im_eps_xy_MatrixVec[i],
        im_eps_yx_MatrixVec[i],
        im_eps_yy_MatrixVec[i],
        im_eps_zz_MatrixVec[i],
        numOfLayer,
        N
      );

    }

    thicknessListVec_ = zeros<RCWAVector>(numOfLayer);
    sourceList_.resize(numOfLayer);
    for(int i = 0; i < numOfLayer; i++){
      thicknessListVec_(i) = (structure_->getLayerByIndex(i))->getThickness();
      sourceList_[i] = (structure_->getLayerByIndex(i))->checkIsSource();
    }

  }

  /*==============================================*/
  // This function computes the flux
  /*==============================================*/
  void Simulation::run(){

    double kxList[numOfKx_], kyList[numOfKy_];
    double scalex[numOfOmega_], scaley[numOfOmega_];
    // here dkx is not normalized
    double dkx = (kxEnd_ - kxStart_) / (numOfKx_ - 1);
    // here kyEnd_ is normalized for 1D case
    double dky = (kyEnd_ - kyStart_) / (numOfKy_ - 1);
    for(int i = 0; i < numOfKx_; i++){
      kxList[i] = kxStart_ + dkx * i;
    }
    for(int i = 0; i < numOfKy_; i++){
      kyList[i] = kyStart_ + dky * i;
    }

    for(int i = 0; i < numOfOmega_; i++){
      switch (dim_) {
        case NO_:{
          scalex[i] = 1;
          scaley[i] = 1;
          break;
        }
        case ONE_:{
          scalex[i] = omegaList_[i] / datum::c_0;
          scaley[i] = 1;
          break;
        }
        case TWO_:{
          scalex[i] = omegaList_[i] / datum::c_0;
          scaley[i] = scalex[i];
          break;
        }
        default: break;
      }
    }

    int totalNum = numOfKx_ * numOfKy_ * numOfOmega_;
    // use dynamic allocation for stack memory problem
    double* resultArray = new double[totalNum];
    MPI_Status status;
    int rank, numProcs, start, end, startPosition, endPosition;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    int offset = 0;
    int chunkSize = totalNum / numProcs;
    int numLeft = totalNum - numProcs * chunkSize;

    /*************************************
    // for master node
    *************************************/

    if(rank == MASTER){
      for(int thread = 1; thread < numProcs; thread++){
        if(numLeft > 0){
          start = thread * chunkSize + offset;
          end = (thread + 1) * chunkSize + offset + 1;
          offset++;
          numLeft--;
        }
        else{
          start = thread * chunkSize + offset;
          end = (thread + 1) * chunkSize + offset;
        }
        if(thread == numProcs - 1) end = totalNum;
        MPI_Send(&start, 1, MPI_INT, thread, SENDTAG, MPI_COMM_WORLD);
        MPI_Send(&end, 1, MPI_INT, thread, SENDTAG, MPI_COMM_WORLD);
        MPI_Send(&resultArray[start], end - start, MPI_DOUBLE, thread, SENDTAG, MPI_COMM_WORLD);
      }

      for(int i = 0; i < chunkSize; i++){
        int omegaIdx = i / (numOfKx_ * numOfKy_);
        int residue = i % (numOfKx_ * numOfKy_);
        int kxIdx = residue / numOfKy_;
        int kyIdx = residue % numOfKy_;
        resultArray[i] = this->getPhiAtKxKy(omegaIdx, kxList[kxIdx] / scalex[omegaIdx], kyList[kyIdx] / scaley[omegaIdx]);
        std::cout << kxList[kxIdx] / scalex[omegaIdx] << "\t" << kyList[kyIdx]  / scaley[omegaIdx] << "\t" << resultArray[i] << std::endl;
      }

      // wait for all the slave process finished
      for(int thread = 1; thread < numProcs; thread++){
        MPI_Recv(&start, 1, MPI_INT, thread, RECVTAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&end, 1, MPI_INT, thread, RECVTAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&resultArray[start], end - start, MPI_DOUBLE, thread, RECVTAG, MPI_COMM_WORLD, &status);
      }
      for(int i = 0; i < numOfOmega_; i++){
        Phi_[i] = 0;
        for(int j = 0; j < numOfKx_ * numOfKy_; j++){
          Phi_[i] += resultArray[i * numOfKx_ * numOfKy_ + j];
        }
        Phi_[i] *= prefactor_ * dkx / scalex[i] * dky / scaley[i];
      }
      this->saveToFile();
    }
    /*************************************
    // for slave nodes
    *************************************/

    else{
      MPI_Recv(&startPosition, 1, MPI_INT, MASTER, SENDTAG, MPI_COMM_WORLD, &status);
      MPI_Recv(&endPosition, 1, MPI_INT, MASTER, SENDTAG, MPI_COMM_WORLD, &status);
      MPI_Recv(&resultArray[startPosition], endPosition - startPosition, MPI_DOUBLE, MASTER, SENDTAG, MPI_COMM_WORLD, &status);

      for(int i = startPosition; i < endPosition; i++){
        int omegaIdx = i / (numOfKx_ * numOfKy_);
        int residue = i % (numOfKx_ * numOfKy_);
        int kxIdx = residue / numOfKy_;
        int kyIdx = residue % numOfKy_;
        resultArray[i] = this->getPhiAtKxKy(omegaIdx, kxList[kxIdx] / scalex[omegaIdx], kyList[kyIdx] / scaley[omegaIdx]);
        std::cout << kxList[kxIdx] / scalex[omegaIdx] << "\t" << kyList[kyIdx]  / scaley[omegaIdx] << "\t" << resultArray[i] << std::endl;
      }

      MPI_Send(&startPosition, 1, MPI_INT, MASTER, RECVTAG, MPI_COMM_WORLD);
      MPI_Send(&endPosition, 1, MPI_INT, MASTER, RECVTAG, MPI_COMM_WORLD);
      MPI_Send(&resultArray[startPosition], endPosition - startPosition, MPI_DOUBLE, MASTER, RECVTAG, MPI_COMM_WORLD);
    }

    delete[] resultArray;
    resultArray = nullptr;
    MPI_Finalize();
  }
  /*==============================================*/
  // Implementaion of the class on planar simulation
  /*==============================================*/
  SimulationPlanar::SimulationPlanar() : Simulation(){
    dim_ = NO_;
    degree_ = DEGREE;
    prefactor_ = 1;
  }
  /*==============================================*/
  // This is a thin wrapper for the usage of smart pointer
  /*==============================================*/
  Ptr<SimulationPlanar> SimulationPlanar::instanceNew(){
    return new SimulationPlanar();
  };
  /*==============================================*/
  // Function setting the integral over kx
  // @args:
  // end: the end of the integration
  /*==============================================*/
  void SimulationPlanar::setKParallelIntegral(const double end){
    kxStart_ = 0;
    numOfKx_ = 0;
    kxEnd_ = end;
  }
  /*==============================================*/
  // Function setting the integral over kx
  // @args:
  // points: number of points considered in the integral
  // end: the end of the integration
  /*==============================================*/
  void SimulationPlanar::setKxIntegral(const int points, const double end){
    kxStart_ = -end;
    numOfKx_ = points;
    kxEnd_ = end;
  }
  /*==============================================*/
  // Function setting the integral over kx (symetric)
  // @args:
  // points: number of points considered in the integral
  // end: the end of the integration
  /*==============================================*/
  void SimulationPlanar::setKxIntegralSym(const int points, const double end){
    kxStart_ = 0;
    numOfKx_ = points;
    kxEnd_ = end;
    prefactor_ *= 2;
  }
  /*==============================================*/
  // Function setting the integral over ky
  // @args:
  // end: the end of the integration
  // points: number of points considered in the integral
  /*==============================================*/
  void SimulationPlanar::setKyIntegral(const int points, const double end){
    kxStart_ = -end;
    numOfKx_ = points;
    kxEnd_ = end;
  }
  /*==============================================*/
  // Function setting the integral over ky (symetric)
  // @args:
  // end: the end of the integration
  // points: number of points considered in the integral
  /*==============================================*/
  void SimulationPlanar::setKyIntegralSym(const int points, const double end){
    kxStart_ = 0;
    numOfKx_ = points;
    kxEnd_ = end;
    prefactor_ *= 2;
  }
  /*==============================================*/
  // Function setting the integral to be the gauss_legendre
  // @args:
  // degree: the degree of gauss_legendre integral, default 512
  /*==============================================*/
  void SimulationPlanar::useQuadgl(const int degree){
    degree_ = degree;
    options_.IntegralMethod = GAUSSLEGENDRE_;
  }
  /*==============================================*/
  // Function setting the integral to be the gauss_kronrod
  // @args:
  // degree: the refinement integral, default 512
  /*==============================================*/
  void SimulationPlanar::useQuadgk(const int degree){
    degree_ = degree;
    options_.IntegralMethod = GAUSSKRONROD_;
  }

  /*==============================================*/
  // This function gets the flux at a given kx
  // @args:
  // omegaIndex: the index of omega
  // KParallel: the KParallel value, normalized
  // NOTE:
  // assuming scalar or
  //  eps_x,  0,   0
  //    0  ,eps_x, 0
  //    0     0,  eps_z
  // make sure you understand your problem whether can be solved by this function
  /*==============================================*/
  double SimulationPlanar::getPhiAtKParallel(const int omegaIdx, const double KParallel){
    if(omegaIdx >= numOfOmega_){
      throw UTILITY::RangeException(std::to_string(omegaIdx) + ": out of range!");
    }
    return POW3(omegaList_[omegaIdx] / datum::c_0) / POW2(datum::pi) * KParallel *
      poyntingFlux(omegaList_[omegaIdx] / datum::c_0, thicknessListVec_, KParallel, 0, EMatricesVec_[omegaIdx],
      grandImaginaryMatrixVec_[omegaIdx], eps_zz_Inv_MatrixVec_[omegaIdx], Gx_mat_, Gy_mat_,
      sourceList_, targetLayer_,1);
  }


  /*==============================================*/
  // This function integrates kx assuming scalar or
  //  eps_x,  0,   0
  //    0  ,eps_x, 0
  //    0     0,  eps_z
  // make sure you understand your problem whether can be solved by this function
  /*==============================================*/
  void SimulationPlanar::runNaive(){
    int rank, numProcs, start, end;
    int offset = 0;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    int displs[numProcs], sendCounts[numProcs];
    int chunkSize = numOfOmega_ / numProcs;
    int numLeft = numOfOmega_ - numProcs * chunkSize;
    double* recvBuf = new double[chunkSize + 1];
    MPI_Datatype dtype;
    // allocating sizes for each processor
    displs[0] = 0;
    sendCounts[0] = chunkSize;
    for(int thread = 1; thread < numProcs; thread++){
      if(numLeft > 0){
        start = thread * chunkSize + offset;
        end = (thread + 1) * chunkSize + offset + 1;
        offset++;
        numLeft--;
      }
      else{
        start = thread * chunkSize + offset;
        end = (thread + 1) * chunkSize + offset;
      }
      if(thread == numProcs - 1) end = numOfOmega_;
      displs[thread] = start;
      sendCounts[thread] = end - start;
    }
    MPI_Type_contiguous(sendCounts[rank], MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    // scattering the Phi to each processor
    MPI_Scatterv(&Phi_[0], sendCounts, displs, MPI_DOUBLE, &recvBuf[0], 1, dtype, MASTER, MPI_COMM_WORLD);
    ArgWrapper wrapper;
    wrapper.thicknessList = thicknessListVec_;
    wrapper.Gx_mat = Gx_mat_;
    wrapper.Gy_mat = Gy_mat_;
    wrapper.sourceList = sourceList_;
    wrapper.targetLayer = targetLayer_;

    for(int i = 0; i < sendCounts[rank]; i++){
      int omegaIdx = i + displs[rank];
      wrapper.omega = omegaList_[omegaIdx] / datum::c_0;
      wrapper.EMatrices = EMatricesVec_[omegaIdx];

      wrapper.grandImaginaryMatrices = grandImaginaryMatrixVec_[omegaIdx];
      wrapper.eps_zz_Inv = eps_zz_Inv_MatrixVec_[omegaIdx];
      switch (options_.IntegralMethod) {
        case GAUSSLEGENDRE_:{
          recvBuf[i] = gauss_legendre(degree_, wrapperFunQuadgl, &wrapper, kxStart_, kxEnd_);
          break;
        }
        case GAUSSKRONROD_:{
          // this part, debug is needed
          Workspace<double> integrationWork(degree_);
          Function<double, ArgWrapper> F(wrapperFunQuadgk, &wrapper);
          double epsabs = integrationWork.get_eps();
          double abserr, epsrel = 0;
          integrationWork.qag(F, kxStart_, kxEnd_, epsabs, epsrel, recvBuf[i], abserr);
          break;
        }
        default:{
          break;
        }
      }
      recvBuf[i] *= POW3(omegaList_[omegaIdx] / datum::c_0) / POW2(datum::pi);
    }

    // gatther the Phi from each processor
    MPI_Gatherv(&recvBuf[0], 1, dtype, &Phi_[0], sendCounts, displs, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == MASTER){
        this->saveToFile();
    }
    delete[] recvBuf;
    recvBuf = nullptr;

    MPI_Finalize();
  }
  /*==============================================*/
  // Implementaion of the class on 1D grating simulation
  /*==============================================*/
  SimulationGrating::SimulationGrating() : Simulation(){
    prefactor_ = 1;
    dim_ = ONE_;
  }
  /*==============================================*/
  // This is a thin wrapper for the usage of smart pointer
  /*==============================================*/
  Ptr<SimulationGrating> SimulationGrating::instanceNew(){
    return new SimulationGrating();
  }
  /*==============================================*/
  // Function setting the integral over kx
  // @args:
  // points: number of points of sampling kx
  /*==============================================*/
  void SimulationGrating::setKxIntegral(const int points){
    numOfKx_ = points;
    if(period_[0] == 0.0){
      throw UTILITY::ValueException("Periodicity not set!");
    }
    kxStart_ = -datum::pi / period_[0];
    kxEnd_ = -kxStart_;
  }

  /*==============================================*/
  // This function set the integral of kx when the system is symmetric in x direction
  // @args:
  // points: number of kx points
  /*==============================================*/
  void SimulationGrating::setKxIntegralSym(const int points){
    numOfKx_ = points;
    kxStart_ = 0;
    if(period_[0] == 0.0){
      throw UTILITY::ValueException("Periodicity not set!");
    }
    kxEnd_ = datum::pi / period_[0];
    prefactor_ *= 2;
  }

  /*==============================================*/
  // This function set the integral of ky
  // @args:
  // points: number of ky points
  // end: the upperbound of the integral
  /*==============================================*/
  void SimulationGrating::setKyIntegral(const int points, const double end){
    kyStart_ = -end;
    numOfKy_ = points;
    kyEnd_ = end;
  }
  /*==============================================*/
  // This function set the integral of ky assume y symmetry
  // @args:
  // points: number of ky points
  // end: the upperbound of the integral
  /*==============================================*/
  void SimulationGrating::setKyIntegralSym(const int points, const double end){
    kyStart_ = 0;
    numOfKy_ = points;
    kyEnd_ = end;
    prefactor_ *= 2;
  }
  /*==============================================*/
  // Implementaion of the class on 2D patterning simulation
  /*==============================================*/
  SimulationPattern::SimulationPattern() : Simulation(){
    prefactor_ = 1;
    dim_ = TWO_;
  }
  /*==============================================*/
  // This is a thin wrapper for the usage of smart pointer
  /*==============================================*/
  Ptr<SimulationPattern> SimulationPattern::instanceNew(){
    return new SimulationPattern();
  }

  /*==============================================*/
  // Function setting the integral over kx
  // @args:
  // points: number of points of sampling kx
  /*==============================================*/
  void SimulationPattern::setKxIntegral(const int points){
    if(period_[0] == 0.0){
      throw UTILITY::ValueException("Periodicity not set!");
    }
    kxStart_ = -datum::pi / period_[0];
    numOfKx_ = points;
    kxEnd_ = -kxStart_;
  }

  /*==============================================*/
  // This function set the integral of kx when the system is symmetric in x direction
  // @args:
  // points: number of kx points
  /*==============================================*/
  void SimulationPattern::setKxIntegralSym(const int points){
    kxStart_ = 0;
    numOfKx_ = points;
    if(period_[0] == 0.0){
      throw UTILITY::ValueException("Periodicity not set!");
    }
    kxEnd_ = datum::pi / period_[0];
    prefactor_ *= 2;
  }

  /*==============================================*/
  // This function set the integral of ky
  // @args:
  // points: number of ky points
  // end: the upperbound of the integral
  /*==============================================*/
  void SimulationPattern::setKyIntegral(const int points){
    if(period_[1] == 0.0){
      throw UTILITY::ValueException("Periodicity not set!");
    }
    kyStart_ = -datum::pi / period_[1];
    numOfKy_ = points;
    kyEnd_ = -kyStart_;
  }

  /*==============================================*/
  // This function set the integral of ky when the system is symmetric in y direction
  // @args:
  // points: number of ky points
  /*==============================================*/
  void SimulationPattern::setKyIntegralSym(const int points){
    kyStart_ = 0;
    numOfKy_ = points;
    if(period_[1] == 0.0){
      throw UTILITY::ValueException("Periodicity not set!");
    }
    kyEnd_ = datum::pi / period_[1];
    prefactor_ *= 2;
  }

}