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
  FileLoader::FileLoader(){
    omegaList_ = nullptr;
    epsilonList_.epsilonVals = nullptr;
  }
  /*==============================================*/
  // This is a thin wrapper for the usage of smart pointer
  /*==============================================*/
  Ptr<FileLoader> FileLoader::instanceNew(){
    return new FileLoader();
  }

  /*==============================================*/
  // Function checking a line contains more than 3 spaces
  /*==============================================*/
  EPSTYPE checkType(std::string line){
    int numOfSpace = 0;
    for(size_t i = 0; i < line.length(); i++){
      if(isspace(line.at(i))) numOfSpace++;
    }
    switch (numOfSpace) {
      case 2: return SCALAR_;
      case 6: return DIAGONAL_;
      case 10: return TENSOR_;
      default:{
        std::cerr << "Input type wrong!" << std::endl;
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
      std::cerr << fileName + " not exists!" << std::endl;
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
      numOfOmega_ = count;
      preSet_ = true;
      omegaList_ = new double[numOfOmega_];
      epsilonList_.epsilonVals = new EpsilonVal[numOfOmega_];
    }
    else{
      if(numOfOmega_ != count){
        std::cerr << "wrong omega points!" << std::endl;
        throw UTILITY::RangeException("wrong omega points!");
      }
    }

    epsilonList_.type_ = type;
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
    if(omegaList_ != nullptr){
      delete [] omegaList_;
      omegaList_ = nullptr;
    }
    if(epsilonList_.epsilonVals != nullptr){
      delete [] epsilonList_.epsilonVals;
      epsilonList_.epsilonVals = nullptr;
    }
  }
  /*==============================================*/
  // This function wraps the data for quad_gaussian_kronrod
  // @args:
  // kx: the kx value (normalized)
  // wrapper: wrapper for all the arguments wrapped in wrapper
  /*==============================================*/
  static void wrapperFunQuadgk(unsigned ndim,
    const double *kx,
    void* data,
    unsigned fdim,
    double *fval
    ){
    ArgWrapper wrapper = *(ArgWrapper*)data;
    fval[0] = kx[0] * poyntingFlux(
      wrapper.omega,
      wrapper.thicknessList,
      kx[0],
      0,
      wrapper.EMatrices,
      wrapper.grandImaginaryMatrices,
      wrapper.eps_zz_Inv,
      wrapper.Gx_mat,
      wrapper.Gy_mat,
      wrapper.sourceList,
      wrapper.targetLayer,
      1,
      wrapper.polar
    );
  }
  /*==============================================*/
  // This function wraps the data for quad_legendre
  // @args:
  // kx: the kx value (normalized)
  // data: wrapper for all the arguments wrapped in wrapper
  /*==============================================*/
  static double wrapperFunQuadgl(const double kx, void* data){
    ArgWrapper wrapper = *(ArgWrapper*) data;
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
      1,
      wrapper.polar
    );
  }
  /*======================================================*/
  // Implementaion of the parent simulation super class
  /*=======================================================*/
  Simulation::Simulation() : nGx_(0), nGy_(0), numOfOmega_(0), Phi_(nullptr), omegaList_(nullptr),
   kxStart_(0), kxEnd_(0), kyStart_(0), kyEnd_(0), numOfKx_(0), numOfKy_(0)
  {
    targetLayer_ = -1;
    dim_ = NO_;
    structure_ = Structure::instanceNew();
    fileLoader_ = FileLoader::instanceNew();
    resultArray_ = nullptr;
  }
  /*==============================================*/
  // Class destructor
  /*==============================================*/
  Simulation::~Simulation(){
    if(Phi_ != nullptr){
      delete[] Phi_;
      Phi_ = nullptr;
    }
    if(resultArray_ != nullptr){
      delete[] resultArray_;
      resultArray_ = nullptr;
    }
  }


  /*==============================================*/
  // This function return the Phi value
  /*==============================================*/
  double* Simulation::getPhi(){
    if(resultArray_ == nullptr){
      std::cerr << "Integration has not done yet!" << std::endl;
      throw UTILITY::MemoryException("Integration has not done yet!");
    }

    if(Phi_ == nullptr){
      std::cerr << "Please do RCWA and integration first!" << std::endl;
      throw UTILITY::MemoryException("Please do RCWA and integration first!");
    }

    double* kxList = new double[numOfKx_];
    double* kyList = new double[numOfKy_];
    double* scalex = new double[numOfOmega_];
    double* scaley = new double[numOfOmega_];
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
          if(options_.kxIntegralPreset) scalex[i] = 1;
          else scalex[i] = omegaList_[i] / datum::c_0;
          scaley[i] = 1;
          break;
        }
        case TWO_:{
          if(options_.kxIntegralPreset) scalex[i] = 1;
          else scalex[i] = omegaList_[i] / datum::c_0;
          if(options_.kyIntegralPreset) scaley[i] = 1;
          else scaley[i] = omegaList_[i] / datum::c_0;
          break;
        }
        default: break;
      }
    }

    for(int i = 0; i < numOfOmega_; i++){
        Phi_[i] = 0;
        for(int j = 0; j < numOfKx_ * numOfKy_; j++){
          Phi_[i] += resultArray_[i * numOfKx_ * numOfKy_ + j];
        }
        Phi_[i] *= prefactor_ * dkx / scalex[i] * dky / scaley[i] * POW2(omegaList_[i] / datum::c_0);
    }

    delete[] kxList;
    kxList = nullptr;
    delete[] kyList;
    kyList = nullptr;
    delete[] scalex;
    scalex = nullptr;
    delete[] scaley;
    scaley = nullptr;
    return Phi_;
  }

  /*==============================================*/
  // This function return the omega value
  /*==============================================*/
  double* Simulation::getOmega(){
    if(omegaList_ == nullptr){
      std::cerr << "omega does not exist!" << std::endl;
      throw UTILITY::MemoryException("omega does not exist!");
    }
    return omegaList_;
  }
  /*==============================================*/
  // This function return reconstructed dielectric at a given point
  /*==============================================*/
  void Simulation::getEpsilon(const int omegaIndex, const double position[3], double* epsilon){
    // TODO
  }
  /*==============================================*/
  // This function return the number of omega
  /*==============================================*/
  int Simulation::getNumOfOmega(){
    return numOfOmega_;
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
    structure_->setPeriodicity(p1, p2);
  }
  /*==============================================*/
  // This function adds material to the system
  // @args:
  // name: the name of the material
  // infile: the input file containing the properties of the material
  /*==============================================*/
  void Simulation::addMaterial(const std::string name, const std::string infile){
    if(materialInstanceMap_.find(name) != materialInstanceMap_.cend()){
      std::cerr << name + ": Material already exist!" << std::endl;
      throw UTILITY::NameInUseException(name + ": Material already exist!");
      return;
    }

    fileLoader_->load(infile);
    Ptr<Material> material = Material::instanceNew(name,
     fileLoader_->getOmegaList(),
     fileLoader_->getEpsilonList(),
     fileLoader_->getNumOfOmega()
    );
    materialInstanceMap_.insert(MaterialMap::value_type(name, material));
    structure_->addMaterial(material);
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
  void Simulation::setMaterial(const std::string name, double** epsilon, const std::string type){
    if(materialInstanceMap_.find(name) == materialInstanceMap_.cend()){
      std::cerr << name + ": Material does not exist!" << std::endl;
      throw UTILITY::IllegalNameException(name + ": Material does not exist!");
      return;
    }
    Ptr<Material> material = materialInstanceMap_.find(name)->second;
    EPSTYPE originalType = material->getType();
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
      std::cerr << "Please choose 'type' from 'scalar', 'diagonal' or 'tensor'!" << std::endl;
      throw UTILITY::AttributeNotSupportedException("Please choose 'type' from 'scalar', 'diagonal' or 'tensor'!");
    }

    material->setEpsilon(newEpsilon, numOfOmega);
    delete [] newEpsilon.epsilonVals;
    newEpsilon.epsilonVals = nullptr;

    for(const_LayerInstanceIter it = layerInstanceMap_.cbegin(); it != layerInstanceMap_.cend(); it++){
      Ptr<Layer> layer = it->second;
      if(layer->hasMaterial(material)){

        if(originalType != TENSOR_ && type == "tensor"){
          layer->containTensor(true);
        }
        if(originalType == TENSOR_ && type != "tensor"){
          layer->containTensor(false);
        }

      }
    }

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
      std::cerr << materialName + ": Material does not exist!" << std::endl;
      throw UTILITY::IllegalNameException(materialName + ": Material does not exist!");
      return;
    }
    if(layerInstanceMap_.find(name) != layerInstanceMap_.cend()){
      std::cerr << name + ": Layer already exists!" << std::endl;
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
      std::cerr << materialName + ": Material does not exist!" << std::endl;
      throw UTILITY::IllegalNameException(materialName + ": Material does not exist!");
      return;
    }
    if(layerInstanceMap_.find(name) == layerInstanceMap_.cend()){
      std::cerr << name + ": Layer does not exist!" << std::endl;
      throw UTILITY::IllegalNameException(name + ": Layer does not exist!");
      return;
    }
    Ptr<Material> material = materialInstanceMap_.find(materialName)->second;
    Ptr<Layer> layer = layerInstanceMap_.find(name)->second;
    if(material->getType() == TENSOR_ && !layer->hasMaterial(material)){
      layer->containTensor(true);
    }
    if((layer->getBackGround())->getType() == TENSOR_ && material->getType() != TENSOR_){
      layer->containTensor(false);
    }
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
      std::cerr << name + ": Layer does not exist!" << std::endl;
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
      std::cerr << originalName + ": Layer does not exist!" << std::endl;
      throw UTILITY::IllegalNameException(originalName + ": Layer does not exist!");
      return;
    }
    if(layerInstanceMap_.find(name) != layerInstanceMap_.cend()){
      std::cerr << name + ": cannot add a layer that already exists!" << std::endl;
      throw UTILITY::NameInUseException(name + ": cannot add a layer that already exists!");
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
      std::cerr << name + ": Layer does not exist!" << std::endl;
      throw UTILITY::NameInUseException(name + ": Layer does not exist!");
      return;
    }
    layerInstanceMap_.erase(name);
    structure_->deleteLayerByName(name);
  }
  /*==============================================*/
  // This function sets a layer as the source
  // @args:
  // name: the name of the source layer
  /*==============================================*/
  void Simulation::setSourceLayer(const std::string name){
    if(layerInstanceMap_.find(name) == layerInstanceMap_.cend()){
      std::cerr << name + ": Layer does not exist!" << std::endl;
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
      std::cerr << name + ": Layer does not exist!" << std::endl;
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
      std::cerr << std::to_string(omegaIdx) + ": out of range!" << std::endl;
      throw UTILITY::RangeException(std::to_string(omegaIdx) + ": out of range!");
    }
    int N = getN(nGx_, nGy_);
    return omegaList_[omegaIdx] / datum::c_0 / POW3(datum::pi) / 2.0 *
      poyntingFlux(omegaList_[omegaIdx] / datum::c_0,
        thicknessListVec_,
        kx,
        ky,
        EMatricesVec_[omegaIdx],
        grandImaginaryMatrixVec_[omegaIdx],
        eps_zz_Inv_MatrixVec_[omegaIdx],
        Gx_mat_,
        Gy_mat_,
        sourceList_,
        targetLayer_,
        N,
        options_.polarization
      );
  }

  /*==============================================*/
  // This function builds up the matrices
  /*==============================================*/
  void Simulation::buildRCWA(){
    // reset simulation first
    this->resetSimulation();
    // essential, get the shared Gx_mat_ and Gy_mat_
    getGMatrices(nGx_, nGy_, period_, Gx_mat_, Gy_mat_, dim_);
    // get constants
    Ptr<Layer> firstLayer = structure_->getLayerByIndex(0);
    Ptr<Material> backGround = firstLayer->getBackGround();
    numOfOmega_ = backGround->getNumOfOmega();
    omegaList_ = backGround->getOmegaList();
    int numOfLayer = structure_->getNumOfLayer();

    thicknessListVec_ = zeros<RCWArVector>(numOfLayer);
    sourceList_.resize(numOfLayer);
    for(int i = 0; i < numOfLayer; i++){
      thicknessListVec_(i) = (structure_->getLayerByIndex(i))->getThickness();
      sourceList_[i] = (structure_->getLayerByIndex(i))->checkIsSource();
      if(sourceList_[i] && i >= targetLayer_){
        std::cerr << "Target Layer needs to be above source layer!" << std::endl;
        throw UTILITY::RangeException("Target Layer needs to be above source layer!");
      }
    }
    // set the first and last layer to have 0 thickness
    thicknessListVec_(0) = 0;
    thicknessListVec_(numOfLayer - 1) = 0;

    if(dim_ != NO_ && period_[0] == 0.0){
      std::cerr << "Periodicity not set!" << std::endl;
      throw UTILITY::ValueException("Periodicity not set!");
    }
    if(dim_ == TWO_ && period_[1] == 0.0){
      std::cerr << "Periodicity not set!" << std::endl;
      throw UTILITY::ValueException("Periodicity not set!");
    }

    // initializing for the output
    Phi_ = new double[numOfOmega_];
    if(dim_ != NO_ && (numOfKx_ == 0 || numOfKy_ == 0)){
      std::cerr << "Set integration range first!" << std::endl;
      throw UTILITY::ValueException("Set integration range first!");
    }
    resultArray_ = new double[numOfKx_ * numOfKy_ * numOfOmega_];
    for(int i = 0; i < numOfKx_ * numOfKy_ * numOfOmega_; i++){
      resultArray_[i] = 0;
    }

    EMatricesVec_.resize(numOfOmega_);
    grandImaginaryMatrixVec_.resize(numOfOmega_);
    eps_zz_Inv_MatrixVec_.resize(numOfOmega_);

    RCWAcMatricesVec eps_xx_MatrixVec(numOfOmega_), eps_xy_MatrixVec(numOfOmega_), eps_yx_MatrixVec(numOfOmega_), eps_yy_MatrixVec(numOfOmega_);
    RCWAcMatricesVec im_eps_xx_MatrixVec(numOfOmega_), im_eps_xy_MatrixVec(numOfOmega_), im_eps_yx_MatrixVec(numOfOmega_), im_eps_yy_MatrixVec(numOfOmega_), im_eps_zz_MatrixVec(numOfOmega_);

    int N = getN(nGx_, nGy_);
    RCWAcMatrix onePadding1N = eye<RCWAcMatrix>(N, N);

    for(int i = 0; i < numOfLayer; i++){
      Ptr<Layer> layer = structure_->getLayerByIndex(i);
      Ptr<Material> backGround = layer->getBackGround();

      for(int j = 0; j < numOfOmega_; j++){
        RCWAcMatrix eps_xx(N, N, fill::zeros), eps_xy(N, N, fill::zeros), eps_yx(N, N, fill::zeros), eps_yy(N, N, fill::zeros), eps_zz_Inv(N, N, fill::zeros);
        RCWAcMatrix im_eps_xx(N, N, fill::zeros), im_eps_xy(N, N, fill::zeros), im_eps_yx(N, N, fill::zeros), im_eps_yy(N, N, fill::zeros), im_eps_zz(N, N, fill::zeros);

        EpsilonVal epsBG = backGround->getEpsilonAtIndex(j);
        EpsilonVal epsBGTensor = FMM::toTensor(epsBG, backGround->getType());

        const_MaterialIter m_it = layer->getMaterialsBegin();
        int count = 0;
        for(const_PatternIter it = layer->getPatternsBegin(); it != layer->getPatternsEnd(); it++){
          Pattern pattern = *it;
          Ptr<Material> material = *(m_it + count);
          EpsilonVal epsilon = material->getEpsilonAtIndex(j);
          count++;

          switch(pattern.type_){
            /*************************************/
            // if the pattern is a grating (1D)
            /************************************/
            case GRATING_:{
              double center = pattern.arg1_.first;
              double width = pattern.arg1_.second;
              FMM::transformGrating(
                eps_xx,
                eps_xy,
                eps_yx,
                eps_yy,
                eps_zz_Inv,
                im_eps_xx,
                im_eps_xy,
                im_eps_yx,
                im_eps_yy,
                im_eps_zz,
                epsBGTensor,
                epsilon,
                material->getType(),
                nGx_,
                center,
                width,
                period_[0],
                layer->hasTensor(),
                options_.FMMRule == INVERSERULE_
              );
              break;
            }

            /*************************************/
            // if the pattern is a rectangle (2D)
            /************************************/
            case RECTANGLE_:{
              double centers[2] = {pattern.arg1_.first, pattern.arg1_.second};
              double widths[2] = {pattern.arg2_.first, pattern.arg2_.second};
              FMM::transformRectangle(
                eps_xx,
                eps_xy,
                eps_yx,
                eps_yy,
                eps_zz_Inv,
                im_eps_xx,
                im_eps_xy,
                im_eps_yx,
                im_eps_yy,
                im_eps_zz,
                epsBGTensor,
                epsilon,
                material->getType(),
                nGx_,
                nGy_,
                centers,
                widths,
                period_,
                layer->hasTensor(),
                options_.FMMRule == INVERSERULE_
              );
              break;
            }
            /*************************************/
            // if the pattern is a circle (2D)
            /************************************/
            case CIRCLE_:{
              double centers[2] = {pattern.arg1_.first, pattern.arg2_.first};
              double radius = pattern.arg1_.second;
              FMM::transformCircle(
                eps_xx,
                eps_xy,
                eps_yx,
                eps_yy,
                eps_zz_Inv,
                im_eps_xx,
                im_eps_xy,
                im_eps_yx,
                im_eps_yy,
                im_eps_zz,
                epsBGTensor,
                epsilon,
                material->getType(),
                nGx_,
                nGy_,
                centers,
                radius,
                period_,
                layer->hasTensor(),
                options_.FMMRule == INVERSERULE_
              );
              break;
            }
            default: break;
          }
        }
        /*************************************/
        // collection information from the background
        /************************************/
        eps_xx += dcomplex(epsBGTensor.tensor[0], epsBGTensor.tensor[1]) * onePadding1N;
        im_eps_xx += epsBGTensor.tensor[1] * onePadding1N;

        eps_yy += dcomplex(epsBGTensor.tensor[6], epsBGTensor.tensor[7]) * onePadding1N;
        im_eps_yy += epsBGTensor.tensor[7] * onePadding1N;

        eps_zz_Inv += dcomplex(epsBGTensor.tensor[8], epsBGTensor.tensor[9]) * onePadding1N;
        eps_zz_Inv = eps_zz_Inv.i();

        im_eps_zz += epsBGTensor.tensor[9] * onePadding1N;

        if(layer->hasTensor()){
          eps_xy += dcomplex(epsBGTensor.tensor[2], epsBGTensor.tensor[3]) * onePadding1N;
          eps_yx += dcomplex(epsBGTensor.tensor[4], epsBGTensor.tensor[5]) * onePadding1N;
          im_eps_xy += epsBGTensor.tensor[3] * onePadding1N;
          im_eps_yx += epsBGTensor.tensor[5] * onePadding1N;
        }

        eps_xx_MatrixVec[j].push_back(eps_xx);
        eps_xy_MatrixVec[j].push_back(eps_xy);
        eps_yx_MatrixVec[j].push_back(eps_yx);
        eps_yy_MatrixVec[j].push_back(eps_yy);
        eps_zz_Inv_MatrixVec_[j].push_back(eps_zz_Inv);
        im_eps_xx_MatrixVec[j].push_back(im_eps_xx);
        im_eps_xy_MatrixVec[j].push_back(im_eps_xy);
        im_eps_yx_MatrixVec[j].push_back(im_eps_yx);
        im_eps_yy_MatrixVec[j].push_back(im_eps_yy);
        im_eps_zz_MatrixVec[j].push_back(im_eps_zz);
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

  }
  /*==============================================*/
  // This function prints out the information of the system
  /*==============================================*/
  void Simulation::outputSysInfo(){
    std::cout << "==================================================" << std::endl;
    std::cout << "The system has in total " << structure_->getNumOfLayer() << " layers." << std::endl;
    if(dim_ == ONE_){
      std::cout << "Periodicity in x direction is " << period_[0] << std::endl;
    }
    else if(dim_ == TWO_){
      std::cout << "Periodicity in x, y directions are " << period_[0] << ", " << period_[1] << std::endl;
    }
    std::cout << "==================================================" << std::endl;
    std::cout << "Printing from bottom to up." << std::endl;
    std::cout << "==================================================" << std::endl;
    for(const_LayerIter it = structure_->getLayersBegin(); it != structure_->getLayersEnd(); it++){
      Ptr<Layer> layer = it->second;
      std::cout << "Layer index " << it->first << ": " << layer->getName() << std::endl;
      std::cout << "Thickness: " << layer->getThickness() << std::endl;

      std::cout << "contains tensor: ";
      if(layer->hasTensor()) std::cout << "YES" << std::endl;
      else std::cout << "NO" << std::endl;

      std::cout << "Is source: ";
      if(layer->checkIsSource()) std::cout << "YES" << std::endl;
      else std::cout << "NO" << std::endl;

      std::cout << "Its background is: " << layer->getBackGround()->getName() << std::endl;
      if(layer->getNumOfMaterial() != 0){
        std::cout << "It has other components:" << std::endl << std::endl;
        int count = 0;
        const_MaterialIter m_it = layer->getMaterialsBegin();
        for(const_PatternIter it = layer->getPatternsBegin(); it != layer->getPatternsEnd(); it++){
          std::cout << "Material for pattern " << count + 1 << ": " << (*(m_it + count))->getName() << std::endl;
          std::cout << "Pattern " << count + 1 << " is: ";
          switch((*it).type_){
            case GRATING_:{
              std::cout << "grating, ";
              std::cout << "(c, w) = (" << (*it).arg1_.first << ", " << (*it).arg1_.second << ")\n";
              break;
            }
            case RECTANGLE_:{
              std::cout << "rectangle, ";
              std::cout << "(c_x, w_x) = (" << (*it).arg1_.first << ", " << (*it).arg2_.first <<"), ";
              std::cout << "(c_y, w_y) = (" << (*it).arg1_.second << ", " << (*it).arg2_.second <<")\n";
              break;
            }
            case CIRCLE_:{
              std::cout << "circle";
              std::cout << "(c, w) = (" << (*it).arg1_.first << ", " << (*it).arg2_.first << "), ";
              std::cout << "r = " << (*it).arg1_.second << std::endl;
              break;
            }
            default: break;
          }
          count++;
        }
      }
      std::cout << "==================================================" << std::endl;
    }
  }
  /*==============================================*/
  // function using naive implementaion
  /*==============================================*/
  void Simulation::optUseNaiveRule(){
    options_.FMMRule = NAIVEFMM_;
  }
  /*==============================================*/
  // function using inverse rule
  /*==============================================*/
  void Simulation::optUseInverseRule(){
    options_.FMMRule = INVERSERULE_;
  }
  /*==============================================*/
  // function print intermediate results
  /*==============================================*/
  void Simulation::optPrintIntermediate(){
    options_.PrintIntermediate = true;
  }
  /*==============================================*/
  // function sets that only TE mode is computed
  /*==============================================*/
  void Simulation::optOnlyComputeTE(){
    options_.polarization = TE_;
  }
  /*==============================================*/
  // function sets that only TM mode is computed
  /*==============================================*/
  void Simulation::optOnlyComputeTM(){
    options_.polarization = TM_;
  }

  /*==============================================*/
  // function print intermediate results
  /*==============================================*/
  void Simulation::setThread(const int thread){
    if(thread <= 0){
      std::cerr << "Number of thread should >= 1!" << std::endl;
      throw UTILITY::RangeException("Number of thread should >= 1!");
    }
    #if defined(_OPENMP)
      numOfThread_ = std::min(thread, omp_get_max_threads());
    #endif
  }
  /*==============================================*/
  // Function setting the integral over kx
  // @args:
  // points: number of points of sampling kx
  // end: the upperbound of the integral
  /*==============================================*/
  void Simulation::setKxIntegral(const int points, const double end){
    if(points < 2){
      std::cerr << "Needs no less than 2 points!" << std::endl;
      throw UTILITY::ValueException("Needs no less than 2 points!");
    }
    numOfKx_ = points;
    if(dim_ != NO_ && period_[0] == 0.0){
      std::cerr << "Periodicity not set!" << std::endl;
      throw UTILITY::ValueException("Periodicity not set!");
    }
    if(dim_ == NO_ && end == 0.0){
      std::cerr << "integral upper bound cannot be zero!" << std::endl;
      throw UTILITY::ValueException("integral upper bound cannot be zero!");
    }
    if(end != 0){
      kxEnd_ = end;
      options_.kxIntegralPreset = true;
    }
    else{
      kxEnd_ = datum::pi / period_[0];
      options_.kxIntegralPreset = false;
    }
    kxStart_ = -kxEnd_;
  }

  /*==============================================*/
  // This function set the integral of kx when the system is symmetric in x direction
  // @args:
  // points: number of kx points
  // end: the upperbound of the integral
  /*==============================================*/
  void Simulation::setKxIntegralSym(const int points, const double end){
    if(points < 2){
      std::cerr << "Needs no less than 2 points!" << std::endl;
      throw UTILITY::ValueException("Needs no less than 2 points!");
    }
    numOfKx_ = points;
    if(dim_ != NO_ && period_[0] == 0.0){
      std::cerr << "Periodicity not set!" << std::endl;
      throw UTILITY::ValueException("Periodicity not set!");
    }
    if(dim_ == NO_ && end == 0.0){
      std::cerr << "integral upper bound cannot be zero!" << std::endl;
      throw UTILITY::ValueException("integral upper bound cannot be zero!");
    }
    if(end != 0){
      kxEnd_ = end;
      options_.kxIntegralPreset = true;
    }
    else{
      kxEnd_ = datum::pi / period_[0];
      options_.kxIntegralPreset = false;
    }
    kxStart_ = 0;
    prefactor_ *= 2;
  }

  /*==============================================*/
  // This function set the integral of ky
  // @args:
  // points: number of ky points
  // end: the upperbound of the integral
  /*==============================================*/
  void Simulation::setKyIntegral(const int points, const double end){
    if(points < 2){
      std::cerr << "Needs no less than 2 points!" << std::endl;
      throw UTILITY::ValueException("Needs no less than 2 points!");
    }
    numOfKy_ = points;
    if(dim_ == TWO_ && period_[1] == 0.0){
      std::cerr << "Periodicity not set!" << std::endl;
      throw UTILITY::ValueException("Periodicity not set!");
    }
    if((dim_ == NO_ || dim_ == ONE_) && end == 0.0){
      std::cerr << "integral upper bound cannot be zero!" << std::endl;
      throw UTILITY::ValueException("integral upper bound cannot be zero!");
    }
    if(end != 0){
      kyEnd_ = end;
      options_.kyIntegralPreset = true;
    }
    else{
      kyEnd_ = datum::pi / period_[1];
      options_.kyIntegralPreset = false;
    }
    kxStart_ = -kxEnd_;
  }
  /*==============================================*/
  // This function set the integral of ky assume y symmetry
  // @args:
  // points: number of ky points
  // end: the upperbound of the integral
  /*==============================================*/
  void Simulation::setKyIntegralSym(const int points, const double end){
    if(points < 2){
      std::cerr << "Needs no less than 2 points!" << std::endl;
      throw UTILITY::ValueException("Needs no less than 2 points!");
    }
    numOfKy_ = points;
    if(dim_ == TWO_ && period_[1] == 0.0){
      std::cerr << "Periodicity not set!" << std::endl;
      throw UTILITY::ValueException("Periodicity not set!");
    }
    if((dim_ == NO_ || dim_ == ONE_) && end == 0.0){
      std::cerr << "integral upper bound cannot be zero!" << std::endl;
      throw UTILITY::ValueException("integral upper bound cannot be zero!");
    }
    if(end != 0){
      kyEnd_ = end;
      options_.kyIntegralPreset = true;
    }
    else{
      kyEnd_ = datum::pi / period_[1];
      options_.kyIntegralPreset = false;
    }
    kyStart_ = 0;
    prefactor_ *= 2;
  }
  /*==============================================*/
  // This function computes the flux for internal usage
  // @args:
  // start: the starting index
  // end: the end index
  /*==============================================*/
  void Simulation::integrateKxKyInternal(const int start, const int end, const bool parallel){
    double* kxList = new double[numOfKx_];
    double* kyList = new double[numOfKy_];
    double* scalex = new double[numOfOmega_];
    double* scaley = new double[numOfOmega_];
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
          if(options_.kxIntegralPreset) scalex[i] = 1;
          else scalex[i] = omegaList_[i] / datum::c_0;
          scaley[i] = 1;
          break;
        }
        case TWO_:{
          if(options_.kxIntegralPreset) scalex[i] = 1;
          else scalex[i] = omegaList_[i] / datum::c_0;
          if(options_.kyIntegralPreset) scaley[i] = 1;
          else scaley[i] = omegaList_[i] / datum::c_0;
          break;
        }
        default: break;
      }
    }

    if(parallel){
      #if defined(_OPENMP)
        #pragma omp parallel for num_threads(numOfThread_)
      #endif
      for(int i = start; i < end; i++){
          int omegaIdx = i / (numOfKx_ * numOfKy_);
          int residue = i % (numOfKx_ * numOfKy_);
          int kxIdx = residue / numOfKy_;
          int kyIdx = residue % numOfKy_;
          resultArray_[i] = this->getPhiAtKxKy(omegaIdx, kxList[kxIdx] / scalex[omegaIdx], kyList[kyIdx] / scaley[omegaIdx]);
          if(options_.PrintIntermediate){
            std::stringstream msg;
            msg << omegaList_[omegaIdx] << "\t" << kxList[kxIdx] / scalex[omegaIdx] << "\t" << kyList[kyIdx]  / scaley[omegaIdx] << "\t" << resultArray_[i] << std::endl;
            std::cout << msg.str();
          }
      }

    }
    else{
      for(int i = start; i < end; i++){
          int omegaIdx = i / (numOfKx_ * numOfKy_);
          int residue = i % (numOfKx_ * numOfKy_);
          int kxIdx = residue / numOfKy_;
          int kyIdx = residue % numOfKy_;
          resultArray_[i] = this->getPhiAtKxKy(omegaIdx, kxList[kxIdx] / scalex[omegaIdx], kyList[kyIdx] / scaley[omegaIdx]);
          if(options_.PrintIntermediate){
            std::stringstream msg;
            msg << omegaList_[omegaIdx] << "\t" << kxList[kxIdx] / scalex[omegaIdx] << "\t" << kyList[kyIdx]  / scaley[omegaIdx] << "\t" << resultArray_[i] << std::endl;
            std::cout << msg.str();
          }
      }
    }

    delete[] kxList;
    kxList = nullptr;
    delete[] kyList;
    kyList = nullptr;
    delete[] scalex;
    scalex = nullptr;
    delete[] scaley;
    scaley = nullptr;
  }

  /*==============================================*/
  // This function computes the flux for general usage
  /*==============================================*/
  void Simulation::integrateKxKy(){
    this->integrateKxKyInternal(0, numOfOmega_ * numOfKx_ * numOfKy_, true);
  }
  /*==============================================*/
  // This function computes the flux for MPI only
  /*==============================================*/
  void Simulation::integrateKxKyMPI(const int rank, const int size){
    int totalNum = numOfOmega_ * numOfKx_ * numOfKy_;
    int chunksize = totalNum / size;
    int left = totalNum % size;
    int start, end;
    if(rank >= left){
      start = left * (chunksize + 1) + (rank - left) * chunksize;
      end = start + chunksize;
    }
    else{
      start = rank * (chunksize + 1);
      end = start + chunksize + 1;
    }
    if(end > totalNum) end = totalNum;
    this->integrateKxKyInternal(start, end, false);

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
    options_.IntegrateKParallel = true;
  }
  /*==============================================*/
  // Function setting the integral to be the gauss_legendre
  // @args:
  // degree: the degree of gauss_legendre integral, default 1024
  /*==============================================*/
  void SimulationPlanar::optUseQuadgl(const int degree){
    degree_ = degree;
    options_.IntegralMethod = GAUSSLEGENDRE_;
  }
  /*==============================================*/
  // Function setting the integral to be the gauss_kronrod
  /*==============================================*/
  void SimulationPlanar::optUseQuadgk(){
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
    if(options_.IntegrateKParallel == false){
      std::cerr << "Cannot use kparallel integral here!" << std::endl;
      throw UTILITY::InternalException("Cannot use kparallel integral here!");
    }
    if(omegaIdx >= numOfOmega_){
      std::cerr << std::to_string(omegaIdx) + ": out of range!" << std::endl;
      throw UTILITY::RangeException(std::to_string(omegaIdx) + ": out of range!");
    }
    return POW2(omegaList_[omegaIdx] / datum::c_0) / POW2(datum::pi) * KParallel *
      poyntingFlux(omegaList_[omegaIdx] / datum::c_0, thicknessListVec_, KParallel, 0, EMatricesVec_[omegaIdx],
      grandImaginaryMatrixVec_[omegaIdx], eps_zz_Inv_MatrixVec_[omegaIdx], Gx_mat_, Gy_mat_,
      sourceList_, targetLayer_,1, options_.polarization);
  }


  /*==============================================*/
  // This function integrates kx assuming scalar or
  //  eps_x,  0,   0
  //    0  ,eps_x, 0
  //    0     0,  eps_z
  // make sure you understand your problem whether can be solved by this function
  /*==============================================*/
  void SimulationPlanar::integrateKParallel(){
    if(options_.IntegrateKParallel == false){
      std::cerr << "Cannot use kparallel integral here!" << std::endl;
      throw UTILITY::InternalException("Cannot use kparallel integral here!");
    }
    if(Phi_ == nullptr){
      std::cerr << "Please do RCWA and integration first!" << std::endl;
      throw UTILITY::MemoryException("Please do RCWA and integration first!");
    }

    #if defined(_OPENMP)
      #pragma omp parallel for num_threads(numOfThread_)
    #endif
    for(int i = 0; i < numOfOmega_; i++){
      ArgWrapper wrapper;
      wrapper.thicknessList = thicknessListVec_;
      wrapper.Gx_mat = Gx_mat_;
      wrapper.Gy_mat = Gy_mat_;
      wrapper.sourceList = sourceList_;
      wrapper.targetLayer = targetLayer_;
      wrapper.omega = omegaList_[i] / datum::c_0;
      wrapper.EMatrices = EMatricesVec_[i];

      wrapper.grandImaginaryMatrices = grandImaginaryMatrixVec_[i];
      wrapper.eps_zz_Inv = eps_zz_Inv_MatrixVec_[i];
      wrapper.polar = options_.polarization;
      switch (options_.IntegralMethod) {
        case GAUSSLEGENDRE_:{
          Phi_[i] = gauss_legendre(degree_, wrapperFunQuadgl, &wrapper, kxStart_, kxEnd_);
          break;
        }
        case GAUSSKRONROD_:{
          double err;
          adapt_integrate(1, wrapperFunQuadgk, &wrapper, 1, &kxStart_, &kxEnd_, 0, ABSERROR, RELERROR, &Phi_[i], &err);
          break;
        }
        default:{
          break;
        }
      }
      Phi_[i] *= POW3(omegaList_[i] / datum::c_0) / POW2(datum::pi);
    }
  }
  /*==============================================*/
  // Implementaion of getPhi to overwrite the base class
  /*==============================================*/
  double* SimulationPlanar::getPhi(){
    return Phi_;
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
  // This function add grating to a layer
  // @args:
  // layerName: the name of the layer
  // materialName: the name of the material
  // center: the center of the grating
  // width: the width of the grating
  /*==============================================*/
  void SimulationGrating::setLayerPatternGrating(
    const std::string layerName,
    const std::string materialName,
    const double center,
    const double width
  ){
    if(materialInstanceMap_.find(materialName) == materialInstanceMap_.cend()){
      std::cerr << materialName + ": Material does not exist!" << std::endl;
      throw UTILITY::IllegalNameException(materialName + ": Material does not exist!");
      return;
    }
    if(layerInstanceMap_.find(layerName) == layerInstanceMap_.cend()){
      std::cerr << layerName + ": Layer does not exist!" << std::endl;
      throw UTILITY::IllegalNameException(layerName + ": Layer does not exist!");
      return;
    }
    Ptr<Material> material = materialInstanceMap_.find(materialName)->second;
    Ptr<Layer> layer = layerInstanceMap_.find(layerName)->second;
    layer->addGratingPattern(material, center, width);
  }
  /*==============================================*/
  // function using adaptive resolution algorithm
  /*==============================================*/
  void SimulationGrating::optUseAdaptive(){
    options_.FMMRule = SPATIALADAPTIVE_;
  }
  /*==============================================*/
  // Implementaion of the class on 2D patterning simulation
  /*==============================================*/
  SimulationPattern::SimulationPattern() : Simulation(){
    prefactor_ = 1;
    dim_ = TWO_;
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
  void SimulationPattern::setLayerPatternRectangle(
    const std::string layerName,
    const std::string materialName,
    const double centerx,
    const double centery,
    const double widthx,
    const double widthy
  ){
    if(materialInstanceMap_.find(materialName) == materialInstanceMap_.cend()){
      std::cerr << materialName + ": Material does not exist!" << std::endl;
      throw UTILITY::IllegalNameException(materialName + ": Material does not exist!");
      return;
    }
    if(layerInstanceMap_.find(layerName) == layerInstanceMap_.cend()){
      std::cerr << layerName + ": Layer does not exist!" << std::endl;
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
  void SimulationPattern::setLayerPatternCircle(
    const std::string layerName,
    const std::string materialName,
    const double centerx,
    const double centery,
    const double radius
  ){
    if(materialInstanceMap_.find(materialName) == materialInstanceMap_.cend()){
      std::cerr << materialName + ": Material does not exist!" << std::endl;
      throw UTILITY::IllegalNameException(materialName + ": Material does not exist!");
      return;
    }
    if(layerInstanceMap_.find(layerName) == layerInstanceMap_.cend()){
      std::cerr << layerName + ": Layer does not exist!" << std::endl;
      throw UTILITY::IllegalNameException(layerName + ": Layer does not exist!");
      return;
    }
    Ptr<Material> material = materialInstanceMap_.find(materialName)->second;
    Ptr<Layer> layer = layerInstanceMap_.find(layerName)->second;
    double arg1[2] = {centerx, centery};
    layer->addCirclePattern(material, arg1, radius);
  }
  /*==============================================*/
  // This is a thin wrapper for the usage of smart pointer
  /*==============================================*/
  Ptr<SimulationPattern> SimulationPattern::instanceNew(){
    return new SimulationPattern();
  }

}