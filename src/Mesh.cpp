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
#include "time.h"
namespace MESH{
  /*==============================================
  This function loads data from disk
  @args:
  fileName: the name of the input file
  omega: the input omega list
  epsilon: the input epsilon list
  size: the size of omega
  ==============================================*/
  void fileLoader(
    const std::string fileName,
    double* omega,
    dcomplex* epsilon,
    const int size
  ){
    std::ifstream inputFile(fileName);
    if(!inputFile.good()){
      throw UTILITY::FileNotExistException(fileName + " not exists!");
    }
    std::string line;
    double real, imag;
    for(size_t i = 0; i < size; i++){
      inputFile >> omega[i] >> real >> imag;
      epsilon[i] = dcomplex(real, -imag);
    }
    inputFile.close();
  }

  /*==============================================
  This function wraps the data for quad_legendre
  @args:
  kx: the kx value (normalized)
  data: wrapper for all the arguments wrapped in wrapper
  ==============================================*/
  double wrapperFun(const double kx, void* data){
    ArgWrapper wrapper = *(ArgWrapper*) data;
    return kx * poyntingFlux(
      wrapper.omega,
      &(wrapper.thicknessList),
      kx,
      0,
      &(wrapper.EMatrices),
      &(wrapper.grandImaginaryMatrices),
      &(wrapper.dielectricMatrixZInv),
      &(wrapper.Gx_mat),
      &(wrapper.Gy_mat),
      &(wrapper.sourceList),
      wrapper.targetLayer,
      1
    );
  }
  /*======================================================
  Implementaion of the parent simulation super class
  =======================================================*/
  Simulation::Simulation() : nGx_(0), nGy_(0), numOfOmega_(0), structure_(nullptr),
  Phi_(nullptr), omegaList_(nullptr), kxStart_(0), kxEnd_(0), kyStart_(0), kyEnd_(0), numOfKx_(0), numOfKy_(0)
  {
    period_ = new double[2];
    period_[0] = 0;
    period_[1] = 0;
    output_ = "output.txt";
    targetLayer_ = -1;
    dim_ = NO_;
  }

  Simulation::~Simulation(){
    delete[] Phi_;
    Phi_ = nullptr;
    delete[] period_;
    period_ = nullptr;
  }

  /*==============================================
  This function saves data to disk
  ==============================================*/
  void Simulation::saveToFile(){
    std::ofstream outputFile(output_);
    for(size_t i = 0; i < numOfOmega_; i++){
      outputFile << omegaList_[i] << "\t" << Phi_[i] << std::endl;
    }
    outputFile.close();
  }

  /*==============================================
  This function adds structure to the simulation
  @args:
  structure: the structure of the simulation
  ==============================================*/
  void Simulation::addStructure(const Ptr<Structure> structure){
    structure_ = structure;
    period_ = structure_->getPeriodicity();
  }
  /*==============================================
  This function set the output file name
  @args:
  name: the output file name
  ==============================================*/
  void Simulation::setOutputFile(std::string name){
    output_ = name;
  }
  /*==============================================
  This function sets number of positive Gx
  @args:
  Gx: number of positive Gx
  ==============================================*/
  void Simulation::setGx(const int Gx){
    nGx_ = Gx;
  }

  /*==============================================
  This function sets number of positive Gy
  @args:
  Gy: number of positive Gy
  ==============================================*/
  void Simulation::setGy(const int Gy){
    nGy_ = Gy;
  }

  /*==============================================
  This function sets the target layer by layer index
  @args:
  index: target layer index
  ==============================================*/
  void Simulation::setTargetLayerByIndex(const int index){
    if(index >= structure_->getNumOfLayer()){
      throw UTILITY::RangeException(std::to_string(index) + ": out of range!");
    }
    targetLayer_ = index;
  }

  /*==============================================
  This function sets the target layer by layer
  @args:
  index: target layer
  ==============================================*/
  void Simulation::setTargetLayerByLayer(const Ptr<Layer> layer){
    for(size_t i = 0; i < structure_->getNumOfLayer(); i++){
      if(structure_->getLayerByIndex(i) == layer){
        targetLayer_ = i;
        return;
      }
    }
  }

  /*==============================================
  This function cleans up the simulation
  ==============================================*/
  void Simulation::resetSimulation(){
    for(size_t i = 0; i < EMatricesVec_.size(); i++){
      EMatricesVec_[i].clear();
      grandImaginaryMatricesVec_.clear();
      dielectricMatrixZInvVec_[i].clear();
    }
    EMatricesVec_.clear();
    grandImaginaryMatricesVec_.clear();
    dielectricMatrixZInvVec_.clear();
  }

  /*==============================================
  This function gets the structure
  ==============================================*/
  Ptr<Structure> Simulation::getStructure(){
    return structure_;
  }

  /*==============================================
  This function gets the omega array
  ==============================================*/
  double* Simulation::getOmegaList(){
    return omegaList_;
  }

  /*==============================================
  This function gets the periodicity
  ==============================================*/
  double* Simulation::getPeriodicity(){
    return period_;
  }

  /*==============================================
  This function gets the Phi at given kx and ky
  @args:
  omegaIndex: the index of omega
  kx: the kx value, normalized
  ky: the ky value, normalized
  @note
  used by grating and patterning
  ==============================================*/
  double Simulation::getPhiAtKxKy(const int omegaIdx, const double kx, const double ky, const int N){
    if(omegaIdx >= numOfOmega_){
      throw UTILITY::RangeException(std::to_string(omegaIdx) + ": out of range!");
    }
    return POW3(omegaList_[omegaIdx] / datum::c_0) / POW3(datum::pi) / 2.0 *
      poyntingFlux(omegaList_[omegaIdx] / datum::c_0, &thicknessListVec_, kx, ky, &(EMatricesVec_[omegaIdx]),
      &(grandImaginaryMatricesVec_[omegaIdx]), &(dielectricMatrixZInvVec_[omegaIdx]), &Gx_mat_, &Gy_mat_,
      &sourceList_, targetLayer_,N);
  }

  /*==============================================
  This function computes the Fourier transform for planar
  @args:
  dielectricMatrixVecTE: the vector containing dielectric matrices
  dielectricImMatrixVec: the vector containing imaginary part matrices
  epsilonBG: the epsilon values for the background
  N: the number of total G
  ==============================================*/
  void Simulation::transformPlanar(
    RCWAMatricesVec* dielectricMatrixVecTE,
    RCWAMatricesVec* dielectricMatrixVecTM,
    RCWAMatricesVec* dielectricImMatrixVec,
    const dcomplex* epsilonBG,
    const int N
  ){
    RCWAMatrix onePadding1N = eye<RCWAMatrix>(N, N);
    for(size_t i = 0; i < numOfOmega_; i++){
      (*dielectricMatrixVecTE)[i].push_back(onePadding1N * epsilonBG[i]);
      (*dielectricMatrixVecTM)[i].push_back(onePadding1N * epsilonBG[i]);
      (*dielectricImMatrixVec)[i].push_back(onePadding1N * (epsilonBG[i]).imag());
      dielectricMatrixZInvVec_[i].push_back(onePadding1N * dcomplex(1, 0) / epsilonBG[i]);
    }
  }
  /*==============================================
  This function computes the Fourier transform for grating
  @args:
  dielectricMatrixVecTE: the vector containing dielectric matrices
  dielectricImMatrixVec: the vector containing imaginary part matrices
  layer: the current layer
  epsilonBG: the epsilon values for the background
  N: the number of total G
  ==============================================*/
  void Simulation::transformGrating(
    RCWAMatricesVec* dielectricMatrixVecTE,
    RCWAMatricesVec* dielectricMatrixVecTM,
    RCWAMatricesVec* dielectricImMatrixVec,
    const Ptr<Layer> layer,
    const dcomplex* epsilonBG,
    const int N
  ){
    dcomplex IMAG_I = dcomplex(0, 1);
    RCWAMatrix G_row(1, N), G_col(N, 1);
    for(int i = 0; i < N; i++){
      G_row(0, i) = -i * 2.0 * datum::pi / period_[0];
      G_col(i, 0) = -G_row(0, i);
    }
    RCWAMatrix G_mat = toeplitz(G_col, G_row);

    RCWAMatrix onePadding1N = eye<RCWAMatrix>(N, N);
    int numOfMaterial = layer->getNumOfMaterial();
    RCWAVector centerVec(numOfMaterial), widthVec(numOfMaterial);
    int count = 0;
    for(const_PatternIter it = layer->getArg1Begin(); it != layer->getArg1End(); it++){
      centerVec(count) = (it->first + it->second) / 2;
      widthVec(count) = it->second - it->first;
      count++;
    }
    for (size_t i = 0; i < numOfOmega_; i++) {
      RCWAMatrix dielectricMatrix(N, N, fill::zeros), dielectricImMatrix(N, N, fill::zeros), dielectricMatrixInv(N, N, fill::zeros);
      count = 0;
      for(const_MaterialIter it = layer->getVecBegin(); it != layer->getVecEnd(); it++){
        dcomplex epsilon = (*it)->getEpsilonAtIndex(i);

        dielectricMatrix += exp(IMAG_I * G_mat * centerVec(count)) * (epsilon - epsilonBG[i])
          * widthVec(count) % sinc(G_mat / 2 * widthVec(count));

        dielectricMatrixInv += exp(IMAG_I * G_mat * centerVec(count)) * (dcomplex(1.0, 0) / epsilon - dcomplex(1.0, 0) / epsilonBG[i])
          * widthVec(count) % sinc(G_mat / 2 * widthVec(count));

        dielectricImMatrix += exp(IMAG_I * G_mat * centerVec(count)) * (epsilon.imag() - epsilonBG[i].imag())
          * widthVec(count) % sinc(G_mat / 2 * widthVec(count));
        count++;
      }
      dielectricMatrix /= period_[0];
      dielectricMatrixInv /= period_[0];
      dielectricImMatrix /= period_[0];
      dielectricMatrix += epsilonBG[i] * eye<RCWAMatrix>(N, N);
      dielectricMatrixInv += dcomplex(1.0,0) / epsilonBG[i] * eye<RCWAMatrix>(N, N);
      dielectricImMatrix += epsilonBG[i].imag() * eye<RCWAMatrix>(N, N);

      (*dielectricMatrixVecTE)[i].push_back(dielectricMatrix);
      (*dielectricMatrixVecTM)[i].push_back(dielectricMatrixInv.i());
      (*dielectricImMatrixVec)[i].push_back(dielectricImMatrix);
      dielectricMatrixZInvVec_[i].push_back(dielectricMatrix.i());
    }

  }

  /*==============================================
  This function computes the Fourier transform for rectangle
  @args:
  dielectricMatrixVecTE: the vector containing dielectric matrices
  dielectricImMatrixVec: the vector containing imaginary part matrices
  layer: the current layer
  epsilonBG: the epsilon values for the background
  N: the number of total G
  ==============================================*/
  void Simulation::transformRectangle(
    RCWAMatricesVec* dielectricMatrixVecTE,
    RCWAMatricesVec* dielectricMatrixVecTM,
    RCWAMatricesVec* dielectricImMatrixVec,
    const Ptr<Layer> layer,
    const dcomplex* epsilonBG,
    const int N
  ){
    dcomplex IMAG_I = dcomplex(0, 1.0);
    RCWAMatrix Gx_r, Gx_l, Gy_r, Gy_l;
    meshGrid(&Gx_mat_, &Gx_mat_, &Gx_r, &Gx_l);
    meshGrid(&Gy_mat_, &Gy_mat_, &Gy_r, &Gy_l);

    RCWAMatrix GxMat = Gx_l - Gx_r;
    RCWAMatrix GyMat = Gy_l - Gy_r;

    int numOfMaterial = layer->getNumOfMaterial();

    RCWAVector centerxVec(numOfMaterial), centeryVec(numOfMaterial), widthxVec(numOfMaterial), widthyVec(numOfMaterial);
    int count = 0;
    for(const_PatternIter it = layer->getArg1Begin(); it != layer->getArg1End(); it++){
      centerxVec(count) = it->first;
      centeryVec(count) = it->second;
      count++;
    }
    count = 0;
    for(const_PatternIter it = layer->getArg2Begin(); it != layer->getArg2End(); it++){
      widthxVec(count) = it->first;
      widthyVec(count) = it->second;
      count++;
    }

    for (size_t i = 0; i < numOfOmega_; i++) {
      count = 0;
      RCWAMatrix dielectricMatrix(N, N, fill::zeros), dielectricImMatrix(N, N, fill::zeros), dielectricMatrixInv(N, N, fill::zeros);
      for(const_MaterialIter it = layer->getVecBegin(); it != layer->getVecEnd(); it++){
        dcomplex epsilon = (*it)->getEpsilonAtIndex(i);

        dielectricMatrix += widthxVec(count) * widthyVec(count) / (period_[0] * period_[1]) * (epsilon - epsilonBG[i])
          * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
          % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);
        dielectricMatrixInv += widthxVec(count) * widthyVec(count) / (period_[0] * period_[1]) * (dcomplex(1.0,0) / epsilon - dcomplex(1.0,0) / epsilonBG[i])
          * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
          % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);

        dielectricMatrix += widthxVec(count) * widthyVec(count) / (period_[0] * period_[1]) * (epsilon.imag() - epsilonBG[i].imag())
          * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
          % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);
        count++;
      }

      dielectricMatrix += epsilonBG[i] * eye<RCWAMatrix>(N, N);
      dielectricMatrixInv += dcomplex(1, 0) / epsilonBG[i] * eye<RCWAMatrix>(N, N);
      dielectricImMatrix += epsilonBG[i].imag() * eye<RCWAMatrix>(N, N);

      (*dielectricMatrixVecTE)[i].push_back(dielectricMatrixInv.i());
      (*dielectricMatrixVecTM)[i].push_back(dielectricMatrixInv.i());
      (*dielectricImMatrixVec)[i].push_back(dielectricImMatrix);
      dielectricMatrixZInvVec_[i].push_back(dielectricMatrix.i());
    }
  }

  /*==============================================
  This function computes the Fourier transform for circle
  @args:
  dielectricMatrixVecTE: the vector containing dielectric matrices
  dielectricImMatrixVec: the vector containing imaginary part matrices
  layer: the current layer
  epsilonBG: the epsilon values for the background
  N: the number of total G
  ==============================================*/
  void Simulation::transformCircle(
      RCWAMatricesVec* dielectricMatrixVecTE,
      RCWAMatricesVec* dielectricMatrixVecTM,
      RCWAMatricesVec* dielectricImMatrixVec,
      const Ptr<Layer> layer,
      const dcomplex* epsilonBG,
      const int N
    ){
    // TODO
  }

  /*==============================================
  This function builds up the matrices
  ==============================================*/
  void Simulation::build(){
    // essential, get the shared Gx_mat_ and Gy_mat_
    getGMatrices(nGx_, nGy_, period_, &Gx_mat_, &Gy_mat_, dim_);
    // get constants
    Ptr<Layer> firstLayer = structure_->getLayerByIndex(0);
    Ptr<Material> backGround = firstLayer->getBackGround();
    numOfOmega_ = backGround->getNumOfOmega();
    omegaList_ = backGround->getOmegaList();
    Phi_ = new double[numOfOmega_];
    EMatricesVec_.resize(numOfOmega_);
    grandImaginaryMatricesVec_.resize(numOfOmega_);
    dielectricMatrixZInvVec_.resize(numOfOmega_);

    RCWAMatricesVec dielectricMatrixVecTE(numOfOmega_), dielectricMatrixVecTM(numOfOmega_), dielectricImMatrixVec(numOfOmega_);
    int numOfLayer = structure_->getNumOfLayer();
    int N = getN(nGx_, nGy_);
    for(size_t i = 0; i < numOfLayer; i++){
      Ptr<Layer> layer = structure_->getLayerByIndex(i);
      dcomplex* epsilonBG = (layer->getBackGround())->getEpsilon();
      switch (layer->getPattern()) {
        /*************************************
        /* if the pattern is a plane */
        case PLANAR_:{
          this->transformPlanar(
            &dielectricMatrixVecTE,
            &dielectricMatrixVecTM,
            &dielectricImMatrixVec,
            epsilonBG,
            N
          );
          break;
        }
        /*************************************
        /* if the pattern is a grating (1D) */
        /************************************/
        case GRATING_:{
          this->transformGrating(
            &dielectricMatrixVecTE,
            &dielectricMatrixVecTM,
            &dielectricImMatrixVec,
            layer,
            epsilonBG,
            N
          );
          break;
        }

        /*************************************
        /* if the pattern is a rectangle (2D) */
        /************************************/
        case RECTANGLE_:{
          this->transformRectangle(
            &dielectricMatrixVecTE,
            &dielectricMatrixVecTM,
            &dielectricImMatrixVec,
            layer,
            epsilonBG,
            N
          );
          break;
        }

      /*************************************
      /* if the pattern is a circle (2D) */
      /************************************/
        case CIRCLE_:{
          this->transformCircle(
            &dielectricMatrixVecTE,
            &dielectricMatrixVecTM,
            &dielectricImMatrixVec,
            layer,
            epsilonBG,
            N
          );
          break;
        }
        default: break;
      }
    }
    for(size_t i = 0; i < numOfOmega_; i++){
      getEMatrices(
        &(EMatricesVec_[i]),
        &(dielectricMatrixVecTE[i]),
        &(dielectricMatrixVecTM[i]),
        numOfLayer,
        N
      );

      getGrandImaginaryMatrices(
        &(grandImaginaryMatricesVec_[i]),
        &(dielectricImMatrixVec[i]),
        numOfLayer,
        N
      );
    }

    thicknessListVec_ = zeros<RCWAVector>(numOfLayer);
    sourceList_.resize(numOfLayer);
    for(size_t i = 0; i < numOfLayer; i++){
      thicknessListVec_(i) = (structure_->getLayerByIndex(i))->getThickness();
      sourceList_[i] = (structure_->getLayerByIndex(i))->checkIsSource();
    }

  }

  /*==============================================
  This function computes the flux
  ==============================================*/
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
          // this function not run
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
    int N = getN(nGx_, nGy_);
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
        resultArray[i] = this->getPhiAtKxKy(omegaIdx, kxList[kxIdx] / scalex[omegaIdx], kyList[kyIdx] / scaley[omegaIdx], N);
        std::cout << kxList[kxIdx] / scalex[omegaIdx] << "\t" << kyList[kyIdx]  / scaley[omegaIdx] << "\t" << resultArray[i] << std::endl;
      }

      // wait for all the slave process finished
      for(size_t thread = 1; thread < numProcs; thread++){
        MPI_Recv(&start, 1, MPI_INT, thread, RECVTAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&end, 1, MPI_INT, thread, RECVTAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&resultArray[start], end - start, MPI_DOUBLE, thread, RECVTAG, MPI_COMM_WORLD, &status);
      }
      for(size_t i = 0; i < numOfOmega_; i++){
        for(size_t j = 0; j < numOfKx_ * numOfKy_; j++){
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
        resultArray[i] = this->getPhiAtKxKy(omegaIdx, kxList[kxIdx] / scalex[omegaIdx], kyList[kyIdx] / scaley[omegaIdx], N);
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
  /*======================================================
  Implementaion of the class on planar simulation
  =======================================================*/
  SimulationPlanar::SimulationPlanar() : Simulation(){
    dim_ = NO_;
    degree_ = DEGREE;
  }
  /*======================================================
  This is a thin wrapper for the usage of smart pointer
  =======================================================*/
  Ptr<SimulationPlanar> SimulationPlanar::instanceNew(){
    return new SimulationPlanar();
  };
  /*======================================================
  Function setting the integral over kx
  @args:
  end: the end of the integration
  degree: the degree of gauss_legendre integral, default 512
  =======================================================*/
  void SimulationPlanar::setKxIntegral(const double end, const int degree){
    kxStart_ = 0;
    numOfKx_ = 0;
    kxEnd_ = end;
    degree_ = degree;
  }

  /*==============================================
  This function gets the flux at a given kx
  @args:
  omegaIndex: the index of omega
  kx: the kx value, normalized
  ==============================================*/
  double SimulationPlanar::getPhiAtKx(const int omegaIdx, const double kx){
    if(omegaIdx >= numOfOmega_){
      throw UTILITY::RangeException(std::to_string(omegaIdx) + ": out of range!");
    }
    return POW3(omegaList_[omegaIdx] / datum::c_0) / POW2(datum::pi) * kx *
      poyntingFlux(omegaList_[omegaIdx] / datum::c_0, &thicknessListVec_, kx, 0, &(EMatricesVec_[omegaIdx]),
      &(grandImaginaryMatricesVec_[omegaIdx]), &(dielectricMatrixZInvVec_[omegaIdx]), &Gx_mat_, &Gy_mat_,
      &sourceList_, targetLayer_,1);
  }


  /*==============================================
  This function integrates kx using gauss_legendre method
  ==============================================*/
  void SimulationPlanar::run(){
    MPI_Status status;
    int rank, numProcs, start, end;
    int offset = 0;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    int displs[numProcs], sendCounts[numProcs];
    int chunkSize = numOfOmega_ / numProcs;
    int numLeft = numOfOmega_ - numProcs * chunkSize;
    double recvBuf[chunkSize + 1];
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
      wrapper.grandImaginaryMatrices = grandImaginaryMatricesVec_[omegaIdx];
      wrapper.dielectricMatrixZInv = dielectricMatrixZInvVec_[omegaIdx];
      recvBuf[i] = POW3(omegaList_[omegaIdx] / datum::c_0) / POW2(datum::pi) *
        gauss_legendre(degree_, wrapperFun, &wrapper, kxStart_, kxEnd_);
    }
    // gatther the Phi from each processor
    MPI_Gatherv(&recvBuf[0], 1, dtype, &Phi_[0], sendCounts, displs, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == MASTER){
        this->saveToFile();
    }
    MPI_Finalize();

  }
  /*======================================================
  Implementaion of the class on 1D grating simulation
  =======================================================*/
  SimulationGrating::SimulationGrating() : Simulation(){
    prefactor_ = 2;
    dim_ = ONE_;
  }
  /*======================================================
  This is a thin wrapper for the usage of smart pointer
  =======================================================*/
  Ptr<SimulationGrating> SimulationGrating::instanceNew(){
    return new SimulationGrating();
  }
  /*======================================================
  Function setting the integral over kx
  @args:
  points: number of points of sampling kx
  =======================================================*/
  void SimulationGrating::setKxIntegral(const int points){
    numOfKx_ = points;
    if(period_[0] == 0.0){
      throw UTILITY::ValueException("Periodicity not set!");
    }
    kxStart_ = -datum::pi / period_[0];
    kxEnd_ = -kxStart_;
  }

  /*==============================================
  This function set the integral of kx when the system is symmetric in x direction
  @args:
  points: number of kx points
  ==============================================*/
  void SimulationGrating::setKxIntegralSym(const int points){
    numOfKx_ = points;
    kxStart_ = 0;
    if(period_[0] == 0.0){
      throw UTILITY::ValueException("Periodicity not set!");
    }
    kxEnd_ = datum::pi / period_[0];
    prefactor_ *= 2;
  }

  /*==============================================
  This function set the integral of ky
  @args:
  points: number of ky points
  end: the upperbound of the integral
  ==============================================*/
  void SimulationGrating::setKyIntegral(const int points, const double end){
    kyStart_ = 0;
    numOfKy_ = points;
    kyEnd_ = end;
  }

  /*======================================================
  Implementaion of the class on 2D patterning simulation
  =======================================================*/
  SimulationPattern::SimulationPattern() : Simulation(){
    prefactor_ = 1;
    dim_ = TWO_;
  }
  /*======================================================
  This is a thin wrapper for the usage of smart pointer
  =======================================================*/
  Ptr<SimulationPattern> SimulationPattern::instanceNew(){
    return new SimulationPattern();
  }

  /*======================================================
  Function setting the integral over kx
  @args:
  points: number of points of sampling kx
  =======================================================*/
  void SimulationPattern::setKxIntegral(const int points){
    if(period_[0] == 0.0){
      throw UTILITY::ValueException("Periodicity not set!");
    }
    kxStart_ = -datum::pi / period_[0];
    numOfKx_ = points;
    kxEnd_ = -kxStart_;
  }

  /*==============================================
  This function set the integral of kx when the system is symmetric in x direction
  @args:
  points: number of kx points
  ==============================================*/
  void SimulationPattern::setKxIntegralSym(const int points){
    kxStart_ = 0;
    numOfKx_ = points;
    if(period_[0] == 0.0){
      throw UTILITY::ValueException("Periodicity not set!");
    }
    kxEnd_ = datum::pi / period_[0];
    prefactor_ *= 2;
  }

  /*==============================================
  This function set the integral of ky
  @args:
  points: number of ky points
  end: the upperbound of the integral
  ==============================================*/
  void SimulationPattern::setKyIntegral(const int points){
    if(period_[1] == 0.0){
      throw UTILITY::ValueException("Periodicity not set!");
    }
    kyStart_ = -datum::pi / period_[1];
    numOfKy_ = points;
    kyEnd_ = -kyStart_;
  }

  /*==============================================
  This function set the integral of ky when the system is symmetric in y direction
  @args:
  points: number of ky points
  ==============================================*/
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