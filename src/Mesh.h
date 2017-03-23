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

#ifndef _MESH_H
#define _MESH_H
//#ifndef ARMA_DONT_USE_WRAPPER
//#define ARMA_DONT_USE_WRAPPER
#include "Rcwa.h"
#include "System.h"
#include "Fmm.h"
#include "Common.h"
#include "config.h"
#include <fstream>
#include "mpi.h"

namespace MESH{
using namespace SYSTEM;
using namespace RCWA;
// using namespace FMM;


enum INTEGRAL {GAUSSLEGENDRE_, GAUSSKRONROD_};
enum METHOD {NAIVEFMM_, INVERSERULE_, SPATIALADAPTIVE_};

typedef struct OPTIONS{
  int FMMRule = NAIVEFMM_;
  int IntegralMethod = GAUSSKRONROD_;
} Options;


typedef struct ARGWEAPPER{
  double omega;
  RCWAVector thicknessList;
  RCWAMatrices EMatrices;
  RCWAMatrices grandImaginaryMatrices;
  RCWAMatrices eps_zz_Inv;
  RCWAMatrix Gx_mat;
  RCWAMatrix Gy_mat;
  SourceList sourceList;
  int targetLayer;
} ArgWrapper;

/*======================================================*/
//  Implementaion of the FileLoader class
/*=======================================================*/
class FileLoader : public PtrInterface{
public:
  static Ptr<FileLoader> instanceNew();
  static Ptr<FileLoader> instanceNew(const int numOfOmega);
  void load(const std::string fileName);
  double* getOmegaList();
  EPSILON getEpsilonList();
  int getNumOfOmega();
  FileLoader(const FileLoader&) = delete;
protected:
  ~FileLoader();
private:
  FileLoader(const int numOfOmega = 0);
  double* omegaList_;
  EPSILON epsilonList_;
  int numOfOmega_;
  bool preSet_ = false;
};

/*======================================================*/
//  Implementaion of the parent simulation super class
/*=======================================================*/
class Simulation : public PtrInterface{
public:
  void addStructure(const Ptr<Structure>& structure);
  void setTargetLayerByIndex(const int index);
  void setTargetLayerByLayer(const Ptr<Layer>& layer);
  void setGx(const int nGx);
  void setGy(const int nGy);
  void setOutputFile(const std::string name);

  double* getOmegaList();
  double* getPeriodicity();

  double getPhiAtKxKy(const int omegaIndex, const double kx, const double ky = 0);
  void build();
  void rebuild();
  void run();

protected:
  Simulation();
  Simulation(const Simulation&) = delete;
  ~Simulation();
  void resetSimulation();
  double* period_;
  double kxStart_;
  double kxEnd_;
  int numOfKx_;

  double kyStart_;
  double kyEnd_;
  int numOfKy_;
  int prefactor_;
  int nGx_;
  int nGy_;
  int numOfOmega_;
  Ptr<Structure> structure_;



  double* Phi_;
  double* omegaList_;
  std::string output_;
  int targetLayer_;

  RCWAMatrix Gx_mat_;
  RCWAMatrix Gy_mat_;

  RCWAMatricesVec EMatricesVec_;
  RCWAMatricesVec grandImaginaryMatrixVec_;
  RCWAMatricesVec eps_zz_Inv_MatrixVec_;

  SourceList sourceList_;
  RCWAVector thicknessListVec_;
  DIMENSION dim_;
  Options options_;

  Ptr<Structure> getStructure();
  void saveToFile();

};


/*======================================================*/
// Implementaion of the class on planar simulation
/*=======================================================*/
class SimulationPlanar : public Simulation{
public:
  static Ptr<SimulationPlanar> instanceNew();
  SimulationPlanar(const SimulationPlanar&) = delete;
  void setGx() = delete;
  void setGy() = delete;

  // this function is used when one knows that the problem is only a kx integral
  void setKParallelIntegral(const double end);

  // this function is used when one knows that the problem is not a simple kx integral
  void setKxIntegral(const int points, const double end);
  void setKxIntegralSym(const int points, const double end);
  void setKyIntegral(const int points, const double end);
  void setKyIntegralSym(const int points, const double end);

  void useQuadgl(int degree = DEGREE);
  void useQuadgk(int degree = DEGREE);

  void runNaive();
  double getPhiAtKParallel(const int omegaIndex, const double KParallel);

protected:
  ~SimulationPlanar(){};

private:

  SimulationPlanar();
  int degree_ = DEGREE;
};

/*======================================================*/
// Implementaion of the class on 1D grating simulation
/*=======================================================*/
class SimulationGrating : public Simulation{
public:

  static Ptr<SimulationGrating> instanceNew();
  SimulationGrating(const SimulationGrating&) = delete;

  void setGy() = delete;
  void setKxIntegral(const int points);
  void setKxIntegralSym(const int points);
  // for ky integral, from 0 to inf
  void setKyIntegral(const int points, const double end);
  void setKyIntegralSym(const int points, const double end);
  void useAdaptive();
  void useNative();
protected:
  ~SimulationGrating(){};
private:
  SimulationGrating();
};

/*======================================================*/
// Implementaion of the class on 2D patterning simulation
/*=======================================================*/
class SimulationPattern : public Simulation{
public:

  static Ptr<SimulationPattern> instanceNew();

  SimulationPattern(const SimulationPattern&) = delete;

  void setKxIntegral(const int points);
  void setKxIntegralSym(const int points);
  void setKyIntegral(const int points);
  void setKyIntegralSym(const int points);
protected:
  ~SimulationPattern(){};

private:
  SimulationPattern();
};
}
#endif
