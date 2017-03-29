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
#include "Cubature.h"
#include "System.h"
#include "Fmm.h"
#include "Common.h"
#include "config.h"
#include <fstream>
#if defined(_OPENMP)
  #include <omp.h>
#endif

namespace MESH{
using namespace SYSTEM;
using namespace RCWA;
// using namespace FMM;


enum INTEGRAL {GAUSSLEGENDRE_, GAUSSKRONROD_};
enum METHOD {NAIVEFMM_, INVERSERULE_, SPATIALADAPTIVE_};

typedef struct OPTIONS{
  int FMMRule = NAIVEFMM_;
  int IntegralMethod = GAUSSKRONROD_;
  bool PrintIntermediate = false;
  bool IntegrateKParallel = true;
  bool kxIntegralPreset = false;
  bool kyIntegralPreset = false;
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
  void load(const std::string fileName);
  double* getOmegaList();
  EPSILON getEpsilonList();
  int getNumOfOmega();
  FileLoader(const FileLoader&) = delete;
protected:
  ~FileLoader();
private:
  FileLoader();
  double* omegaList_;
  EPSILON epsilonList_;
  int numOfOmega_;
  bool preSet_ = false;
};
/*======================================================*/
//  definition of maps used in the simulation
/*=======================================================*/
typedef std::map< std::string, Ptr<Layer> > LayerInstanceMap;
typedef LayerInstanceMap::const_iterator const_LayerInstanceIter;
typedef std::map< std::string, Ptr<Material> > MaterialInstanceMap;

/*======================================================*/
//  Implementaion of the parent simulation super class
/*=======================================================*/
class Simulation : public PtrInterface{
public:

  // adding new function for the same interface as S4

  void setPeriodicity(const double p1, const double p2 = 0);
  void addMaterial(const std::string name, const std::string infile);
  void setMaterial(const std::string name, double** epsilon, const std::string type);

  void addLayer(const std::string name, const double thick, const std::string materialName);
  void setLayer(const std::string name, const double thick, const std::string materialName);
  void setLayerThickness(const std::string name, const double thick);
  void addLayerCopy(const std::string name, const std::string originalName);
  void deleteLayer(const std::string name);

  void setLayerPatternGrating(
    const std::string layerName,
    const std::string materialName,
    const double center,
    const double width
  );

  void setLayerPatternRectangle(
    const std::string layerName,
    const std::string materialName,
    const double centerx,
    const double centery,
    const double widthx,
    const double widthy
  );

  void setLayerPatternCircle(
    const std::string layerName,
    const std::string materialName,
    const double centerx,
    const double centery,
    const double radius
  );

  void setSourceLayer(const std::string name);
  void setProbeLayer(const std::string name);

  void setGx(const int nGx);
  void setGy(const int nGy);
  void saveToFile(const std::string fileName);
  double* getPhi();
  double* getOmega();
  int getNumOfOmega();

  double getPhiAtKxKy(const int omegaIndex, const double kx, const double ky = 0);
  void buildRCWA();

  void getSysInfo();

  void optUseInverseRule();
  void optUseNaiveRule();
  void optPrintIntermediate();
  void setThread(const int numThread);

  void setKxIntegral(const int points, const double end = 0);
  void setKxIntegralSym(const int points, const double end = 0);
  void setKyIntegral(const int points, const double end = 0);
  void setKyIntegralSym(const int points, const double end = 0);

  void integrateKxKy();
  ~Simulation();
protected:
  Simulation();
  Simulation(const Simulation&) = delete;

  void resetSimulation();
  void setTargetLayerByLayer(const Ptr<Layer>& layer);
  Ptr<Structure> getStructure();


  int nGx_;
  int nGy_;
  int numOfOmega_;
  double* Phi_;
  double* omegaList_;
  double kxStart_;
  double kxEnd_;
  double kyStart_;
  double kyEnd_;

  int numOfKx_;
  int numOfKy_;

  double period_[2];
  int prefactor_;

  LayerInstanceMap layerInstanceMap_;
  MaterialInstanceMap materialInstanceMap_;
  Ptr<Structure> structure_;
  Ptr<FileLoader> fileLoader_;



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

  int MPI_Initialized_, MPI_Finalized_;
  int numOfThread_ = 1;
};


/*======================================================*/
// Implementaion of the class on planar simulation
/*=======================================================*/
class SimulationPlanar : public Simulation{
public:
  static Ptr<SimulationPlanar> instanceNew();
  SimulationPlanar(const SimulationPlanar&) = delete;
  void useInverseRule() = delete;
  void setGx() = delete;
  void setGy() = delete;

  // this function is used when one knows that the problem is only a kx integral
  void setKParallelIntegral(const double end);

  void optUseQuadgl(int degree = DEGREE);
  void optUseQuadgk();

  void integrateKParallel();
  double getPhiAtKParallel(const int omegaIndex, const double KParallel);

  SimulationPlanar();
protected:

private:

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

  void optUseAdaptive();
  SimulationGrating();
protected:

private:
};

/*======================================================*/
// Implementaion of the class on 2D patterning simulation
/*=======================================================*/
class SimulationPattern : public Simulation{
public:

  static Ptr<SimulationPattern> instanceNew();
  SimulationPattern();
  SimulationPattern(const SimulationPattern&) = delete;

protected:

private:
};
}
#endif
