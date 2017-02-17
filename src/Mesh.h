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
#define ARMA_DONT_USE_WRAPPER
#include "Rcwa.h"
#include "System.h"
#include "Common.h"
#include "config.h"
#include "gauss_legendre.h"
#include <fstream>
#include "mpi.h"

namespace MESH{
using namespace SYSTEM;
using namespace RCWA;
typedef std::vector<RCWAMatrices> RCWAMatricesVec;
typedef RCWAMatricesVec::iterator RCWAMatricesIter;
typedef RCWAMatricesVec::const_iterator const_RCWAMatricesIter;

typedef struct ARGWEAPPER{
  double omega;
  RCWAVector thicknessList;
  RCWAMatrices EMatrices;
  RCWAMatrices grandImaginaryMatrices;
  RCWAMatrices dielectricMatrixInverse;
  RCWAMatrix Gx_mat;
  RCWAMatrix Gy_mat;
  SourceList sourceList;
  int targetLayer;
} ArgWrapper;

void fileLoader(std::string fileName, double* omega, dcomplex* epsilon, int size);
/*======================================================
Implementaion of the parent simulation super class
=======================================================*/
class Simulation{
public:
  Simulation();
  ~Simulation();
  void addStructure(Structure* structure);
  void enableMPI(int numOfCore = 1);
  void resetSimulation();
  void setTargetLayerByIndex(int index);
  void setTargetLayerByLayer(Layer* layer);
  void setGx(int nGx);
  void setGy(int nGy);
  void setOutputFile(std::string name);

  double* getOmegaList();
  double* getPeriodicity();

  double getPhiAtKxKy(int omegaIndex, double kx, double ky = 0);

  void build();
  void run();

protected:

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
  Structure* structure_;
  double* Phi_;
  double* omegaList_;
  std::string output_;
  int targetLayer_;
  RCWAMatricesVec EMatricesVec_;
  RCWAMatricesVec grandImaginaryMatricesVec_;
  RCWAMatricesVec dielectricMatrixInverseVec_;
  RCWAMatrix Gx_mat_;
  RCWAMatrix Gy_mat_;

  SourceList sourceList_;
  RCWAVector thicknessListVec_;

  Structure* getStructure();
  void saveToFile();
  void transformPlanar(
    RCWAMatricesVec* dielectricMatrixVecTE,
    RCWAMatricesVec* dielectricMatrixVecTM,
    RCWAMatricesVec* dielectricImMatrixVec,
    const dcomplex* epsilon,
    const int N
  );

  void transformGrating(
    RCWAMatricesVec* dielectricMatrixVecTE,
    RCWAMatricesVec* dielectricMatrixVecTM,
    RCWAMatricesVec* dielectricImMatrixVec,
    Layer* layer,
    const dcomplex* epsilonBG,
    const int N
  );

  void transformRectangle(
    RCWAMatricesVec* dielectricMatrixVecTE,
    RCWAMatricesVec* dielectricMatrixVecTM,
    RCWAMatricesVec* dielectricImMatrixVec,
    Layer* layer,
    const dcomplex* epsilonBG,
    const int N
  );

  void transformCircle(
    RCWAMatricesVec* dielectricMatrixVecTE,
    RCWAMatricesVec* dielectricMatrixVecTM,
    RCWAMatricesVec* dielectricImMatrixVec,
    Layer* layer,
    const dcomplex* epsilonBG,
    const int N
  );
};


/*======================================================
Implementaion of the class on planar simulation
=======================================================*/
class SimulationPlanar : public Simulation{
public:
  SimulationPlanar() : Simulation(){};
  ~SimulationPlanar();

  void setGx() = delete;
  void setGy() = delete;
  double getPhiAtKxKy(int omegaIndex, double kx, double ky) = delete;

  void setKxIntegral(double end);
  void run();
  double getPhiAtKx(int omegaIndex, double kx);

private:
};

/*======================================================
Implementaion of the class on 1D grating simulation
=======================================================*/
class SimulationGrating : public Simulation{
public:
  SimulationGrating() : Simulation(){
    prefactor_ = 1;
  };
  ~SimulationGrating();
  void setGy() = delete;
  void setKxIntegral(int points);
  void setKxIntegralSym(int points);
  // for ky integral, from 0 to inf
  void setKyIntegral(int points, double end);


private:
};

/*======================================================
Implementaion of the class on 2D patterning simulation
=======================================================*/
class SimulationPattern : public Simulation{
public:
  SimulationPattern() : Simulation(){
    prefactor_ = 1;
  };
  ~SimulationPattern();

  void setKxIntegral(int points);
  void setKxIntegralSym(int points);
  void setKyIntegral(int points);
  void setKyIntegralSym(int points);


private:
};
}
#endif
