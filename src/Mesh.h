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
}ArgWrapper;

void filerLoader(std::string fileName, double* omega, dcomplex* epsilon, int size);
void saveData(std::string fileName, double* omega, double* fluxSpectrum, int size);
double wrapperFun(double kx, ArgWrapper* wrapper);
/*======================================================
Implementaion of the parent simulation super class
=======================================================*/
class Simulation{
public:
  Simulation();
  ~Simulation();
  void addStructure(Structure* structure);
  void enableMPI(int numOfCore = 1);
  void setPeriod(double p1, double p2 = 0);
  void resetSimulation();
  void setTargetLayerByIndex(int index);
  void setTargetLayerByLayer(Layer* layer);

  Structure* getStructure();
  double* getOmegaList();
  double* getPhi();
  double* getPeriod();

  double getPhiAtKxKy(int omegaIndex, double kx, double ky = 0);

  void build();
  void run();

protected:
  int nGx_;
  int nGy_;
  int numOfCore_;
  int numOfOmega_;
  double period_[2];
  Structure* structure_;
  double* Phi_;
  double* omegaList_;

  double kxStart_;
  double kxEnd_;
  double numOfKx_;

  double kyStart_;
  double kyEnd_;
  double numOfKy_;
  int targetLayer_;

  RCWAMatricesVec EMatricesVec_;
  RCWAMatricesVec grandImaginaryMatricesVec_;
  RCWAMatricesVec dielectricMatrixInverseVec_;
  RCWAMatrix Gx_mat_;
  RCWAMatrix Gy_mat_;

  SourceList sourceList_;
  RCWAVector thicknessListVec_;
};


/*======================================================
Implementaion of the class on planar simulation
=======================================================*/
class SimulationPlanar : public Simulation{
public:
  SimulationPlanar() : Simulation(){};
  ~SimulationPlanar();

  void setKxIntegral(double points, double end);
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
  void setKxIntegral(double points);
  void setKxIntegralSym(double points);
  // for ky integral, from 0 to inf
  void setKyIntegral(double points, double end);


private:
  int prefactor_;
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

  void setKxIntegral(double points);
  void setKxIntegralSym(double points);
  void setKyIntegral(double points);
  void setKyIntegralSym(double points);


private:
  int prefactor_;
};
}
#endif
