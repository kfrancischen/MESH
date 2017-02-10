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

namespace MESH{
using namespace SYSTEM;
using namespace RCWA;
typedef std::vector<RCWAMatrices> RCWAMatricesVec;

typedef struct ARGWEAPPER{
  double ky;
  RCWAMatrices* EMatrices;
  RCWAMatrices* grandImaginaryMatrices;
  RCWAMatrices* dielectricMatrixInverse;
  RCWAMatrix* Gx_mat;
  RCWAMatrix* Gy_mat;
  SourceList* sourceList;
  int targetLayer;
  int N = 0;
}ArgWrapper;

void filerLoader(std::string fileName, double* omega, dcomplex* epsilon, int size);
void saveData(std::string fileName, double* omega, double* fluxSpectrum, int size);
/*======================================================
Implementaion of the parent simulation super class
=======================================================*/
class Simulation{
public:
  Simulation();
  ~Simulation();
  void addStructure(Structure* structure);
  void setNumOfPoints(int numOfPoints);
  void setOmega(double* omegaList);
  void enableMPI(int numOfCore = 1);
  void setPeriod(double p1 = 0, double p2 = 0);
  void resetSimulation();

  Structure* getStructure();
  double* getOmegaList();
  double* getFluxSpectrum();
  double* getPeriod();

  virtual void initMatrices() = 0;
  virtual void setKxIntegral(double start, double points, double end) = 0;
  virtual void setKxIntegral(double start, double points) = 0;
  virtual void setKyIntegral(double start, double points) = 0;
  virtual void run() = 0;

private:
  int nGx_;
  int nGy_;
  int numOfCore_;
  int numOfPoints_;
  double period_[2];
  Structure* structure_;
  double* fluxSpectrum_;
  double* omegaList_;

  RCWAMatricesVec EMatricesVec;
  RCWAMatricesVec dielectricImMatrixVec;
  RCWAMatricesVec dielectricMatrixVec;
  RCWAMatricesVec dielectricMatrixInverseVec;
  RCWAMatrix Gx_mat;
  RCWAMatrix Gy_mat;
};


/*======================================================
Implementaion of the class on planar simulation
=======================================================*/
class SimulationPlanar : public Simulation{
public:
  SimulationPlanar() : Simulation(){};
  ~SimulationPlanar();

  void initMatrices();
  void setKxIntegral(double start, double points, double end);
  void setKxIntegral(double start, double points);
  void setKyIntegral(double start, double points);
  void run();
private:

};

/*======================================================
Implementaion of the class on 1D grating simulation
=======================================================*/
class SimulationGrating : public Simulation{
public:
  SimulationGrating() : Simulation(){};
  ~SimulationGrating();

  void initMatrices();
  void setKxIntegral(double start, double points, double end);
  void setKxIntegral(double start, double points);
  void setKyIntegral(double start, double points);
  void run();
};

/*======================================================
Implementaion of the class on 2D patterning simulation
=======================================================*/
class SimulationPattern : public Simulation{
public:
  SimulationPattern() : Simulation(){};
  ~SimulationPattern();

  void initMatrices();
  void setKxIntegral(double start, double points, double end);
  void setKxIntegral(double start, double points);
  void setKyIntegral(double start, double points);
  void run();
};
}
#endif
