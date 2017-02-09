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
#ifndef _RCWA_H
#define _RCWA_H
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <complex>
#include <vector>
#include "Common.h"

namespace RCWA{

  using namespace arma;
  typedef cx_mat RCWAMatrix;
  typedef std::vector< RCWAMatrix > RCWAMatrices;
  typedef fvec RCWAVector;
  typedef std::vector< LOSSY > SourceList;



  /*============================================================
  * Function initializing S matrix
  @arg:
   n: size of S matrix
   S: the S matrix, initialized to be an identity matrix with size nxn
  ==============================================================*/
  void initSMatrix(const int n, RCWAMatrix* S);

  /*============================================================
  * Function computing S matrix for each layers
  @arg:
   startLayer: the starting layer for the propogation
   Nx: number of G in x direction
   Ny: number of G in y direction
   thicknessList: the thickness for each layer
   MMatrices: matrix corresponding to the propogation in each layer
   FMatrices: matrix corresponding to the phase in each layer
  ==============================================================*/
  RCWAMatrices getSMatrix(
    const int startLayer,
    const int Nx,
    const int Ny,
    const RCWAVector thicknessList,
    RCWAMatrices MMatrices,
    RCWAMatrices FMatrices
  );

  /*============================================================
  * Function similar to meshgrid in matlab
  @arg:
   vL, vR: the input two vectors (here everything is a matrix)
   qL, qR: the input two vectors
    [qL, qR] = meshgrid(vl, vR)
  ==============================================================*/
  void populateQ(
    RCWAMatrix vL,
    RCWAMatrix vR,
    RCWAMatrix* qL,
    RCWAMatrix* qR
  );

  /*============================================================
  * Function computing the poynting vector at given (kx, ky)
  @arg:
   omega: the angular frequency (normalized to c)
   thicknessList: the thickness for each layer
   kx: the k vector at x direction (absolute value)
   ky: the y vector at x direction (absolute value)
   dielectricMatrixInverse: the inverse of dielectric matrix
   dielectricMatrix: the dielectric after Fourier transform
   dielectricImMatrix: the imaginary part of the dielectric function after Fourier transform
   sourceList: list of 0 or 1 with the same size of thicknessList
   targetLayer: the targetLayer for the flux measurement
   nGx: number of positive G in x direction
   nGy: number of positive G in y direction
   period: the periodicity of the system (varies for different dimensions)
   d: the dimension of the grating
  ==============================================================*/
  double poyntingFlux(
    const double omega,
    const RCWAVector thicknessList,
    double kx,
    double ky,
    const RCWAMatrices dielectricMatrixInverse,
    const RCWAMatrices dielectricMatrix,
    const RCWAMatrices dielectricImMatrix,
    const SourceList sourceList,
    const int targetLayer,
    const int nGx,
    const int nGy,
    const double *period,
    const DIMENSION d
  );
}
#endif
