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
//#ifndef ARMA_DONT_USE_WRAPPER
//#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <vector>
#include "Common.h"

namespace RCWA{
  enum DIRECTION {UP_, DOWN_, ALL_};
  using namespace arma;
  typedef cx_mat RCWAMatrix;
  typedef std::vector< RCWAMatrix > RCWAMatrices;
  typedef std::vector< RCWAMatrices > RCWAMatricesVec;
  typedef vec RCWAVector;

  /*============================================================
  * Function similar to meshgrid in matlab
  @arg:
   vL, vR: the input two vectors (here everything is a matrix)
   qL, qR: the input two vectors
    [qL, qR] = meshgrid(vl, vR)
  ==============================================================*/
  void meshGrid(
    const RCWAMatrix& vL,
    const RCWAMatrix& vR,
    RCWAMatrix& qL,
    RCWAMatrix& qR
  );
  /*============================================================
  * Function computing S matrix for each layers
  @arg:
   startLayer: the starting layer for the propogation
   N: number of G total G
   numOfLayer: the number of layers
   MMatrices: matrix corresponding to the propogation in each layer
   FMatrices: matrix corresponding to the phase in each layer
   SMatrices: the S matrix for each layer
  ==============================================================*/
  void getSMatrices(
    const int startLayer,
    const int N,
    const int numOfLayer,
    const RCWAMatrices& MMatrices,
    const RCWAMatrices& FMatrices,
    RCWAMatrices& SMatrices,
    const DIRECTION direction
  );

  /*============================================================
  * Function computing numbef of G for the system
  @arg:
   nGx: positive G along x direction
   nGy: positive G along y direction
  ==============================================================*/
 int getN(
   const int nGx,
   const int nGy
 );

 /*============================================================
 * Function computing the sinc function (sin(x) / x)
 @arg:
  x: the input argument
 ==============================================================*/
 RCWAMatrix sinc(const RCWAMatrix x);
 /*============================================================
 * Function computing the jinc function (J1(x) / x)
 @arg:
  x: the input argument
 ==============================================================*/
 double jinc(const double x);
  /*============================================================
  * Function computing G matrix for the system
  @arg:
   startLayer: the starting layer for the propogation
   nGx: positive G along x direction
   nGy: positive G along y direction
   period: periodicity along two directions
   Gx_mat: the Gx matrix
   Gy_mat: the Gy matrix
   d: the dimension of the structure
  ==============================================================*/
  void getGMatrices(
    const int nGx,
    const int nGy,
    const double period[2],
    RCWAMatrix& Gx_mat,
    RCWAMatrix& Gy_mat,
    const DIMENSION d
  );

  /*============================================================
  * Function computing imaginary dielectric matrix for the system
  @arg:
  grandImaginaryMatrices: the matrices containing the imaginary part matrix
  im_eps_xx: the imaginary part of eps_xx
  im_eps_xy: the imaginary part of eps_xy
  im_eps_yx: the imaginary part of eps_yx
  im_eps_yy: the imaginary part of eps_yy
  im_eps_zz: the imaginary part of eps_zz
  numOfLayer: the number of layer in the system
  N: the number of G
  ==============================================================*/
  void getGrandImaginaryMatrices(
    RCWAMatrices& grandImaginaryMatrices,
    const RCWAMatrices& im_eps_xx,
    const RCWAMatrices& im_eps_xy,
    const RCWAMatrices& im_eps_yx,
    const RCWAMatrices& im_eps_yy,
    const RCWAMatrices& im_eps_zz,
    int numOfLayer,
    int N
  );

  /*============================================================
  * Function computing E matrices for the system
  @arg:
  EMatrices: the Ematrices
  eps_xx: the epsilon in xx direction
  eps_xy: the epsilon in xy direction
  eps_yx: the epsilon in yx direction
  eps_yy: the epsilon in yy direction
  numOfLayer: the number of layer in the system
  N: the number of G
  ==============================================================*/
  void getEMatrices(
    RCWAMatrices& EMatrices,
    const RCWAMatrices& eps_xx,
    const RCWAMatrices& eps_xy,
    const RCWAMatrices& eps_yx,
    const RCWAMatrices& eps_yy,
    const int numOfLayer,
    const int N
  );

  /*============================================================
  * Function computing the poynting vector at given (kx, ky)
  @arg:
   omega: the angular frequency (normalized to c)
   thicknessList: the thickness for each layer
   kx: the k vector at x direction (normalized value)
   ky: the y vector at x direction (normalized value)
   EMatrices:  the E matrices for all layers
   grandImaginaryMatrices: collection of all imaginary matrices in all layers
   eps_zz_inv: the inverse of eps_zz
   Gx_mat: the Gx matrix
   Gy_mat: the Gy matrix
   sourceList: list of 0 or 1 with the same size of thicknessList
   targetLayer: the targetLayer for the flux measurement
   N: total number of G
  ==============================================================*/
  double poyntingFlux(
    const double omega,
    const RCWAVector& thicknessList,
    double kx,
    double ky,
    const RCWAMatrices& EMatrices,
    const RCWAMatrices& grandImaginaryMatrices,
    const RCWAMatrices& eps_zz_inv,
    const RCWAMatrix& Gx_mat,
    const RCWAMatrix& Gy_mat,
    const SourceList& sourceList,
    const int targetLayer,
    const int N
  );
}
#endif
