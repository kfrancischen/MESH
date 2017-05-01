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
#include <cmath>
#include "Common.h"

namespace RCWA{
  enum DIRECTION {UP_, DOWN_, ALL_};
  using namespace arma;
  typedef cx_mat RCWAcMatrix;
  typedef mat RCWArMatrix;
  typedef std::vector< RCWAcMatrix > RCWAcMatrices;
  typedef std::vector< RCWAcMatrices > RCWAcMatricesVec;
  typedef vec RCWArVector;

  /*============================================================
  * Function similar to meshgrid in matlab for real numbers
  @arg:
   vL, vR: the input two vectors (here everything is a matrix)
   qL, qR: the input two vectors
    [qL, qR] = meshgrid(vl, vR)
  @note:
    both vL and vR should be a colum vector
  ==============================================================*/
  void meshGrid(
    const RCWArMatrix& vL,
    const RCWArMatrix& vR,
    RCWArMatrix& qL,
    RCWArMatrix& qR
  );
  /*============================================================
  * Function similar to meshgrid in matlab for complex numbers
  @arg:
   vL, vR: the input two vectors (here everything is a matrix)
   qL, qR: the input two vectors
    [qL, qR] = meshgrid(vl, vR)
  @note:
    both vL and vR should be a colum vector
  ==============================================================*/
  void meshGrid(
    const RCWAcMatrix& vL,
    const RCWAcMatrix& vR,
    RCWAcMatrix& qL,
    RCWAcMatrix& qR
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
    const RCWAcMatrices& MMatrices,
    const RCWAcMatrices& FMatrices,
    RCWAcMatrices& SMatrices,
    const DIRECTION direction
  );


 /*============================================================
 * Function computing the sinc function (sin(x) / x) for matrix x
 @arg:
  x: the input argument
 ==============================================================*/
 RCWArMatrix sinc(const RCWArMatrix x);
 /*============================================================
 * Function computing the sinc function (sin(x) / x) for double input x
 @arg:
  x: the input argument
 ==============================================================*/
 double sinc(const double x);
 /*============================================================
 * Function computing the jinc function (J1(x) / x)
 @arg:
  x: the input argument
 ==============================================================*/
 double jinc(const double x);

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
    RCWAcMatrices& grandImaginaryMatrices,
    const RCWAcMatrices& im_eps_xx,
    const RCWAcMatrices& im_eps_xy,
    const RCWAcMatrices& im_eps_yx,
    const RCWAcMatrices& im_eps_yy,
    const RCWAcMatrices& im_eps_zz,
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
    RCWAcMatrices& EMatrices,
    const RCWAcMatrices& eps_xx,
    const RCWAcMatrices& eps_xy,
    const RCWAcMatrices& eps_yx,
    const RCWAcMatrices& eps_yy,
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
    const RCWArVector& thicknessList,
    double kx,
    double ky,
    const RCWAcMatrices& EMatrices,
    const RCWAcMatrices& grandImaginaryMatrices,
    const RCWAcMatrices& eps_zz_inv,
    const RCWArMatrix& Gx_mat,
    const RCWArMatrix& Gy_mat,
    const SourceList& sourceList,
    const int targetLayer,
    const int N,
    const POLARIZATION polar
  );

}
#endif
