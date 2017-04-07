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

#ifndef _FMM_H
#define _FMM_H
#include <armadillo>
#include "Common.h"
#include "Rcwa.h"

namespace FMM{
  using namespace arma;
  using RCWA::RCWAcMatrix;
  using RCWA::RCWArMatrix;
  using RCWA::RCWAVector;
  using RCWA::sinc;
  using RCWA::jinc;
  using RCWA::meshGrid;

  /*==============================================*/
  // helper function to change a random dielectric to a tensor
  // @args:
  // epsilon: field
  // type: the type of the dielectric
  /*==============================================*/
  EpsilonVal toTensor(const EpsilonVal epsilon, const EPSTYPE type);

  /*==============================================*/
  // This function computes the Fourier transform for grating geometry
  // @args:
  // eps_xx: the Fourier trainsform for eps_xx
  // eps_xy: the Fourier trainsform for eps_xy
  // eps_zx: the Fourier trainsform for eps_yx
  // eps_yy: the Fourier trainsform for eps_yy
  // eps_zz_Inv: the inverse of Fourier transform of eps_zz
  // im_eps_xx: the Fourier trainsform for imaginary part
  // im_eps_xy: the Fourier trainsform for imaginary part
  // im_eps_yx: the Fourier trainsform for imaginary part
  // im_eps_yy: the Fourier trainsform for imaginary part
  // im_eps_zz: the Fourier trainsform for imaginary part
  // epsilonBGTensor: the epsilon of bacground (transformed to tensor already)
  // nGx: the total number of G
  // center: the center of the grating
  // width: the width of the grating
  // period: the periodicity
  // hasTensor: whether this layer contains tensor
  // userInverse: whether to use inverse rule
  /*==============================================*/
  void transformGrating(
    RCWAcMatrix& eps_xx,
    RCWAcMatrix& eps_xy,
    RCWAcMatrix& eps_yx,
    RCWAcMatrix& eps_yy,
    RCWAcMatrix& eps_zz_Inv,
    RCWAcMatrix& im_eps_xx,
    RCWAcMatrix& im_eps_xy,
    RCWAcMatrix& im_eps_yx,
    RCWAcMatrix& im_eps_yy,
    RCWAcMatrix& im_eps_zz,
    const EpsilonVal& epsilonBGTensor,
    const EpsilonVal& epsilon,
    const EPSTYPE epsilonType,
    const int nGx,
    const double center,
    const double width,
    const double period,
    const bool hasTensor,
    const bool useInverse
  );

  /*==============================================*/
  // This function computes the Fourier transform for rectangle geometry
  // @args:
  // eps_xx: the Fourier trainsform for eps_xx
  // eps_xy: the Fourier trainsform for eps_xy
  // eps_zx: the Fourier trainsform for eps_yx
  // eps_yy: the Fourier trainsform for eps_yy
  // eps_zz_Inv: the inverse of Fourier transform of eps_zz
  // im_eps_xx: the Fourier trainsform for imaginary part
  // im_eps_xy: the Fourier trainsform for imaginary part
  // im_eps_yx: the Fourier trainsform for imaginary part
  // im_eps_yy: the Fourier trainsform for imaginary part
  // im_eps_zz: the Fourier trainsform for imaginary part
  // epsilonBGTensor: the epsilon of bacground (transformed to tensor already)
  // nG_x: the total number of G in x direction
  // nG_y: the total number of G in y direction
  // centers: the centers of the rectangle
  // widths: the widths of the rectangle
  // period: the periodicity
  // hasTensor: whether this layer contains tensor
  // userInverse: whether to use inverse rule
  /*==============================================*/
  void transformRectangle(
    RCWAcMatrix& eps_xx,
    RCWAcMatrix& eps_xy,
    RCWAcMatrix& eps_yx,
    RCWAcMatrix& eps_yy,
    RCWAcMatrix& eps_zz_Inv,
    RCWAcMatrix& im_eps_xx,
    RCWAcMatrix& im_eps_xy,
    RCWAcMatrix& im_eps_yx,
    RCWAcMatrix& im_eps_yy,
    RCWAcMatrix& im_eps_zz,
    const EpsilonVal& epsBGTensor,
    const EpsilonVal& epsilon,
    const EPSTYPE epsilonType,
    const int nGx,
    const int nGy,
    const double centers[2],
    const double widths[2],
    const double period[2],
    const bool hasTensor,
    const bool useInverse
  );

  /*==============================================*/
  // This function computes the Fourier transform for circle geometry
  // @args:
  // eps_xx: the Fourier trainsform for eps_xx
  // eps_xy: the Fourier trainsform for eps_xy
  // eps_zx: the Fourier trainsform for eps_yx
  // eps_yy: the Fourier trainsform for eps_yy
  // eps_zz_Inv: the inverse of Fourier transform of eps_zz
  // im_eps_xx: the Fourier trainsform for imaginary part
  // im_eps_xy: the Fourier trainsform for imaginary part
  // im_eps_yx: the Fourier trainsform for imaginary part
  // im_eps_yy: the Fourier trainsform for imaginary part
  // im_eps_zz: the Fourier trainsform for imaginary part
  // epsilonBGTensor: the epsilon of bacground (transformed to tensor already)
  // nG_x: the total number of G in x direction
  // nG_y: the total number of G in y direction
  // centers: the centers of the circle
  // radius: the radius of the circle
  // period: the periodicity
  // hasTensor: whether this layer contains tensor
  // userInverse: whether to use inverse rule
  /*==============================================*/
  void transformCircle(
    RCWAcMatrix& eps_xx,
    RCWAcMatrix& eps_xy,
    RCWAcMatrix& eps_yx,
    RCWAcMatrix& eps_yy,
    RCWAcMatrix& eps_zz_Inv,
    RCWAcMatrix& im_eps_xx,
    RCWAcMatrix& im_eps_xy,
    RCWAcMatrix& im_eps_yx,
    RCWAcMatrix& im_eps_yy,
    RCWAcMatrix& im_eps_zz,
    const EpsilonVal& epsBGTensor,
    const EpsilonVal& epsilon,
    const EPSTYPE epsilonType,
    const int nGx,
    const int nGy,
    const double centers[2],
    const double radius,
    const double period[2],
    const bool hasTensor,
    const bool useInverse
  );


}

#endif