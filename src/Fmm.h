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
#include "Rcwa.h"
#include "System.h"

namespace FMM{
  using RCWA::RCWAMatrix;
  using RCWA::RCWAMatrices;
  using RCWA::RCWAVector;
  using SYSTEM::Layer;
  using SYSTEM::Material;
  typedef std::vector<RCWAMatrices> RCWAMatricesVec;


  /*==============================================*/
  // This function computes the Fourier transform for planar geometry for tensor case
  // @args:
  // eps_xx_MatrixVec: the Fourier trainsform for eps_xx for all omega
  // eps_xy_MatrixVec: the Fourier trainsform for eps_xy for all omega
  // eps_zx_MatrixVec: the Fourier trainsform for eps_yx for all omega
  // eps_yy_MatrixVec: the Fourier trainsform for eps_yy for all omega
  // im_eps_xx_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // im_eps_xy_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // im_eps_yx_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // im_eps_yy_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // im_eps_zz_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // eps_zz_Inv_MatrixVec: the inverse of Fourier transform of eps_zz
  // Layer: the layer considered
  // N: the total number of G
  /*==============================================*/
  void transformPlanar(
    RCWAMatricesVec& eps_xx_MatrixVec,
    RCWAMatricesVec& eps_xy_MatrixVec,
    RCWAMatricesVec& eps_yx_MatrixVec,
    RCWAMatricesVec& eps_yy_MatrixVec,
    RCWAMatricesVec& eps_zz_Inv_MatrixVec,
    RCWAMatricesVec& im_eps_xx_MatrixVec,
    RCWAMatricesVec& im_eps_xy_MatrixVec,
    RCWAMatricesVec& im_eps_yx_MatrixVec,
    RCWAMatricesVec& im_eps_yy_MatrixVec,
    RCWAMatricesVec& im_eps_zz_MatrixVec,
    Ptr<Layer>& layer,
    const int N
  );


  /*==============================================*/
  // This function computes the Fourier transform for planar geometry for tensor case
  // @args:
  // eps_xx_MatrixVec: the Fourier trainsform for eps_xx for all omega
  // eps_xy_MatrixVec: the Fourier trainsform for eps_xy for all omega
  // eps_zx_MatrixVec: the Fourier trainsform for eps_yx for all omega
  // eps_yy_MatrixVec: the Fourier trainsform for eps_yy for all omega
  // im_eps_xx_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // im_eps_xy_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // im_eps_yx_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // im_eps_yy_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // im_eps_zz_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // eps_zz_Inv_MatrixVec: the inverse of Fourier transform of eps_zz
  // Layer: the layer considered
  // N: the total number of G
  // period: the periodicity
  /*==============================================*/
  void transformGratingNaive(
    RCWAMatricesVec& eps_xx_MatrixVec,
    RCWAMatricesVec& eps_xy_MatrixVec,
    RCWAMatricesVec& eps_yx_MatrixVec,
    RCWAMatricesVec& eps_yy_MatrixVec,
    RCWAMatricesVec& eps_zz_Inv_MatrixVec,
    RCWAMatricesVec& im_eps_xx_MatrixVec,
    RCWAMatricesVec& im_eps_xy_MatrixVec,
    RCWAMatricesVec& im_eps_yx_MatrixVec,
    RCWAMatricesVec& im_eps_yy_MatrixVec,
    RCWAMatricesVec& im_eps_zz_MatrixVec,
    Ptr<Layer>& layer,
    const int N,
    const double period,
    bool useInverseRule = true
  );

  void transformGratingDiagonalAdaptive();

  void transformGratingTensorAdaptive();

  /*==============================================*/
  // This function computes the Fourier transform for planar geometry for tensor case
  // @args:
  // eps_xx_MatrixVec: the Fourier trainsform for eps_xx for all omega
  // eps_xy_MatrixVec: the Fourier trainsform for eps_xy for all omega
  // eps_zx_MatrixVec: the Fourier trainsform for eps_yx for all omega
  // eps_yy_MatrixVec: the Fourier trainsform for eps_yy for all omega
  // im_eps_xx_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // im_eps_xy_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // im_eps_yx_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // im_eps_yy_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // im_eps_zz_MatrixVec: the Fourier trainsform for imaginary part for all omega
  // eps_zz_Inv_MatrixVec: the inverse of Fourier transform of eps_zz
  // Layer: the layer considered
  // N: the total number of G
  // period: the periodicity
  // useInverseRule: whether use inverse rule
  /*==============================================*/
  void transformRectangle(
    RCWAMatricesVec& eps_xx_MatrixVec,
    RCWAMatricesVec& eps_xy_MatrixVec,
    RCWAMatricesVec& eps_yx_MatrixVec,
    RCWAMatricesVec& eps_yy_MatrixVec,
    RCWAMatricesVec& eps_zz_Inv_MatrixVec,
    RCWAMatricesVec& im_eps_xx_MatrixVec,
    RCWAMatricesVec& im_eps_xy_MatrixVec,
    RCWAMatricesVec& im_eps_yx_MatrixVec,
    RCWAMatricesVec& im_eps_yy_MatrixVec,
    RCWAMatricesVec& im_eps_zz_MatrixVec,
    Ptr<Layer>& layer,
    const int N,
    const double* period,
    bool useInverseRule = true
  );

  void transformCircleDiagonal(
    RCWAMatricesVec& eps_xx_MatrixVec,
    RCWAMatricesVec& eps_yy_MatrixVec,
    RCWAMatricesVec& eps_zz_Inv_MatrixVec,
    RCWAMatricesVec& im_eps_xx_MatrixVec,
    RCWAMatricesVec& im_eps_yy_MatrixVec,
    RCWAMatricesVec& im_eps_zz_MatrixVec,
    Ptr<Layer>& layer,
    const int N,
    const double* period
  );

  void transformCircleTensor(
    RCWAMatricesVec& eps_xx_MatrixVec,
    RCWAMatricesVec& eps_xy_MatrixVec,
    RCWAMatricesVec& eps_yx_MatrixVec,
    RCWAMatricesVec& eps_yy_MatrixVec,
    RCWAMatricesVec& eps_zz_Inv_MatrixVec,
    RCWAMatricesVec& im_eps_xx_MatrixVec,
    RCWAMatricesVec& im_eps_xy_MatrixVec,
    RCWAMatricesVec& im_eps_yx_MatrixVec,
    RCWAMatricesVec& im_eps_yy_MatrixVec,
    RCWAMatricesVec& im_eps_zz_MatrixVec,
    Ptr<Layer>& layer,
    const int N,
    const double* period
  );

}

#endif