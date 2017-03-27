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

#ifndef _COMMON_H
#define _COMMON_H

#include <complex>
#include <vector>
#include "utility/Utility.h"
#include "mathlib/Mathlib.h"
#include "config.h"
#define DEGREE 1024

using UTILITY::Ptr;
using UTILITY::PtrInterface;
using UTILITY::NamedInterface;

enum DIMENSION { NO_, ONE_, TWO_ };
enum PATTERN {GRATING_, RECTANGLE_, CIRCLE_};
enum EPSTYPE {SCALAR_, DIAGONAL_, TENSOR_};
#define POW2(x) pow(x, 2)
#define POW3(x) pow(x, 3)
#define SENDTAG 2
#define RECVTAG 1
#define MASTER 0
#define ABSERROR 1e-10
#define RELERROR 1e-10

typedef std::complex<double> dcomplex;
const dcomplex IMAG_I = dcomplex(0.0, 1.0);
const dcomplex REAL_I = dcomplex(1.0, 0.0);
union EpsilonVal{
  double scalar[2];
  double diagonal[6];
  double tensor[10];
};

struct EPSILON
{
  EpsilonVal* epsilonVals;
  ~EPSILON(){
    //delete [] epsilonVals;
  }
  EPSTYPE type_;
};

typedef std::vector<bool> SourceList;

#endif