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

#include <complex>
#include "utility/Utility.h"
#ifndef _COMMON_H
#define _COMMON_H

#define DEGREE 512

using UTILITY::Ptr;
using UTILITY::PtrInterface;
enum DIMENSION { NO_, ONE_, TWO_ };
enum SOURCE {ISSOURCE_, ISNOTSOURCE_};
enum PATTEN {PLANAR_, GRATING_, RECTANGLE_, CIRCLE_};
typedef std::complex<double> dcomplex;
#define POW2(x) pow(x, 2)
#define POW3(x) pow(x, 3)
#define SENDTAG 2
#define RECVTAG 1
#define MASTER 0

typedef std::vector< SOURCE > SourceList;
#endif