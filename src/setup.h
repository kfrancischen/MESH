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
#ifndef _SETUP_H
#define _SETUP_H

#include "Rcwa.h"
#include "Gsel.h"
#include "Cubature.h"
#include "System.h"
#include "Fmm.h"
#include "Mesh.h"

using namespace MESH;

typedef struct CONSTANT{
  double pi = datum::pi;
  double k_B = datum::k;
  double eps_0 = datum::eps_0;
  double m_e = datum::m_e;
  double eV = datum::eV;
  double mu_0 = datum::mu_0;
  double h = datum::h;
  double h_bar = datum::h_bar;
  double c_0 = datum::c_0;
  double q = datum::ec;
  double sigma = datum::sigma;
} Constant;

const Constant constants;
#endif