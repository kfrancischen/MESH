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

#ifndef _INTERPOLATOR_H
#define _INTERPOLATOR_H
#include <string>
#include <vector>
#include <algorithm>
#include "Common.h"

namespace MESH{
class Interpolator{
public:
  Interpolator(const std::vector<double>& x, const std::vector< std::vector<double> >& y);
  std::vector<double> getVal(const double& x);
  Interpolator(const Interpolator&) = delete;
  ~Interpolator();

protected:
private:
  double* x_vals_;
  double** y_vals_;
  int size_x_, size_y_;
};
}

#endif