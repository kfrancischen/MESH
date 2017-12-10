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

#include "Interpolator.h"

// constructor
MESH::Interpolator::Interpolator(const std::vector<double>& x, const std::vector< std::vector<double> >& y){
  if(!std::is_sorted(x.begin(), x.end())){
    std::cerr << "input value should be sorted in accending order!" << std::endl;
    throw UTILITY::InternalException("input value should be sorted in accending order!");
  }
  
  size_x_ = x.size();
  if(size_x_ <= 1){
    std::cerr << "Need at least two points!" << std::endl;
    throw UTILITY::ValueException("Need at least two points!");
  }
  size_y_ = y[0].size();

  x_vals_ = new double[size_x_];
  y_vals_ = new double*[size_x_];
  for(int i = 0; i < size_x_; i++){
    x_vals_[i] = x[i];
    y_vals_[i] = new double[size_y_];
    for (int j = 0; j < size_y_; j++){
      y_vals_[i][j] = y[i][j];
    }
  }
}

// destructor
MESH::Interpolator::~Interpolator(){
  if(x_vals_ != nullptr){
    delete [] x_vals_;
    x_vals_ = nullptr;
  }

  if(y_vals_!= nullptr){
    for (int i = 0; i < size_x_; i++){
      delete [] y_vals_[i];
    }
    delete [] y_vals_;
    y_vals_ = nullptr;
  }
}


// do the interpolation
std::vector<double> MESH::Interpolator::getVal(const double& x){
  if(x < x_vals_[0] || x > x_vals_[size_x_ - 1]){
    std::cerr << "Input value out of range!" << std::endl;
    throw UTILITY::RangeException("Input value out of range!");
  }
  std::vector<double> result;

  for(int i = 1; i < size_x_; i++){
    if(x <= x_vals_[i]){
      double ratio = 0;
      if(x_vals_[i] != x_vals_[i-1]){
        ratio = (x - x_vals_[i-1]) / (x_vals_[i] - x_vals_[i-1]);
      }
      for(int j = 0; j < size_y_; j++){
        result.push_back( ratio * y_vals_[i][j] + (1-ratio) * y_vals_[i-1][j] );
      }
      break;
    }
  }

  return result;
}