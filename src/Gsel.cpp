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


#include "Gsel.h"
#ifndef DBL_EPSILON
  #define DBL_EPSILON 2.2204460492503131e-16
#endif
namespace GSEL{
  /*============================================================
  // this part of code is copied from S4
  ============================================================*/

  static double Gcmp_d(const void *i, const void *j, void *arg){
  	// Comparing lengths of two vectors:
  	//  Vector u is defined by
  	//   u[0] = G[2*i+0] * Lk[0] + G[2*i+1] * Lk[2]
  	//   u[1] = G[2*i+0] * Lk[1] + G[2*i+1] * Lk[3]
  	//  Similarly for v and index j.
  	//  The length is then:
  	//   |u|^2 = w^T [ G[2*i+0]*G[2*i+0] ]
  	//               [ G[2*i+0]*G[2*i+1] ]
  	//               [ G[2*i+1]*G[2*i+1] ]
  	//  where
  	//   w = [    Lk[0]*Lk[0]+Lk[1]*Lk[1]  ]
  	//       [ 2*(Lk[0]*Lk[2]+Lk[1]*Lk[3]) ]
  	//       [    Lk[2]*Lk[2]+Lk[3]*Lk[3]  ]
  	const int *Gi = (int*)i;
  	const int *Gj = (int*)j;

  	const double *w = (const double*)arg;
  	int dv[3] = {
  		Gi[0]*Gi[0] - Gj[0]*Gj[0],
  		Gi[0]*Gi[1] - Gj[0]*Gj[1],
  		Gi[1]*Gi[1] - Gj[1]*Gj[1]
  	};
  	return dv[0]*w[0] + dv[1]*w[1] + dv[2]*w[2];
  }

  static int Gcmp(const void *i, const void *j, void *arg){
  	double d = Gcmp_d(i, j, arg);
  	if(d > 0){ return 1; }
  	if(d < 0){ return -1; }
  	return 0;
  }

  static bool G_same(const int Gi[2], const int Gj[2], const double Lk[4]){
  	double Lkprod[3] = {
  		Lk[0]*Lk[0]+Lk[1]*Lk[1],
  		2.*(Lk[0]*Lk[2]+Lk[1]*Lk[3]),
  		Lk[2]*Lk[2]+Lk[3]*Lk[3]
  	};
  	const double ilen = hypot(
  		Gi[0]*Lk[0]+Gi[1]*Lk[2],
  		Gi[0]*Lk[1]+Gi[1]*Lk[3]
  	);
  	const double jlen = hypot(
  		Gj[0]*Lk[0]+Gj[1]*Lk[2],
  		Gj[0]*Lk[1]+Gj[1]*Lk[3]
  	);
  	const double maxlen = (ilen > jlen) ? ilen : jlen;
  	double d = fabs(Gcmp_d(&Gi[0], &Gj[0], &Lkprod[0]));
  	return d < 2*DBL_EPSILON*maxlen;
  }

  /*============================================================
  * Function computing G matrix for the system using circular truncation
  @args:
  nG: the user given nG, will be changed
  reciprocalLattice: the reciprocal lattice
  Gx_mat: the G matrix in x direction
  Gy_mat: the G matrix in y direction
  ==============================================================*/
  static void GSelCircular(
    int& nG,
    const Lattice& reciprocalLattice,
    RCWArMatrix& Gx_mat,
    RCWArMatrix& Gy_mat
  ){

  }
  /*============================================================
  * Function computing G matrix for the system using parallelogramic truncation
  @args:
  nG: the user given nG, will be changed
  reciprocalLattice: the reciprocal lattice
  Gx_mat: the G matrix in x direction
  Gy_mat: the G matrix in y direction
  option: flag differentiate 1D and 2D
  ==============================================================*/
  static void GSelParallelogramic(
    int& nG,
    const Lattice& reciprocalLattice,
    RCWArMatrix& Gx_mat,
    RCWArMatrix& Gy_mat,
    int option
  ){
    // case of 2D
    if(option == 1){
      int NRoot = (int)std::sqrt(nG);
      if(NRoot % 2 == 0 && NRoot > 0) NRoot--;
      int M = NRoot / 2;

      RCWArMatrix G_list(NRoot, 1);
      for(int i = -M; i <= M; i++){
        G_list(i+M, 0) = i;
      }
      RCWA::meshGrid(G_list, G_list, Gx_mat, Gy_mat);
      nG = POW2(NRoot);
    }
    // case of 1D
    else{
      int M = nG / 2;
      nG = 2 * M + 1;
      RCWArMatrix Gx_list(nG, 1), Gy_list(1, 1);
      for(int i = -M; i <= M; i++){
        Gx_list(i+M, 0) = i;
      }
      Gy_list(0, 0) = 0;
      RCWA::meshGrid(Gx_list, Gy_list, Gx_mat, Gy_mat);
    }

    RCWArMatrix GxMat_temp = Gx_mat * reciprocalLattice.bx[0] + Gy_mat * reciprocalLattice.by[0];
    Gy_mat = Gx_mat * reciprocalLattice.bx[1] + Gy_mat * reciprocalLattice.by[1];
    Gx_mat = GxMat_temp;
  }

  /*============================================================
  * Function computing G matrix for the system
  @args:
  nG: the user given nG, will be changed
  reciprocalLattice: the reciprocal lattice
  Gx_mat: the G matrix in x direction
  Gy_mat: the G matrix in y direction
  d: dimension of the problem
  truncation: the option of truncation
  ==============================================================*/
  void getGMatrices(
    int& nG,
    const Lattice& reciprocalLattice,
    RCWArMatrix& Gx_mat,
    RCWArMatrix& Gy_mat,
    const DIMENSION d,
    const TRUNCATION truncation
  ){
    if(nG <= 0){
      std::cerr << "Need G value more than 1!" << std::endl;
      throw UTILITY::ValueException("Need G value more than 1!");
    }
    if(d == NO_){
      nG = 1;
      GSelParallelogramic(nG, reciprocalLattice, Gx_mat, Gy_mat, 1);
    }
    else if(d == ONE_){
      GSelParallelogramic(nG, reciprocalLattice, Gx_mat, Gy_mat, 0);
    }
    else{
      switch (truncation) {
        case CIRCULAR_:{
          GSelCircular(nG, reciprocalLattice, Gx_mat, Gy_mat);
          break;
        }
        case PARALLELOGRAMIC_:{
          GSelParallelogramic(nG, reciprocalLattice, Gx_mat, Gy_mat, 1);
          break;
        }
        default: break;
      }
    }

  }
}


