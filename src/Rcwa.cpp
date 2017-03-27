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
#include "Rcwa.h"

/*============================================================
* Function similar to meshgrid in matlab
@arg:
 vL, vR: the input two vectors (here everything is a matrix)
 qL, qR: the input two vectors
 [qL, qR] = meshgrid(vl, vR)
==============================================================*/
void RCWA::meshGrid(
  const RCWAMatrix& vL,
  const RCWAMatrix& vR,
  RCWAMatrix& qL,
  RCWAMatrix& qR
){

  qL = repmat(vL.st(), vR.n_rows, 1);
  qR = repmat(vR, 1, vL.n_rows);
}

/*============================================================
* Function computing S matrix for each layers
@arg:
 startLayer: the starting layer for the propogation
 N: number of G total G
 numOfLayer: the number of layers
 MMatrices: matrix corresponding to the propogation in each layer
 FMatrices: matrix corresponding to the phase in each layer
 SMatrices: the S matrix for each layer
 DIRECTION: the direction of propogation
==============================================================*/
void RCWA::getSMatrices(
  const int startLayer,
  const int N,
  const int numOfLayer,
  const RCWAMatrices& MMatrices,
  const RCWAMatrices& FMatrices,
  RCWAMatrices& SMatrices,
  const DIRECTION direction
){

  int r1 = 0, r2 = 2*N -1, r3 = 2*N, r4 = 4*N -1;
// propogating down
  if(direction == ALL_ || direction == DOWN_){
    // propogating down
    for(int i = startLayer; i >=1; i--){
      RCWAMatrix I = solve(MMatrices[i], MMatrices[i-1], solve_opts::fast);

      RCWAMatrix leftTop = I(span(r3, r4), span(r3, r4));
      RCWAMatrix rightTop = I(span(r3, r4), span(r1, r2));
      RCWAMatrix leftBottom = I(span(r1, r2), span(r3, r4));
      RCWAMatrix rightBottom = I(span(r1, r2), span(r1, r2));

      SMatrices[i-1](span(r1, r2), span(r1, r2)) = solve(
        leftTop - FMatrices[i] * SMatrices[i](span(r1, r2), span(r3, r4)) * leftBottom,
        FMatrices[i],
        solve_opts::fast
      ) * SMatrices[i](span(r1, r2), span(r1,r2));

      RCWAMatrix test = leftTop - FMatrices[i] * SMatrices[i](span(r1, r2), span(r3, r4)) * leftBottom;

      SMatrices[i-1](span(r1, r2), span(r3, r4)) = solve(
        leftTop - FMatrices[i] * SMatrices[i](span(r1, r2), span(r3, r4)) * leftBottom,
        FMatrices[i] * SMatrices[i](span(r1, r2), span(r3,r4)) * rightBottom - rightTop,
        solve_opts::fast
      ) * FMatrices[i-1];

      SMatrices[i-1](span(r3,r4), span(r1,r2)) = SMatrices[i](span(r3,r4), span(r1,r2)) +
        SMatrices[i](span(r3,r4), span(r3,r4)) * leftBottom * SMatrices[i-1](span(r1,r2), span(r1,r2));

      SMatrices[i-1](span(r3,r4), span(r3,r4)) = SMatrices[i](span(r3,r4), span(r3,r4)) *
        (leftBottom * SMatrices[i-1](span(r1,r2), span(r3,r4)) + rightBottom * FMatrices[i-1]);
    }
  }
  if(direction == ALL_ || direction == UP_){
    // propogating up
    for(int i = startLayer; i < numOfLayer - 1; i++){
      RCWAMatrix I = solve(MMatrices[i], MMatrices[i+1], solve_opts::fast);
      RCWAMatrix leftTop = I(span(r1, r2), span(r1, r2));
      RCWAMatrix rightTop = I(span(r1, r2), span(r3, r4));
      RCWAMatrix leftBottom = I(span(r3, r4), span(r1, r2));
      RCWAMatrix rightBottom = I(span(r3, r4), span(r3, r4));

      SMatrices[i+1](span(r1, r2), span(r1, r2)) = solve(
        leftTop - FMatrices[i] * SMatrices[i](span(r1, r2), span(r3, r4)) * leftBottom,
        FMatrices[i],
        solve_opts::fast
      ) * SMatrices[i](span(r1, r2), span(r1,r2));

      SMatrices[i+1](span(r1, r2), span(r3, r4)) = solve(
        leftTop - FMatrices[i] * SMatrices[i](span(r1, r2), span(r3, r4)) * leftBottom,
        FMatrices[i] * SMatrices[i](span(r1, r2), span(r3,r4)) * rightBottom - rightTop,
        solve_opts::fast
      ) * FMatrices[i+1];

      SMatrices[i+1](span(r3,r4), span(r1,r2)) = SMatrices[i](span(r3,r4), span(r1,r2)) +
        SMatrices[i](span(r3,r4), span(r3,r4)) * leftBottom * SMatrices[i+1](span(r1,r2), span(r1,r2));

      SMatrices[i+1](span(r3,r4), span(r3,r4)) = SMatrices[i](span(r3,r4), span(r3,r4)) *
        (leftBottom * SMatrices[i+1](span(r1,r2), span(r3,r4)) + rightBottom * FMatrices[i+1]);
    }
  }

}

/*============================================================
* Function computing numbef of G for the system
@arg:
 nGx: positive G along x direction
 nGy: positive G along y direction
==============================================================*/
int RCWA::getN(
 const int nGx,
 const int nGy
){
  return (2*nGx + 1) * (2*nGy + 1);
}
/*============================================================
* Function computing the sinc function (sin(x) / x)
@arg:
 x: the input argument
==============================================================*/
RCWA::RCWAMatrix RCWA::sinc(const RCWAMatrix x){
  RCWAMatrix output = sin(x) / x;
  output.elem( find(x == 0.0) ).ones();
  return output;
}
/*============================================================
* Function computing the jinc function (J1(x) / x)
@arg:
 x: the input argument
==============================================================*/
double jinc(const double x){
  if(x == 0.0) return 1;
  double j1 = 0;
  double j0, y0, y1, j0p, j1p, y0p, y1p;
  bessjy01a(x, j0, j1, y0, y1, j0p, j1p, y0p, y1p);
  return j1 / x;
}
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
void RCWA::getGMatrices(
  const int nGx,
  const int nGy,
  const double period[2],
  RCWAMatrix& Gx_mat,
  RCWAMatrix& Gy_mat,
  const DIMENSION d
)
{

  int Nx = 2 * nGx + 1;
  int Ny = 2 * nGy + 1;
  int N = Nx * Ny;
  int Gx = 0, Gy = 0;
  switch (d) {
    case NO_:
      break;
    case ONE_:{
      Gx = 2 * datum::pi / period[0];
      break;
    }
    case TWO_:{
      Gx = 2 * datum::pi / period[0];
      Gy = 2 * datum::pi / period[1];
      break;
    }
  }
  RCWAMatrix Gx_list(1, Nx), Gy_list(1, Ny);
  for(int i = -nGx; i <= nGx; i++){
    Gx_list(0, i+nGx) = i * Gx;
  }
  for(int i = -nGy; i <= nGy; i++){
    Gy_list(0, i+nGy) = i * Gy;
  }

  meshGrid(Gx_list, Gy_list, Gx_mat, Gy_mat);

  Gx_mat.reshape(N, 1);
  Gy_mat.reshape(N, 1);
}

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
// IMPORTANT: this functoin need to be changed to be compatible with tensor interface
void RCWA::getGrandImaginaryMatrices(
  RCWAMatrices& grandImaginaryMatrices,
  const RCWAMatrices& im_eps_xx,
  const RCWAMatrices& im_eps_xy,
  const RCWAMatrices& im_eps_yx,
  const RCWAMatrices& im_eps_yy,
  const RCWAMatrices& im_eps_zz,
  int numOfLayer,
  int N
)
{
  for(int i = 0; i < numOfLayer; i++){
    RCWAMatrix grandImaginaryMatrix = zeros<RCWAMatrix>(3*N, 3*N);
    grandImaginaryMatrix(span(0, N-1), span(0, N-1)) = im_eps_xx[i];
    grandImaginaryMatrix(span(0, N-1), span(N, 2*N-1)) = im_eps_xy[i];
    grandImaginaryMatrix(span(N, 2*N-1), span(0, N-1)) = im_eps_yx[i];
    grandImaginaryMatrix(span(N, 2*N-1), span(N, 2*N-1)) = im_eps_yy[i];
    grandImaginaryMatrix(span(2*N, 3*N-1), span(2*N, 3*N-1)) = im_eps_zz[i];
    grandImaginaryMatrices.push_back(grandImaginaryMatrix);
  }
}

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
// IMPORTANT: this functoin need to be changed to be compatible with tensor interface
void RCWA::getEMatrices(
  RCWAMatrices& EMatrices,
  const RCWAMatrices& eps_xx,
  const RCWAMatrices& eps_xy,
  const RCWAMatrices& eps_yx,
  const RCWAMatrices& eps_yy,
  const int numOfLayer,
  const int N
){

  for(int i = 0; i < numOfLayer; i++){
    RCWAMatrix EMatrix = join_vert(
      join_horiz(eps_yy[i], -eps_yx[i]),
      join_horiz(-eps_xy[i], eps_xx[i])
    );
    EMatrices.push_back(EMatrix);
  }
}


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
// IMPORTANT: there is no change in this function even for a tensor
double RCWA::poyntingFlux(
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
){

  /*======================================================
  this part initializes parameters
  =======================================================*/
  kx = kx * omega;
  ky = ky * omega;
  int r1 = 0, r2 = 2 * N -1, r3 = 2 * N, r4 = 4 * N -1;
  // dcomplex IMAG_I = dcomplex(0, 1);
  RCWAMatrix onePadding4N(4*N, 4*N, fill::eye);
  RCWAMatrix onePadding2N(2*N, 2*N, fill::eye);
  RCWAMatrix onePadding1N(N, N, fill::eye);
  RCWAMatrix zeroPadding2N(2*N, 2*N, fill::zeros);
  RCWAMatrix zeroPadding4N(4*N, 4*N, fill::zeros);
  int numOfLayer = thicknessList.n_elem;

  // populate Gx and Gy matrices
  RCWAMatrix kxMat = diagmat(kx + Gx_mat);
  RCWAMatrix kyMat = diagmat(ky + Gy_mat);
  /*======================================================
  this part initializes structure matrices
  =======================================================*/
  RCWAMatrices TMatrices(numOfLayer), MMatrices(numOfLayer);
  RCWAMatrices EigenValMatrices(numOfLayer), EigenVecMatrices(numOfLayer), FMatrices(numOfLayer);


  // initialize K matrix
  RCWAMatrix KMatrix = join_vert(
    join_horiz(kxMat * kxMat, kxMat * kyMat),
    join_horiz(kyMat * kxMat, kyMat * kyMat)
  );
  /*======================================================
  This part solves RCWA
  e.g initialize M and F matrices, and compute the Eigen value problem
  =======================================================*/
  for(int i = 0; i < numOfLayer; i++){

    TMatrices[i] = join_vert(
      join_horiz(kyMat * eps_zz_inv[i] * kyMat, -kyMat * eps_zz_inv[i] * kxMat),
      join_horiz(-kxMat * eps_zz_inv[i] * kyMat, kxMat * eps_zz_inv[i] * kxMat)
    );
    RCWAMatrix eigMatrix = EMatrices[i] * (POW2(omega) * onePadding2N - TMatrices[i]) - KMatrix;
    cx_vec eigVal;
    // here is the problem
    eig_gen(eigVal, EigenVecMatrices[i], eigMatrix);

    eigVal = sqrt(eigVal);

    eigVal = -eigVal % sign(imag(eigVal));
    EigenValMatrices[i] = diagmat(eigVal);
    if(i == 0 || i == numOfLayer - 1){
      FMatrices[i] = onePadding2N;
    }
    else{
      FMatrices[i] = diagmat(exp(-IMAG_I * dcomplex(thicknessList(i),0) * eigVal));
    }
    MMatrices[i] = zeroPadding4N;
    MMatrices[i](span(r1, r2), span(r1, r2)) = (omega * onePadding2N - TMatrices[i] / omega) *
      EigenVecMatrices[i] * (EigenValMatrices[i]).i();

    MMatrices[i](span(r1, r2), span(r3, r4)) = -MMatrices[i](span(r1, r2), span(r1, r2));
    MMatrices[i](span(r3, r4), span(r1, r2)) = EigenVecMatrices[i];
    MMatrices[i](span(r3, r4), span(r3, r4)) = MMatrices[i](span(r3, r4), span(r1, r2));
    // normalization
    // MMatrices[i] = normalise(MMatrices[i], 2, 0);
    // MMatrices[i] = MMatrices[i] * (diagmat(sqrt(diagvec(MMatrices[i].t() * MMatrices[i])))).i();
  }

  /*======================================================
  This part initialize matrix for flux computation
  =======================================================*/

  double flux = 0;
  RCWAMatrices S_matrices_target;


  for(int i = 0; i < numOfLayer; i++){
    S_matrices_target.push_back(onePadding4N);
  }

  getSMatrices(targetLayer, N, numOfLayer,
      MMatrices, FMatrices, S_matrices_target, UP_);
  RCWAMatrix q_R, q_L, targetFields, P1, P2, Q1, Q2, R;
  RCWAMatrix integralSelf, integralMutual, integral, poyntingMat;
  RCWAMatrices S_matrices(numOfLayer), NewFMatrices(numOfLayer);

  RCWAMatrix source = zeros<RCWAMatrix>(4*N, 3*N);
  /*======================================================
  This part compute flux by collecting emission from source layers
  =======================================================*/

  for(int layerIdx = 0; layerIdx < targetLayer; layerIdx++){

    // if is not source layer, then continue
    if(sourceList[layerIdx] == false) continue;

    // initial steps, propogate S matrix
    RCWAMatrix q(diagvec(EigenValMatrices[layerIdx]));
    meshGrid(q, q, q_R, q_L);

    // defining source
    source(span(0,N-1), span(2*N, 3*N-1)) = -kyMat * eps_zz_inv[layerIdx] / omega;
    source(span(N, 2*N-1), span(2*N, 3*N-1)) = kxMat * eps_zz_inv[layerIdx] / omega;
    source(span(2*N, 3*N-1), span(N, 2*N-1)) = onePadding1N;
    source(span(3*N, 4*N-1), span(0, N-1)) = -onePadding1N;

    // treat as if the source layer has no thickness
    NewFMatrices = FMatrices;
    NewFMatrices[layerIdx] = onePadding2N;

    for(int i = 0; i < numOfLayer; i++){
      S_matrices[i] = onePadding4N;
    }

    getSMatrices(layerIdx, N, numOfLayer,
        MMatrices, NewFMatrices, S_matrices, ALL_);


    // solve the source
    targetFields = solve(MMatrices[layerIdx], source, solve_opts::fast);

    // calculating the P1 and P2
    P1 = solve(
      onePadding2N - S_matrices[targetLayer](span(r1, r2), span(r3, r4)) * S_matrices_target[numOfLayer-1](span(r3, r4), span(r1, r2)),
      S_matrices[targetLayer](span(r1, r2), span(r1, r2)),
      solve_opts::fast
    );

    P2 = S_matrices_target[numOfLayer-1](span(r3, r4), span(r1, r2)) * P1;

    // calculating the Q1 and Q2
    Q1 = onePadding2N - FMatrices[layerIdx] * S_matrices[0](span(r3, r4), span(r1, r2)) *
      FMatrices[layerIdx] * S_matrices[numOfLayer-1](span(r3, r4), span(r1, r2));

    Q2 = -FMatrices[layerIdx] * S_matrices[0](span(r3, r4), span(r1, r2));

    // calculating R
    R = MMatrices[targetLayer] * join_vert(FMatrices[targetLayer] * P1, P2) *
      Q1.i() * join_horiz(onePadding2N, Q2);

    // calculating integrands
    if(layerIdx == 0 || layerIdx == numOfLayer - 1){
      integralSelf = 1 / (IMAG_I * (q_L - conj(q_R)));
      integralMutual = zeroPadding2N;
    }
    else{
      integralSelf = (1 - exp(-IMAG_I * dcomplex(thicknessList[layerIdx], 0) * (q_L - conj(q_R)))) /
        (IMAG_I * (q_L - conj(q_R)));
      integralMutual = (exp(IMAG_I * dcomplex(thicknessList[layerIdx], 0) * conj(q_R)) - exp(-IMAG_I * dcomplex(thicknessList[layerIdx], 0) * q_L)) /
        (IMAG_I * (q_L + conj(q_R)));
    }

    integral = join_vert(
      join_horiz(integralSelf, integralMutual),
      join_horiz(integralMutual, integralSelf)
    );

    // calculating kernel
    poyntingMat = (targetFields * grandImaginaryMatrices[layerIdx] * targetFields.t()) % integral;

    poyntingMat = -R * poyntingMat * R.t();
    flux += real(trace(poyntingMat(span(r1, r2), span(r3, r4))));
  }
  return flux;

}

