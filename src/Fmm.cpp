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

 #include "Fmm.h"
 #include "Common.h"

 namespace FMM{
   /*==============================================*/
   // helper function to change a scalar dielectric to a tensor
   // @args:
   // epsilon: scalar field
   /*==============================================*/
   static EpsilonVal fromScalarToDiagonal(const EpsilonVal epsilon){
     EpsilonVal result;
     for(int i = 0; i < 3; i++){
       result.diagonal[2*i] = epsilon.scalar[0];
       result.diagonal[2*i + 1] = epsilon.scalar[1];
     }
     return result;
   }
   /*==============================================*/
   // helper function to change a diagonal dielectric to a tensor
   // @args:
   // epsilon: diagonal field
   /*==============================================*/
   static EpsilonVal fromDiagonalToTensor(const EpsilonVal epsilon){
     EpsilonVal result;
     for(int i = 0; i < 10; i++){
       result.tensor[i] = 0;
     }
     result.tensor[0] = epsilon.diagonal[0];
     result.tensor[1] = epsilon.diagonal[1];
     result.tensor[6] = epsilon.diagonal[2];
     result.tensor[7] = epsilon.diagonal[3];
     result.tensor[8] = epsilon.diagonal[4];
     result.tensor[9] = epsilon.diagonal[5];

     return result;
   }
   /*==============================================*/
   // helper function to change a scalar dielectric to a tensor
   // @args:
   // epsilon: scalar field
   /*==============================================*/
   static EpsilonVal fromScalarToTensor(const EpsilonVal epsilon){
     EpsilonVal diag = fromScalarToDiagonal(epsilon);
     return fromDiagonalToTensor(diag);
   }

   /*==============================================*/
   // helper function to change a random dielectric to a tensor
   // @args:
   // epsilon: field
   // type: the type of the dielectric
   /*==============================================*/
  EpsilonVal toTensor(const EpsilonVal epsilon, const EPSTYPE type){
     switch (type) {
       case SCALAR_: return fromScalarToTensor(epsilon);
       case DIAGONAL_: return fromDiagonalToTensor(epsilon);
       default: return epsilon;
     }
   }

  /*==============================================*/
   // helper function to do fourier transform for one value
   // @args:
   // epsVal: the value need to be transformed
   // epsBG: the background value need to be transformed
   // G_mat: the G_mat for grating
   // center: the center position
   // width: the width of the grating
   /*==============================================*/
   static RCWAMatrix transformGratingElement(
     const dcomplex epsVal,
     const dcomplex epsBG,
     const RCWAMatrix& G_mat,
     const double center,
     const double width
   ){
    return exp(IMAG_I * G_mat * center) * (epsVal - epsBG)
             * width % sinc(G_mat / 2 * width);
   }

   /*==============================================*/
   // helper function to do fourier transform for one value
   // @args:
   // epsVal: the value need to be transformed
   // epsBG: the background value need to be transformed
   // G_mat: the G_mat for grating
   // centerx: the center position in x direction
   // centery: the center position in y direction
   // widthx: the width of the rectangle in x direction
   // widthy: the width of the rectangle in y direction
   /*==============================================*/
   static RCWAMatrix transformRectangleElement(
     const dcomplex epsVal,
     const dcomplex epsBG,
     const RCWAMatrix& GxMat,
     const RCWAMatrix& GyMat,
     const double centerx,
     const double centery,
     const double widthx,
     const double widthy
   ){
     return widthx * widthy * (epsVal - epsBG)
             * exp(IMAG_I * (GxMat * centerx + GyMat * centery))
             % sinc(GxMat * widthx / 2) % sinc(GyMat * widthy / 2);
   }

   /*==============================================*/
   // helper function to do fourier transform for one value
   // @args:
   // epsVal: the value need to be transformed
   // epsBG: the background value need to be transformed
   // G_mat: the G_mat for grating
   // centerx: the center position in x direction
   // centery: the center position in y direction
   // radius: the radius of the circle
   /*==============================================*/
   static RCWAMatrix transformCircleElement(
     const dcomplex epsVal,
     const dcomplex epsBG,
     const RCWAMatrix& GxMat,
     const RCWAMatrix& GyMat,
     const double centerx,
     const double centery,
     const double radius
   ){
     arma::mat rho = sqrt(square(real(GxMat)) + square(real(GyMat))) * 2 * datum::pi * radius;
     arma::mat jincMat = zeros<arma::mat>( size(rho) );

     mat::iterator jincMat_it = jincMat.begin();
     int count = 0;
     for(mat::const_iterator it = rho.begin(); it != rho.end(); it++){
       if(*it == 0.0){
         *(jincMat_it + count) = 1;
       }
       else{
         *(jincMat_it + count) = jinc(*it);
       }
       count++;
     }

     return (epsVal - epsBG) * exp(IMAG_I * (GxMat * centerx + GyMat * centery))
            * 2 * datum::pi * POW2(radius) % jincMat;

   }
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
   // N: the total number of G
   // center: the center of the grating
   // width: the width of the grating
   // period: the periodicity
   // hasTensor: whether this layer contains tensor
   // useInverse: whether to use inverse rule
   /*==============================================*/
   void transformGrating(
    RCWAMatrix& eps_xx,
    RCWAMatrix& eps_xy,
    RCWAMatrix& eps_yx,
    RCWAMatrix& eps_yy,
    RCWAMatrix& eps_zz_Inv,
    RCWAMatrix& im_eps_xx,
    RCWAMatrix& im_eps_xy,
    RCWAMatrix& im_eps_yx,
    RCWAMatrix& im_eps_yy,
    RCWAMatrix& im_eps_zz,
    const EpsilonVal& epsBGTensor,
    const EpsilonVal& epsilon,
    const EPSTYPE epsilonType,
    const int N,
    const double center,
    const double width,
    const double period,
    const bool hasTensor,
    const bool useInverse
   ){

     RCWAMatrix Gx_mat, Gy_mat;
     double periods[2] = {period, 0};
     RCWA::getGMatrices((N-1)/2, 0, periods, Gx_mat, Gy_mat, ONE_);
     // dcomplex IMAG_I = dcomplex(0, 1.0);
     RCWAMatrix Gx_r, Gx_l, Gy_r, Gy_l;
     meshGrid(Gx_mat, Gx_mat, Gx_r, Gx_l);

     RCWAMatrix G_mat(real(Gx_l) - real(Gx_r), zeros<arma::mat>(N, N));


     dcomplex eps_BG_xx = dcomplex(epsBGTensor.tensor[0], epsBGTensor.tensor[1]);
     dcomplex eps_BG_xy = dcomplex(epsBGTensor.tensor[2], epsBGTensor.tensor[3]);
     dcomplex eps_BG_yx = dcomplex(epsBGTensor.tensor[4], epsBGTensor.tensor[5]);
     dcomplex eps_BG_yy = dcomplex(epsBGTensor.tensor[6], epsBGTensor.tensor[7]);
     dcomplex eps_BG_zz = dcomplex(epsBGTensor.tensor[8], epsBGTensor.tensor[9]);
     dcomplex im_eps_BG_xx = dcomplex(epsBGTensor.tensor[1], 0);
     dcomplex im_eps_BG_xy = dcomplex(epsBGTensor.tensor[3], 0);
     dcomplex im_eps_BG_yx = dcomplex(epsBGTensor.tensor[5], 0);
     dcomplex im_eps_BG_yy = dcomplex(epsBGTensor.tensor[7], 0);
     dcomplex im_eps_BG_zz = dcomplex(epsBGTensor.tensor[9], 0);

     EpsilonVal epsTensor = toTensor(epsilon, epsilonType);

     eps_xx += transformGratingElement(dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]),
         eps_BG_xx, G_mat, center, width) / period;
     im_eps_xx += transformGratingElement(dcomplex(epsTensor.tensor[1], 0),
         im_eps_BG_xx, G_mat, center, width) / period;

     eps_yy += transformGratingElement(dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]),
         eps_BG_yy, G_mat, center, width) / period;
     im_eps_yy += transformGratingElement(dcomplex(epsTensor.tensor[7], 0),
         im_eps_BG_yy, G_mat, center, width) / period;

     eps_zz_Inv += transformGratingElement(dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]),
         eps_BG_zz, G_mat, center, width) / period;
     im_eps_zz += transformGratingElement(dcomplex(epsTensor.tensor[9], 0),
         im_eps_BG_zz, G_mat, center, width) / period;

     if(hasTensor){
         eps_xy += transformGratingElement(dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]),
           eps_BG_xy, G_mat, center, width) / period;
         eps_yx += transformGratingElement(dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]),
           eps_BG_yx, G_mat, center, width) / period;
         im_eps_xy += transformGratingElement(dcomplex(epsTensor.tensor[3], 0),
           im_eps_BG_xy, G_mat, center, width) / period;
         im_eps_yx += transformGratingElement(dcomplex(epsTensor.tensor[5], 0),
           im_eps_BG_yx, G_mat, center, width) / period;
     }
   }

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
   // nGx: the total number of G in x direction
   // nGy: the total number of G in y direction
   // centers: the centers of the rectangle
   // widths: the widths of the rectangle
   // period: the periodicity
   // hasTensor: whether this layer contains tensor
   // useInverse: whether to use inverse rule
   /*==============================================*/
   void transformRectangle(
    RCWAMatrix& eps_xx,
    RCWAMatrix& eps_xy,
    RCWAMatrix& eps_yx,
    RCWAMatrix& eps_yy,
    RCWAMatrix& eps_zz_Inv,
    RCWAMatrix& im_eps_xx,
    RCWAMatrix& im_eps_xy,
    RCWAMatrix& im_eps_yx,
    RCWAMatrix& im_eps_yy,
    RCWAMatrix& im_eps_zz,
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
   ){

     int N = RCWA::getN(nGx, nGy);
     double area = period[0] * period[1];

     RCWAMatrix Gx_mat, Gy_mat;
     RCWA::getGMatrices(nGx, nGy, period, Gx_mat, Gy_mat, TWO_);
     // dcomplex IMAG_I = dcomplex(0, 1.0);
     RCWAMatrix Gx_r, Gx_l, Gy_r, Gy_l;
     meshGrid(Gx_mat, Gx_mat, Gx_r, Gx_l);
     meshGrid(Gy_mat, Gy_mat, Gy_r, Gy_l);

     RCWAMatrix GxMat(real(Gx_l) - real(Gx_r), zeros<arma::mat>(N, N));
     RCWAMatrix GyMat(real(Gy_l) - real(Gy_r), zeros<arma::mat>(N, N));


     dcomplex eps_BG_xx = dcomplex(epsBGTensor.tensor[0], epsBGTensor.tensor[1]);
     dcomplex eps_BG_xy = dcomplex(epsBGTensor.tensor[2], epsBGTensor.tensor[3]);
     dcomplex eps_BG_yx = dcomplex(epsBGTensor.tensor[4], epsBGTensor.tensor[5]);
     dcomplex eps_BG_yy = dcomplex(epsBGTensor.tensor[6], epsBGTensor.tensor[7]);
     dcomplex eps_BG_zz = dcomplex(epsBGTensor.tensor[8], epsBGTensor.tensor[9]);
     dcomplex im_eps_BG_xx = dcomplex(epsBGTensor.tensor[1], 0);
     dcomplex im_eps_BG_xy = dcomplex(epsBGTensor.tensor[3], 0);
     dcomplex im_eps_BG_yx = dcomplex(epsBGTensor.tensor[5], 0);
     dcomplex im_eps_BG_yy = dcomplex(epsBGTensor.tensor[7], 0);
     dcomplex im_eps_BG_zz = dcomplex(epsBGTensor.tensor[9], 0);

     EpsilonVal epsTensor = toTensor(epsilon, epsilonType);

     eps_xx += transformRectangleElement(dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]), eps_BG_xx,
         GxMat, GyMat, centers[0], centers[1], widths[0], widths[1]) / area;
     im_eps_xx += transformRectangleElement(dcomplex(epsTensor.tensor[1], 0), im_eps_BG_xx,
         GxMat, GyMat, centers[0], centers[1], widths[0], widths[1]) / area;

     eps_yy += transformRectangleElement(dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]), eps_BG_yy,
         GxMat, GyMat, centers[0], centers[1], widths[0], widths[1]) / area;
     im_eps_yy += transformRectangleElement(dcomplex(epsTensor.tensor[7], 0), im_eps_BG_yy,
         GxMat, GyMat, centers[0], centers[1], widths[0], widths[1]) / area;

     eps_zz_Inv += transformRectangleElement(dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]), eps_BG_zz,
         GxMat, GyMat, centers[0], centers[1], widths[0], widths[1]) / area;
     im_eps_zz += transformRectangleElement(dcomplex(epsTensor.tensor[9], 0), im_eps_BG_zz,
         GxMat, GyMat, centers[0], centers[1], widths[0], widths[1]) / area;


     if(hasTensor){
        eps_xy += transformRectangleElement(dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]), eps_BG_xy,
           GxMat, GyMat, centers[0], centers[1], widths[0], widths[1]) / area;
        eps_yx += transformRectangleElement(dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]), eps_BG_yx,
           GxMat, GyMat, centers[0], centers[1], widths[0], widths[1]) / area;

        im_eps_xy += transformRectangleElement(dcomplex(epsTensor.tensor[3], 0), im_eps_BG_xy,
           GxMat, GyMat, centers[0], centers[1], widths[0], widths[1]) / area;
        im_eps_yx += transformRectangleElement(dcomplex(epsTensor.tensor[5], 0), im_eps_BG_yx,
           GxMat, GyMat, centers[0], centers[1], widths[0], widths[1]) / area;
     }
   }

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
   // useInverse: whether to use inverse rule
   /*==============================================*/
   void transformCircle(
    RCWAMatrix& eps_xx,
    RCWAMatrix& eps_xy,
    RCWAMatrix& eps_yx,
    RCWAMatrix& eps_yy,
    RCWAMatrix& eps_zz_Inv,
    RCWAMatrix& im_eps_xx,
    RCWAMatrix& im_eps_xy,
    RCWAMatrix& im_eps_yx,
    RCWAMatrix& im_eps_yy,
    RCWAMatrix& im_eps_zz,
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
   ){
     int N = RCWA::getN(nGx, nGy);
     double area = period[0] * period[1];

     RCWAMatrix Gx_mat, Gy_mat;
     RCWA::getGMatrices(nGx, nGy, period, Gx_mat, Gy_mat, TWO_);
     // dcomplex IMAG_I = dcomplex(0, 1.0);
     RCWAMatrix Gx_r, Gx_l, Gy_r, Gy_l;
     meshGrid(Gx_mat, Gx_mat, Gx_r, Gx_l);
     meshGrid(Gy_mat, Gy_mat, Gy_r, Gy_l);

     RCWAMatrix GxMat = Gx_l - Gx_r;
     RCWAMatrix GyMat = Gy_l - Gy_r;

     RCWAMatrix onePadding1N = eye<RCWAMatrix>(N, N);

     dcomplex eps_BG_xx = dcomplex(epsBGTensor.tensor[0], epsBGTensor.tensor[1]);
     dcomplex eps_BG_xy = dcomplex(epsBGTensor.tensor[2], epsBGTensor.tensor[3]);
     dcomplex eps_BG_yx = dcomplex(epsBGTensor.tensor[4], epsBGTensor.tensor[5]);
     dcomplex eps_BG_yy = dcomplex(epsBGTensor.tensor[6], epsBGTensor.tensor[7]);
     dcomplex eps_BG_zz = dcomplex(epsBGTensor.tensor[8], epsBGTensor.tensor[9]);
     dcomplex im_eps_BG_xx = dcomplex(epsBGTensor.tensor[1], 0);
     dcomplex im_eps_BG_xy = dcomplex(epsBGTensor.tensor[3], 0);
     dcomplex im_eps_BG_yx = dcomplex(epsBGTensor.tensor[5], 0);
     dcomplex im_eps_BG_yy = dcomplex(epsBGTensor.tensor[7], 0);
     dcomplex im_eps_BG_zz = dcomplex(epsBGTensor.tensor[9], 0);

     EpsilonVal epsTensor = toTensor(epsilon, epsilonType);

     eps_xx += transformCircleElement(dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]), eps_BG_xx,
         GxMat, GyMat, centers[0], centers[1], radius) / area;
     im_eps_xx += transformCircleElement(dcomplex(epsTensor.tensor[1], 0), im_eps_BG_xx,
         GxMat, GyMat, centers[0], centers[1], radius) / area;

     eps_yy += transformCircleElement(dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]), eps_BG_yy,
         GxMat, GyMat, centers[0], centers[1], radius) / area;
     im_eps_yy += transformCircleElement(dcomplex(epsTensor.tensor[7], 0), im_eps_BG_yy,
         GxMat, GyMat, centers[0], centers[1], radius) / area;

     eps_zz_Inv += transformCircleElement(dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]), eps_BG_zz,
         GxMat, GyMat, centers[0], centers[1], radius) / area;
     im_eps_zz += transformCircleElement(dcomplex(epsTensor.tensor[9], 0), im_eps_BG_zz,
         GxMat, GyMat, centers[0], centers[1], radius) / area;


     if(hasTensor){
        eps_xy += transformCircleElement(dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]), eps_BG_xy,
           GxMat, GyMat, centers[0], centers[1], radius) / area;
        eps_yx += transformCircleElement(dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]), eps_BG_yx,
           GxMat, GyMat, centers[0], centers[1], radius) / area;

        im_eps_xy += transformCircleElement(dcomplex(epsTensor.tensor[3], 0), im_eps_BG_xy,
           GxMat, GyMat, centers[0], centers[1], radius) / area;
        im_eps_yx += transformCircleElement(dcomplex(epsTensor.tensor[5], 0), im_eps_BG_yx,
           GxMat, GyMat, centers[0], centers[1], radius) / area;
     }
    }
 }