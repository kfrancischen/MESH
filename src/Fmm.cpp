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
   static RCWAcMatrix transformGratingElement(
     const dcomplex epsVal,
     const dcomplex epsBG,
     const RCWArMatrix& G_mat,
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
   static RCWAcMatrix transformRectangleElement(
     const dcomplex epsVal,
     const dcomplex epsBG,
     const RCWArMatrix& GxMat,
     const RCWArMatrix& GyMat,
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
   static RCWAcMatrix transformCircleElement(
     const dcomplex epsVal,
     const dcomplex epsBG,
     const RCWArMatrix& GxMat,
     const RCWArMatrix& GyMat,
     const double centerx,
     const double centery,
     const double radius
   ){
     RCWArMatrix rho = sqrt(square(GxMat) + square(GyMat)) * radius;
     RCWArMatrix jincMat = zeros<RCWArMatrix>( size(rho) );

     RCWArMatrix::iterator jincMat_it = jincMat.begin();
     int count = 0;
     for(mat::const_iterator it = rho.begin(); it != rho.end(); it++){
       *(jincMat_it + count) = jinc(*it);
       count++;
     }

     return (epsVal - epsBG) * exp(IMAG_I * (GxMat * centerx + GyMat * centery))
            * 2 * datum::pi * POW2(radius) % jincMat;

   }
   /*==============================================*/
   // helper function to do fourier transform for one value
   // @args:
   // epsVal: the value need to be transformed
   // epsBG: the background value need to be transformed
   // G_mat: the G_mat for grating
   // centerx: the center position in x direction
   // centery: the center position in y direction
   // a: the half width of the rectangle in x direction
   // b: the half width of the rectangle in y direction
   /*==============================================*/
   static RCWAcMatrix transformEllipseElement(
     const dcomplex epsVal,
     const dcomplex epsBG,
     const RCWArMatrix& GxMat,
     const RCWArMatrix& GyMat,
     const double centerx,
     const double centery,
     const double a,
     const double b
   ){
     RCWArMatrix rho = sqrt(square(GxMat * a) + square(GyMat * b));
     RCWArMatrix jincMat = zeros<RCWArMatrix>( size(rho) );

     RCWArMatrix::iterator jincMat_it = jincMat.begin();
     int count = 0;
     for(RCWArMatrix::const_iterator it = rho.begin(); it != rho.end(); it++){
       *(jincMat_it + count) = jinc(*it);
       count++;
     }

     return (epsVal - epsBG) * exp(IMAG_I * (GxMat * centerx + GyMat * centery))
            * 2 * datum::pi * a * b % jincMat;
   }
   /*==============================================*/
   // helper function to do fourier transform for one value
   // @args:
   // epsVal: the value need to be transformed
   // epsBG: the background value need to be transformed
   // G_mat: the G_mat for grating
   // centerx: the center position in x direction
   // centery: the center position in y direction
   // edgeList: the list of edges
   // area: the area of the polygon
   /*==============================================*/
   static RCWAcMatrix transformPolygonElement(
     const dcomplex epsVal,
     const dcomplex epsBG,
     const RCWArMatrix& GxMat,
     const RCWArMatrix& GyMat,
     const double centerx,
     const double centery,
     const EdgeList& edgeList,
     const double area
   ){
     RCWAcMatrix result = zeros<RCWAcMatrix>( size(GxMat) );
     for(size_t i = 0; i < edgeList.size(); i++){
       double x_cur = edgeList[i].first, y_cur = edgeList[i+1].second;
       double x_pre, y_pre, x_next, y_next;
       if(i == 0){
         x_pre = edgeList[edgeList.size()-1].first;
         y_pre = edgeList[edgeList.size()-1].second;
       }
       else{
         x_pre = edgeList[i-1].first;
         y_pre = edgeList[i-1].second;
       }
       if(i == edgeList.size() - 1){
         x_next = edgeList[0].first;
         y_next = edgeList[0].second;
       }
       else{
         x_next = edgeList[i+1].first;
         y_next = edgeList[i+1].second;
       }

       RCWArMatrix::const_iterator GyMat_it = GyMat.cbegin();
       RCWAcMatrix::iterator result_it = result.begin();
       int count = 0;
       for(RCWArMatrix::const_iterator GxMat_it = GxMat.cbegin(); GxMat_it != GxMat.cend(); GxMat_it++){
         if( (*GxMat_it) == 0 && *(GyMat_it + count) == 0 ){
           *(result_it + count) = area; // this should be only computed once per iteration
         }
         else{
           *(result_it + count) += exp(IMAG_I * ((*GxMat_it) * x_cur +  (*(GyMat_it + count)) * y_cur))
              * ((y_cur - y_pre) * (x_next - x_cur) - (y_next - y_cur) * (x_cur - x_pre))
              / ( ( (*GxMat_it) * (x_cur - x_pre) + (*(GyMat_it + count)) * (y_cur - y_pre))  * ( (*GxMat_it) * (x_next - x_cur) +  (*(GyMat_it + count)) * (y_next - y_cur)));
         }
         count++;
       }
     }
     return (epsVal - epsBG) * exp(IMAG_I * (GxMat * centerx + GyMat * centery)) % result;
   }
   /*==============================================*/
   // This function computes the Fourier transform for grating geometry
   // @args:
   // eps_xx: the Fourier transform for eps_xx
   // eps_xy: the Fourier transform for eps_xy
   // eps_zx: the Fourier transform for eps_yx
   // eps_yy: the Fourier transform for eps_yy
   // eps_zz: the Fourier transform for eps_zz
   // im_eps_xx: the Fourier transform for imaginary part
   // im_eps_xy: the Fourier transform for imaginary part
   // im_eps_yx: the Fourier transform for imaginary part
   // im_eps_yy: the Fourier transform for imaginary part
   // im_eps_zz: the Fourier transform for imaginary part
   // epsilonBGTensor: the epsilon of bacground (transformed to tensor already)
   // nGx: the total number of G
   // center: the center of the grating
   // width: the width of the grating
   // period: the periodicity
   // hasTensor: whether this layer contains tensor
   /*==============================================*/
   void transformGrating(
    RCWAcMatrix& eps_xx,
    RCWAcMatrix& eps_xy,
    RCWAcMatrix& eps_yx,
    RCWAcMatrix& eps_yy,
    RCWAcMatrix& eps_zz,
    RCWAcMatrix& im_eps_xx,
    RCWAcMatrix& im_eps_xy,
    RCWAcMatrix& im_eps_yx,
    RCWAcMatrix& im_eps_yy,
    RCWAcMatrix& im_eps_zz,
    const EpsilonVal& epsBGTensor,
    const EpsilonVal& epsilon,
    const EPSTYPE epsilonType,
    const int nGx,
    const double center,
    const double width,
    const double period,
    const bool hasTensor
   ){

     RCWArMatrix Gx_mat, Gy_mat;
     double periods[2] = {period, 0};
     RCWA::getGMatrices(nGx, 0, periods, Gx_mat, Gy_mat, ONE_);
     // dcomplex IMAG_I = dcomplex(0, 1.0);
     RCWArMatrix Gx_r, Gx_l, Gy_r, Gy_l;
     meshGrid(Gx_mat, Gx_mat, Gx_r, Gx_l);

     RCWArMatrix G_mat = Gx_l - Gx_r;


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

     eps_zz += transformGratingElement(dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]),
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
   // eps_xx: the Fourier transform for eps_xx
   // eps_xy: the Fourier transform for eps_xy
   // eps_zx: the Fourier transform for eps_yx
   // eps_yy: the Fourier transform for eps_yy
   // eps_zz: the Fourier transform for eps_zz
   // im_eps_xx: the Fourier transform for imaginary part
   // im_eps_xy: the Fourier transform for imaginary part
   // im_eps_yx: the Fourier transform for imaginary part
   // im_eps_yy: the Fourier transform for imaginary part
   // im_eps_zz: the Fourier transform for imaginary part
   // epsilonBGTensor: the epsilon of bacground (transformed to tensor already)
   // nGx: the total number of G in x direction
   // nGy: the total number of G in y direction
   // centers: the centers of the rectangle
   // widths: the widths of the rectangle
   // period: the periodicity
   // hasTensor: whether this layer contains tensor
   /*==============================================*/
   void transformRectangle(
    RCWAcMatrix& eps_xx,
    RCWAcMatrix& eps_xy,
    RCWAcMatrix& eps_yx,
    RCWAcMatrix& eps_yy,
    RCWAcMatrix& eps_zz,
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
    const bool hasTensor
   ){

     double area = period[0] * period[1];

     RCWArMatrix Gx_mat, Gy_mat;
     RCWA::getGMatrices(nGx, nGy, period, Gx_mat, Gy_mat, TWO_);
     // dcomplex IMAG_I = dcomplex(0, 1.0);
     RCWArMatrix Gx_r, Gx_l, Gy_r, Gy_l;
     meshGrid(Gx_mat, Gx_mat, Gx_r, Gx_l);
     meshGrid(Gy_mat, Gy_mat, Gy_r, Gy_l);

     RCWArMatrix GxMat = Gx_l - Gx_r;
     RCWArMatrix GyMat = Gy_l - Gy_r;

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

     eps_zz += transformRectangleElement(dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]), eps_BG_zz,
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
   // eps_xx: the Fourier transform for eps_xx
   // eps_xy: the Fourier transform for eps_xy
   // eps_zx: the Fourier transform for eps_yx
   // eps_yy: the Fourier transform for eps_yy
   // eps_zz: the Fourier transform for eps_zz
   // im_eps_xx: the Fourier transform for imaginary part
   // im_eps_xy: the Fourier transform for imaginary part
   // im_eps_yx: the Fourier transform for imaginary part
   // im_eps_yy: the Fourier transform for imaginary part
   // im_eps_zz: the Fourier transform for imaginary part
   // epsilonBGTensor: the epsilon of bacground (transformed to tensor already)
   // nG_x: the total number of G in x direction
   // nG_y: the total number of G in y direction
   // centers: the centers of the circle
   // radius: the radius of the circle
   // period: the periodicity
   // hasTensor: whether this layer contains tensor
   /*==============================================*/
   void transformCircle(
    RCWAcMatrix& eps_xx,
    RCWAcMatrix& eps_xy,
    RCWAcMatrix& eps_yx,
    RCWAcMatrix& eps_yy,
    RCWAcMatrix& eps_zz,
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
    const bool hasTensor
   ){
     int N = RCWA::getN(nGx, nGy);
     double area = period[0] * period[1];

     RCWArMatrix Gx_mat, Gy_mat;
     RCWA::getGMatrices(nGx, nGy, period, Gx_mat, Gy_mat, TWO_);
     // dcomplex IMAG_I = dcomplex(0, 1.0);
     RCWArMatrix Gx_r, Gx_l, Gy_r, Gy_l;
     meshGrid(Gx_mat, Gx_mat, Gx_r, Gx_l);
     meshGrid(Gy_mat, Gy_mat, Gy_r, Gy_l);

     RCWArMatrix GxMat = Gx_l - Gx_r;
     RCWArMatrix GyMat = Gy_l - Gy_r;
     //cout << GxMat << std::endl << GyMat << std::endl;
     RCWAcMatrix onePadding1N = eye<RCWAcMatrix>(N, N);

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

     eps_zz += transformCircleElement(dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]), eps_BG_zz,
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

    /*==============================================*/
    // This function computes the Fourier transform for ellipse geometry
    // @args:
    // eps_xx: the Fourier transform for eps_xx
    // eps_xy: the Fourier transform for eps_xy
    // eps_zx: the Fourier transform for eps_yx
    // eps_yy: the Fourier transform for eps_yy
    // eps_zz: the Fourier transform for eps_zz
    // im_eps_xx: the Fourier transform for imaginary part
    // im_eps_xy: the Fourier transform for imaginary part
    // im_eps_yx: the Fourier transform for imaginary part
    // im_eps_yy: the Fourier transform for imaginary part
    // im_eps_zz: the Fourier transform for imaginary part
    // epsilonBGTensor: the epsilon of bacground (transformed to tensor already)
    // nG_x: the total number of G in x direction
    // nG_y: the total number of G in y direction
    // centers: the centers of the ellipse
    // halfwidths: the halfwidths of the ellipse
    // period: the periodicity
    // hasTensor: whether this layer contains tensor
    /*==============================================*/
    void transformEllipse(
     RCWAcMatrix& eps_xx,
     RCWAcMatrix& eps_xy,
     RCWAcMatrix& eps_yx,
     RCWAcMatrix& eps_yy,
     RCWAcMatrix& eps_zz,
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
     const double halfwidths[2],
     const double period[2],
     const bool hasTensor
    ){
      int N = RCWA::getN(nGx, nGy);
      double area = period[0] * period[1];

      RCWArMatrix Gx_mat, Gy_mat;
      RCWA::getGMatrices(nGx, nGy, period, Gx_mat, Gy_mat, TWO_);
      // dcomplex IMAG_I = dcomplex(0, 1.0);
      RCWArMatrix Gx_r, Gx_l, Gy_r, Gy_l;
      meshGrid(Gx_mat, Gx_mat, Gx_r, Gx_l);
      meshGrid(Gy_mat, Gy_mat, Gy_r, Gy_l);

      RCWArMatrix GxMat = Gx_l - Gx_r;
      RCWArMatrix GyMat = Gy_l - Gy_r;

      RCWAcMatrix onePadding1N = eye<RCWAcMatrix>(N, N);

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

      eps_xx += transformEllipseElement(dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]), eps_BG_xx,
          GxMat, GyMat, centers[0], centers[1], halfwidths[0], halfwidths[1]) / area;
      im_eps_xx += transformEllipseElement(dcomplex(epsTensor.tensor[1], 0), im_eps_BG_xx,
          GxMat, GyMat, centers[0], centers[1], halfwidths[0], halfwidths[1]) / area;

      eps_yy += transformEllipseElement(dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]), eps_BG_yy,
          GxMat, GyMat, centers[0], centers[1], halfwidths[0], halfwidths[1]) / area;
      im_eps_yy += transformEllipseElement(dcomplex(epsTensor.tensor[7], 0), im_eps_BG_yy,
          GxMat, GyMat, centers[0], centers[1], halfwidths[0], halfwidths[1]) / area;

      eps_zz += transformEllipseElement(dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]), eps_BG_zz,
          GxMat, GyMat, centers[0], centers[1], halfwidths[0], halfwidths[1]) / area;
      im_eps_zz += transformEllipseElement(dcomplex(epsTensor.tensor[9], 0), im_eps_BG_zz,
          GxMat, GyMat, centers[0], centers[1], halfwidths[0], halfwidths[1]) / area;


      if(hasTensor){
         eps_xy += transformEllipseElement(dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]), eps_BG_xy,
            GxMat, GyMat, centers[0], centers[1], halfwidths[0], halfwidths[1]) / area;
         eps_yx += transformEllipseElement(dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]), eps_BG_yx,
            GxMat, GyMat, centers[0], centers[1], halfwidths[0], halfwidths[1]) / area;

         im_eps_xy += transformEllipseElement(dcomplex(epsTensor.tensor[3], 0), im_eps_BG_xy,
            GxMat, GyMat, centers[0], centers[1], halfwidths[0], halfwidths[1]) / area;
         im_eps_yx += transformEllipseElement(dcomplex(epsTensor.tensor[5], 0), im_eps_BG_yx,
            GxMat, GyMat, centers[0], centers[1], halfwidths[0], halfwidths[1]) / area;
      }
     }
     /*==============================================*/
     // This function computes the Fourier transform for polygon geometry
     // @args:
     // eps_xx: the Fourier transform for eps_xx
     // eps_xy: the Fourier transform for eps_xy
     // eps_zx: the Fourier transform for eps_yx
     // eps_yy: the Fourier transform for eps_yy
     // eps_zz: the Fourier transform for eps_zz
     // im_eps_xx: the Fourier transform for imaginary part
     // im_eps_xy: the Fourier transform for imaginary part
     // im_eps_yx: the Fourier transform for imaginary part
     // im_eps_yy: the Fourier transform for imaginary part
     // im_eps_zz: the Fourier transform for imaginary part
     // epsilonBGTensor: the epsilon of bacground (transformed to tensor already)
     // nG_x: the total number of G in x direction
     // nG_y: the total number of G in y direction
     // centers: the centers of the polygon
     // edgeList: the edges of the polygon
     // period: the periodicity
     // hasTensor: whether this layer contains tensor
     /*==============================================*/
     void transformPolygon(
      RCWAcMatrix& eps_xx,
      RCWAcMatrix& eps_xy,
      RCWAcMatrix& eps_yx,
      RCWAcMatrix& eps_yy,
      RCWAcMatrix& eps_zz,
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
      const EdgeList& edgeList,
      const double period[2],
      const bool hasTensor
    ){
      int N = RCWA::getN(nGx, nGy);
      double area = period[0] * period[1];
      double polygonArea = getPolygonArea(edgeList);
      RCWArMatrix Gx_mat, Gy_mat;
      RCWA::getGMatrices(nGx, nGy, period, Gx_mat, Gy_mat, TWO_);
      // dcomplex IMAG_I = dcomplex(0, 1.0);
      RCWArMatrix Gx_r, Gx_l, Gy_r, Gy_l;
      meshGrid(Gx_mat, Gx_mat, Gx_r, Gx_l);
      meshGrid(Gy_mat, Gy_mat, Gy_r, Gy_l);

      RCWArMatrix GxMat = Gx_l - Gx_r;
      RCWArMatrix GyMat = Gy_l - Gy_r;

      RCWAcMatrix onePadding1N = eye<RCWAcMatrix>(N, N);

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

      eps_xx += transformPolygonElement(dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]), eps_BG_xx,
          GxMat, GyMat, centers[0], centers[1], edgeList, polygonArea) / area;
      im_eps_xx += transformPolygonElement(dcomplex(epsTensor.tensor[1], 0), im_eps_BG_xx,
          GxMat, GyMat, centers[0], centers[1], edgeList, polygonArea) / area;

      eps_yy += transformPolygonElement(dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]), eps_BG_yy,
          GxMat, GyMat, centers[0], centers[1], edgeList, polygonArea) / area;
      im_eps_yy += transformPolygonElement(dcomplex(epsTensor.tensor[7], 0), im_eps_BG_yy,
          GxMat, GyMat, centers[0], centers[1], edgeList, polygonArea) / area;

      eps_zz += transformPolygonElement(dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]), eps_BG_zz,
          GxMat, GyMat, centers[0], centers[1], edgeList, polygonArea) / area;
      im_eps_zz += transformPolygonElement(dcomplex(epsTensor.tensor[9], 0), im_eps_BG_zz,
          GxMat, GyMat, centers[0], centers[1], edgeList, polygonArea) / area;


      if(hasTensor){
         eps_xy += transformPolygonElement(dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]), eps_BG_xy,
            GxMat, GyMat, centers[0], centers[1], edgeList, polygonArea) / area;
         eps_yx += transformPolygonElement(dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]), eps_BG_yx,
            GxMat, GyMat, centers[0], centers[1], edgeList, polygonArea) / area;

         im_eps_xy += transformPolygonElement(dcomplex(epsTensor.tensor[3], 0), im_eps_BG_xy,
            GxMat, GyMat, centers[0], centers[1], edgeList, polygonArea) / area;
         im_eps_yx += transformPolygonElement(dcomplex(epsTensor.tensor[5], 0), im_eps_BG_yx,
            GxMat, GyMat, centers[0], centers[1], edgeList, polygonArea) / area;
      }
    }
 }