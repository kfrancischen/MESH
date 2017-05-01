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
   // G_mat: the G_mat
   // width: the width of the grating
   /*==============================================*/
   static RCWArMatrix transformGratingElement(
     const RCWArMatrix& G_mat,
     const double width
   ){
    return width * RCWA::sinc(G_mat / 2 * width);
   }

   /*==============================================*/
   // helper function to do fourier transform for one value
   // @args:
   // Gx_mat: the Gx_mat
   // Gy_mat: the Gy_mat
   // widthx: the width of the rectangle in x direction
   // widthy: the width of the rectangle in y direction
   /*==============================================*/
   static RCWArMatrix transformRectangleElement(
     const RCWArMatrix& GxMat,
     const RCWArMatrix& GyMat,
     const double widthx,
     const double widthy
   ){
     return widthx * widthy
             * RCWA::sinc( GxMat * widthx / 2 )
             % RCWA::sinc( GyMat * widthy / 2 );
   }

   /*==============================================*/
   // helper function to do fourier transform for one value
   // @args:
   // Gx_mat: the Gx_mat
   // Gy_mat: the Gy_mat
   // radius: the radius of the circle
   /*==============================================*/
   static RCWArMatrix transformCircleElement(
     const RCWArMatrix& GxMat,
     const RCWArMatrix& GyMat,
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

     return 2 * datum::pi * POW2(radius) * jincMat;

   }
   /*==============================================*/
   // helper function to do fourier transform for one value
   // @args:
   // Gx_mat: the Gx_mat
   // Gy_mat: the Gy_mat
   // a: the half width of the rectangle in x direction
   // b: the half width of the rectangle in y direction
   /*==============================================*/
   static RCWArMatrix transformEllipseElement(
     const RCWArMatrix& GxMat,
     const RCWArMatrix& GyMat,
     const double a,
     const double b
   ){
     RCWArMatrix rho = sqrt( square(a * GxMat) + square(b * GyMat) );
     RCWArMatrix jincMat = zeros<RCWArMatrix>( size(rho) );

     RCWArMatrix::iterator jincMat_it = jincMat.begin();
     int count = 0;
     for(RCWArMatrix::const_iterator it = rho.begin(); it != rho.end(); it++){
       *(jincMat_it + count) = jinc(*it);
       count++;
     }

     return 2 * datum::pi * a * b * jincMat;
   }
   /*==============================================*/
   // helper function to do fourier transform for one value
   // @args:
   // Gx_mat: the Gx_mat
   // Gy_mat: the Gy_mat
   // edgeList: the list of edges
   // area: the area of the polygon
   /*==============================================*/
   static RCWAcMatrix transformPolygonElement(
     const RCWArMatrix& GxMat,
     const RCWArMatrix& GyMat,
     const EdgeList& edgeList,
     const double area
   ){
     RCWAcMatrix result = zeros<RCWAcMatrix>( size(GxMat) );
     RCWArMatrix::const_iterator GyMat_it = GyMat.cbegin();
     RCWAcMatrix::iterator result_it = result.begin();
     int count = 0;
     for(RCWArMatrix::const_iterator GxMat_it = GxMat.cbegin(); GxMat_it != GxMat.cend(); GxMat_it++){
       double u = (*GxMat_it), v = *(GyMat_it + count);
       if(u == 0.0 && v == 0.0){
         *(result_it + count) = area;
       }
       else{
         for(size_t i = 0; i < edgeList.size(); i++){
           double x_cur = edgeList[i].first, y_cur = edgeList[i].second;
           double x_next, y_next;
           if(i == edgeList.size() - 1){
             x_next = edgeList[0].first;
             y_next = edgeList[0].second;
           }
           else{
             x_next = edgeList[i+1].first;
             y_next = edgeList[i+1].second;
           }

           // case when u = 0
           if( u == 0.0 ){
             *(result_it + count) += IMAG_I/v * (x_next - x_cur) * exp(IMAG_I * (u*(x_next+x_cur)/2 + v*(y_next+y_cur)/2))
                * RCWA::sinc((x_next-x_cur)*u/2 + (y_next-y_cur)*v/2);
           }
           else{
             *(result_it + count) += -IMAG_I/u * (y_next - y_cur) * exp(IMAG_I * (u*(x_next+x_cur)/2 + v*(y_next+y_cur)/2))
                * RCWA::sinc((x_next-x_cur)*u/2 + (y_next-y_cur)*v/2);
           }
         }
       }
       count++;
     }
     return result;
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
   // Gx_mat: the Gx matrix
   // center: the center of the grating
   // width: the width of the grating
   // area: the area of one periodicity
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
    const RCWArMatrix& Gx_Mat,
    const double center,
    const double width,
    const double area,
    const bool hasTensor
   ){

     RCWArMatrix Gx_r, Gx_l, Gy_r, Gy_l;
     meshGrid(Gx_Mat, Gx_Mat, Gx_r, Gx_l);

     RCWArMatrix G_mat = Gx_l - Gx_r;
     RCWAcMatrix phase = exp(IMAG_I * G_mat * center);

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
     RCWArMatrix geoMat = transformGratingElement(G_mat, width);
     eps_xx += (dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]) - eps_BG_xx) / area * phase % geoMat;

     im_eps_xx += (dcomplex(epsTensor.tensor[1], 0) - im_eps_BG_xx) / area * phase % geoMat;

     eps_yy += (dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]) - eps_BG_yy) / area * phase % geoMat;

     im_eps_yy += (dcomplex(epsTensor.tensor[7], 0) - im_eps_BG_yy) / area * phase % geoMat;

     eps_zz += (dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]) - eps_BG_zz) / area * phase % geoMat;

     im_eps_zz +=(dcomplex(epsTensor.tensor[9], 0) - im_eps_BG_zz) / area * phase % geoMat;

     if(hasTensor){
       eps_xy += (dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]) - eps_BG_xy) / area * phase % geoMat;

       eps_yx += (dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]) - eps_BG_yx) / area * phase % geoMat;

       im_eps_xy += (dcomplex(epsTensor.tensor[3], 0) - im_eps_BG_xy) / area * phase % geoMat;

       im_eps_yx += (dcomplex(epsTensor.tensor[5], 0) - im_eps_BG_yx) / area * phase % geoMat;
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
   // Gx_mat: the Gx matrix
   // Gy_mat: the Gy matrix
   // centers: the centers of the rectangle
   // angle: the rotated angle with respect to x axis
   // widths: the widths of the rectangle
   // area: the area of one periodicity
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
    const RCWArMatrix& Gx_Mat,
    const RCWArMatrix& Gy_Mat,
    const double centers[2],
    const double angle,
    const double widths[2],
    const double area,
    const bool hasTensor
   ){


     RCWArMatrix Gx_r, Gx_l, Gy_r, Gy_l;
     meshGrid(Gx_Mat, Gx_Mat, Gx_r, Gx_l);
     meshGrid(Gy_Mat, Gy_Mat, Gy_r, Gy_l);

     RCWArMatrix GxMat = Gx_l - Gx_r;
     RCWArMatrix GyMat = Gy_l - Gy_r;

     RCWAcMatrix phase = exp(IMAG_I * (GxMat * centers[0] + GyMat * centers[1]));
     RCWArMatrix G_temp = GxMat * cos(angle) + GyMat * sin(angle);
     GyMat = -GxMat * sin(angle) + GyMat * cos(angle);
     GxMat = G_temp;

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
     RCWArMatrix geoMat = transformRectangleElement(GxMat, GyMat, widths[0], widths[1]);

     eps_xx += (dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]) - eps_BG_xx) / area * phase % geoMat;

     im_eps_xx += (dcomplex(epsTensor.tensor[1], 0) - im_eps_BG_xx) / area * phase % geoMat;

     eps_yy += (dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]) - eps_BG_yy) / area * phase % geoMat;

     im_eps_yy += (dcomplex(epsTensor.tensor[7], 0) - im_eps_BG_yy) / area * phase % geoMat;

     eps_zz += (dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]) - eps_BG_zz) / area * phase % geoMat;

     im_eps_zz += (dcomplex(epsTensor.tensor[9], 0) - im_eps_BG_zz) / area * phase % geoMat;


     if(hasTensor){
       eps_xy += (dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]) - eps_BG_xy) / area * phase % geoMat;

       eps_yx += (dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]) - eps_BG_yx) / area * phase % geoMat;

       im_eps_xy += (dcomplex(epsTensor.tensor[3], 0) - im_eps_BG_xy) / area * phase % geoMat;

       im_eps_yx += (dcomplex(epsTensor.tensor[5], 0) - im_eps_BG_yx) / area * phase % geoMat;
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
   // Gx_mat: the Gx matrix
   // Gy_mat: the Gy matrix
   // centers: the centers of the circle
   // angle: the rotated angle with respect to x axis
   // radius: the radius of the circle
   // area: the area of one periodicity
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
    const RCWArMatrix& Gx_Mat,
    const RCWArMatrix& Gy_Mat,
    const double centers[2],
    const double radius,
    const double area,
    const bool hasTensor
   ){

     // dcomplex IMAG_I = dcomplex(0, 1.0);
     RCWArMatrix Gx_r, Gx_l, Gy_r, Gy_l;
     meshGrid(Gx_Mat, Gx_Mat, Gx_r, Gx_l);
     meshGrid(Gy_Mat, Gy_Mat, Gy_r, Gy_l);

     RCWArMatrix GxMat = Gx_l - Gx_r;
     RCWArMatrix GyMat = Gy_l - Gy_r;
     RCWAcMatrix phase = exp(IMAG_I * (GxMat * centers[0] + GyMat * centers[1]));


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
     RCWArMatrix geoMat = transformCircleElement(GxMat, GyMat, radius);

     eps_xx += (dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]) - eps_BG_xx) / area * phase % geoMat;

     im_eps_xx += (dcomplex(epsTensor.tensor[1], 0) - im_eps_BG_xx) / area * phase % geoMat;

     eps_yy += (dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]) - eps_BG_yy) / area * phase % geoMat;

     im_eps_yy += (dcomplex(epsTensor.tensor[7], 0) - im_eps_BG_yy) / area * phase % geoMat;

     eps_zz += (dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]) - eps_BG_zz) / area * phase % geoMat;

     im_eps_zz += (dcomplex(epsTensor.tensor[9], 0) - im_eps_BG_zz) / area * phase % geoMat;

     if(hasTensor){
       eps_xy += (dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]) - eps_BG_xy) / area * phase % geoMat;

       eps_yx += (dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]) - eps_BG_yx) / area * phase % geoMat;

       im_eps_xy += (dcomplex(epsTensor.tensor[3], 0) - im_eps_BG_xy) / area * phase % geoMat;

       im_eps_yx += (dcomplex(epsTensor.tensor[5], 0) - im_eps_BG_yx) / area * phase % geoMat;
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
    // Gx_mat: the Gx matrix
    // Gy_mat: the Gy matrix
    // centers: the centers of the ellipse
    // angle: the rotated angle with respect to x axis
    // halfwidths: the halfwidths of the ellipse
    // area: the area of one periodicity
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
     const RCWArMatrix& Gx_Mat,
     const RCWArMatrix& Gy_Mat,
     const double centers[2],
     const double angle,
     const double halfwidths[2],
     const double area,
     const bool hasTensor
    ){
      // dcomplex IMAG_I = dcomplex(0, 1.0);
      RCWArMatrix Gx_r, Gx_l, Gy_r, Gy_l;
      meshGrid(Gx_Mat, Gx_Mat, Gx_r, Gx_l);
      meshGrid(Gy_Mat, Gy_Mat, Gy_r, Gy_l);

      RCWArMatrix GxMat = Gx_l - Gx_r;
      RCWArMatrix GyMat = Gy_l - Gy_r;
      RCWAcMatrix phase = exp(IMAG_I * (GxMat * centers[0] + GyMat * centers[1]));
      RCWArMatrix G_temp = GxMat * cos(angle) + GyMat * sin(angle);
      GyMat = -GxMat * sin(angle) + GyMat * cos(angle);
      GxMat = G_temp;


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
      RCWArMatrix geoMat = transformEllipseElement(GxMat, GyMat, halfwidths[0], halfwidths[1]);

      eps_xx += (dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]) - eps_BG_xx) / area * phase % geoMat;

      im_eps_xx += (dcomplex(epsTensor.tensor[1], 0) - im_eps_BG_xx) / area * phase % geoMat;

      eps_yy += (dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]) - eps_BG_yy) / area * phase % geoMat;

      im_eps_yy += (dcomplex(epsTensor.tensor[7], 0) - im_eps_BG_yy) / area * phase % geoMat;

      eps_zz += (dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]) - eps_BG_zz) / area * phase % geoMat;

      im_eps_zz += (dcomplex(epsTensor.tensor[9], 0) - im_eps_BG_zz) / area * phase % geoMat;

      if(hasTensor){
        eps_xy += (dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]) - eps_BG_xy) / area * phase % geoMat;
        eps_yx += (dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]) - eps_BG_yx) / area * phase % geoMat;

        im_eps_xy += (dcomplex(epsTensor.tensor[3], 0) - im_eps_BG_xy) / area * phase % geoMat;

        im_eps_yx += (dcomplex(epsTensor.tensor[5], 0) - im_eps_BG_yx) / area * phase % geoMat;
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
     // Gx_mat: the Gx matrix
     // Gy_mat: the Gy matrix
     // centers: the centers of the polygon
     // angle: the rotated angle with respect to x axis
     // edgeList: the edges of the polygon
     // area: the area of one periodicity
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
      const RCWArMatrix& Gx_Mat,
      const RCWArMatrix& Gy_Mat,
      const double centers[2],
      const double angle,
      const EdgeList& edgeList,
      const double area,
      const bool hasTensor
    ){
      double polygonArea = getPolygonArea(edgeList);
      // dcomplex IMAG_I = dcomplex(0, 1.0);
      RCWArMatrix Gx_r, Gx_l, Gy_r, Gy_l;
      meshGrid(Gx_Mat, Gx_Mat, Gx_r, Gx_l);
      meshGrid(Gy_Mat, Gy_Mat, Gy_r, Gy_l);

      RCWArMatrix GxMat = Gx_l - Gx_r;
      RCWArMatrix GyMat = Gy_l - Gy_r;
      RCWAcMatrix phase = exp(IMAG_I * (GxMat * centers[0] + GyMat * centers[1]));
      RCWArMatrix G_temp = GxMat * cos(angle) + GyMat * sin(angle);
      GyMat = -GxMat * sin(angle) + GyMat * cos(angle);
      GxMat = G_temp;


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
      RCWAcMatrix geoMat = transformPolygonElement(GxMat, GyMat, edgeList, polygonArea);

      eps_xx += (dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]) - eps_BG_xx) / area * phase % geoMat;

      im_eps_xx += (dcomplex(epsTensor.tensor[1], 0) - im_eps_BG_xx) / area * phase % geoMat;

      eps_yy += (dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]) - eps_BG_yy) / area * phase % geoMat;

      im_eps_yy += (dcomplex(epsTensor.tensor[7], 0) - im_eps_BG_yy) / area * phase % geoMat;

      eps_zz += (dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]) - eps_BG_zz) / area * phase % geoMat;

      im_eps_zz += (dcomplex(epsTensor.tensor[9], 0) - im_eps_BG_zz) / area * phase % geoMat;

      if(hasTensor){
        eps_xy += (dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]) - eps_BG_xy) / area * phase % geoMat;

        eps_yx += (dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]) - eps_BG_yx) / area * phase % geoMat;

        im_eps_xy += (dcomplex(epsTensor.tensor[3], 0) - im_eps_BG_xy) / area * phase % geoMat;

        im_eps_yx += (dcomplex(epsTensor.tensor[5], 0) - im_eps_BG_yx) / area * phase % geoMat;
      }
    }
 }