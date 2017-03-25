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
   // helper function to change a random dielectric to a diagonal
   // @args:
   // epsilon: field
   // type: the type of the dielectric
   /*==============================================*/
   static EpsilonVal toDiagonal(const EpsilonVal epsilon, const EPSTYPE type){
    switch(type){
      case SCALAR_: return fromScalarToDiagonal(epsilon);
      case DIAGONAL_: return epsilon;
      default: throw UTILITY::AttributeNotSupportedException("Cannot convert tensor to diagonal");
    }
   }
   /*==============================================*/
   // helper function to change a random dielectric to a tensor
   // @args:
   // epsilon: field
   // type: the type of the dielectric
   /*==============================================*/
   static EpsilonVal toTensor(const EpsilonVal epsilon, const EPSTYPE type){
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
   // N: the total number of G
   /*==============================================*/
   static RCWAMatrix transformPlanarElement(const dcomplex epsVal, const int N){
     return epsVal * eye<RCWAMatrix>(N, N);
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
     const double width, 
     const int N
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
     const double widthy,
     const int N
   ){
     return widthx * widthy * (epsVal - epsBG)
             * exp(IMAG_I * (GxMat * centerx + GyMat * centery))
             % sinc(GxMat * widthx / 2) % sinc(GyMat * widthy / 2);
   }
   /*==============================================*/
   // This function computes the Fourier transform for planar geometry for tensor case
   // @args:
   // eps_xx_MatrixVec: the Fourier trainsform for eps_xx for all omega
   // eps_xy_MatrixVec: the Fourier trainsform for eps_xy for all omega
   // eps_zx_MatrixVec: the Fourier trainsform for eps_yx for all omega
   // eps_yy_MatrixVec: the Fourier trainsform for eps_yy for all omega
   // im_eps_xx_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_xy_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_yx_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_yy_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_zz_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // eps_zz_Inv_MatrixVec: the inverse of Fourier transform of eps_zz
   // Layer: the layer considered
   // N: the total number of G
   /*==============================================*/
   void transformPlanar(
     RCWAMatricesVec& eps_xx_MatrixVec,
     RCWAMatricesVec& eps_xy_MatrixVec,
     RCWAMatricesVec& eps_yx_MatrixVec,
     RCWAMatricesVec& eps_yy_MatrixVec,
     RCWAMatricesVec& eps_zz_Inv_MatrixVec,
     RCWAMatricesVec& im_eps_xx_MatrixVec,
     RCWAMatricesVec& im_eps_xy_MatrixVec,
     RCWAMatricesVec& im_eps_yx_MatrixVec,
     RCWAMatricesVec& im_eps_yy_MatrixVec,
     RCWAMatricesVec& im_eps_zz_MatrixVec,
     const Ptr<Layer>& layer,
     const int N
   ){
     int numOfOmega = eps_xx_MatrixVec.size();
     Ptr<Material> backGround = layer->getBackGround();
     RCWAMatrix onePadding1N = eye<RCWAMatrix>(N, N);

     for(int i = 0; i < numOfOmega; i++){
       EpsilonVal epsilonBG = backGround->getEpsilonAtIndex(i);
       EpsilonVal epsilonBGTensor = toTensor(epsilonBG, backGround->getType());

       eps_xx_MatrixVec[i].push_back( transformPlanarElement(dcomplex(epsilonBGTensor.tensor[0], epsilonBGTensor.tensor[1]), N) );
       im_eps_xx_MatrixVec[i].push_back( transformPlanarElement(dcomplex(epsilonBGTensor.tensor[1], 0), N) );
       eps_yy_MatrixVec[i].push_back( transformPlanarElement(dcomplex(epsilonBGTensor.tensor[6], epsilonBGTensor.tensor[7]), N) );
       im_eps_yy_MatrixVec[i].push_back( transformPlanarElement(dcomplex(epsilonBGTensor.tensor[7], 0), N) );
       eps_zz_Inv_MatrixVec[i].push_back( transformPlanarElement(REAL_I / dcomplex(epsilonBGTensor.tensor[8], epsilonBGTensor.tensor[9]), N) );
       im_eps_zz_MatrixVec[i].push_back( transformPlanarElement(dcomplex(epsilonBGTensor.tensor[9], 0), N) );
       
       if(layer->hasTensor()){
         eps_xy_MatrixVec[i].push_back( transformPlanarElement(dcomplex(epsilonBGTensor.tensor[2], epsilonBGTensor.tensor[3]), N) );
         im_eps_xy_MatrixVec[i].push_back( transformPlanarElement(dcomplex(epsilonBGTensor.tensor[3], 0), N) );
         eps_yx_MatrixVec[i].push_back( transformPlanarElement(dcomplex(epsilonBGTensor.tensor[4], epsilonBGTensor.tensor[5]), N) );
         im_eps_yx_MatrixVec[i].push_back( transformPlanarElement(dcomplex(epsilonBGTensor.tensor[5], 0), N) );
       }
       else{
         eps_xy_MatrixVec[i].push_back(zeros<RCWAMatrix>(N, N));
         im_eps_xy_MatrixVec[i].push_back(zeros<RCWAMatrix>(N, N));
         eps_yx_MatrixVec[i].push_back(zeros<RCWAMatrix>(N, N));
         im_eps_yx_MatrixVec[i].push_back(zeros<RCWAMatrix>(N, N));
       }
     }
   }

   /*==============================================*/
   // This function computes the Fourier transform for grating geometry
   // @args:
   // eps_xx_MatrixVec: the Fourier trainsform for eps_xx for all omega
   // eps_xy_MatrixVec: the Fourier trainsform for eps_xy for all omega
   // eps_zx_MatrixVec: the Fourier trainsform for eps_yx for all omega
   // eps_yy_MatrixVec: the Fourier trainsform for eps_yy for all omega
   // im_eps_xx_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_xy_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_yx_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_yy_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_zz_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // eps_zz_Inv_MatrixVec: the inverse of Fourier transform of eps_zz
   // Layer: the layer considered
   // N: the total number of G
   // period: the periodicity
   /*==============================================*/
   void transformGratingNaive(
     RCWAMatricesVec& eps_xx_MatrixVec,
     RCWAMatricesVec& eps_xy_MatrixVec,
     RCWAMatricesVec& eps_yx_MatrixVec,
     RCWAMatricesVec& eps_yy_MatrixVec,
     RCWAMatricesVec& eps_zz_Inv_MatrixVec,
     RCWAMatricesVec& im_eps_xx_MatrixVec,
     RCWAMatricesVec& im_eps_xy_MatrixVec,
     RCWAMatricesVec& im_eps_yx_MatrixVec,
     RCWAMatricesVec& im_eps_yy_MatrixVec,
     RCWAMatricesVec& im_eps_zz_MatrixVec,
     const Ptr<Layer>& layer,
     const int N,
     const double period,
     bool useInverseRule
   ){
     // if exists a tensor, no inverser rule
     if(layer->hasTensor()) useInverseRule = false;
            
     int numOfOmega = eps_xx_MatrixVec.size();
     // dcomplex IMAG_I = dcomplex(0, 1);
     RCWAMatrix G_row(1, N), G_col(N, 1);
     for(int i = 0; i < N; i++){
       G_row(0, i) = -i * 2.0 * datum::pi / period;
       G_col(i, 0) = -G_row(0, i);
     }
     RCWAMatrix G_mat = toeplitz(G_col, G_row);

     RCWAMatrix onePadding1N = eye<RCWAMatrix>(N, N);
     int numOfMaterial = layer->getNumOfMaterial();
     RCWAVector centerVec(numOfMaterial), widthVec(numOfMaterial);
     int count = 0;
     for(const_PatternIter it = layer->getArg1Begin(); it != layer->getArg1End(); it++){
       centerVec(count) = it->first;
       widthVec(count) = it->second;
       count++;
     }
     
     Ptr<Material> backGround = layer->getBackGround();
     for(int i = 0; i < numOfOmega; i++){
       RCWAMatrix eps_xx(N, N, fill::zeros), eps_xy(N, N, fill::zeros), eps_yx(N, N, fill::zeros), eps_yy(N, N, fill::zeros), eps_zz_Inv(N, N, fill::zeros);
       RCWAMatrix im_eps_xx(N, N, fill::zeros), im_eps_xy(N, N, fill::zeros), im_eps_yx(N, N, fill::zeros), im_eps_yy(N, N, fill::zeros), im_eps_zz(N, N, fill::zeros);
       count = 0; // reset count
       EpsilonVal epsilonBG = backGround->getEpsilonAtIndex(i);
       EpsilonVal epsBGTensor = toTensor(epsilonBG, backGround->getType());

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

       for(const_MaterialIter it = layer->getVecBegin(); it != layer->getVecEnd(); it++){
         EpsilonVal epsilon = (*it)->getEpsilonAtIndex(i);
         EpsilonVal epsTensor = toTensor(epsilon, (*it)->getType());
   
         eps_xx += transformGratingElement(dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]),
            eps_BG_xx, G_mat, centerVec(count), widthVec(count), N);
         im_eps_xx += transformGratingElement(dcomplex(epsTensor.tensor[1], 0),
            im_eps_BG_xx, G_mat, centerVec(count), widthVec(count), N);

         eps_yy += transformGratingElement(dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]),
            eps_BG_yy, G_mat, centerVec(count), widthVec(count), N);
         im_eps_yy += transformGratingElement(dcomplex(epsTensor.tensor[7], 0),
            im_eps_BG_yy, G_mat, centerVec(count), widthVec(count), N);

         eps_zz_Inv += transformGratingElement(dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]),
            eps_BG_zz, G_mat, centerVec(count), widthVec(count), N);
         im_eps_zz += transformGratingElement(dcomplex(epsTensor.tensor[9], 0),
            im_eps_BG_zz, G_mat, centerVec(count), widthVec(count), N);

        if(layer->hasTensor()){
          eps_xy += transformGratingElement(dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]),
            eps_BG_xy, G_mat, centerVec(count), widthVec(count), N);
          eps_yx += transformGratingElement(dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]),
            eps_BG_yx, G_mat, centerVec(count), widthVec(count), N);
          im_eps_xy += transformGratingElement(dcomplex(epsTensor.tensor[3], 0),
            im_eps_BG_xy, G_mat, centerVec(count), widthVec(count), N);
          im_eps_yx += transformGratingElement(dcomplex(epsTensor.tensor[5], 0),
            im_eps_BG_yx, G_mat, centerVec(count), widthVec(count), N);
        }
         count += 1;
       }

       eps_xx = eps_xx / period + eps_BG_xx * onePadding1N;
       im_eps_xx = im_eps_xx / period + im_eps_BG_xx * onePadding1N;

       eps_yy = eps_yy / period + eps_BG_yy * onePadding1N;
       im_eps_yy = im_eps_yy / period + im_eps_BG_yy * onePadding1N;

       eps_zz_Inv = eps_zz_Inv / period + eps_BG_zz * onePadding1N;
       eps_zz_Inv = eps_zz_Inv.i();
       im_eps_zz = im_eps_zz / period + im_eps_BG_zz * onePadding1N;

       if(layer->hasTensor()){
        eps_xy = eps_xy / period + eps_BG_xy * onePadding1N;
        eps_yx = eps_yx / period + eps_BG_yx * onePadding1N;
        im_eps_xy = im_eps_xy / period + im_eps_BG_xy * onePadding1N;
        im_eps_yx = im_eps_yx / period + im_eps_BG_yx * onePadding1N;
       }

       eps_xx_MatrixVec[i].push_back(eps_xx);
       eps_xy_MatrixVec[i].push_back(eps_xy);
       eps_yx_MatrixVec[i].push_back(eps_yx);
       eps_yy_MatrixVec[i].push_back(eps_yy);
       eps_zz_Inv_MatrixVec[i].push_back(eps_zz_Inv);
       im_eps_xx_MatrixVec[i].push_back(im_eps_xx);
       im_eps_xy_MatrixVec[i].push_back(im_eps_xy);
       im_eps_yx_MatrixVec[i].push_back(im_eps_yx);
       im_eps_yy_MatrixVec[i].push_back(im_eps_yy);
       im_eps_zz_MatrixVec[i].push_back(im_eps_zz);
     }
   }


   void transformGratingDiagonalAdaptive(){

   }

   /*==============================================*/
   // This function computes the Fourier transform for rectangle case
   // @args:
   // eps_xx_MatrixVec: the Fourier trainsform for eps_xx for all omega
   // eps_xy_MatrixVec: the Fourier trainsform for eps_xy for all omega
   // eps_zx_MatrixVec: the Fourier trainsform for eps_yx for all omega
   // eps_yy_MatrixVec: the Fourier trainsform for eps_yy for all omega
   // im_eps_xx_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_xy_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_yx_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_yy_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_zz_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // eps_zz_Inv_MatrixVec: the inverse of Fourier transform of eps_zz
   // Layer: the layer considered
   // nGx: the total number of G in x direction
   // nGy: the total number of G in y direction
   // period: the periodicity
   // useInverseRule: whether use inverse rule
   /*==============================================*/
   void transformRectangle(
        RCWAMatricesVec& eps_xx_MatrixVec,
        RCWAMatricesVec& eps_xy_MatrixVec,
        RCWAMatricesVec& eps_yx_MatrixVec,
        RCWAMatricesVec& eps_yy_MatrixVec,
        RCWAMatricesVec& eps_zz_Inv_MatrixVec,
        RCWAMatricesVec& im_eps_xx_MatrixVec,
        RCWAMatricesVec& im_eps_xy_MatrixVec,
        RCWAMatricesVec& im_eps_yx_MatrixVec,
        RCWAMatricesVec& im_eps_yy_MatrixVec,
        RCWAMatricesVec& im_eps_zz_MatrixVec,
        const Ptr<Layer>& layer,
        const int nGx,
        const int nGy,
        const double* period,
        bool useInverseRule
  ){
     if(layer->hasTensor()) useInverseRule = false;
     int numOfOmega = eps_xx_MatrixVec.size();
     int N = RCWA::getN(nGx, nGy);

     RCWAMatrix Gx_mat, Gy_mat;
     RCWA::getGMatrices(nGx, nGy, period, Gx_mat, Gy_mat, TWO_);

     // dcomplex IMAG_I = dcomplex(0, 1.0);
     RCWAMatrix Gx_r, Gx_l, Gy_r, Gy_l;
     meshGrid(Gx_mat, Gx_mat, Gx_r, Gx_l);
     meshGrid(Gy_mat, Gy_mat, Gy_r, Gy_l);

     RCWAMatrix GxMat = Gx_l - Gx_r;
     RCWAMatrix GyMat = Gy_l - Gy_r;

     RCWAMatrix onePadding1N = eye<RCWAMatrix>(N, N);
     int numOfMaterial = layer->getNumOfMaterial();

     RCWAVector centerxVec(numOfMaterial), centeryVec(numOfMaterial), widthxVec(numOfMaterial), widthyVec(numOfMaterial);
     int count = 0;// reset count
     for(const_PatternIter it = layer->getArg1Begin(); it != layer->getArg1End(); it++){
       centerxVec(count) = it->first;
       centeryVec(count) = it->second;
       count++;
     }
     count = 0;
     for(const_PatternIter it = layer->getArg2Begin(); it != layer->getArg2End(); it++){
       widthxVec(count) = it->first;
       widthyVec(count) = it->second;
       count++;
     }

     Ptr<Material> backGround = layer->getBackGround();
     for(int i = 0; i < numOfOmega; i++){
       RCWAMatrix eps_xx(N, N, fill::zeros), eps_xy(N, N, fill::zeros), eps_yx(N, N, fill::zeros), eps_yy(N, N, fill::zeros), eps_zz_Inv(N, N, fill::zeros);
       RCWAMatrix im_eps_xx(N, N, fill::zeros), im_eps_xy(N, N, fill::zeros), im_eps_yx(N, N, fill::zeros), im_eps_yy(N, N, fill::zeros), im_eps_zz(N, N, fill::zeros);
       count = 0;
       EpsilonVal epsilonBG = backGround->getEpsilonAtIndex(i);
       EpsilonVal epsBGTensor = toTensor(epsilonBG, backGround->getType());

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


       for(const_MaterialIter it = layer->getVecBegin(); it != layer->getVecEnd(); it++){
         EpsilonVal epsilon = (*it)->getEpsilonAtIndex(i);
         EpsilonVal epsTensor = toTensor(epsilon, (*it)->getType());

         eps_xx += transformRectangleElement(dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]), eps_BG_xx,
            GxMat, GyMat, centerxVec(count), centeryVec(count), widthxVec(count), widthyVec(count), N);
         im_eps_xx += transformRectangleElement(dcomplex(epsTensor.tensor[1], 0), im_eps_BG_xx,
            GxMat, GyMat, centerxVec(count), centeryVec(count), widthxVec(count), widthyVec(count), N);

         eps_yy += transformRectangleElement(dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]), eps_BG_yy,
            GxMat, GyMat, centerxVec(count), centeryVec(count), widthxVec(count), widthyVec(count), N);
         im_eps_yy += transformRectangleElement(dcomplex(epsTensor.tensor[7], 0), im_eps_BG_yy,
            GxMat, GyMat, centerxVec(count), centeryVec(count), widthxVec(count), widthyVec(count), N);

         eps_zz_Inv += transformRectangleElement(dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]), eps_BG_zz,
            GxMat, GyMat, centerxVec(count), centeryVec(count), widthxVec(count), widthyVec(count), N);
         im_eps_zz += transformRectangleElement(dcomplex(epsTensor.tensor[9], 0), im_eps_BG_zz,
            GxMat, GyMat, centerxVec(count), centeryVec(count), widthxVec(count), widthyVec(count), N); 
        

         if(layer->hasTensor()){
            eps_xy += transformRectangleElement(dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]), eps_BG_xy,
            GxMat, GyMat, centerxVec(count), centeryVec(count), widthxVec(count), widthyVec(count), N);
            eps_yx += transformRectangleElement(dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]), eps_BG_yx,
            GxMat, GyMat, centerxVec(count), centeryVec(count), widthxVec(count), widthyVec(count), N);

            im_eps_xy += transformRectangleElement(dcomplex(epsTensor.tensor[3], 0), im_eps_BG_xy,
            GxMat, GyMat, centerxVec(count), centeryVec(count), widthxVec(count), widthyVec(count), N);
            im_eps_yx += transformRectangleElement(dcomplex(epsTensor.tensor[5], 0), im_eps_BG_yx,
            GxMat, GyMat, centerxVec(count), centeryVec(count), widthxVec(count), widthyVec(count), N);
         }
         count += 1;
       }

       eps_xx = eps_xx / (period[0] * period[1]) + eps_BG_xx * onePadding1N;
       im_eps_xx = im_eps_xx / (period[0] * period[1]) + im_eps_BG_xx * onePadding1N;

       eps_yy = eps_yy / (period[0] * period[1]) + eps_BG_yy * onePadding1N;
       im_eps_yy = im_eps_yy / (period[0] * period[1]) + im_eps_BG_yy * onePadding1N;

       eps_zz_Inv = eps_zz_Inv / (period[0] * period[1]) + eps_BG_zz * onePadding1N;
       eps_zz_Inv = eps_zz_Inv.i();
       im_eps_zz = im_eps_zz / (period[0] * period[1]) + im_eps_BG_zz * onePadding1N;

       if(layer->hasTensor()){
          eps_xy = eps_xy / (period[0] * period[1]) + eps_BG_xy * onePadding1N;
          eps_yx = eps_yx / (period[0] * period[1]) + eps_BG_yx * onePadding1N;
          im_eps_xy = im_eps_xy / (period[0] * period[1]) + im_eps_BG_xy * onePadding1N;
          im_eps_yx = im_eps_yx / (period[0] * period[1]) + im_eps_BG_yx * onePadding1N;
       }

       eps_xx_MatrixVec[i].push_back(eps_xx);
       eps_xy_MatrixVec[i].push_back(eps_xy);
       eps_yx_MatrixVec[i].push_back(eps_yx);
       eps_yy_MatrixVec[i].push_back(eps_yy);
       eps_zz_Inv_MatrixVec[i].push_back(eps_zz_Inv);
       im_eps_xx_MatrixVec[i].push_back(im_eps_xx);
       im_eps_xy_MatrixVec[i].push_back(im_eps_xy);
       im_eps_yx_MatrixVec[i].push_back(im_eps_yx);
       im_eps_yy_MatrixVec[i].push_back(im_eps_yy);
       im_eps_zz_MatrixVec[i].push_back(im_eps_zz);
    }
  }

   /*==============================================*/
   // This function computes the Fourier transform for circular case
   // @args:
   // eps_xx_MatrixVec: the Fourier trainsform for eps_xx for all omega
   // eps_xy_MatrixVec: the Fourier trainsform for eps_xy for all omega
   // eps_zx_MatrixVec: the Fourier trainsform for eps_yx for all omega
   // eps_yy_MatrixVec: the Fourier trainsform for eps_yy for all omega
   // im_eps_xx_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_xy_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_yx_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_yy_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // im_eps_zz_MatrixVec: the Fourier trainsform for imaginary part for all omega
   // eps_zz_Inv_MatrixVec: the inverse of Fourier transform of eps_zz
   // Layer: the layer considered
   // nGx: the total number of G in x direction
   // nGy: the total number of G in y direction
   // period: the periodicity
   /*==============================================*/
   void transformCircle(
    RCWAMatricesVec& eps_xx_MatrixVec,
    RCWAMatricesVec& eps_xy_MatrixVec,
    RCWAMatricesVec& eps_yx_MatrixVec,
    RCWAMatricesVec& eps_yy_MatrixVec,
    RCWAMatricesVec& eps_zz_Inv_MatrixVec,
    RCWAMatricesVec& im_eps_xx_MatrixVec,
    RCWAMatricesVec& im_eps_xy_MatrixVec,
    RCWAMatricesVec& im_eps_yx_MatrixVec,
    RCWAMatricesVec& im_eps_yy_MatrixVec,
    RCWAMatricesVec& im_eps_zz_MatrixVec,
    const Ptr<Layer>& layer,
    const int nGx,
    const int nGy,
    const double* period
   ){
    // TODO
   }
 }