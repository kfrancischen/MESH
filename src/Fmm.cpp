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
   static EpsilonVal fromScalarToDiagonal(EpsilonVal epsilon){
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
   static EpsilonVal fromDiagonalToTensor(EpsilonVal epsilon){
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
   }
   /*==============================================*/
   // helper function to change a scalar dielectric to a tensor
   // @args:
   // epsilon: scalar field
   /*==============================================*/
   static EpsilonVal fromScalarToTensor(EpsilonVal epsilon){
     EpsilonVal diag = fromScalarToTensor(epsilon);
     return fromDiagonalToTensor(diag);
   }

   /*==============================================*/
   // helper function to change a random dielectric to a tensor
   // @args:
   // epsilon: field
   // type: the type of the dielectric
   /*==============================================*/
   static EpsilonVal toTensor(EpsilonVal epsilon, EPSTYPE type){
     switch (type) {
       case SCALAR_: return fromScalarToTensor(epsilon);
       case DIAGONAL_: return fromDiagonalToTensor(epsilon);
       default: return epsilon;
     }
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
     Ptr<Layer>& layer,
     const int N
   ){
     int numOfOmega = eps_xx_MatrixVec.size();
     Ptr<Material> backGround = layer->getBackGround();
     RCWAMatrix onePadding1N = eye<RCWAMatrix>(N, N);

     for(int i = 0; i < numOfOmega; i++){
       EpsilonVal epsilonBG = backGround->getEpsilonAtIndex(i);

       EpsilonVal epsilonBGTensor = toTensor(epsilonBG, backGround->getType());
       eps_xx_MatrixVec[i].push_back(onePadding1N * dcomplex(epsilonBG.tensor[0], epsilonBG.tensor[1]));
       im_eps_xx_MatrixVec[i].push_back(onePadding1N * epsilonBG.tensor[1]);
       eps_xy_MatrixVec[i].push_back(onePadding1N * dcomplex(epsilonBG.tensor[2], epsilonBG.tensor[3]));
       im_eps_xy_MatrixVec[i].push_back(onePadding1N * epsilonBG.tensor[3]);
       eps_yx_MatrixVec[i].push_back(onePadding1N * dcomplex(epsilonBG.tensor[4], epsilonBG.tensor[5]));
       im_eps_yx_MatrixVec[i].push_back(onePadding1N * epsilonBG.tensor[5]);
       eps_yy_MatrixVec[i].push_back(onePadding1N * dcomplex(epsilonBG.tensor[6], epsilonBG.tensor[7]));
       im_eps_yy_MatrixVec[i].push_back(onePadding1N * epsilonBG.tensor[7]);
       eps_zz_MatrixVec[i].push_back(onePadding1N * dcomplex(1.0, 0) / dcomplex(epsilonBG.tensor[8], epsilonBG.tensor[9]));
       im_eps_zz_MatrixVec[i].push_back(onePadding1N * epsilonBG.tensor[9]);
     }
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
     Ptr<Layer>& layer,
     const int N,
     const double period,
     bool useInverseRule
   ){
     // if exists a tensor, no inverser rule
     if(layer->hasTensor()) useInverseRule = false;

     int numOfOmega = eps_xx_MatrixVec.size();
     dcomplex IMAG_I = dcomplex(0, 1);
     RCWAMatrix G_row(1, N), G_col(N, 1);
     for(int i = 0; i < N; i++){
       G_row(0, i) = -i * 2.0 * datum::pi / period_[0];
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
       RCWAMatrix im_eps_xx(N, N, fill:zeros), im_eps_xy(N, N, fill:zeros), im_eps_yx(N, N, fill:zeros), im_eps_yy(N, N, fill::zeros), im_eps_zz(N, N, fill::zeros);
       count = 0; // reset count
       EpsilonVal epsilonBG = backGround->getEpsilonAtIndex(i);

       EpsilonVal epsBGTensor = toTensor(epsilonBG, backGround->getType());

       dcomplex eps_BG_xx = dcomplex(epsBGTensor.tensor[0], epsBGTensor.tensor[1]);
       dcomplex eps_BG_xy = dcomplex(epsBGTensor.tensor[2], epsBGTensor.tensor[3]);
       dcomplex eps_BG_yx = dcomplex(epsBGTensor.tensor[4], epsBGTensor.tensor[5]);
       dcomplex eps_BG_yy = dcomplex(epsBGTensor.tensor[6], epsBGTensor.tensor[7]);
       dcomplex eps_BG_zz = dcomplex(epsBGTensor.tensor[8], epsBGTensor.tensor[9]);
       dcomplex im_eps_BG_xx = dcomplex(0, epsBGTensor.tensor[1]);
       dcomplex im_eps_BG_xy = dcomplex(0, epsBGTensor.tensor[3]);
       dcomplex im_eps_BG_yx = dcomplex(0, epsBGTensor.tensor[5]);
       dcomplex im_eps_BG_yy = dcomplex(0, epsBGTensor.tensor[7]);
       dcomplex im_eps_BG_zz = dcomplex(0, epsBGTensor.tensor[9]);


       for(const_MaterialIter it = layer->getVecBegin(); it != layer->getVecEnd(); it++){
         EpsilonVal epsilon = (*it)->getEpsilonAtIndex(i);
         EpsilonVal epsTensor = toTensor(epsTensor. (*it)->getType());

         if(useInverseRule){
           eps_xx += exp(IMAG_I * G_mat * centerVec(count)) * (dcomplex(1.0, 0) / dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]) - dcomplex(1.0, 0) / eps_BG_xx)
             * widthVec(count) % sinc(G_mat / 2 * widthVec(count));
         }
         else{
           eps_xx += exp(IMAG_I * G_mat * centerVec(count)) * (dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]) - eps_BG_xx)
             * widthVec(count) % sinc(G_mat / 2 * widthVec(count));
         }

         eps_xy += exp(IMAG_I * G_mat * centerVec(count)) * (dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]) - eps_BG_xy)
           * widthVec(count) % sinc(G_mat / 2 * widthVec(count));
         eps_yx += exp(IMAG_I * G_mat * centerVec(count)) * (dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]) - eps_BG_yx)
           * widthVec(count) % sinc(G_mat / 2 * widthVec(count));
         eps_yy += exp(IMAG_I * G_mat * centerVec(count)) * (dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]) - eps_BG_yy)
           * widthVec(count) % sinc(G_mat / 2 * widthVec(count));
         eps_zz_Inv += exp(IMAG_I * G_mat * centerVec(count)) * (dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]) - eps_BG_zz)
           * widthVec(count) % sinc(G_mat / 2 * widthVec(count));


         im_eps_xx += exp(IMAG_I * G_mat * centerVec(count)) * (epsTensor.tensor[1] - im_eps_BG_xx)
           * widthVec(count) % sinc(G_mat / 2 * widthVec(count));
         im_eps_xy += exp(IMAG_I * G_mat * centerVec(count)) * (epsTensor.tensor[3] - im_eps_BG_xy)
           * widthVec(count) % sinc(G_mat / 2 * widthVec(count));
         im_eps_yx += exp(IMAG_I * G_mat * centerVec(count)) * (epsTensor.tensor[5] - im_eps_BG_yx)
           * widthVec(count) % sinc(G_mat / 2 * widthVec(count));
         im_eps_yy += exp(IMAG_I * G_mat * centerVec(count)) * (epsTensor.tensor[7] - im_eps_BG_yy)
           * widthVec(count) % sinc(G_mat / 2 * widthVec(count));
         im_eps_zz += exp(IMAG_I * G_mat * centerVec(count)) * (epsTensor.tensor[9] - im_eps_BG_zz)
           * widthVec(count) % sinc(G_mat / 2 * widthVec(count));

         count += 1;
       }

       if(useInverseRule){
         eps_xx = eps_xx / period + dcomplex(1.0, 0) / eps_BG_xx * onePadding1N;
         eps_xx = eps_xx.i();
       }
       else{
         eps_xx = eps_xx / period + eps_BG_xx * onePadding1N;
       }
       eps_xy = eps_xy / period + eps_BG_xy * onePadding1N;
       eps_yx = eps_yx / period + eps_BG_yx * onePadding1N;
       eps_yy = eps_yy / period + eps_BG_yy * onePadding1N;
       eps_zz_Inv = eps_zz_Inv / period + eps_BG_zz * onePadding1N;
       eps_zz_Inv = eps_zz_Inv.i();

       im_eps_xx = im_eps_xx / period + im_eps_BG_xx * onePadding1N;
       im_eps_xy = im_eps_xy / period + im_eps_BG_xy * onePadding1N;
       im_eps_yx = im_eps_yx / period + im_eps_BG_yx * onePadding1N;
       im_eps_yy = im_eps_yy / period + im_eps_BG_yy * onePadding1N;
       im_eps_zz = im_eps_zz / period + im_eps_BG_zz * onePadding1N;

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

   void ransformGratingTensorAdaptive(){

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
        Ptr<Layer>& layer,
        const int N,
        const double* period,
        bool useInverseRule
  ){
     if(laye->hasTensor()) useInverseRule = false;

     dcomplex IMAG_I = dcomplex(0, 1.0);
     RCWAMatrix Gx_r, Gx_l, Gy_r, Gy_l;
     meshGrid(Gx_mat_, Gx_mat_, Gx_r, Gx_l);
     meshGrid(Gy_mat_, Gy_mat_, Gy_r, Gy_l);

     RCWAMatrix GxMat = Gx_l - Gx_r;
     RCWAMatrix GyMat = Gy_l - Gy_r;

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
       RCWAMatrix im_eps_xx(N, N, fill:zeros), im_eps_xy(N, N, fill:zeros), im_eps_yx(N, N, fill:zeros), im_eps_yy(N, N, fill::zeros), im_eps_zz(N, N, fill::zeros);
       count = 0;
       EpsilonVal epsilonBG = backGround->getEpsilonAtIndex(i);
       EpsilonVal epsBGTensor = toTensor(epsilonBG, backGround->getType());

       dcomplex eps_BG_xx = dcomplex(epsBGTensor.tensor[0], epsBGTensor.tensor[1]);
       dcomplex eps_BG_xy = dcomplex(epsBGTensor.tensor[2], epsBGTensor.tensor[3]);
       dcomplex eps_BG_yx = dcomplex(epsBGTensor.tensor[4], epsBGTensor.tensor[5]);
       dcomplex eps_BG_yy = dcomplex(epsBGTensor.tensor[6], epsBGTensor.tensor[7]);
       dcomplex eps_BG_zz = dcomplex(epsBGTensor.tensor[8], epsBGTensor.tensor[9]);
       dcomplex im_eps_BG_xx = dcomplex(0, epsBGTensor.tensor[1]);
       dcomplex im_eps_BG_xy = dcomplex(0, epsBGTensor.tensor[3]);
       dcomplex im_eps_BG_yx = dcomplex(0, epsBGTensor.tensor[5]);
       dcomplex im_eps_BG_yy = dcomplex(0, epsBGTensor.tensor[7]);
       dcomplex im_eps_BG_zz = dcomplex(0, epsBGTensor.tensor[9]);


       for(const_MaterialIter it = layer->getVecBegin(); it != layer->getVecEnd(); it++){
         EpsilonVal epsilon = (*it)->getEpsilonAtIndex(i);
         EpsilonVal epsTensor = toTensor(epsTensor. (*it)->getType());

         if(useInverseRule){
           eps_xx += widthxVec(count) * widthyVec(count) / (period[0] * period[1]) * (dcomplex(1.0, 0) / dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]) - dcomplex(1.0, 0) / eps_BG_xx)
             * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
             % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);
         }
         else{
           eps_xx += widthxVec(count) * widthyVec(count) / (period[0] * period[1]) * (dcomplex(epsTensor.tensor[0], epsTensor.tensor[1]) - eps_BG_xx)
             * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
             % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);
         }

         eps_xy += widthxVec(count) * widthyVec(count) / (period[0] * period[1]) * (dcomplex(epsTensor.tensor[2], epsTensor.tensor[3]) - eps_BG_xy)
           * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
           % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);
         eps_yx += widthxVec(count) * widthyVec(count) / (period[0] * period[1]) * (dcomplex(epsTensor.tensor[4], epsTensor.tensor[5]) - eps_BG_yx)
           * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
           % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);
         eps_yy += widthxVec(count) * widthyVec(count) / (period[0] * period[1]) * (dcomplex(epsTensor.tensor[6], epsTensor.tensor[7]) - eps_BG_yy)
           * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
           % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);
         eps_zz_Inv += widthxVec(count) * widthyVec(count) / (period[0] * period[1]) * (dcomplex(epsTensor.tensor[8], epsTensor.tensor[9]) - eps_BG_zz)
           * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
           % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);

         im_eps_xx += widthxVec(count) * widthyVec(count) / (period[0] * period[1]) * (epsTensor.tensor[1] - im_eps_BG_xx)
           * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
           % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);
         im_eps_xy += widthxVec(count) * widthyVec(count) / (period[0] * period[1]) * (epsTensor.tensor[3] - im_eps_BG_xy)
           * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
           % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);
         im_eps_yx += widthxVec(count) * widthyVec(count) / (period[0] * period[1]) * (epsTensor.tensor[5] - im_eps_BG_yx)
           * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
           % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);
         im_eps_yy += widthxVec(count) * widthyVec(count) / (period[0] * period[1]) * (epsTensor.tensor[7] - im_eps_BG_yy)
           * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
           % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);
         im_eps_zz += widthxVec(count) * widthyVec(count) / (period[0] * period[1]) * (epsTensor.tensor[9] - im_eps_BG_zz)
           * exp(IMAG_I * (GxMat * centerxVec(count) + GyMat * centeryVec(count)))
           % sinc(GxMat * widthxVec(count) / 2) % sinc(GyMat * widthyVec(count) / 2);

         count += 1;
       }
       if(useInverseRule){
         eps_xx += dcomplex(1.0, 0) / eps_BG_xx * onePadding1N;
         eps_xx = eps_xx.i();
       }
       else{
         eps_xx += eps_BG_xx * onePadding1N;
       }
       eps_xy += eps_BG_xy * onePadding1N;
       eps_yx += eps_BG_yx * onePadding1N;
       eps_yy += eps_BG_yy * onePadding1N;
       eps_zz_Inv += eps_BG_zz * onePadding1N;
       eps_zz_Inv = eps_zz_Inv.i();

       im_eps_xx += im_eps_BG_xx * onePadding1N;
       im_eps_xy += im_eps_BG_xy * onePadding1N;
       im_eps_yx += im_eps_BG_yx * onePadding1N;
       im_eps_yy += im_eps_BG_yy * onePadding1N;
       im_eps_zz += im_eps_BG_zz * onePadding1N;

       eps_xx_MatrixVec[i].push_back(eps_xx);
       eps_xy_MatrixVec[i].push_back(eps_xy);
       eps_yx_MatrixVec[i].push_back(eps_yx);
       eps_yy_MatrixVec[i].push_back(eps_yy);
       eps_zz_Inv_MatrixVec[i].push_back(eps_zz_Inv);
       im_eps_xx_MatrixVec[i].push_back(im_eps_xx);
       im_eps_xy_MatrixVec[i].push_back(im_eps_xy);
       im_eps_yx_MatrixVec[i].push_back(im_eps_yx);
       im_eps_yy_MatrixVec.push_back(im_eps_yy);
       im_eps_zz_MatrixVec[i].push_back(im_eps_zz);

   }

   void transformCircleDiagonal(
     RCWAMatricesVec& eps_xx_MatrixVec,
     RCWAMatricesVec& eps_yy_MatrixVec,
     RCWAMatricesVec& eps_zz_Inv_MatrixVec,
     RCWAMatricesVec& im_eps_xx_MatrixVec,
     RCWAMatricesVec& im_eps_yy_MatrixVec,
     RCWAMatricesVec& im_eps_zz_MatrixVec,
     Ptr<Layer>& layer,
     const int N
   ){
     // TODO
   }

   void transformCircleTensor(
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
     Ptr<Layer>& layer,
     const int N
   ){
    // TODO
   }
 }