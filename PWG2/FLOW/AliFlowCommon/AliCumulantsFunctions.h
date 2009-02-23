/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/********************************** 
 * functions and equations needed * 
 * for calculation of cumulants   *
 * and final flow estimates       *
 *                                *   
 * author: Ante Bilandzic         * 
 *          (anteb@nikhef.nl)     *
 *********************************/ 

#ifndef ALICUMULANTSFUNCTIONS_H
#define ALICUMULANTSFUNCTIONS_H

#include "AliFlowCommonConstants.h"
#include "AliFlowCumuConstants.h"

class TH1;
class TProfile;
class TProfile2D;
class TProfile3D;

class TObjArray;
class TList;
class TFile;

class AliFlowCommonHistResults;

//================================================================================================================

class AliCumulantsFunctions{
 public:
  AliCumulantsFunctions();
  virtual ~AliCumulantsFunctions();
  //AliCumulantsFunctions(TProfile2D *intGenFun, TProfile3D *diffGenFunRe, TProfile3D *diffGenFunIm, TH1D *ifr, TH1D *dfr2, TH1D *dfr4, TH1D *dfr6, TH1D *dfr8, TProfile *avMult, TProfile *fQVector, TH1D *fQDistrib, TProfile2D *fdRe0, TProfile2D *fdRe1, TProfile2D *fdRe2, TProfile2D *fdRe3, TProfile2D *fdRe4, TProfile2D *fdRe5, TProfile2D *fdRe6, TProfile2D *fdRe7, TProfile2D *fdIm0, TProfile2D *fdIm1, TProfile2D *fdIm2, TProfile2D *fdIm3, TProfile2D *fdIm4, TProfile2D *fdIm5, TProfile2D *fdIm6, TProfile2D *fdIm7);
  
  AliCumulantsFunctions(TProfile2D *intGenFun, TProfile2D *intGenFun4, TProfile2D *intGenFun6, TProfile2D *intGenFun8, TProfile2D *intGenFun16, TProfile *avMult4, TProfile *avMult6, TProfile *avMult8, TProfile *avMult16, TProfile3D *diffPtRPGenFunRe, TProfile3D *diffPtRPGenFunIm, TProfile *ptBinRPNoOfParticles, TProfile3D *diffEtaRPGenFunRe, TProfile3D *diffEtaRPGenFunIm, TProfile *etaBinRPNoOfParticles, TProfile3D *diffPtPOIGenFunRe, TProfile3D *diffPtPOIGenFunIm, TProfile *ptBinPOINoOfParticles, TProfile3D *diffEtaPOIGenFunRe, TProfile3D *diffEtaPOIGenFunIm, TProfile *etaBinPOINoOfParticles, TH1D *ifr, TH1D *dfr2, TH1D *dfr4, TH1D *dfr6, TH1D *dfr8, TProfile *avMult, TProfile *qVector, TProfile *averageOfSquaredWeight, AliFlowCommonHistResults *chr2nd, AliFlowCommonHistResults *chr4th, AliFlowCommonHistResults *chr6th, AliFlowCommonHistResults *chr8th, AliFlowCommonHist *ch);
 
  void Calculate();

 private:
  AliCumulantsFunctions(const AliCumulantsFunctions& fun);
  AliCumulantsFunctions& operator=(const AliCumulantsFunctions& fun);
  
  TProfile2D *fIntGenFun;    //average value of generating function for int. flow
  
  TProfile2D *fIntGenFun4;   //average value of generating function for int. flow (only for other system of Eq.)
  TProfile2D *fIntGenFun6;   //average value of generating function for int. flow (only for other system of Eq.)
  TProfile2D *fIntGenFun8;   //average value of generating function for int. flow (only for other system of Eq.)
  TProfile2D *fIntGenFun16;  //average value of generating function for int. flow (only for other system of Eq.)
  
  TProfile *fAvMult4;        //average selected multiplicity for int. flow (only for other system of Eq.)
  TProfile *fAvMult6;        //average selected multiplicity for int. flow (only for other system of Eq.)
  TProfile *fAvMult8;        //average selected multiplicity for int. flow (only for other system of Eq.)
  TProfile *fAvMult16;       //average selected multiplicity for int. flow (only for other system of Eq.) 
  
  TProfile3D *fDiffPtRPGenFunRe;   //average value of generating function for diff. flow in pt (real part)
  TProfile3D *fDiffPtRPGenFunIm;   //average value of generating function for diff. flow in pt (imaginary part)
  TProfile *fPtBinRPNoOfParticles; //number of particles per pt bin
  
  TProfile3D *fDiffEtaRPGenFunRe;   //average value of generating function for diff. flow in eta (real part)
  TProfile3D *fDiffEtaRPGenFunIm;   //average value of generating function for diff. flow in eta (imaginary part)
  TProfile *fEtaBinRPNoOfParticles; //number of particles per eta bin
 
  TProfile3D *fDiffPtPOIGenFunRe;   //average value of generating function for diff. flow in pt (real part)
  TProfile3D *fDiffPtPOIGenFunIm;   //average value of generating function for diff. flow in pt (imaginary part)
  TProfile *fPtBinPOINoOfParticles; //number of particles per pt bin
  
  TProfile3D *fDiffEtaPOIGenFunRe;   //average value of generating function for diff. flow in eta (real part)
  TProfile3D *fDiffEtaPOIGenFunIm;   //average value of generating function for diff. flow in eta (imaginary part)
  TProfile *fEtaBinPOINoOfParticles; //number of particles per eta bin
    
  TH1D *fifr;                //integrated flow final results 
  TH1D *fdfr2;               //differential flow final results
  TH1D *fdfr4;               //differential flow final results
  TH1D *fdfr6;               //differential flow final results
  TH1D *fdfr8;               //differential flow final results 
   
  TProfile *fAvMult;         //average selected multiplicity for int. flow
  TProfile *fQVector;        //average values of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, 3rd bin: <(Q_x)^2>, 4th bin: <(Q_y)^2>) 
  
  TProfile *fAverageOfSquaredWeight;    //<w^2>
  
  AliFlowCommonHistResults *fchr2nd; //final results for 2nd order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults *fchr4th; //final results for 4th order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults *fchr6th; //final results for 6th order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults *fchr8th; //final results for 8th order int. and diff. flow stored in the common histograms
  
  AliFlowCommonHist *fch; //common control histogram  
  
  /*
  TProfile2D *fdRe0,*fdRe1,*fdRe2,*fdRe3,*fdRe4,*fdRe5,*fdRe6,*fdRe7;//differential flow 
  TProfile2D *fdIm0,*fdIm1,*fdIm2,*fdIm3,*fdIm4,*fdIm5,*fdIm6,*fdIm7;//differential flow
  */
  
  ClassDef(AliCumulantsFunctions, 0);
};

//================================================================================================================

#endif





