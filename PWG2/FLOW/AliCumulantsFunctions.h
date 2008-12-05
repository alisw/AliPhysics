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

#ifndef AliCumulantsFunctions_H
#define AliCumulantsFunctions_H

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
  //AliCumulantsFunctions(TProfile2D *IntGenFun, TProfile3D *DiffGenFunRe, TProfile3D *DiffGenFunIm, TH1D *ifr, TH1D *dfr2, TH1D *dfr4, TH1D *dfr6, TH1D *dfr8, TProfile *AvMult, TProfile *fQVector, TH1D *fQDistrib, TProfile2D *fdRe0, TProfile2D *fdRe1, TProfile2D *fdRe2, TProfile2D *fdRe3, TProfile2D *fdRe4, TProfile2D *fdRe5, TProfile2D *fdRe6, TProfile2D *fdRe7, TProfile2D *fdIm0, TProfile2D *fdIm1, TProfile2D *fdIm2, TProfile2D *fdIm3, TProfile2D *fdIm4, TProfile2D *fdIm5, TProfile2D *fdIm6, TProfile2D *fdIm7);
  
  AliCumulantsFunctions(TProfile2D *IntGenFun, TProfile3D *DiffGenFunRe, TProfile3D *DiffGenFunIm, TProfile *BinNoOfParticles, TH1D *ifr, TH1D *dfr2, TH1D *dfr4, TH1D *dfr6, TH1D *dfr8, TProfile *AvMult, TProfile *QVector, AliFlowCommonHistResults *chr2nd, AliFlowCommonHistResults *chr4th, AliFlowCommonHistResults *chr6th, AliFlowCommonHistResults *chr8th);
 
  void Calculate();

 private:
  AliCumulantsFunctions(const AliCumulantsFunctions& fun);
  AliCumulantsFunctions& operator=(const AliCumulantsFunctions& fun);
  
  TProfile2D *fIntGenFun;    //average value of generating function for int. flow
  TProfile3D *fDiffGenFunRe; //average value of generating function for diff. flow (real part)
  TProfile3D *fDiffGenFunIm; //average value of generating function for diff. flow (imaginary part)
  
  TProfile *fBinNoOfParticles; //number of particles per pt bin
    
  TH1D *fifr;                //integrated flow final results 
  TH1D *fdfr2;               //differential flow final results
  TH1D *fdfr4;               //differential flow final results
  TH1D *fdfr6;               //differential flow final results
  TH1D *fdfr8;               //differential flow final results 
   
  TProfile *fAvMult;         //average selected multiplicity for int. flow
  TProfile *fQVector;        //average values of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, 3rd bin: <(Q_x)^2>, 4th bin: <(Q_y)^2>) 
  
  AliFlowCommonHistResults *fchr2nd; //final results for 2nd order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults *fchr4th; //final results for 4th order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults *fchr6th; //final results for 6th order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults *fchr8th; //final results for 8th order int. and diff. flow stored in the common histograms
  
  /*
  TProfile2D *fdRe0,*fdRe1,*fdRe2,*fdRe3,*fdRe4,*fdRe5,*fdRe6,*fdRe7;//differential flow 
  TProfile2D *fdIm0,*fdIm1,*fdIm2,*fdIm3,*fdIm4,*fdIm5,*fdIm6,*fdIm7;//differential flow
  */
  
  ClassDef(AliCumulantsFunctions, 0);
};

//================================================================================================================

#endif





