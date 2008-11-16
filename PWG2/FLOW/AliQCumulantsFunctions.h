/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/********************************** 
 * functions and equations needed * 
 * for calculation of Q-cumulants *
 * and final flow estimates       *
 *                                *   
 * author:  Ante Bilandzic        * 
 *           (anteb@nikhef.nl)    *
 *********************************/ 

#ifndef AliQCumulantsFunctions_H
#define AliQCumulantsFunctions_H

#include "AliFlowCommonConstants.h"

class TH1;
class TProfile;
class TProfile2D;
class TProfile3D;

class TObjArray;
class TList;
class TFile;

//================================================================================================================

class AliQCumulantsFunctions{
 public:
  AliQCumulantsFunctions();
  virtual ~AliQCumulantsFunctions();
  AliQCumulantsFunctions(TH1D *intRes, TH1D *diffRes2nd, TH1D *diffRes4th, TProfile *AvMult, TProfile *QVector, TProfile *QCorr, TProfile *QCovar, TProfile *QCorrPerBin, TProfile *Direct, TProfile *bin2_1n1n, TProfile *bin2_2n2n, TProfile *bin3_2n1n1n, TProfile *bin3_1n1n2n, TProfile *bin4_1n1n1n1n);
 
  void Calculate();

 private:
  AliQCumulantsFunctions(const AliQCumulantsFunctions& Qfun);
  AliQCumulantsFunctions& operator=(const AliQCumulantsFunctions& Qfun);
  
  TH1D     *fIntRes;         //results for integrated flow
  TH1D     *fDiffRes2nd;     //results for differential flow (2nd order)
  TH1D     *fDiffRes4th;     //results for differential flow (4th order)
          
  TProfile *fAvMult;         //avarage selected multiplicity for int. flow
  TProfile *fQVector;        //avarage values of Q-vector components
  TProfile *fQCorr;          //QCorrelations
  TProfile *fQCovar;         //QCorrelations

  TProfile *fQCorrPerBin;    //QCorrelationsPerBin
  TProfile *fDirect;         //direct correlations
  
  TProfile *fbin2_1n1n;
  TProfile *fbin2_2n2n;
  TProfile *fbin3_2n1n1n;
  TProfile *fbin3_1n1n2n;
  TProfile *fbin4_1n1n1n1n;
        
  ClassDef(AliQCumulantsFunctions, 0);
};

//================================================================================================================

#endif





