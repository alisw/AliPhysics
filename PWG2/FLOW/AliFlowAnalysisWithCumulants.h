/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/******************************** 
 * flow analysis with cumulants * 
 *                              * 
 * author: Ante Bilandzic       * 
 *          (anteb@nikhef.nl)   *
 *******************************/ 

#ifndef AliFlowAnalysisWithCumulants_H
#define AliFlowAnalysisWithCumulants_H

#include "AliFlowCommonConstants.h"
#include "AliFlowCumuConstants.h"

class TObjArray;
class TList;
class TFile;

class TH1;
class TProfile;
class TProfile2D;
class TProfile3D;

class AliFlowEventSimple;
class AliFlowTrackSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;
class AliFlowVector;

//================================================================================================================

class AliFlowAnalysisWithCumulants{
 public:
  AliFlowAnalysisWithCumulants();
  virtual ~AliFlowAnalysisWithCumulants(); 
  
  virtual void CreateOutputObjects();
  virtual void Make(AliFlowEventSimple* anEvent);
  virtual void Finish();
  virtual void WriteHistograms(TString* outputFileName);
  
  TList* GetHistList() const {return this->fHistList;}      //output histogram list
 
 private:
  AliFlowAnalysisWithCumulants(const AliFlowAnalysisWithCumulants& afawc);
  AliFlowAnalysisWithCumulants& operator=(const AliFlowAnalysisWithCumulants& afawc);
  AliFlowTrackSimple* fTrack;                               //track
  static const Int_t fgkQmax=AliFlowCumuConstants::kQmax;   //needed for numerics
  static const Int_t fgkPmax=AliFlowCumuConstants::kPmax;   //needed for numerics  
  static const Int_t fgkFlow=AliFlowCumuConstants::kFlow;   //integrated flow coefficient to be calculated
  static const Int_t fgkMltpl=AliFlowCumuConstants::kMltpl; //the multiple in p=m*n (diff. flow) 
  static const Int_t fgknBins=100;                          //number of pt bins
  TList* fHistList;                                         //list to hold all output histograms

  Double_t fR0;       //needed for numerics
  Double_t fPtMax;    //maximum pt
  Double_t fPtMin;    //minimum pt
  Double_t fBinWidth; //width of pt bin (in GeV)
      
  Double_t fAvQx;  //<Q_x>
  Double_t fAvQy;  //<Q_y>
  Double_t fAvQ2x; //<(Q_x)^2>
  Double_t fAvQ2y; //<(Q_y)^2>

  TProfile*          fAvMultIntFlow;     //avarage selected multiplicity
 
  TProfile*          fQVectorComponents; //avarages of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, 3rd bin: <(Q_x)^2>, 4th bin: <(Q_y)^2>)
    
  TH1D*              fIntFlowResults;   //integrated flow final results
  
  TH1D*              fDiffFlowResults2; //differential flow final results (2nd order estimate) 
  TH1D*              fDiffFlowResults4; //differential flow final results (4th order estimate)
  TH1D*              fDiffFlowResults6; //differential flow final results (6th order estimate)
  TH1D*              fDiffFlowResults8; //differential flow final results (8th order estimate)
  
  TProfile2D*        fIntFlowGenFun;    //avarage of the generating function for integrated flow 
  TProfile3D*        fDiffFlowGenFunRe; //avarage of the generating function for differential flow (real part)
  TProfile3D*        fDiffFlowGenFunIm; //avarage of the generating function for differential flow (imaginary part)
  
  TProfile2D *fDiffFlowGenFunRe0,*fDiffFlowGenFunRe1,*fDiffFlowGenFunRe2,*fDiffFlowGenFunRe3;//differential flow 
  TProfile2D *fDiffFlowGenFunRe4,*fDiffFlowGenFunRe5,*fDiffFlowGenFunRe6,*fDiffFlowGenFunRe7;//differential flow
  TProfile2D *fDiffFlowGenFunIm0,*fDiffFlowGenFunIm1,*fDiffFlowGenFunIm2,*fDiffFlowGenFunIm3;//differential flow
  TProfile2D *fDiffFlowGenFunIm4,*fDiffFlowGenFunIm5,*fDiffFlowGenFunIm6,*fDiffFlowGenFunIm7;//differential flow
 
  AliFlowCommonHist* fCommonHists;      //common control histograms
  
  TH1D*              fQDist;            //q-distribution
      
  ClassDef(AliFlowAnalysisWithCumulants, 0);
};

//================================================================================================================

#endif





