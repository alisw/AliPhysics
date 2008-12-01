/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/********************************** 
 * flow analysis with Q-cumulants * 
 *                                * 
 * author:  Ante Bilandzic        * 
 *           (anteb@nikhef.nl)    *
 *********************************/ 

#ifndef AliFlowAnalysisWithQCumulants_H
#define AliFlowAnalysisWithQCumulants_H

#include "AliFlowCommonConstants.h"

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

class AliFlowAnalysisWithQCumulants{
 public:
  AliFlowAnalysisWithQCumulants();
  virtual ~AliFlowAnalysisWithQCumulants(); 
  
  virtual void CreateOutputObjects();
  virtual void Make(AliFlowEventSimple* anEvent);
  virtual void Finish();
  virtual void WriteHistograms(TString* outputFileName);
  
  TList* GetHistList() const {return this->fHistList;} //output histogram list
 
 private:
  AliFlowAnalysisWithQCumulants(const AliFlowAnalysisWithQCumulants& afawQc);
  AliFlowAnalysisWithQCumulants& operator=(const AliFlowAnalysisWithQCumulants& afawQc);
  AliFlowTrackSimple* fTrack;             //track
  TList*              fHistList;          //list to hold all output histograms
  TProfile*           fAvMultIntFlowQC;   //avarage selected multiplicity
 
  TProfile*           fQvectorComponents; //avarages of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, 3rd bin: <(Q_x)^2>, 4th bin: <(Q_y)^2>, 5th bin: <|Q|^2>)
            
  TH1D*               fIntFlowResultsQC;               //integrated flow results from Q-cumulants
  TH1D*               fDiffFlowResults2ndOrderQC;      //differential flow results from 2nd order Q-cumulant
  TH1D*               fDiffFlowResults4thOrderQC;      //differential flow results from 4th order Q-cumulant
  TH1D*               fCovariances;                    //final results for covariances: 1st bin: <2*4>-<2>*<4>, 2nd bin: <2*6>-<2>*<6>, ...
  
  TProfile*           fQCorrelations;                  //multi-particle correlations calculated from Q-vectors 
  TProfile*           fQCovariance;                    //covariance of multi-particle correlations 
  TProfile*           fDirectCorrelations;             //multi-particle correlations calculated with nested loop  
  
  TProfile*           fReD;
  TProfile*           fImD;
  
  TProfile*           fReq1n;                           //real part of q-vector evaluated in harmonic n for each pt-bin
  TProfile*           fImq1n;                           //imaginary part of q-vector evaluated in harmonic n for each pt-bin
  TProfile*           fReq2n;                           //real part of q-vector evaluated in harmonic 2n for each pt-bin
  TProfile*           fImq2n;                           //imaginary part of q-vector evaluated in harmonic 2n for each pt-bin

  TProfile*           f2_1n1n;                          //<<2'>>_{n|n} per pt-bin
  TProfile*           f2_2n2n;                          //<<2'>>_{2n|2n} per pt-bin
  TProfile*           f3_2n1n1n;                        //<<3'>>_{2n|n,n} per pt-bin
  TProfile*           f3_1n1n2n;                        //<<3'>>_{n,n|2n} per pt-bin
  TProfile*           f4_1n1n1n1n;                      //<<4'>>_{n,n|n,n} per pt-bin
  
  TProfile*           fQCorrelationsPerBin;
 
  AliFlowCommonHist*         fCommonHists;              //common control histograms
  
  /*
  AliFlowCommonHistResults*  fCommonHistsResults2nd;    //final results for 2nd order int. and diff. flow stored in the common histograms 
  AliFlowCommonHistResults*  fCommonHistsResults4th;    //final results for 4th order int. and diff. flow stored in the common histograms 
  AliFlowCommonHistResults*  fCommonHistsResults6th;    //final results for 6th order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults*  fCommonHistsResults8th;    //final results for 8th order int. and diff. flow stored in the common histograms
  */
      
  TH1D*               f2Distribution;                   //distribution of <2>_{n|n}
  TH1D*               f4Distribution;                   //distribution of <4>_{n,n|n,n}
  TH1D*               f6Distribution;                   //distribution of <6>_{n,n,n|n,n,n}        
                        
  ClassDef(AliFlowAnalysisWithQCumulants, 0);
};

//================================================================================================================

#endif





