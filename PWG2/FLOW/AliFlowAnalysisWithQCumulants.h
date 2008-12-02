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
 
//----------------------------------------------------------------------------------------------------------------
//                                            setters and getters                                                 
//----------------------------------------------------------------------------------------------------------------
  TList* GetHistList() const {return this->fHistList;} //output histogram list
 
  void SetIntFlowResults(TH1D* ifr)  {this->fIntFlowResultsQC = ifr;};
  TH1D* GetIntFlowResults() const    {return this->fIntFlowResultsQC;};
  
  void SetDiffFlowResults2nd(TH1D* diff2nd)  {this->fDiffFlowResults2ndOrderQC = diff2nd;};
  TH1D* GetDiffFlowResults2nd() const        {return this->fDiffFlowResults2ndOrderQC;};
  
  void SetDiffFlowResults4th(TH1D* diff4th)  {this->fDiffFlowResults4thOrderQC = diff4th;};
  TH1D* GetDiffFlowResults4th() const        {return this->fDiffFlowResults4thOrderQC;};
  
  void SetCovariances(TH1D* cov)  {this->fCovariances = cov;};
  TH1D* GetCovariances() const    {return this->fCovariances;};
  
  void SetCommonHistsResults2nd(AliFlowCommonHistResults* chr2nd)  {this->fCommonHistsResults2nd = chr2nd;};
  AliFlowCommonHistResults* GetCommonHistsResults2nd() const       {return this->fCommonHistsResults2nd;};
  
  void SetCommonHistsResults4th(AliFlowCommonHistResults* chr4th)  {this->fCommonHistsResults4th = chr4th;};
  AliFlowCommonHistResults* GetCommonHistsResults4th() const       {return this->fCommonHistsResults4th;};
  
  void SetCommonHistsResults6th(AliFlowCommonHistResults* chr6th)  {this->fCommonHistsResults6th = chr6th;};
  AliFlowCommonHistResults* GetCommonHistsResults6th() const       {return this->fCommonHistsResults6th;};
  
  void SetCommonHistsResults8th(AliFlowCommonHistResults* chr8th)  {this->fCommonHistsResults8th = chr8th;};
  AliFlowCommonHistResults* GetCommonHistsResults8th() const       {return this->fCommonHistsResults8th;};
  
  void SetAverageMultiplicity(TProfile* am)        {this->fAvMultIntFlowQC = am;};
  TProfile* GetAverageMultiplicity() const       {return this->fAvMultIntFlowQC;};
  
  void SetQCorrelations(TProfile* QCorr)   {this->fQCorrelations = QCorr;};
  TProfile* GetQCorrelations() const       {return this->fQCorrelations;};
  
  void SetQProduct(TProfile* qp)      {this->fQProduct = qp;};
  TProfile* GetQProduct() const       {return this->fQProduct;};
  
  void SetQVectorComponents(TProfile* qvc)      {this->fQvectorComponents = qvc;};
  TProfile* GetQVectorComponents() const        {return this->fQvectorComponents;};
  
  void SetTwo_1n1nPerBin(TProfile* pb2_1n1n)  {this->f2_1n1n = pb2_1n1n;};
  TProfile* GetTwo_1n1nPerBin() const       {return this->f2_1n1n;};
  
  void SetTwo_2n2nPerBin(TProfile* pb2_2n2n)  {this->f2_2n2n = pb2_2n2n;};
  TProfile* GetTwo_2n2nPerBin() const       {return this->f2_2n2n;};
  
  void SetThree_2n1n1nPerBin(TProfile* pb3_2n1n1n)  {this->f3_2n1n1n = pb3_2n1n1n;};
  TProfile* GetThree_2n1n1nPerBin() const         {return this->f3_2n1n1n;};
  
  void SetThree_1n1n2nPerBin(TProfile* pb3_1n1n2n)  {this->f3_1n1n2n = pb3_1n1n2n;};
  TProfile* GetThree_1n1n2nPerBin() const         {return this->f3_1n1n2n;};
  
  void SetFour_1n1n1n1nPerBin(TProfile* pb4_1n1n1n1n)  {this->f4_1n1n1n1n = pb4_1n1n1n1n;};
  TProfile* GetFour_1n1n1n1nPerBin() const           {return this->f4_1n1n1n1n;}; 
  
  void SetDirectCorrelations(TProfile* dc)     {this->fDirectCorrelations = dc;};
  TProfile* GetDirectCorrelations() const      {return this->fDirectCorrelations;};
//----------------------------------------------------------------------------------------------------------------
 
 private:
  AliFlowAnalysisWithQCumulants(const AliFlowAnalysisWithQCumulants& afawQc);
  AliFlowAnalysisWithQCumulants& operator=(const AliFlowAnalysisWithQCumulants& afawQc);
  
  AliFlowTrackSimple* fTrack;                           //track
  TList*              fHistList;                        //list to hold all output histograms
  TProfile*           fAvMultIntFlowQC;                 //average selected multiplicity (for int. flow)
 
  TProfile*           fQvectorComponents;               //averages of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, ...)
            
  TH1D*               fIntFlowResultsQC;                //integrated flow results from Q-cumulants
  TH1D*               fDiffFlowResults2ndOrderQC;       //differential flow results from 2nd order Q-cumulant
  TH1D*               fDiffFlowResults4thOrderQC;       //differential flow results from 4th order Q-cumulant
  TH1D*               fCovariances;                     //final results for covariances: 1st bin: <2*4>-<2>*<4>, 2nd bin: <2*6>-<2>*<6>, ...
  
  TProfile*                  fQCorrelations;            //multi-particle correlations calculated from Q-vectors 
  TProfile*                  fQProduct;                 //average of products: 1st bin: <2*4>, 2nd bin: <2*6>, ...
  
  TProfile*                  fDirectCorrelations;       //multi-particle correlations calculated with nested loop  
  
  TProfile*                  fReq1n;                    //real part of q-vector evaluated in harmonic n for each pt-bin
  TProfile*                  fImq1n;                    //imaginary part of q-vector evaluated in harmonic n for each pt-bin
  TProfile*                  fReq2n;                    //real part of q-vector evaluated in harmonic 2n for each pt-bin
  TProfile*                  fImq2n;                    //imaginary part of q-vector evaluated in harmonic 2n for each pt-bin

  TProfile*                  f2_1n1n;                   //<<2'>>_{n|n} per pt-bin
  TProfile*                  f2_2n2n;                   //<<2'>>_{2n|2n} per pt-bin
  TProfile*                  f3_2n1n1n;                 //<<3'>>_{2n|n,n} per pt-bin
  TProfile*                  f3_1n1n2n;                 //<<3'>>_{n,n|2n} per pt-bin
  TProfile*                  f4_1n1n1n1n;               //<<4'>>_{n,n|n,n} per pt-bin
 
  AliFlowCommonHist*         fCommonHists;              //common control histograms
  
  AliFlowCommonHistResults*  fCommonHistsResults2nd;    //final results for 2nd order int. and diff. flow stored in the common histograms 
  AliFlowCommonHistResults*  fCommonHistsResults4th;    //final results for 4th order int. and diff. flow stored in the common histograms 
  AliFlowCommonHistResults*  fCommonHistsResults6th;    //final results for 6th order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults*  fCommonHistsResults8th;    //final results for 8th order int. and diff. flow stored in the common histograms
      
  TH1D*                      f2Distribution;            //distribution of <2>_{n|n}
  TH1D*                      f4Distribution;            //distribution of <4>_{n,n|n,n}
  TH1D*                      f6Distribution;            //distribution of <6>_{n,n,n|n,n,n} 
  
  Int_t                      fnBinsPt;                  //number of pt bins

  Double_t                   fPtMin;                    //minimum pt                
  Double_t                   fPtMax;                    //maximum pt
  
  ClassDef(AliFlowAnalysisWithQCumulants, 0);
};

//================================================================================================================

#endif





