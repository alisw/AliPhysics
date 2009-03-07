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

#ifndef ALIQCUMULANTSFUNCTIONS_H
#define ALIQCUMULANTSFUNCTIONS_H

#include "AliFlowCommonConstants.h"

class TH1;
class TProfile;
class TProfile2D;
class TProfile3D;

class TObjArray;
class TList;
class TFile;

class AliFlowCommonHistResults;

//================================================================================================================

class AliQCumulantsFunctions{
 public:
  AliQCumulantsFunctions();
  virtual ~AliQCumulantsFunctions();
  AliQCumulantsFunctions(TH1D *intRes, TH1D *diffRes2nd, TH1D *diffRes4th, TH1D *covar, TProfile *AvMult, TProfile *QVector, TProfile *QCorr, TProfile *QCorrW, TProfile *QProd, TProfile *Direct, TProfile *binPt2p1n1nRP, TProfile *binPt2p2n2nRP, TProfile *binPt3p2n1n1nRP, TProfile *binPt3p1n1n2nRP, TProfile *binPt4p1n1n1n1nRP, TProfile *binEta2p1n1nRP, TProfile *binEta2p2n2nRP, TProfile *binEta3p2n1n1nRP, TProfile *binEta3p1n1n2nRP, TProfile *binEta4p1n1n1n1nRP, TProfile *binPt2p1n1nPOI, TProfile *binPt2p2n2nPOI, TProfile *binPt3p2n1n1nPOI, TProfile *binPt3p1n1n2nPOI, TProfile *binPt4p1n1n1n1nPOI, TProfile *binEta2p1n1nPOI, TProfile *binEta2p2n2nPOI, TProfile *binEta3p2n1n1nPOI, TProfile *binEta3p1n1n2nPOI, TProfile *binEta4p1n1n1n1nPOI, TProfile *binWPt2p1n1nPOI, TProfile *binWPt2p2n2nPOI, TProfile *binWPt3p2n1n1nPOI, TProfile *binWPt3p1n1n2nPOI, TProfile *binWPt4p1n1n1n1nPOI, TProfile *binWEta2p1n1nPOI, TProfile *binWEta4p1n1n1n1nPOI,  TProfile *binWPt2p1n1nRP, TProfile *binWPt4p1n1n1n1nRP,  TProfile *binWEta2p1n1nRP, TProfile *binWEta4p1n1n1n1nRP, AliFlowCommonHist *ch2nd, AliFlowCommonHist *ch4th, AliFlowCommonHist *ch6th, AliFlowCommonHist *ch8th, AliFlowCommonHistResults *chr2nd, AliFlowCommonHistResults *chr4th, AliFlowCommonHistResults *chr6th, AliFlowCommonHistResults *chr8th);
 
  void Calculate();

 private:
  AliQCumulantsFunctions(const AliQCumulantsFunctions& Qfun);
  AliQCumulantsFunctions& operator=(const AliQCumulantsFunctions& Qfun);
  
  TH1D     *fIntRes;          //results for integrated flow
  TH1D     *fDiffRes2nd;      //results for differential flow (2nd order)
  TH1D     *fDiffRes4th;      //results for differential flow (4th order)
  TH1D     *fCovar;           //results for covariances (1st bin: <2*4>-<2>*<4>, 2nd bin: <2*6>-<2>*<6>, ...)  
          
  TProfile *fAvMult;          //avarage selected multiplicity for int. flow
  TProfile *fQVector;         //avarage values of Q-vector components
  TProfile *fQCorr;           //multi-particle correlations calculated from Q-vectors 
  TProfile *fQCorrW;   //weighted multi-particle correlations calculated from Q-vectors 
  TProfile *fQProd;           //average of products: 1st bin: <2*4>, 2nd bin: <2*6>, ...

  TProfile *fDirect;          //direct correlations (correlations calculated with nested loopps)
  
  TProfile *fbinPt2p1n1nRP;     //<<2'>>_{n|n} per pt-bin
  TProfile *fbinPt2p2n2nRP;     //<<2'>>_{2n|2n} per pt-bin
  TProfile *fbinPt3p2n1n1nRP;   //<<3'>>_{2n,n|n} per pt-bin
  TProfile *fbinPt3p1n1n2nRP;   //<<3'>>_{n,n|2n} per pt-bin
  TProfile *fbinPt4p1n1n1n1nRP; //<<4'>>_{n,n|n,n} per pt-bin
  
  TProfile *fbinEta2p1n1nRP;     //<<2'>>_{n|n} per eta-bin
  TProfile *fbinEta2p2n2nRP;     //<<2'>>_{2n|2n} per eta-bin
  TProfile *fbinEta3p2n1n1nRP;   //<<3'>>_{2n,n|n} per eta-bin
  TProfile *fbinEta3p1n1n2nRP;   //<<3'>>_{n,n|2n} per eta-bin
  TProfile *fbinEta4p1n1n1n1nRP; //<<4'>>_{n,n|n,n} per eta-bin  
  
  TProfile *fbinPt2p1n1nPOI;     //<<2'>>_{n|n} per pt-bin
  TProfile *fbinPt2p2n2nPOI;     //<<2'>>_{2n|2n} per pt-bin
  TProfile *fbinPt3p2n1n1nPOI;   //<<3'>>_{2n,n|n} per pt-bin
  TProfile *fbinPt3p1n1n2nPOI;   //<<3'>>_{n,n|2n} per pt-bin
  TProfile *fbinPt4p1n1n1n1nPOI; //<<4'>>_{n,n|n,n} per pt-bin
  
  TProfile *fbinEta2p1n1nPOI;     //<<2'>>_{n|n} per eta-bin
  TProfile *fbinEta2p2n2nPOI;     //<<2'>>_{2n|2n} per eta-bin
  TProfile *fbinEta3p2n1n1nPOI;   //<<3'>>_{2n,n|n} per eta-bin
  TProfile *fbinEta3p1n1n2nPOI;   //<<3'>>_{n,n|2n} per eta-bin
  TProfile *fbinEta4p1n1n1n1nPOI; //<<4'>>_{n,n|n,n} per eta-bin  
  
  
  
  
  
  
  TProfile *fbinWPt2p1n1nPOI;     //<<2'>>_{n|n} per pt-bin
  TProfile *fbinWPt2p2n2nPOI;     //<<2'>>_{2n|2n} per pt-bin
  TProfile *fbinWPt3p2n1n1nPOI;   //<<3'>>_{2n,n|n} per pt-bin
  TProfile *fbinWPt3p1n1n2nPOI;   //<<3'>>_{n,n|2n} per pt-bin
  TProfile *fbinWPt4p1n1n1n1nPOI; //<<4'>>_{n,n|n,n} per pt-bin
  
  TProfile *fbinWEta2p1n1nPOI;     //<<2'>>_{n|n} per pt-bin
  TProfile *fbinWEta4p1n1n1n1nPOI; //<<4'>>_{n,n|n,n} per pt-bin
  
  TProfile *fbinWPt2p1n1nRP;     //<<2'>>_{n|n} per pt-bin
  TProfile *fbinWPt4p1n1n1n1nRP; //<<4'>>_{n,n|n,n} per pt-bin
  
  TProfile *fbinWEta2p1n1nRP;     //<<2'>>_{n|n} per pt-bin
  TProfile *fbinWEta4p1n1n1n1nRP; //<<4'>>_{n,n|n,n} per pt-bin
  
  
  
  
  
  
  AliFlowCommonHist *fch2nd; //common control histograms (taking into account only the events with 2 and more particles) 
  AliFlowCommonHist *fch4th; //common control histograms (taking into account only the events with 4 and more particles) 
  AliFlowCommonHist *fch6th; //common control histograms (taking into account only the events with 6 and more particles) 
  AliFlowCommonHist *fch8th; //common control histograms (taking into account only the events with 8 and more particles) 
  
  AliFlowCommonHistResults *fchr2nd; //final results for 2nd order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults *fchr4th; //final results for 4th order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults *fchr6th; //final results for 6th order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults *fchr8th; //final results for 8th order int. and diff. flow stored in the common histograms
              
  ClassDef(AliQCumulantsFunctions, 0);
};

//================================================================================================================

#endif





