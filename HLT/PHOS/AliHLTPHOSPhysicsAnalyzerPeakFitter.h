#ifndef ALIHLTPHOSPHYSICSANALYZERPEAKFITTER_H
#define ALIHLTPHOSPHYSICSANALYZERPEAKFITTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "Rtypes.h"
#include "AliHLTPHOSClusterDataStruct.h"
#include "AliHLTDataTypes.h"
#include "TH1F.h"
#include "TMath.h"


class AliHLTPHOSPhysicsAnalyzerPeakFitter
{

 public:

  AliHLTPHOSPhysicsAnalyzerPeakFitter();
  virtual ~AliHLTPHOSPhysicsAnalyzerPeakFitter();
  AliHLTPHOSPhysicsAnalyzerPeakFitter(const AliHLTPHOSPhysicsAnalyzerPeakFitter &);
  AliHLTPHOSPhysicsAnalyzerPeakFitter & operator = (const AliHLTPHOSPhysicsAnalyzerPeakFitter &) {return *this;}

  void    SetHistogram(TH1F* histPtr)            { fRootHistPtr = histPtr; }

  Int_t   FitGaussian();
  Int_t   FitLorentzian();
  

 private:
  
  Float_t fGainLow;
  Float_t fGainHigh;
  TH1F* fRootHistPtr;
 
  ClassDef(AliHLTPHOSPhysicsAnalyzerPeakFitter, 1);

};






#endif
