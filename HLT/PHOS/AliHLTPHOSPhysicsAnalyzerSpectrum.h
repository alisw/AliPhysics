#ifndef ALIHLTPHOSPHYSICSANALYZERSPECTRUM_H
#define ALIHLTPHOSPHYSICSANALYZERSPECTRUM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



#include "AliHLTPHOSPhysicsAnalyzer.h"


const Int_t USE_DISTANCE_CUT = 0;

class AliHLTPHOSPhysicsAnalyzerSpectrum : public AliHLTPHOSPhysicsAnalyzer
{
 public:
  AliHLTPHOSPhysicsAnalyzerSpectrum();
  AliHLTPHOSPhysicsAnalyzerSpectrum(const AliHLTPHOSPhysicsAnalyzerSpectrum &);
  AliHLTPHOSPhysicsAnalyzerSpectrum & operator = (const AliHLTPHOSPhysicsAnalyzerSpectrum)
    {
      return *this; 
    }

  virtual ~AliHLTPHOSPhysicsAnalyzerSpectrum();

  Int_t SetThreshold(Float_t photonEnergy0, Float_t photonEnergy1);
  Float_t EvalDistance();
  Float_t EvalCutDistance(Float_t cutMass);
 
  virtual void Analyze(AliHLTPHOSClusterDataStruct* clustersPtr[10000], Int_t nClusters);

 private:
  Float_t* fPos0Ptr;
  Float_t* fPos1Ptr;
  Float_t* fThresholdPtr;
  Float_t* fEnergyPtr;
  Bool_t fUseDistanceCut;

  ClassDef(AliHLTPHOSPhysicsAnalyzerSpectrum, 1);
  
};


#endif
