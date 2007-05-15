/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#ifndef ALIHLTPHOSPHYSICSANALYZERSPECTRUM_H
#define ALIHLTPHOSPHYSICSANALYZERSPECTRUM_H

#include "AliHLTPHOSPhysicsAnalyzer.h"
#include "Rtypes.h"


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
 
  virtual void Analyze(AliHLTPHOSClusterDataStruct* clustersPtr[10000], Int_t nClusters);

 private:
  Float_t* fPos0Ptr;                        //! /**<Position of the first cluster*/
  Float_t* fPos1Ptr;                        //! /**</Position of the second cluster*/
  Float_t* fThresholdPtr;                   //! /**<Cut thresholds*/
  Float_t* fEnergyPtr;                      //! /**<Energy of the clusters*/

  ClassDef(AliHLTPHOSPhysicsAnalyzerSpectrum, 1);
  
};

#endif
 
