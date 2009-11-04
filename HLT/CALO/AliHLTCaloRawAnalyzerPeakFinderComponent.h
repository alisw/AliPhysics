//-*- Mode: C++ -*-
// $Id: AliHLTPHOSRawAnalyzerPeakFinderComponent.h 29824 2008-11-10 13:43:55Z richterm $

#ifndef ALIHLTCALORAWANALYZERPEAKFINDERCOMPONENT_H
#define ALIHLTCALORAWANALYZERPEAKFINDERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 



#include "AliHLTCaloRawAnalyzerComponentv3.h"

class AliHLTCaloRawAnalyzerPeakFinderComponent: public AliHLTCaloRawAnalyzerComponentv3
{
 public:
  AliHLTCaloRawAnalyzerPeakFinderComponent(TString det);
  virtual ~AliHLTCaloRawAnalyzerPeakFinderComponent();

  virtual int Deinit();
  virtual const char* GetComponentID() = 0;
  virtual AliHLTComponent* Spawn() = 0;
 private:
  //  Bool_t LoadPFVector(); 
  //  Bool_t LoadPFVector(int startindex, int Nsamples, int tau, int fs);
  
  // virtual const Bool_t LoadPFVector() = 0; 
  // virtual const Bool_t LoadPFVector(const int startindex, const int Nsamples, const int tau, const int fs) = 0;
  
  virtual Bool_t LoadPFVector() { return true; }; 
  virtual Bool_t LoadPFVector(const int /*startindex*/, const int /*Nsamples*/, const int /*tau*/, const int /*fs*/) {  return true; };

  AliHLTCaloRawAnalyzerPeakFinderComponent(const AliHLTCaloRawAnalyzerPeakFinderComponent &, TString det ); 
  
  /*
  AliHLTCaloRawAnalyzerPeakFinderComponent & operator = (const AliHLTCaloRawAnalyzerPeakFinderComponent)
    {
      return *this;
    };
  */

private:
    AliHLTCaloRawAnalyzerPeakFinderComponent();
};



#endif
