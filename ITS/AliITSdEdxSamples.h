#ifndef ALIITSDEDXSAMPLES_H
#define ALIITSDEDXSAMPLES_H
/* Copyright(c) 2009-2012, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to store information for PID with ITS                   //
// and truncated mean computation methods                        //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TObject.h>
#include "AliPID.h"

class AliITSPidParams;

class AliITSdEdxSamples : public TObject {

 public:
  AliITSdEdxSamples();
  AliITSdEdxSamples(Int_t nSamples, Double_t* samples, Double_t mom, Int_t specie=0);
  virtual ~AliITSdEdxSamples(){};

  void SetSamples(Int_t nSamples, Double_t* samples);
  void SetSamplesAndMomenta(Int_t nSamples, Double_t* samples, Double_t* mom);
  void SetMomentum(Double_t mom){
    fP=mom;
  }
  void SetParticleSpecieMC(Int_t specie){
    fParticleSpecie=specie;
  }

  Int_t GetNumberOfSamples() const {
    return fNSamples;
  }
  Double_t GetdEdxSample(Int_t i) const {
    if(i<fNSamples) return fdEdxSamples[i];
    else return 0.;
  }
  Double_t GetMomentum() const {
    return fP;
  }
  Double_t GetMomentumAtSample(Int_t i) const{
    if(i<fNSamples) return fPAtSample[i];
    else return 0.;
  } 
  Int_t GetParticleSpecieMC() const {
    return fParticleSpecie;
  }

  Double_t GetTruncatedMean(Double_t frac=0.5, Double_t mindedx=0.) const;
  Double_t GetWeightedMean(Double_t mindedx=0.) const;
  void     GetConditionalProbabilities(AliITSPidParams* pars, Double_t condprob[AliPID::kSPECIES], Double_t mindedx=0.) const;

 protected:
  static const UShort_t fgkMaxSamples=10;  // max. n. of layers with dE/dx info
  
  Int_t    fNSamples;                   // number of samples
  Double_t fdEdxSamples[fgkMaxSamples]; // dE/dx samples
  Double_t fP;                          // track momentum
  Int_t    fParticleSpecie;             // MC generated particle
  Double_t fPAtSample[fgkMaxSamples];   // track momentum at specific samples
  
  ClassDef(AliITSdEdxSamples,1);

};
#endif 
