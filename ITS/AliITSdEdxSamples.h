#ifndef ALIITSDEDXSAMPLES_H
#define ALIITSDEDXSAMPLES_H
/* Copyright(c) 2009-2012, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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
  AliITSdEdxSamples(Int_t nSamples, Double_t* esamples, Double_t* xsamples, Double_t mom, Int_t specie=0);
  AliITSdEdxSamples(const AliITSdEdxSamples& source);
  AliITSdEdxSamples& operator=(const AliITSdEdxSamples& source); 
  virtual ~AliITSdEdxSamples(){};

  void SetdESamples(Int_t nSamples, Double_t* samples);
  void SetdxSamples(Int_t nSamples, Double_t* samples);
  void SetSamplesAndMomenta(Int_t nSamples, Double_t* esamples, Double_t* xsamples, Double_t* mom);
  void SetNSamples(Int_t nSamples){
    fNSamples=nSamples;
  }
  void SetLayerSample(Int_t iLayer, Bool_t haspoint, Double_t dE=0., Double_t dx=0., Double_t p=0.);

  void SetMomentum(Double_t mom){
    fP=mom;
  }
  void SetParticleSpecieMC(Int_t specie){
    fParticleSpecie=specie;
  }
  void SetClusterMap(UInt_t map){
    fClusterMap=map;
  }
  void SetPointOnLayer(Int_t iLay){
    fClusterMap|=(1<<iLay);
  }
  void SetLayersForPID(UInt_t laypid){
    fLayersForPid=laypid;
  }
  void SetUseLayerForPid(Int_t iLay){
    fLayersForPid|=(1<<iLay);
  }
  void SetNotUseLayerForPid(Int_t iLay){
    if(fLayersForPid&(1<<iLay)) fLayersForPid-=(1<<iLay);
  }

  Int_t GetNumberOfSamples() const {
    return fNSamples;
  }
  Int_t GetNumberOfEffectiveSamples() const{
    Int_t nS=0;
    for (Int_t il=0; il<fNSamples; il++)  if(HasPointOnLayer(il)) nS++;
    return nS;
  }

  Double_t GetdESample(Int_t i) const {
    if(i<fNSamples) return fdESamples[i];
    else return 0.;
  }
  Double_t GetNehPairs(Int_t i) const {
    if(i<fNSamples) return fdESamples[i]*1000./3.63;
    else return 0.;
  }
  Double_t GetQfC(Int_t i) const{
    return GetNehPairs(i)*1.6E-4;
  }
  Double_t GetdxSample(Int_t i) const {
    if(i<fNSamples) return fdxSamples[i];
    else return 0.;
  }
  Double_t GetdEdxSample(Int_t i) const { // keV/100um
    if(i<fNSamples && fdxSamples[i]>0.) 
      return fdESamples[i]/(fdxSamples[i]*100.);
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
  UInt_t GetClusterMap() const{
    return fClusterMap;
  }
  Bool_t HasPointOnLayer(Int_t iLay) const{
    return fClusterMap&(1<<iLay);
  }
  UInt_t GetLayersForPid() const{
    return fLayersForPid;
  }
  Bool_t UseLayerForPid(Int_t iLay) const{
    return fLayersForPid&(1<<iLay);
  }

  void PrintAll() const;
  void PrintClusterMap() const;

  Double_t GetTruncatedMean(Double_t frac=0.5, Double_t mindedx=0.) const;
  Double_t GetWeightedMean(Double_t mindedx=0.) const;
  void     GetConditionalProbabilities(AliITSPidParams* pars, Double_t condprob[AliPID::kSPECIES], Double_t mindedx=0.) const;


 protected:

  void SetClusterMapFromdE(){
    fClusterMap=0;
    for(Int_t i=0; i<fNSamples; i++)
      if(fdESamples[i]>0.) SetPointOnLayer(i);
  }

  enum{kMaxSamples=10}; // max. n. of layers with dE/dx info

  Int_t    fNSamples;                // number of samples
  UInt_t   fClusterMap;              // map of clusters in layers
  Double_t fdESamples[kMaxSamples];  // dE samples (keV)
  Double_t fdxSamples[kMaxSamples];  // dx samples (cm)
  Double_t fP;                       // track momentum
  Int_t    fParticleSpecie;          // MC generated particle
  Double_t fPAtSample[kMaxSamples];  // track momentum at specific samples
  UInt_t   fLayersForPid;            // bit-map to enable/disable layers in PID

  ClassDef(AliITSdEdxSamples,3);

};
#endif 
