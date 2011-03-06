#ifndef ALIJETHADRONCORRECTIONV1_H
#define ALIJETHADRONCORRECTIONV1_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//                  
//*-- Author: Mark Horner (LBL/UCT)
//
#include "AliJetHadronCorrection.h"

#define HCPARAMETERS    8 
#define HCPARAMETERSETS 3 

class AliJetDummyGeo;

class AliJetHadronCorrectionv1 : public AliJetHadronCorrection 
{
 public:
  static AliJetHadronCorrectionv1* Instance();
  virtual ~AliJetHadronCorrectionv1() {}

  virtual Double_t GetEnergy(Double_t pmom, Double_t eta, Int_t gid); 
  Double_t GetEnergy(Double_t pmom, Double_t eta){return GetEnergy(pmom,eta,7);}
  
  void SetGeometry(TString name, Double_t fs = 1.); 
  void SetGeometry2(const AliJetDummyGeo *geometry);
  void TrackPositionEMCal(const AliAODTrack* track,Double_t &eta, Double_t &phi);

 protected:
  AliJetHadronCorrectionv1():fSamplingFraction(0) {for (Int_t i = 0; i < 8; i++) fPar[i] = 0.;}
  AliJetHadronCorrectionv1(const char *name, const char *title);

 private:
  void SetParameters(TString name = "") {Warning("SetParameter","Dummy method with argument %s",name.Data());}

  static AliJetHadronCorrectionv1* fgHadrCorr;  // Pointer to global instance (singleton)
  static Double_t fgParLookup[HCPARAMETERS][HCPARAMETERSETS]; // Global array with parameters for hadronic response
  Double_t fPar[8];            // Parameters
  Float_t  fSamplingFraction;  // Sampling fraction
    
  ClassDef(AliJetHadronCorrectionv1,2) // Hadron correction for EMC (version for MDC)
};
	
#endif // ALIJETHADRONCORRECTIONV1_H
	
