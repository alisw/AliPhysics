#ifndef ALIGENMUONUNCORR_H
#define ALIGENMUONUNCORR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


// Uncorrelated muon pairs generator
// author: alessandro.de.falco@ca.infn.it

#include "AliGenerator.h"
#include "TF1.h"
class AliGenMuonUncorr : public AliGenerator
{
 public:

  AliGenMuonUncorr();
  virtual ~AliGenMuonUncorr() {}
  virtual void Generate();
  virtual void Init();
  virtual void SetPart(Int_t part) {fIpart=part;}
  virtual void SetPtRange(Float_t ptmin, Float_t ptmax) {fPtMin = ptmin; fPtMax=ptmax;}
  virtual void SetParticleType(Int_t part) {SetPart(part);}
  virtual void SetSeed(UInt_t /*seed*/) {;}
  virtual void SetPtDistributionParameter(Int_t ipar, Double_t par) { fPt->SetParameter(ipar,par);}
  virtual TF1* GetPtDistribution() { return fPt;}
  virtual TF1* GetYDistribution() { return fY;}
  virtual void SetYDistributionParameter(Int_t ipar, Double_t par) { fY->SetParameter(ipar,par);}
protected:

  Int_t fIpart; // Particle type
  TF1 *fY;
  TF1 *fPt; 
  ClassDef(AliGenMuonUncorr,1) // uncorrelated muon pairs random generator
};

#endif
