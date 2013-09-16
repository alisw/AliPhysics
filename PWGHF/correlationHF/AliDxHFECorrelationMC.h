//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliDxHFECorrelationMC.h
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-11-02
/// @brief  Worker class for DxHFE correlation on MC
///

#ifndef ALIDXHFECORRELATIONMC_H
#define ALIDXHFECORRELATIONMC_H

#include "AliDxHFECorrelation.h"

class AliHFCorrelator;
class AliVEvent;
class AliVParticle;
class TObjArray;

class AliDxHFECorrelationMC : public AliDxHFECorrelation {
 public:
  /// default constructor
  AliDxHFECorrelationMC(const char* name=NULL);
  /// destructor
  virtual ~AliDxHFECorrelationMC();

  /// fill histograms from particles
  virtual int Fill(const TObjArray* candidatesD0, TObjArray* candidatesElectron, const AliVEvent* pEvent);

  /// histogram event properties
  virtual THnSparse* DefineTHnSparse();
  virtual int FillParticleProperties(AliVParticle* tr, AliVParticle *as, Double_t* data, int dimension) const;
  // Tests the particle
  virtual Bool_t TestParticle(AliVParticle* p, Int_t id);
  // Get the D0 efficiency
  virtual double GetD0Eff(AliVParticle* tr);
  virtual void SetD0EffMap(TH1* eff, int whichMap=AliDxHFECorrelation::kPrompt);
  //    if(whichMap==AliDxHFECorrelation::kPrompt) { fD0EffMapP=eff; }
  //  if(whichMap==AliDxHFECorrelation::kFeedDown) { fD0EffMapFD=eff;}
  //};

  // parse argument string
  virtual int ParseArguments(const char* arguments);

  virtual void SetEventType(int type){fMCEventType=type;}
  //virtual void SetStoreOrigin(int D,int el){fStoreOriginD=D; fStoreOriginEl=el;}

  enum{
    kAll = 0,
    kC   = 1,
    kB   = 2,
    kHF  = 3,
    kNonHF= 4,
    kHadrons=5
  };

 protected:

 private:
  /// copy constructor
  AliDxHFECorrelationMC(const AliDxHFECorrelationMC& other);
  /// assignment operator
  AliDxHFECorrelationMC& operator=(const AliDxHFECorrelationMC& other);

  int fMCEventType;  // Holds MC Event type, retrieved from MCHeader
  int fStoreOriginEl;// Which origin to store for electrons
  int fStoreOriginD; // Which origin to store for Ds
  TH1* fD0EffMapP;   // Eff map for Prompt D0 
  TH1* fD0EffMapFD;  // Eff map for Feeddown D0 

  ClassDef(AliDxHFECorrelationMC, 2)
};
#endif
