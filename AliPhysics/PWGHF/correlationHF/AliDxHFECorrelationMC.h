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
  virtual double GetD0Eff(AliVParticle* tr, Double_t evMult);

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
  enum{
    kOriginC = 0,
    kOriginB = 1,
    kOriginNonHF = 2,
    kOriginHadron = 3
  };


  enum{
    kReducedMode=0,
    kFullMode=1,
    kReducedModeFullMCInfo=2
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
  Int_t fRunMode;    // Which mode to run in (bigger thnsparse)
  Short_t fSystem;               // Which system pp/PbPb
  Bool_t fUseReducedOrigin; // Flag to store full or reduced origin info. Default false (full info)
  int fReducedOriginEl;    // Reduced origin info for electrons
  int fReducedOriginD0;   // Reduced origin info for D0s
  Bool_t fStorePoolbin; // Flag to store poolbin information (only valid for pPb reduced mode for now)
  ClassDef(AliDxHFECorrelationMC, 4)
};
#endif
