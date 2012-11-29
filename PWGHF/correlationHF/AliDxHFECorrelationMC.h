//-*- Mode: C++ -*-

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

  virtual void SetEventType(int type){fMCEventType=type;}

 protected:

 private:
  /// copy constructor
  AliDxHFECorrelationMC(const AliDxHFECorrelationMC& other);
  /// assignment operator
  AliDxHFECorrelationMC& operator=(const AliDxHFECorrelationMC& other);

  int  fMCEventType;  // Holds MC Event type, retrieved from MCHeader

  ClassDef(AliDxHFECorrelationMC, 1)
};
#endif
