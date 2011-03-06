#ifndef ALIJETHADRONCORRECTIONV0_H
#define ALIJETHADRONCORRECTIONV0_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//                  
//*-- Author: Aleksei Pavlinov (WSU)
#include "AliJetHadronCorrection.h"

class AliJetHadronCorrectionv0: public AliJetHadronCorrection {

  public:
  static  AliJetHadronCorrectionv0* Instance();
  virtual Double_t GetEnergy(Double_t pmom, Double_t eta, Int_t gid); 
  Double_t GetEnergy(Double_t pmom, Double_t eta){return GetEnergy(pmom,eta,7);}
  virtual ~AliJetHadronCorrectionv0() {}

  protected:
  AliJetHadronCorrectionv0() {;}
  AliJetHadronCorrectionv0(const char *name, const char *title);


  private:
    static AliJetHadronCorrectionv0* fHadrCorr; // pointer to instance

  ClassDef(AliJetHadronCorrectionv0,1) // Hadron correction for EMC (version for MDC)
};

#endif // ALIJETHADRONCORRECTIONV0_H
