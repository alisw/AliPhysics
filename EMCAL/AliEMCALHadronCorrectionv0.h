#ifndef ALIEMCALHADRONCORRECTIONV0_H
#define ALIEMCALHADRONCORRECTIONV0_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//                  
//*-- Author: Aleksei Pavlinov (WSU)
#include "AliEMCALHadronCorrection.h"

class AliEMCALHadronCorrectionv0: public AliEMCALHadronCorrection {

  private:
  static AliEMCALHadronCorrectionv0* fHadrCorr;

  protected:
  AliEMCALHadronCorrectionv0(const char *name="HadronCorrectionv0", const char *title="title");

  public:
  static  AliEMCALHadronCorrectionv0* Instance();
  virtual Double_t GetEnergy(Double_t pmom, Double_t eta, Int_t gid); 
  Double_t GetEnergy(Double_t pmom, Double_t eta) 
  {return GetEnergy(pmom,eta,7);}

  virtual ~AliEMCALHadronCorrectionv0() {}

  ClassDef(AliEMCALHadronCorrectionv0,1) // Hadron correction for EMC (version for MDC)
};

#endif // ALIEMCALHADRONCORRECTIONV0_H
