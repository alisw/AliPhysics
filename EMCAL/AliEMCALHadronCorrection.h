#ifndef ALIEMCALHADRONCORRECTION_H
#define ALIEMCALHADRONCORRECTION_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//                  
//*-- Author: Aleksei Pavlinov (WSU)
// This pure abstract class which defines only interface
// 
#include "TNamed.h"

class AliEMCALHadronCorrection : public TNamed {

  public:
  AliEMCALHadronCorrection(const char *name="name", const char *title="title");
  virtual ~AliEMCALHadronCorrection() {} 

  // Add for particle
  virtual Double_t GetEnergy(Double_t pmom, Double_t eta, Int_t gid)=0; 

  ClassDef(AliEMCALHadronCorrection,1) // Hadron correction for EMC (abstract class)
};

#endif // ALIEMCALHADRONCORRECTION_H
