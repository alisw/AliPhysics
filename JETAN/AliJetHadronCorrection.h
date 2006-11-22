#ifndef ALIJETHADRONCORRECTION_H
#define ALIJETHADRONCORRECTION_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//                  
//*-- Author: Aleksei Pavlinov (WSU)
// This pure abstract class which defines only interface
// 
#include "TNamed.h"

class AliJetHadronCorrection : public TNamed {

  public:
  AliJetHadronCorrection(const char *name="name", const char *title="title");
  virtual ~AliJetHadronCorrection() {} 

  // Add for particle
  virtual Double_t GetEnergy(Double_t pmom, Double_t eta, Int_t gid)=0; 

  ClassDef(AliJetHadronCorrection,1) // Hadron correction for EMC (abstract class)
};

#endif // ALIJETHADRONCORRECTION_H
