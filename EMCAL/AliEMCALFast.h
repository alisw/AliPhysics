#ifndef ALIEMCALFAST_H
#define ALIEMCALFAST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */
#include <TObject.h>
//*-- Author: Andreas Morsch (CERN)

class AliEMCALFast : public TObject {
 public:
    virtual ~AliEMCALFast(){;}
    static Float_t SmearMomentum  (Int_t ind,   Float_t p);
    static Float_t Efficiency     (Int_t ind,   Float_t p);
    static Bool_t  EmcalAcceptance(Float_t eta, Float_t phi);
    static Bool_t  RandomReject   (Float_t eff);    
 protected:
  ClassDef(AliEMCALFast,1) // Jet for EMCAL
};
      

#endif // ALIEMCALJet_H
