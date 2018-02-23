#ifndef AliAODDIJet_H
#define AliAODDIJet_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD di-jet class
//     The present version is for test purposes only
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TLorentzVector.h>
#include <TRefArray.h>
#include <TRef.h>
#include "AliAODJet.h"


class AliAODDiJet : public AliAODJet {

 public:
    AliAODDiJet();
    AliAODDiJet(Double_t px, Double_t py, Double_t pz, Double_t e);
    AliAODDiJet(TLorentzVector & p);
    virtual ~AliAODDiJet();

    void SetJetRefs(AliAODJet* jet1, AliAODJet* jet2);
    AliAODJet* Jet(Int_t i) {return ((AliAODJet*) (fJetR->At(i)));}
    Float_t    DeltaPhi();
    Float_t    PhiImbalance();

 private:
    AliAODDiJet(const AliAODDiJet& jet);
    AliAODDiJet& operator=(const AliAODDiJet& jet);

 private:
    TRefArray*  fJetR;  // References to jets
    TRef        fJet1;  // Reference to Jet 1
    TRef        fJet2;  // Reference to Jet 2
    ClassDef(AliAODDiJet, 1);
};
#endif
