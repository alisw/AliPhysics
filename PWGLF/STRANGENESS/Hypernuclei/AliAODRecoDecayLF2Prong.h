#ifndef ALIAODRECODECAYLF2PRONG_H
#define ALIAODRECODECAYLF2PRONG_H
/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliAODRecoDecayLF2Prong
// base class for AOD reconstructed 2-prong decays
// (3LH->3Hepi, ...)
//  strongly based on AliAODRecoDecayHF2Prong 
//  by A.Dainese and G.E.Bruno
//
//  Author: Ramona Lea (ramona.lea@cern.ch) 
// 
//***********************************************************

#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayLF.h"

class AliAODRecoDecayLF2Prong : public AliAODRecoDecayLF {

 public:

  AliAODRecoDecayLF2Prong();
  AliAODRecoDecayLF2Prong(AliAODVertex *vtx2,
			  Double_t *px,Double_t *py,Double_t *pz,
			  Double_t *d0,Double_t *d0err,Float_t dca);
  AliAODRecoDecayLF2Prong(AliAODVertex *vtx2,
			  Double_t *d0,Double_t *d0err,Float_t dca);    
  AliAODRecoDecayLF2Prong(const AliAODRecoDecayLF2Prong& source);
  AliAODRecoDecayLF2Prong& operator=(const AliAODRecoDecayLF2Prong& source); 

  virtual ~AliAODRecoDecayLF2Prong() {}  
 
  Double_t Prodd0d0() const {return AliAODRecoDecay::Prodd0d0(0,1);} 

  

 private:

  ClassDef(AliAODRecoDecayLF2Prong,1)  // base class for AOD reconstructed 
                                       // heavy-flavour 2-prong decays
};

#endif
