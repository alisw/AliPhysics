#ifndef ALIPHOSCPVHIT_H
#define ALIPHOSCPVHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Hit class for CPV                         //
//                                            //
//  Author: Yuri Kharlov, IHEP, Protvino      //
//  e-mail: Yuri.Kharlov@cern.ch              //
//  Last modified: 23 March 2000              //
////////////////////////////////////////////////
 
// --- ROOT system ---
#include <TLorentzVector.h>

// --- galice header files ---
#include <AliHit.h>

//==============================================================================
//                              AliPHOSCPVHit
//==============================================================================

class AliPHOSCPVHit : public AliHit {
  
public:
  virtual ~AliPHOSCPVHit() {}
           AliPHOSCPVHit() {}
           AliPHOSCPVHit(Int_t shunt, Int_t track, TLorentzVector p, Float_t *xy, Int_t ipart);
  
  TLorentzVector GetMomentum()  { return  fMomentum; }
  Int_t          GetIpart()     { return  fIpart;    }
  void           Print();

private:
  TLorentzVector fMomentum;   // 4-momentum of the particle
  Int_t          fIpart;      // Hit's particle type
  
  ClassDef(AliPHOSCPVHit,1)  // Hit object in one CPV module
};
 
#endif
