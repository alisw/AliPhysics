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

//==============================================================================
//                              AliPHOSCPVHit
//==============================================================================

class AliPHOSCPVHit : public TObject {
  
public:
  virtual ~AliPHOSCPVHit() {}
           AliPHOSCPVHit() {}
           AliPHOSCPVHit(TLorentzVector p, Float_t *xy, Int_t ipart);
  
  TLorentzVector GetMomentum()  { return  fMomentum; }
  Float_t        GetX()         { return  fXhit;     }
  Float_t        GetY()         { return  fYhit;     }
  Int_t          GetIpart()     { return  fIpart;    }
  void           Print();

private:
  TLorentzVector fMomentum;   // 4-momentum of the particle
  Float_t        fXhit;       // Hit's X-coordinates
  Float_t        fYhit;       // Hit's Y-coordinates
  Int_t          fIpart;      // Hit's particle type
  
  ClassDef(AliPHOSCPVHit,1)  // Hit object in one CPV module
};
 
#endif
