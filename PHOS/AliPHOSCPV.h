#ifndef ALIPHOSCPV_H
#define ALIPHOSCPV_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:CPV      //
//                                            //
//  Author: Yuri Kharlov, IHEP, Protvino      //
//  e-mail: Yuri.Kharlov@cern.ch              //
//  Last modified: 23 March 2000              //
////////////////////////////////////////////////
 
// --- ROOT system ---
#include <TClonesArray.h> 
#include <TLorentzVector.h>

// --- galice header files ---
#include "AliHit.h"

//==============================================================================
//                              CPVModule
//==============================================================================

class CPVModule : public TObject {

public:

  virtual ~CPVModule(void);
           CPVModule(void);
  
  void     Clear(Option_t *opt="");
  void     Print(Option_t *opt="");
  void     AddHit(TLorentzVector p, Float_t *xy, Int_t ipart);
  void     MakeBranch(Int_t i);
  void     SetTreeAddress(Int_t i);
  
  TClonesArray *Hits    (void) {return fHits;}

private:
  
  TClonesArray *fHits;              // List of hits in the Module

  ClassDef(CPVModule,1)             // CPV Module
};

//==============================================================================
//                              CPVHit
//==============================================================================

class CPVHit : public TObject {
  
private:
  TLorentzVector fMomentum;   // 4-momentum of the particle
  Float_t        fXhit;       // Hit's X-coordinates
  Float_t        fYhit;       // Hit's Y-coordinates
  Int_t          fIpart;      // Hit's particle type
  
public:
  virtual ~CPVHit() {}
           CPVHit() {}
           CPVHit(TLorentzVector p, Float_t *xy, Int_t ipart);
  
  TLorentzVector GetMomentum()  { return  fMomentum; }
  Float_t        GetX()         { return  fXhit;     }
  Float_t        GetY()         { return  fYhit;     }
  Int_t          GetIpart()     { return  fIpart;    }
  void           Print();

  ClassDef(CPVHit,1)  // Hit object in one CPV module
};
 
//==============================================================================
//                              CPVDigit
//==============================================================================

class CPVDigit : public TObject {
  
private:
  Int_t    fXpad;       // Digit's pad number in Phi
  Int_t    fYpad;       // Digit's pad number in Z
  Float_t  fQpad;       // Digit's pad amplitude
  
public:
  virtual ~CPVDigit() {}
           CPVDigit() {}
           CPVDigit(Int_t x, Int_t y, Float_t q);
  
  void     SetQpad(Float_t q) { fQpad = q;     }
  Int_t    GetXpad()          { return  fXpad; }
  Int_t    GetYpad()          { return  fYpad; }
  Float_t  GetQpad()          { return  fQpad; }

  ClassDef(CPVDigit,1)  // Digit object in one CPV pad
};
 
#endif
