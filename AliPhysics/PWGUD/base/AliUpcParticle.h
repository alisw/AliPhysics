#ifndef AliUpcParticle_h
#define AliUpcParticle_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class to store basic track parameters in the upc tree
// evgeny.kryshen@cern.ch


#include "AliVParticle.h"
#include "AliLog.h"
#include "TArrayI.h"
#include "TArrayF.h"

class AliUpcParticle : public AliVParticle {
 public:
  AliUpcParticle():AliVParticle(),fArrayI(),fArrayF(){}
  AliUpcParticle(Int_t ni, Int_t nf):AliVParticle(),fArrayI(ni),fArrayF(nf){}

  void SetI(Int_t i, Int_t var  ) { fArrayI.SetAt(var,i); }
  void SetF(Int_t i, Float_t var) { fArrayF.SetAt(var,i); }
  Int_t   GetI(Int_t i) { return fArrayI.GetAt(i); }
  Float_t GetF(Int_t i) { return fArrayF.GetAt(i); }

  virtual ~AliUpcParticle(){}

  virtual Double_t Pt()    const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Phi()   const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Eta()   const { AliFatal("Not implemented"); return 0; }
  virtual Short_t Charge() const { AliFatal("Not implemented"); return 0; }
  virtual UInt_t Mask()    const { AliFatal("Not implemented"); return 0; }
  
  virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t P()  const { AliFatal("Not implemented"); return 0; }
  virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

  virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
  virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

  virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }
  virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
  virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }
  virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
  
  virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
  virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }
 protected:
  TArrayI fArrayI;
  TArrayF fArrayF;
  ClassDef(AliUpcParticle,2);
};

#endif

