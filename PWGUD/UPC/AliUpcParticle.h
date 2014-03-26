#ifndef AliUpcParticle_h
#define AliUpcParticle_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class to store basic track parameters in the upc tree
// evgeny.kryshen@cern.ch


#include "AliVParticle.h"
#include "AliLog.h"
#include "TArrayF.h"

class AliUpcParticle : public AliVParticle, public TArrayF {
 public:
  AliUpcParticle():AliVParticle(),TArrayF(),fPt(0),fEta(0),fPhi(0),fCharge(0),fMask(0)
  {
  }
  AliUpcParticle(Float_t pt, Float_t eta, Float_t phi, Short_t charge, UInt_t mask,UInt_t nData=0)
  :AliVParticle(),TArrayF(nData),fPt(pt),fEta(eta),fPhi(phi),fCharge(charge),fMask(mask)
  {
  }
  virtual ~AliUpcParticle(){}

  virtual Double_t Pt()    const { return fPt;      }
  virtual Double_t Phi()   const { return fPhi;     }
  virtual Double_t Eta()   const { return fEta;     }
  virtual Short_t Charge() const { return fCharge;  }
  virtual UInt_t Mask()    const { return fMask;    }
  
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

  virtual Bool_t IsEqual(const TObject* obj) const { return (obj->GetUniqueID() == GetUniqueID()); }

  virtual void SetPhi(Double_t phi) { fPhi = phi; }
 protected:
  Float_t fPt;
  Float_t fEta;
  Float_t fPhi;
  Short_t fCharge;
  UInt_t  fMask;     // Filter bit mask
  ClassDef(AliUpcParticle,1);
};

#endif

