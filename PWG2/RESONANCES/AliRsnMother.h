#ifndef ALIRSNMOTHER_H
#define ALIRSNMOTHER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Interface to candidate resonance decaying into 2 bodies.
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRsnDaughter.h"

class AliRsnMother : public TObject {
public:

   AliRsnMother();
   AliRsnMother(const AliRsnMother &obj);
   AliRsnMother& operator=(const AliRsnMother &obj);
   virtual ~AliRsnMother();

   AliRsnDaughter*   GetDaughter(const Int_t &index)    const {if (index < 1) return   fDaughter[0] ; return   fDaughter[1] ;}
   AliRsnDaughter&   GetDaughterRef(const Int_t &index) const {if (index < 1) return (*fDaughter[0]); return (*fDaughter[1]);}
   TLorentzVector&   Sum()                                    {return fSum;}
   TLorentzVector&   SumMC()                                  {return fSumMC;}
   
   Double_t          AngleTo(AliRsnDaughter track, Bool_t mc = kFALSE) const {return fSum.Angle(track.P(mc).Vect());}
   Double_t          CosThetaStar(Bool_t first = kTRUE, Bool_t useMC = kFALSE);

   Bool_t            IsLabelEqual() const {return abs(fDaughter[0]->GetLabel()) == abs(fDaughter[1]->GetLabel());}
   Bool_t            IsIndexEqual() const {return (fDaughter[0]->GetID() == fDaughter[1]->GetID());}
   Int_t             CommonMother(Int_t &m0, Int_t &m1) const;
   Int_t             CommonMother() const {Int_t d0, d1; return CommonMother(d0, d1);}

   void              SetDaughters(AliRsnDaughter * const daughter1, Double_t mass1, AliRsnDaughter * const daughter2, Double_t mass2);
   void              ResetPair();
   void              PrintInfo(const Option_t *option = "ALL") const;
   Bool_t            CheckPair(Bool_t checkMC = kFALSE) const;

private:

   AliRsnDaughter  *fDaughter[2];      // elements of the pair
   TLorentzVector   fSum;              // sum computed from the two daughters
   TLorentzVector   fSumMC;            // sum computed from the two daughters

   ClassDef(AliRsnMother, 1)
};

#endif
