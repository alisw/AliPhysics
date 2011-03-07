#ifndef ALIRSNMOTHER_H
#define ALIRSNMOTHER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Interface to candidate resonance decaying into 2 bodies.
//
////////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include "AliRsnDaughter.h"

class AliRsnEvent;

class AliRsnMother : public TObject {
public:

   AliRsnMother() : fRefEvent(0), fSum(), fSumMC() {fDaughter[0] = fDaughter[1] = 0;}
   AliRsnMother(const AliRsnMother &obj);
   AliRsnMother& operator=(const AliRsnMother &obj);
   virtual ~AliRsnMother();

   // setters (4-vectors cannot be set)
   void  SetDaughter(Int_t i, AliRsnDaughter *d) {fDaughter[CkID(i)] = d;}
   void  SetRefEvent(AliRsnEvent *event)         {fRefEvent = event;}
   
   // getters
   AliRsnEvent*      GetRefEvent()               {return fRefEvent;}
   AliRsnDaughter*   GetDaughter(const Int_t &i) {return fDaughter[CkID(i)];}
   TLorentzVector&   Sum()                       {return fSum;}
   TLorentzVector&   SumMC()                     {return fSumMC;}
   
   // checks
   Bool_t    IsLabelEqual()  const {return TMath::Abs(fDaughter[0]->GetLabel()) == TMath::Abs(fDaughter[1]->GetLabel());}
   Bool_t    IsIndexEqual()  const {return (fDaughter[0]->GetID() == fDaughter[1]->GetID());}
   Bool_t    IsOwnerEqual()  const {return (fDaughter[0]->GetOwnerEvent() == fDaughter[1]->GetOwnerEvent());}
   Int_t     CommonMother()  const {Int_t d0, d1; return CommonMother(d0, d1);}
   Int_t     CommonMother(Int_t &m0, Int_t &m1) const;
   
   // useful computations/operations
   void      ComputeSum(Double_t mass1, Double_t mass2);
   Double_t  AngleTo(AliRsnDaughter track, Bool_t mc = kFALSE);
   Double_t  AngleToLeading(Bool_t &success);
   Double_t  CosThetaStar(Bool_t first = kTRUE, Bool_t useMC = kFALSE);
   void      ResetPair();
   void      PrintInfo(const Option_t *option = "ALL") const;
   Bool_t    CheckPair(Bool_t checkMC = kFALSE) const;

private:

   Int_t CkID(Int_t i) {if (i < 1) return 0; else return 1;}

   AliRsnDaughter  *fDaughter[2];      // elements of the pair
   AliRsnEvent     *fRefEvent;         // reference event
   TLorentzVector   fSum;              // sum computed from the two daughters
   TLorentzVector   fSumMC;            // sum computed from the two daughters

   ClassDef(AliRsnMother, 1)
};

#endif
