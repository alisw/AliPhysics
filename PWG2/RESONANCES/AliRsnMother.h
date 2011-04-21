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

   AliRsnMother() : fRefEvent(0), fSum(), fRef() {fDaughter[0] = fDaughter[1] = 0;}
   AliRsnMother(const AliRsnMother &obj);
   AliRsnMother& operator=(const AliRsnMother &obj);
   virtual ~AliRsnMother();

   // setters (4-vectors cannot be set)
   void  Reset();
   void  SetDaughter(Int_t i, AliRsnDaughter *d) {fDaughter[CkID(i)] = d;}
   void  SetRefEvent(AliRsnEvent *event)         {fRefEvent = event;}
   
   // getters
   AliRsnEvent*      GetRefEvent()               {return fRefEvent;}
   AliRsnDaughter*   GetDaughter(const Int_t &i) {return fDaughter[CkID(i)];}
   TLorentzVector&   Sum(Bool_t mc)              {return (mc ? fSumMC : fSum);}
   TLorentzVector&   Ref(Bool_t mc)              {return (mc ? fRefMC : fRef);}
   Bool_t            GetResolution(Double_t &value);
   
   // checks
   Bool_t    IsLabelEqual()  const {return TMath::Abs(fDaughter[0]->GetLabel()) == TMath::Abs(fDaughter[1]->GetLabel());}
   Bool_t    IsIndexEqual()  const {return (fDaughter[0]->GetID() == fDaughter[1]->GetID());}
   Bool_t    IsOwnerEqual()  const {return (fDaughter[0]->GetOwnerEvent() == fDaughter[1]->GetOwnerEvent());}
   Int_t     CommonMother()  const;
   
   // angles
   Double_t  AngleTo(AliRsnDaughter *track, Bool_t mc = kFALSE) {return track->P(mc).Angle(Sum(mc).Vect());}
   Double_t  AngleToLeading(Bool_t &success);
   
   // computations
   void      ComputeSum(Double_t mass1, Double_t mass2, Double_t motherMass);
   Double_t  CosThetaStar(Bool_t first = kTRUE, Bool_t useMC = kFALSE);   
   void      PrintInfo(const Option_t *option = "ALL") const;
   Bool_t    CheckPair(Bool_t checkMC = kFALSE) const;

private:

   Int_t CkID(Int_t i) {if (i < 1) return 0; else return 1;}

   AliRsnDaughter  *fDaughter[2];      // elements of the pair
   AliRsnEvent     *fRefEvent;         // reference event
   TLorentzVector   fSum;              // sum computed from the two daughters (rec)
   TLorentzVector   fSumMC;            // sum computed from the two daughters (sim)
   TLorentzVector   fRef;              // same to sum, but with fixed mass hypothesis (rec)
   TLorentzVector   fRefMC;            // same to sum, but with fixed mass hypothesis (sim)

   ClassDef(AliRsnMother, 1)
};

#endif
