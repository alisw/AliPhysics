#ifndef ALIRSNEVENT_H
#define ALIRSNEVENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Interface to full event.
//
////////////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>

#include "AliLog.h"

#include "AliStack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

#include "AliRsnDaughter.h"

class AliRsnCutSet;
class AliPIDResponse;

class AliRsnEvent : public TObject {
public:

   AliRsnEvent(AliVEvent *ref = 0, AliVEvent *refMC = 0);
   AliRsnEvent(const AliRsnEvent& copy);
   AliRsnEvent& operator= (const AliRsnEvent& copy);
   virtual ~AliRsnEvent();

   // basic setters/getters
   void             SetRef(AliVEvent *ref)              {fRef = ref;}
   void             SetRefMC(AliVEvent *refmc);
   void             SetPIDResponse(AliPIDResponse *pid) {fPID = pid;}
   AliVEvent*       GetRef()                            {return fRef;}
   AliVEvent*       GetRefMC()                          {return fRefMC;}
   Int_t            GetLeadingIndex() const             {return fLeading;}
   AliPIDResponse*  GetPIDResponse()                    {return fPID;}

   // getters which convert into allowed input types
   Bool_t           Match(AliVEvent *ev, TClass *ref) {if (ev) return (ev->InheritsFrom(ref)); return kFALSE;}
   Bool_t           IsESD()                           {return (Match(fRef, AliESDEvent::Class()));}
   Bool_t           IsAOD()                           {return (Match(fRef, AliAODEvent::Class()));}
   Bool_t           InputOK();                        
   AliESDEvent*     GetRefESD()                       {if (IsESD()) return (AliESDEvent*)fRef;   return 0x0;}
   AliMCEvent*      GetRefMCESD()                     {if (IsESD()) return (AliMCEvent *)fRefMC; return 0x0;}
   AliAODEvent*     GetRefAOD()                       {if (IsAOD()) return (AliAODEvent*)fRef;   return 0x0;}
   AliAODEvent*     GetRefMCAOD()                     {if (IsAOD()) return (AliAODEvent*)fRefMC; return 0x0;}
   TClonesArray*    GetAODList()                      {return fAODList;}
   Bool_t           HasMC()                           {if (IsESD()) return (fRefMC != 0x0); else if (IsAOD()) return (fRefMC != 0x0 && fAODList != 0x0); return kFALSE;}

   // setters for a daughter
   void             SetDaughter          (AliRsnDaughter &daughter, Int_t index, Bool_t fromMC = kFALSE);
   AliRsnDaughter   GetDaughter          (Int_t index, Bool_t fromMC);
   void             SetDaughterESDtrack  (AliRsnDaughter &target, Int_t index);
   void             SetDaughterESDv0     (AliRsnDaughter &target, Int_t index);
   void             SetDaughterESDcascade(AliRsnDaughter &target, Int_t index);
   void             SetDaughterESDMCtrack(AliRsnDaughter &target, Int_t index);
   void             SetDaughterAODtrack  (AliRsnDaughter &target, Int_t index);
   void             SetDaughterAODv0     (AliRsnDaughter &target, Int_t index);
   void             SetDaughterAODcascade(AliRsnDaughter &target, Int_t index);
   void             SetDaughterAODMCtrack(AliRsnDaughter &target, Int_t index);
   Bool_t           SetMCInfoESD         (AliRsnDaughter &target);
   Bool_t           SetMCInfoAOD         (AliRsnDaughter &target);
   
   // counters/converters of candidates
   Int_t            GetAbsoluteSum() {if (fRef) return (fRef->GetNumberOfTracks() + fRef->GetNumberOfV0s() + fRef->GetNumberOfCascades()); return 0;}
   Bool_t           ConvertAbsoluteIndex(Int_t index, Int_t &realIndex, AliRsnDaughter::ERefType &type);
   Int_t            ConvertRealIndex(Int_t index, AliRsnDaughter::ERefType type);

   // leading particle stuff
   void             SetLeadingParticle(AliRsnDaughter &leading) {if (fLeading >= 0) SetDaughter(leading, fLeading, kFALSE);}
   Int_t            SelectLeadingParticle(AliRsnCutSet *cuts = 0x0);
   
private:

   AliVEvent      *fRef;            //  pointer to input event
   AliVEvent      *fRefMC;          //  pointer to reference MC event (if any)
   Int_t           fLeading;        //  index of leading track
   AliPIDResponse *fPID;            //  pointer to PID response
   TClonesArray   *fAODList;        //  pointer to AOD list of particles (if any)

   ClassDef(AliRsnEvent, 6);
};

inline Bool_t AliRsnEvent::InputOK()
{
//
// Check that input is ESD or AOD
//

   if (IsESD()) {
      AliDebugClass(1, "Input is ESD");
      return kTRUE;
   } else if (IsAOD()) {
      AliDebugClass(1, "Input is AOD");
      return kTRUE;
   } else {
      AliError("Need to process ESD or AOD input");
      return kFALSE;
   }
}

inline void AliRsnEvent::SetRefMC(AliVEvent *mc) 
{
//
// Assign pointer to MC event.
// If it is an AOD, retrieve the list of MC particles
//

   if (!mc) {
      fRefMC = 0x0;
      fAODList = 0x0;
   }

   fRefMC = mc;
   fAODList = 0x0;
   if (fRefMC->InheritsFrom(AliAODEvent::Class())) {
      AliAODEvent *aod = (AliAODEvent*)mc;
      fAODList = (TClonesArray*)(aod->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
      if (!fAODList) fRefMC = 0x0;
   }
}

#endif
