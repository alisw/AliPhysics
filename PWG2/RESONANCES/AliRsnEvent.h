//
// *** Class AliRsnEvent ***
//
// A container for a collection of AliRsnDaughter objects from an event.
// Contains also the primary vertex, useful for some cuts.
// In order to retrieve easily the tracks which have been identified
// as a specific type and charge, there is an array of indexes which
// allows to avoid to loop on all tracks and have only the neede ones.
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#ifndef ALIRSNEVENT_H
#define ALIRSNEVENT_H

#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliRsnDaughter.h"

class AliVEvent;
class AliRsnCutPID;
class AliESDtrackCuts;

class AliRsnEvent : public TObject
{
  public:

    AliRsnEvent(AliVEvent *ref = 0, AliVEvent *refMC = 0);
    AliRsnEvent(const AliRsnEvent& copy);
    AliRsnEvent& operator= (const AliRsnEvent& copy);
    virtual ~AliRsnEvent() { /*nothing*/ }
    
    // basic setters/getters
    void       SetRef  (AliVEvent *ref)     {fRef = ref;}
    void       SetRefMC(AliVEvent *refmc)   {fRefMC = refmc;}
    void       SetLeadingIndex(Int_t i)     {fLeading = i;}
    AliVEvent* GetRef()                     {return fRef;}
    AliVEvent* GetRefMC()                   {return fRefMC;}
    Int_t      GetLeadingIndex() const      {return fLeading;}
    Int_t      GetLeadingParticleID() const {return fLeading;}
    
    // getters which convert into allowed input types
    AliESDEvent* GetRefESD()   {return dynamic_cast<AliESDEvent*>(fRef);}
    AliAODEvent* GetRefAOD()   {return dynamic_cast<AliAODEvent*>(fRef);}
    AliMCEvent*  GetRefMCESD() {return dynamic_cast<AliMCEvent*> (fRefMC);}
    AliAODEvent* GetRefMCAOD() {return dynamic_cast<AliAODEvent*>(fRefMC);}
    Bool_t       IsESD()       {return (GetRefESD() != 0x0);}
    Bool_t       IsAOD()       {return (GetRefAOD() != 0x0);}
    
    // advanced getters
    Double_t         GetVz();
    Int_t            GetMultiplicity(AliESDtrackCuts *cuts = 0x0);
    
    // setters for a daughter
    Bool_t           SetDaughter  (AliRsnDaughter &daughter, Int_t index, AliRsnDaughter::ERefType type = AliRsnDaughter::kTrack);
    Bool_t           SetDaughterMC(AliRsnDaughter &daughter, Int_t index);
    AliRsnDaughter   GetDaughter  (Int_t i, AliRsnDaughter::ERefType type = AliRsnDaughter::kTrack);
    AliRsnDaughter   GetDaughterMC(Int_t i);
    Int_t            GetAbsoluteSum();
    Bool_t           ConvertAbsoluteIndex(Int_t index, Int_t &realIndex, AliRsnDaughter::ERefType &type);
    
    // leading particle stuff
    void             SetLeadingParticle(AliRsnDaughter &leading) {if (fLeading >= 0) SetDaughter(leading, fLeading);}
    Int_t            SelectLeadingParticle(Double_t ptMin = 0.0, AliRsnCutPID *cutPID = 0x0);
    Double_t         GetAverageMomentum(Int_t &count, AliRsnCutPID *cutPID = 0x0);
    Bool_t           GetAngleDistr(Double_t &angleMean, Double_t &angleRMS, AliRsnDaughter reference);

  private:
  
    Bool_t SetDaughterESDtrack  (AliRsnDaughter &target, Int_t index);
    Bool_t SetDaughterAODtrack  (AliRsnDaughter &target, Int_t index);
    Bool_t SetDaughterESDv0     (AliRsnDaughter &target, Int_t index);
    Bool_t SetDaughterAODv0     (AliRsnDaughter &target, Int_t index);
    Bool_t SetDaughterESDcascade(AliRsnDaughter &target, Int_t index);
    Bool_t SetDaughterAODcascade(AliRsnDaughter &target, Int_t index);
    Bool_t SetMCInfoESD         (AliRsnDaughter &target);
    Bool_t SetMCInfoAOD         (AliRsnDaughter &target);

    AliVEvent       *fRef;      // pointer to input event
    AliVEvent       *fRefMC;    // pointer to reference MC event (if any)
    Int_t            fLeading;  // index of leading track

    ClassDef(AliRsnEvent, 4);
};

#endif
