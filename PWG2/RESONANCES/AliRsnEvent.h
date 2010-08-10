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

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliRsnDaughter.h"

class AliVEvent;
class AliMCEvent;
class AliRsnCutPID;

class AliRsnEvent : public TObject
{
  public:

    AliRsnEvent(AliVEvent *ref = 0, AliMCEvent *refMC = 0);
    AliRsnEvent(const AliRsnEvent& copy);
    AliRsnEvent& operator= (const AliRsnEvent& copy);
    virtual ~AliRsnEvent();

    void             SetRef(AliVEvent * const event, AliMCEvent * const mc = 0) {fRef = event; fRefMC = mc;}
    AliVEvent*       GetRef() const {return fRef;}
    AliMCEvent*      GetRefMC() const {return fRefMC;}
    AliESDEvent*     GetRefESD() const {return dynamic_cast<AliESDEvent*>(fRef);}
    AliAODEvent*     GetRefAOD() const {return dynamic_cast<AliAODEvent*>(fRef);}
    Bool_t           IsESD() const {return (GetRefESD() != 0x0);}
    Bool_t           IsAOD() const {return (GetRefAOD() != 0x0);}
    
    Double_t         GetVz();
    Int_t            GetMultiplicity();
    
    Bool_t           SetDaughter(AliRsnDaughter &daughter, Int_t index, AliRsnDaughter::ERefType type = AliRsnDaughter::kTrack);
    Bool_t           SetDaughterMC(AliRsnDaughter &daughter, Int_t index);
    AliRsnDaughter   GetDaughter(Int_t i, AliRsnDaughter::ERefType type = AliRsnDaughter::kTrack);
    AliRsnDaughter   GetDaughterMC(Int_t i);
    
    Int_t            SelectLeadingParticle(Double_t ptMin = 0.0, AliRsnCutPID *cutPID = 0x0);
    Int_t            GetLeadingParticleID() {return fLeading;}
    void             SetLeadingParticle(AliRsnDaughter &leading) {if (fLeading >= 0) SetDaughter(leading, fLeading);}
    Double_t         GetAverageMomentum(Int_t &count, AliRsnCutPID *cutPID = 0x0);
    Bool_t           GetAngleDistr(Double_t &angleMean, Double_t &angleRMS, AliRsnDaughter reference);

  private:

    AliVEvent       *fRef;         // pointer to input event
    AliMCEvent      *fRefMC;       // pointer to reference MC event (if any)
    Int_t            fLeading;     // index of leading track

    ClassDef(AliRsnEvent, 3);
};

#endif
