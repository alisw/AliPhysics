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

#include "AliPID.h"
#include "AliRsnPIDDefESD.h"
#include "AliRsnDaughter.h"

class AliVEvent;
class AliMCEvent;

class AliRsnEvent : public TObject
{
  public:

    AliRsnEvent(AliVEvent *ref = 0, AliMCEvent *refMC = 0);
    AliRsnEvent(const AliRsnEvent& copy);
    AliRsnEvent& operator= (const AliRsnEvent& copy);
    virtual ~AliRsnEvent();

    void             SetRef(AliVEvent *const event, AliMCEvent *const mc = 0) {fRef = event; fRefMC = mc;}
    void             SetRefMC(AliMCEvent * const mc) {fRefMC = mc;}
    AliVEvent*       GetRef() {return fRef;}
    AliMCEvent*      GetRefMC() {return fRefMC;}
    AliRsnPIDDefESD* GetPIDDefESD() {return &fPIDDefESD;}

    void             SetDaughter(AliRsnDaughter &daughter, Int_t index);
    AliRsnDaughter   GetDaughter(Int_t i);
    Int_t            GetMultiplicity();
    Double_t         GetVz();
    AliRsnDaughter   GetLeadingParticle(Double_t ptMin = 0.0, AliPID::EParticleType type = AliPID::kUnknown);
    Double_t         GetAverageMomentum(Int_t &count, AliPID::EParticleType type = AliPID::kUnknown);
    Bool_t           GetAngleDistr(Double_t &angleMean, Double_t &angleRMS, AliRsnDaughter d);

    void             SetPriorProbability(AliPID::EParticleType type, Double_t p) {if (type>=0&&type<(Int_t)AliPID::kSPECIES)fPrior[type]=p;}
    void             SetPriorProbability(Double_t* const out);
    void             DumpPriors();
    void             GetPriorProbability(Double_t *out) const;

  private:

    Bool_t AcceptTrackPID(AliRsnDaughter*const d, AliPID::EParticleType type = AliPID::kUnknown);

    Double_t         fPrior[AliPID::kSPECIES]; // prior probabilities 
    AliVEvent       *fRef;                     // pointer to input event
    AliMCEvent      *fRefMC;                   // pointer to reference MC event (if any)
    AliRsnPIDDefESD  fPIDDefESD;               // (optional) customization of PID weights for ESD

    ClassDef(AliRsnEvent, 3);
};

#endif
