//
// Class AliRsnDaughter
//
// Interface to candidate daughters of a resonance (tracks).
// Points to the source of information, which is generally an AliVParticle-derived object
// and contains few internal data-members to store "on fly" some important information
// for the computations required during resonance analysis.
// ---
// Since the package revision, this object is not supposed to be stacked in memory
// but created "on fly" during analysis and used just for computations, as an interface.
//
// authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          M. Vala (martin.vala@cern.ch)
//

#ifndef ALIRSNDAUGHTER_H
#define ALIRSNDAUGHTER_H

#include <TLorentzVector.h>

#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDv0.h"
#include "AliAODv0.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"

typedef AliPID::EParticleType EPARTYPE;

class AliRsnDaughter : public TObject
{
  public:
    
    enum ERefType
    {
      kTrack = 0,
      kV0,
      kNoType
    };

    AliRsnDaughter();
    AliRsnDaughter(const AliRsnDaughter &copy);
    virtual ~AliRsnDaughter();
    AliRsnDaughter& operator= (const AliRsnDaughter& copy);
    
    // utilities
    void    Reset();
    void    Print(Option_t* const option = "") const;
    
    // flags and labels
    Bool_t  IsOK() const {return fOK;}
    Int_t   GetLabel() const {return fLabel;}
    Int_t   GetPDG(Bool_t abs = kTRUE) const;
    Int_t   GetMotherPDG() const {return fMotherPDG;}
    Int_t   GetID() const;
    Bool_t  HasFlag(ULong_t flag) const;
    void    SetBad() {fOK = kFALSE;}
    void    SetGood() {fOK = kTRUE;}
    void    SetLabel(Int_t label) {fLabel = label;}
    void    SetMotherPDG(Int_t value) {fMotherPDG = value;}

    // 4-momentum
    TLorentzVector& P(Bool_t mc = kFALSE) {return (mc ? fPMC : fP);}
    TLorentzVector& Prec() {return fP;}
    TLorentzVector& Psim() {return fPMC;}
    Bool_t          SetMass(Double_t mass);
    
    // charge
    Bool_t  IsPos()             const {return (fRef->Charge() > 0);}
    Bool_t  IsNeg()             const {return (fRef->Charge() < 0);}
    Bool_t  IsNeutral()         const {return (!IsPos() && !IsNeg());}
    Bool_t  IsSign(Char_t sign) const {if (sign=='+') return IsPos(); else if (sign=='-') return IsNeg(); else return IsNeutral();}
    Short_t Charge()            const {if (IsPos()) return  1 ; else if (IsNeg()) return -1 ; else return  0 ;}
    Char_t  ChargeChar()        const {if (IsPos()) return '+'; else if (IsNeg()) return '-'; else return '0';}

    // getters which automatically convert refs into allowed types
    AliVParticle*     GetRef()         const {return fRef;}
    AliVParticle*     GetRefMC()       const {return fRefMC;}
    AliMCParticle*    GetRefMCtrack()  const {return dynamic_cast<AliMCParticle*>(fRef);}
    AliESDtrack*      GetRefESDtrack() const {return dynamic_cast<AliESDtrack*>(fRef);}
    AliAODTrack*      GetRefAODtrack() const {return dynamic_cast<AliAODTrack*>(fRef);}
    AliESDv0*         GetRefESDv0()    const {return dynamic_cast<AliESDv0*>(fRef);}
    AliAODv0*         GetRefAODv0()    const {return dynamic_cast<AliAODv0*>(fRef);}
    AliMCParticle*    GetRefMCESD()    const {return dynamic_cast<AliMCParticle*>(fRefMC);}
    AliAODMCParticle* GetRefMCAOD()    const {return dynamic_cast<AliAODMCParticle*>(fRefMC);}
    Bool_t            IsMC()           const {if (GetRefMCtrack()) return kTRUE; return kFALSE;}
    Bool_t            IsAOD()          const {if (GetRefAODtrack() || GetRefAODv0()) return kTRUE; return kFALSE;}
    Bool_t            IsESD()          const {if (GetRefESDtrack() || GetRefESDv0()) return kTRUE; return kFALSE;}
    Bool_t            IsTrack()        const {if (GetRefESDtrack() || GetRefAODtrack() || GetRefMCtrack()) return kTRUE; return kFALSE;}
    Bool_t            IsV0()           const {if (GetRefESDv0() || GetRefAODv0()) return kTRUE; return kFALSE;}
    ERefType          RefType()        const {if (IsTrack()) return kTrack; if (IsV0()) return kV0; return kNoType;}
    void              SetRef(AliVParticle *ref) {fRef = ref;}
    void              SetRefMC(AliVParticle *refMC) {fRefMC = refMC;}

  private:

    Bool_t         fOK;                // internal utility flag which is kFALSE when this object should not be used
    Int_t          fLabel;             // GEANT label of corresponding MC (not trivial for V0s)
    Int_t          fMotherPDG;         // PDG code of mother (makes sense only if fRefMC is defined)
    
    TLorentzVector fP;                 // 4-vector filled with track info from default ref (if present)
    TLorentzVector fPMC;               // 4-vector filled with track info from MC ref (if present)
    
    AliVParticle  *fRef;               // reference to reconstructed track in ESD/AOD
    AliVParticle  *fRefMC;             // reference to corresponding MC particle

    ClassDef(AliRsnDaughter, 8)
};

#endif
