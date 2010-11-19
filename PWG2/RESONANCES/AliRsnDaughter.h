//
// Class AliRsnDaughter
//
// Interface to candidate daughters of a resonance (tracks or V0s).
// It contains two pointers to two AliVParticle-derived objects.
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
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliMCParticle.h"

#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODMCParticle.h"

typedef AliPID::EParticleType EPARTYPE;

class AliRsnDaughter : public TObject
{
  public:
    
    enum ERefType
    {
      kTrack = 0,
      kV0,
      kCascade,
      kNoType
    };

    AliRsnDaughter();
    AliRsnDaughter(const AliRsnDaughter &copy);
    AliRsnDaughter& operator= (const AliRsnDaughter& copy);
    virtual ~AliRsnDaughter() { /*nothing*/ }
    
    // basic getters (for data members)
    Bool_t          IsOK() const {return fOK;}
    Int_t           GetLabel() const {return fLabel;}
    Int_t           GetMotherPDG() const {return fMotherPDG;}
    TLorentzVector& Prec() {return fPrec;}
    TLorentzVector& Psim() {return fPsim;}
    TLorentzVector& P(Bool_t mc = kFALSE) {return (mc ? fPsim : fPrec);}
    AliVParticle*   GetRef() {return fRef;}
    AliVParticle*   GetRefMC() {return fRefMC;}
    
    // basic setters (for data members)
    void    SetBad() {fOK = kFALSE;}
    void    SetGood() {fOK = kTRUE;}
    void    SetLabel(Int_t label) {fLabel = label;}
    void    SetMotherPDG(Int_t value) {fMotherPDG = value;}
    void    SetRef(AliVParticle *p) {fRef = p;}
    void    SetRefMC(AliVParticle *p = 0x0) {fRefMC = p;}
    
    // charge checkers
    Bool_t  IsPos()             const {return (fRef->Charge() > 0);}
    Bool_t  IsNeg()             const {return (fRef->Charge() < 0);}
    Bool_t  IsNeutral()         const {return (!IsPos() && !IsNeg());}
    Bool_t  IsSign(Char_t sign) const {if (sign=='+') return IsPos(); else if (sign=='-') return IsNeg(); else return IsNeutral();}
    Short_t Charge()            const {if (IsPos()) return  1 ; else if (IsNeg()) return -1 ; else return  0 ;}
    Char_t  ChargeChar()        const {if (IsPos()) return '+'; else if (IsNeg()) return '-'; else return '0';}
    
    // utilities
    void     Reset();
    Int_t    GetPDG(Bool_t abs = kTRUE);
    Int_t    GetID();
    Bool_t   HasFlag(ULong_t flag);
    Bool_t   SetMass(Double_t mass);
    Bool_t   IsKinkDaughter();

    // getters which automatically convert refs into allowed types
    AliESDtrack*      GetRefESDtrack()   {return dynamic_cast<AliESDtrack*>(fRef);}
    AliESDv0*         GetRefESDv0()      {return dynamic_cast<AliESDv0*>(fRef);}
    AliESDcascade*    GetRefESDcascade() {return dynamic_cast<AliESDcascade*>(fRef);}
    AliAODTrack*      GetRefAODtrack()   {return dynamic_cast<AliAODTrack*>(fRef);}
    AliAODv0*         GetRefAODv0()      {return dynamic_cast<AliAODv0*>(fRef);}
    AliAODcascade*    GetRefAODcascade() {return dynamic_cast<AliAODcascade*>(fRef);}
    AliMCParticle*    GetRefMCtrack()    {return dynamic_cast<AliMCParticle*>(fRef);}
    AliMCParticle*    GetRefMCESD()      {return dynamic_cast<AliMCParticle*>(fRefMC);}
    AliAODMCParticle* GetRefMCAOD()      {return dynamic_cast<AliAODMCParticle*>(fRefMC);}
    
    // check the input type
    Bool_t    IsMC()       {if (GetRefMCtrack()) return kTRUE; return kFALSE;}
    Bool_t    IsAOD()      {if (GetRefAODtrack() || GetRefAODv0()) return kTRUE; return kFALSE;}
    Bool_t    IsESD()      {if (GetRefESDtrack() || GetRefESDv0()) return kTRUE; return kFALSE;}
    Bool_t    IsTrack()    {if (GetRefESDtrack() || GetRefAODtrack() || GetRefMCtrack()) return kTRUE; return kFALSE;}
    Bool_t    IsV0()       {if (GetRefESDv0() || GetRefAODv0()) return kTRUE; return kFALSE;}
    Bool_t    IsCascade()  {if (GetRefESDcascade() || GetRefAODcascade()) return kTRUE; return kFALSE;}
    ERefType  RefType()    {if (IsTrack()) return kTrack; if (IsV0()) return kV0; if (IsCascade()) return kCascade; return kNoType;}

  private:

    Bool_t         fOK;          // internal utility flag which is kFALSE when this object should not be used
    Int_t          fLabel;       // GEANT label of corresponding MC
    Int_t          fMotherPDG;   // PDG code of mother (makes sense only if fRefMC is defined)
    
    TLorentzVector fPrec;        // 4-vector filled with track info from default ref (if present)
    TLorentzVector fPsim;        // 4-vector filled with track info from MC ref (if present)
    
    AliVParticle  *fRef;         // reference to reconstructed track in ESD/AOD
    AliVParticle  *fRefMC;       // reference to corresponding MC particle

    ClassDef(AliRsnDaughter, 8)
};

#endif
