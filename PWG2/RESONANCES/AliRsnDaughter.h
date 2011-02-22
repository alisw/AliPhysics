//
// Class AliRsnDaughter
//
// Interface to candidate daughters of a resonance (tracks).
// Points to the source of information, which is generally an AliVParticle-derived object
// and contains few internal data-members to store "on fly" some important information
// for the computations required during resonance analysis.
// It contains a TLorentzVector data-member which, provided that a meaningful mass was assigned,
// eases a lot the computation of invariant masses from summing up several of these objects.
//
// authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          M. Vala (martin.vala@cern.ch)
//

#ifndef ALIRSNDAUGHTER_H
#define ALIRSNDAUGHTER_H

#include <TLorentzVector.h>

#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODMCParticle.h"

typedef AliPID::EParticleType EPARTYPE;

class AliRsnDaughter : public TObject {
public:

   enum ERefType {
      kTrack = 0,
      kV0,
      kCascade,
      kNoType
   };

   AliRsnDaughter();
   AliRsnDaughter(const AliRsnDaughter &copy);
   AliRsnDaughter& operator= (const AliRsnDaughter& copy);
   virtual ~AliRsnDaughter() { /*empty, since pointers must not be deleted*/ }

   // basic getters (for data members)
   Bool_t          IsOK()         const {return fOK;}
   Int_t           GetLabel()     const {return fLabel;}
   Int_t           GetMotherPDG() const {return fMotherPDG;}
   Int_t           GetRsnID()     const {return fRsnID;}
   TLorentzVector& Prec()               {return fPrec;}
   TLorentzVector& Psim()               {return fPsim;}
   TLorentzVector& P(Bool_t mc)         {return (mc ? fPsim : fPrec);}
   AliVParticle*   GetRef()             {return fRef;}
   AliVParticle*   GetRefMC()           {return fRefMC;}

   // basic setters (for data members)
   void    SetBad()                  {fOK = kFALSE;}
   void    SetGood()                 {fOK = kTRUE;}
   void    SetLabel(Int_t label)     {fLabel = label;}
   void    SetRsnID(Int_t id)        {fRsnID = id;}
   void    SetMotherPDG(Int_t value) {fMotherPDG = value;}
   void    SetRef(AliVParticle *p);
   void    SetRefMC(AliVParticle *p);

   // charge checkers
   Bool_t  IsPos()             const {return (fRef->Charge() > 0);}
   Bool_t  IsNeg()             const {return (fRef->Charge() < 0);}
   Bool_t  IsNeutral()         const {return (!IsPos() && !IsNeg());}
   Bool_t  IsSign(Char_t sign) const {if (sign == '+') return IsPos(); else if (sign == '-') return IsNeg(); else return IsNeutral();}
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
   AliVTrack*        GetRefVtrack()     {if (classMatchRef  (AliVTrack       ::Class())) return static_cast<AliVTrack*>       (fRef)  ; return 0x0;}
   AliESDtrack*      GetRefESDtrack()   {if (classMatchRef  (AliESDtrack     ::Class())) return static_cast<AliESDtrack*>     (fRef)  ; return 0x0;}
   AliESDv0*         GetRefESDv0()      {if (classMatchRef  (AliESDv0        ::Class())) return static_cast<AliESDv0*>        (fRef)  ; return 0x0;}
   AliESDcascade*    GetRefESDcascade() {if (classMatchRef  (AliESDcascade   ::Class())) return static_cast<AliESDcascade*>   (fRef)  ; return 0x0;}
   AliAODTrack*      GetRefAODtrack()   {if (classMatchRef  (AliAODTrack     ::Class())) return static_cast<AliAODTrack*>     (fRef)  ; return 0x0;}
   AliAODv0*         GetRefAODv0()      {if (classMatchRef  (AliAODv0        ::Class())) return static_cast<AliAODv0*>        (fRef)  ; return 0x0;}
   AliAODcascade*    GetRefAODcascade() {if (classMatchRef  (AliAODcascade   ::Class())) return static_cast<AliAODcascade*>   (fRef)  ; return 0x0;}
   AliMCParticle*    GetRefMCtrack()    {if (classMatchRef  (AliMCParticle   ::Class())) return static_cast<AliMCParticle*>   (fRef)  ; return 0x0;}
   AliMCParticle*    GetRefMCESD()      {if (classMatchRefMC(AliMCParticle   ::Class())) return static_cast<AliMCParticle*>   (fRefMC); return 0x0;}
   AliAODMCParticle* GetRefMCAOD()      {if (classMatchRefMC(AliAODMCParticle::Class())) return static_cast<AliAODMCParticle*>(fRefMC); return 0x0;}

   // check the input type
   Bool_t    IsMC()      {if (GetRefMCtrack()) return kTRUE; return kFALSE;}
   Bool_t    IsAOD()     {if (GetRefAODtrack() || GetRefAODv0()) return kTRUE; return kFALSE;}
   Bool_t    IsESD()     {if (GetRefESDtrack() || GetRefESDv0()) return kTRUE; return kFALSE;}
   Bool_t    IsTrack()   {if (GetRefESDtrack() || GetRefAODtrack() || GetRefMCtrack()) return kTRUE; return kFALSE;}
   Bool_t    IsV0()      {if (GetRefESDv0() || GetRefAODv0()) return kTRUE; return kFALSE;}
   Bool_t    IsCascade() {if (GetRefESDcascade() || GetRefAODcascade()) return kTRUE; return kFALSE;}
   ERefType  RefType()   {if (IsTrack()) return kTrack; if (IsV0()) return kV0; if (IsCascade()) return kCascade; return kNoType;}

private:

   Bool_t classMatchRef  (TClass *ref) {if (fRef  ) return (fRef  ->IsA() == ref); return kFALSE;}
   Bool_t classMatchRefMC(TClass *ref) {if (fRefMC) return (fRefMC->IsA() == ref); return kFALSE;}

   Bool_t         fOK;          // internal utility flag which is kFALSE when this object should not be used
   Int_t          fLabel;       // GEANT label of corresponding MC
   Int_t          fMotherPDG;   // PDG code of mother (makes sense only if fRefMC is defined)
   Int_t          fRsnID;       // id in rsn manager cout

   TLorentzVector fPrec;        // 4-vector filled with track info from default ref (if present)
   TLorentzVector fPsim;        // 4-vector filled with track info from MC ref (if present)

   AliVParticle  *fRef;         // reference to reconstructed track in ESD/AOD
   AliVParticle  *fRefMC;       // reference to corresponding MC particle

   ClassDef(AliRsnDaughter, 9)
};



#endif
