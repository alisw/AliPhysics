#ifndef ALIRSNDAUGHTER_H
#define ALIRSNDAUGHTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Interface to single daughter candidate.
//
////////////////////////////////////////////////////////////////////////////////

#include <TLorentzVector.h>

#include "AliPID.h"
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
   
   enum ESpecies {
      kElectron,
      kMuon,
      kPion,
      kKaon,
      kProton,
      kKaon0,
      kLambda,
      kXi,
      kOmega,
      kUnknown
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
   Short_t ChargeS()           const {if (IsPos()) return  1 ; else if (IsNeg()) return -1 ; else return  0 ;}
   Char_t  ChargeC()           const {if (IsPos()) return '+'; else if (IsNeg()) return '-'; else return '0';}

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
   Bool_t    IsAOD()     {if (GetRefAODtrack() || GetRefAODv0() || GetRefESDcascade()) return kTRUE; return kFALSE;}
   Bool_t    IsESD()     {if (GetRefESDtrack() || GetRefESDv0() || GetRefAODcascade()) return kTRUE; return kFALSE;}
   Bool_t    IsTrack()   {if (GetRefESDtrack() || GetRefAODtrack() || GetRefMCtrack()) return kTRUE; return kFALSE;}
   Bool_t    IsV0()      {if (GetRefESDv0() || GetRefAODv0()) return kTRUE; return kFALSE;}
   Bool_t    IsCascade() {if (GetRefESDcascade() || GetRefAODcascade()) return kTRUE; return kFALSE;}
   ERefType  RefType()   {if (IsTrack()) return kTrack; if (IsV0()) return kV0; if (IsCascade()) return kCascade; return kNoType;}
   
   // utilities
   void     Reset();
   Int_t    GetPDG(Bool_t abs = kTRUE);
   Int_t    GetID();
   Bool_t   SetMass(Double_t mass);
   Bool_t   IsKinkDaughter();
   
   // static utilities
   static ERefType    RefType(ESpecies species);
   static Bool_t      IsCharged(ESpecies species) {return (species <= kProton);}
   static const char* SpeciesName(ESpecies species);
   static Int_t       SpeciesPDG(ESpecies species);
   static Double_t    SpeciesMass(ESpecies species);
   static EPARTYPE    ToAliPID(ESpecies species);
   static ESpecies    FromAliPID(EPARTYPE species);

private:

   Bool_t classMatchRef  (TClass *ref) {if (fRef  ) return (fRef  ->InheritsFrom(ref)); return kFALSE;}
   Bool_t classMatchRefMC(TClass *ref) {if (fRefMC) return (fRefMC->InheritsFrom(ref)); return kFALSE;}

   Bool_t         fOK;          // internal utility flag which is kFALSE when this object should not be used
   Int_t          fLabel;       // GEANT label of corresponding MC particle
   Int_t          fMotherPDG;   // PDG code of mother (makes sense only if fRefMC is defined)
   Int_t          fRsnID;       // internal ID for monitoring purposes

   TLorentzVector fPrec;        // 4-vector filled with track info from default ref (if present)
   TLorentzVector fPsim;        // 4-vector filled with track info from MC ref (if present)

   AliVParticle  *fRef;         // reference to reconstructed object
   AliVParticle  *fRefMC;       // reference to corresponding MC particle

   ClassDef(AliRsnDaughter, 9)
};



#endif
