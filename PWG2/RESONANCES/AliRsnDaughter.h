#ifndef ALIRSNDAUGHTER_H
#define ALIRSNDAUGHTER_H

//
//  Interface to single daughter candidate.
//

#include <TMath.h>
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

class AliRsnEvent;

typedef AliPID::EParticleType EPARTYPE;

class AliRsnDaughter : public TObject {
public:

   enum ERefType {
      kTrack,
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

   AliRsnDaughter() : fOK(kFALSE), fLabel(-1), fMotherPDG(0), fRsnID(-1),
                      fPrec(), fPsim(), fRef(0x0), fRefMC(0x0), fOwnerEvent(0x0) { }
   AliRsnDaughter(const AliRsnDaughter &copy);
   AliRsnDaughter& operator= (const AliRsnDaughter& copy);
   virtual ~AliRsnDaughter() { /*empty, since pointers must not be deleted*/ }

   // basic getters
   Bool_t          IsOK()          const {return fOK;}
   Int_t           GetLabel()      const {return fLabel;}
   Int_t           GetMotherPDG()  const {return fMotherPDG;}
   Int_t           GetRsnID()      const {return fRsnID;}
   AliVParticle*   GetRef()              {return fRef;}
   AliVParticle*   GetRefMC()            {return fRefMC;}
   AliRsnEvent*    GetOwnerEvent()       {return fOwnerEvent;}
   TLorentzVector& P(Bool_t mc)          {return (mc ? fPsim : fPrec);}
   TLorentzVector& Prec()                {return fPrec;}
   TLorentzVector& Psim()                {return fPsim;}
   Int_t           GetPDG();
   Int_t           GetPDGAbs()           {return TMath::Abs(GetPDG());}
   Int_t           GetID();
   Int_t           GetMother();

   // basic setters (for data members)
   void  Reset();
   void  SetOK(Bool_t ok)              {fOK = ok;}
   void  SetBad()                      {fOK = kFALSE;}
   void  SetGood()                     {fOK = kTRUE;}
   void  SetLabel(Int_t label)         {fLabel = label;}
   void  SetMotherPDG(Int_t value)     {fMotherPDG = value;}
   void  SetRsnID(Int_t id)            {fRsnID = id;}
   void  SetRef(AliVParticle *p)       {fRef = p;}
   void  SetRefMC(AliVParticle *p)     {fRefMC = p;}
   void  SetOwnerEvent(AliRsnEvent *e) {fOwnerEvent = e;}
   
   // additional functions
   void  FillP(Double_t mass);
   void  Print(Option_t *o = "") const;
   
   // getters related to charge
   Bool_t   IsPos()      const {if (fRef) return (fRef->Charge() > 0); return kFALSE;}
   Bool_t   IsNeg()      const {if (fRef) return (fRef->Charge() < 0); return kFALSE;}
   Bool_t   IsCharged()  const {if (fRef) return (IsPos() || IsNeg()); return kFALSE;}
   Bool_t   IsNeutral()  const {if (fRef) return (!IsCharged());       return kFALSE;}
   Short_t  ChargeS()    const {if (IsPos()) return  1 ; else if (IsNeg()) return -1 ; else return  0 ;}
   Char_t   ChargeC()    const {if (IsPos()) return '+'; else if (IsNeg()) return '-'; else return '0';}
   
    // getters which automatically convert refs into allowed types
   static Bool_t     Match(AliVParticle *p, TClass *ref) {if (p) return (p->InheritsFrom(ref)); return kFALSE;}
   Bool_t            MatchRef(TClass *ref)               {return Match(fRef, ref);}
   Bool_t            MatchRefMC(TClass *ref)             {return Match(fRefMC, ref);}
   ERefType          RefType();
   Bool_t            IsESD();
   Bool_t            IsAOD();

   AliVTrack*        Ref2Vtrack()        {if (Match(fRef, AliVTrack    ::Class())) return static_cast<AliVTrack*>       (fRef); return 0x0;}
   AliESDtrack*      Ref2ESDtrack()      {if (Match(fRef, AliESDtrack  ::Class())) return static_cast<AliESDtrack*>     (fRef); return 0x0;}
   AliAODTrack*      Ref2AODtrack()      {if (Match(fRef, AliAODTrack  ::Class())) return static_cast<AliAODTrack*>     (fRef); return 0x0;}
   AliMCParticle*    Ref2MCparticle()    {if (Match(fRef, AliMCParticle::Class())) return static_cast<AliMCParticle*>   (fRef); return 0x0;}
   AliAODMCParticle* Ref2AODMCparticle() {if (Match(fRef, AliMCParticle::Class())) return static_cast<AliAODMCParticle*>(fRef); return 0x0;}
      
   AliESDv0*         Ref2ESDv0()         {if (Match(fRef, AliESDv0     ::Class())) return static_cast<AliESDv0*>(fRef); return 0x0;}
   AliAODv0*         Ref2AODv0()         {if (Match(fRef, AliAODv0     ::Class())) return static_cast<AliAODv0*>(fRef); return 0x0;}
   
   AliESDcascade*    Ref2ESDcascade()    {if (Match(fRef, AliESDcascade::Class())) return static_cast<AliESDcascade*>(fRef); return 0x0;}
   AliAODcascade*    Ref2AODcascade()    {if (Match(fRef, AliAODcascade::Class())) return static_cast<AliAODcascade*>(fRef); return 0x0;}
   
   AliMCParticle*    RefMC2ESD()         {if (Match(fRefMC, AliMCParticle   ::Class())) return static_cast<AliMCParticle*>   (fRef)  ; return 0x0;}
   AliAODMCParticle* RefMC2AOD()         {if (Match(fRefMC, AliAODMCParticle::Class())) return static_cast<AliAODMCParticle*>(fRefMC); return 0x0;}
   
   // static functions related to internal ESpecies enum
   static ERefType    RefType(ESpecies species);
   static Bool_t      IsCharged(ESpecies species) {return (species <= kProton);}
   static const char* SpeciesName(ESpecies species);
   static Int_t       SpeciesPDG(ESpecies species);
   static Double_t    SpeciesMass(ESpecies species);
   static EPARTYPE    ToAliPID(ESpecies species);
   static ESpecies    FromAliPID(EPARTYPE species);

private:

   Bool_t         fOK;          // internal utility flag which is kFALSE when this object should not be used
   Int_t          fLabel;       // index of MC particle
   Int_t          fMotherPDG;   // PDG code of mother (makes sense only if fRefMC is defined)
   Int_t          fRsnID;       // internal ID for monitoring purposes
   
   TLorentzVector fPrec;        // 4-momentum for rec
   TLorentzVector fPsim;        // 4-momentum for MC

   AliVParticle  *fRef;         // reference to reconstructed object
   AliVParticle  *fRefMC;       // reference to corresponding MC particle
   AliRsnEvent   *fOwnerEvent;  // pointer to owner event

   ClassDef(AliRsnDaughter, 12)
};

//__________________________________________________________________________________________________
inline AliRsnDaughter::ERefType AliRsnDaughter::RefType()
{
//
// Returns the enum value corresponding to the real nature
// of the object pointed by the fRef data member
//

   if (Match(fRef, AliESDtrack  ::Class())) return kTrack;
   if (Match(fRef, AliESDv0     ::Class())) return kV0;
   if (Match(fRef, AliESDcascade::Class())) return kCascade;
   if (Match(fRef, AliAODTrack  ::Class())) return kTrack;
   if (Match(fRef, AliAODv0     ::Class())) return kV0;
   if (Match(fRef, AliAODcascade::Class())) return kCascade;
   if (Match(fRef, AliMCParticle::Class())) return kTrack;
   
   return kNoType;
}

//__________________________________________________________________________________________________
inline Bool_t AliRsnDaughter::IsESD()
{
//
// Tells if the object pointed by fRef data member is ESD
// NOTE: it is true even when fRef is the MC corresponding 
//       object (= AliMCParticle)
//

   if (Match(fRef, AliESDtrack  ::Class())) return kTRUE;
   if (Match(fRef, AliESDv0     ::Class())) return kTRUE;
   if (Match(fRef, AliESDcascade::Class())) return kTRUE;
   if (Match(fRef, AliMCParticle::Class())) return kTRUE;
   
   return kFALSE;
}

//__________________________________________________________________________________________________
inline Bool_t AliRsnDaughter::IsAOD()
{
//
// Tells if the object pointed by fRef data member is AOD
// NOTE: it is true even when fRef is the MC corresponding 
//       object (= AliAODMCParticle)
//

   if (Match(fRef, AliAODTrack     ::Class())) return kTRUE;
   if (Match(fRef, AliAODv0        ::Class())) return kTRUE;
   if (Match(fRef, AliAODcascade   ::Class())) return kTRUE;
   if (Match(fRef, AliAODMCParticle::Class())) return kTRUE;
   
   return kFALSE;
}

//__________________________________________________________________________________________________
inline AliRsnDaughter::ERefType AliRsnDaughter::RefType(ESpecies species)
{
//
// Returns the expected object type for a candidate daughter
// of the given species.
//

   switch (species) {
      case kElectron:
      case kMuon:
      case kPion:
      case kKaon:
      case kProton:
         return kTrack;
      case kKaon0:
      case kLambda:
         return kV0;
      case kXi:
      case kOmega:
         return kCascade;
      default:
         return kNoType;
   }
}

//__________________________________________________________________________________________________
inline void AliRsnDaughter::FillP(Double_t mass)
{
//
// Fills the 4-momentum data member using the values in the
// AliVParticle pointer data members, choosing according to arguments
//

   if (fRefMC) fPsim.SetXYZM(fRefMC->Px(), fRefMC->Py(), fRefMC->Pz(), mass);
   if (fRef)   fPrec.SetXYZM(fRef  ->Px(), fRef  ->Py(), fRef  ->Pz(), mass);
}

#endif
