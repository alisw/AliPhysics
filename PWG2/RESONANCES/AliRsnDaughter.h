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

#include <TMath.h>

#include "AliPID.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliMCParticle.h"

class TParticle;
class AliStack;
class AliMCParticle;
// class AliESDtrack;
// class AliESDtrack;
// class AliAODTrack;
class AliAODEvent;
class AliRsnPIDDefESD;

class AliRsnDaughter : public TObject
{
  public:

    enum EPIDMethod {
      kNoPID = 0,
      kRealistic,
      kPerfect,
      kMethods
    };

    AliRsnDaughter(AliVParticle *ref = 0, TParticle *refMC = 0);
    AliRsnDaughter(const AliRsnDaughter &copy);
    virtual ~AliRsnDaughter();
    AliRsnDaughter& operator= (const AliRsnDaughter& copy);

    // momentum
    Double_t Px() const {return fRef->Px();}
    Double_t Py() const {return fRef->Py();}
    Double_t Pz() const {return fRef->Pz();}
    Double_t Pt() const {return fRef->Pt();}
    Double_t P2() const {return Pt()*Pt() + Pz()*Pz();}
    Double_t P()  const {return TMath::Sqrt(P2());}
    Double_t Eta() const {return fRef->Eta();}
    Double_t Y() const {return fRef->Y();}
    Double_t Y(Double_t mass) const {return 0.5*TMath::Log((E(mass) + Pz()) / (E(mass) - Pz()));}
    Double_t E() const {return fRef->E();}
    Double_t E(Double_t mass) const {return TMath::Sqrt(mass*mass + P2());}
    Double_t Phi() const {return fRef->Phi();}
    Double_t Theta() const {return fRef->Theta();}
    Double_t PhiDeg() const {return TMath::RadToDeg() * Phi();}
    Double_t ThetaDeg() const {return TMath::RadToDeg() * Theta();}
    void     RotateP(Double_t angle, Double_t &x, Double_t &y, Bool_t isDegrees = kTRUE);
    Double_t AngleTo(AliRsnDaughter d, Bool_t outInDegrees = kTRUE);

    // DCA vertex
    Double_t Xv() const {return fRef->Xv();}
    Double_t Yv() const {return fRef->Yv();}
    Double_t Zv() const {return fRef->Zv();}
    Double_t Dr() const {return fDr;}
    Double_t Dz() const {return fDz;}
    void     SetDr(Double_t value) {fDr = value;}
    void     SetDz(Double_t value) {fDz = value;}

    // PID
    const Double_t        *PID() const {return fRef->PID();}
    const Double_t        *ComputedWeights() const {return fPID;}
    Bool_t                 CombineWithPriors(const Double_t* priors, AliRsnPIDDefESD* pidDef = 0x0);
    AliPID::EParticleType  RealisticPID() const;
    AliPID::EParticleType  PerfectPID() const;
    Bool_t                 IsPerfectPID(const AliPID::EParticleType type) const  {return (PerfectPID() == type);}
    Bool_t                 IsRealisticPID(const AliPID::EParticleType type) const  {return (RealisticPID() == type);}
    Double_t               PIDProb(AliPID::EParticleType type) const {return PID()[(Int_t)type];}
    AliPID::EParticleType  PIDType(Double_t &prob) const;
    AliPID::EParticleType  AssignedPID() const;
    AliPID::EParticleType  RequiredPID() const {return fReqPID;}
    void                   SetRequiredPID(AliPID::EParticleType type) {fReqPID = type;}

    // integer parameters
    Bool_t  IsPos() const {return (fRef->Charge() > 0);}
    Bool_t  IsNeg() const {return (fRef->Charge() < 0);}
    Bool_t  IsNeutral() const {return (!IsPos() && !IsNeg());}
    Bool_t  IsSign(Char_t sign) { if (sign=='+') return IsPos(); else if (sign=='-') return IsNeg(); else return IsNeutral();}
    Short_t Charge() const {if (IsPos()) return 1; else if (IsNeg()) return -1; else return 0;}
    Char_t  ChargeC() const {if (IsPos()) return '+'; else if (IsNeg()) return '-'; else return '0';}
    Int_t   GetLabel() const {return fRef->GetLabel();}
    Int_t   GetID() const;
    void    SetStatus(ULong_t value) {fStatus = value;}
    ULong_t GetStatus() const {return fStatus;}
    Bool_t  CheckFlag(ULong_t flag) const {return ((fStatus & flag) > 0);}
    void    SetGood() {fOK = kTRUE;}
    void    SetBad() {fOK = kFALSE;}
    Bool_t  IsOK() const { return fOK; }

    // Kinkness
    Char_t  KinkIndex() const {return fKinkIndex;}
    Bool_t  IsKinkMother() const {return (fKinkIndex < 0);}
    Bool_t  IsKinkDaughter() const {return (fKinkIndex > 0);}
    Bool_t  IsKink() const {return (IsKinkMother() || IsKinkDaughter());}
    void    SetKink(Char_t kink) {fKinkIndex = kink;}
    void    SetKinkMother() {fKinkIndex = -1;}
    void    SetKinkDaughter() {fKinkIndex = 1;}
    void    SetNoKink() {fKinkIndex = 0;}
    void    FindKinkIndex(const AliESDtrack* track);
    void    FindKinkIndex(AliAODEvent*const event);

    // MC info & references
    AliVParticle*  GetRef() const {return fRef;}
    AliESDtrack*   GetRefESD() {return dynamic_cast<AliESDtrack*>(fRef);}
    AliAODTrack*   GetRefAOD() {return dynamic_cast<AliAODTrack*>(fRef);}
    AliMCParticle* GetRefMC()  {return dynamic_cast<AliMCParticle*>(fRef);}
    TParticle*     GetParticle() const {return fParticle;}
    Int_t          GetMotherPDG() const {return fMotherPDG;}
    void           SetRef(AliVParticle * const ref) {fRef = ref;}
    void           SetParticle(TParticle * const p) {fParticle = p;}
    void           SetMotherPDG(Int_t value) {fMotherPDG = value;}
    void           FindMotherPDG(AliStack *const stack);
    Double_t       GetMCEnergy(Double_t mass);

    // utilities
    void          Reset();
    void          Print(Option_t*const option = "ALL") const;

    // static functions
    static void                   SetPIDMethod(EPIDMethod method) {fgPIDMethod = method;}
    static AliPID::EParticleType  InternalType(Int_t pdgCode);
    static const char*            MethodName(EPIDMethod method);

  private:

    Bool_t         fOK;                     // status flag for usability
    Int_t          fKinkIndex;              // indicator of kinkness of the track
    TParticle     *fParticle;               // pointer to (eventual) MC information
    Int_t          fMotherPDG;              // PDG code of mother (if any)
    ULong_t        fStatus;                 // track status (if available)

    Double_t       fDr;                     // transverse impact parameter
    Double_t       fDz;                     // longitudinal impact parameter
    Double_t       fPID[AliPID::kSPECIES];  // PID probabilities
    AliPID::EParticleType fReqPID;          // PID assigned by pairdef

    AliVParticle  *fRef;                    // reference to read object

    static EPIDMethod fgPIDMethod;          // PID method used for analysis

    ClassDef(AliRsnDaughter, 7)
};

#endif
