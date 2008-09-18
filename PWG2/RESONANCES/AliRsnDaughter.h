//
// Class AliRsnDaughter
//
// Light-weight AOD object which contains all required track details
// which are used for resonance analysis.
// Provides converters from all kinds of input track type: ESD, AOD and MC.
//
// author: A. Pulvirenti --- email: alberto.pulvirenti@ct.infn.it
//

#ifndef ALIRSNDAUGHTER_H
#define ALIRSNDAUGHTER_H

#include <TMath.h>

#include "AliVParticle.h"
#include "AliRsnPID.h"

class TParticle;

class AliESDtrack;
class AliAODTrack;
class AliMCParticle;

class AliRsnMCInfo;

class AliRsnDaughter : public AliVParticle
{
  public:

    enum EPIDMethod
    {
      kNoPID = 0,
      kRealistic,
      kWeighted,
      kPerfect,
      kMethods
    };

    enum EPIDType
    {
      kEsd=0,
      kITS,
      kTPC,
      kTOF,
      kITS_TPC,
      kITS_TOF,
      kTPC_TOF,
      kITS_TPC_TOF,
      kITS_TPC_TOF_SP,
      kLastPIDType
    };


    AliRsnDaughter();
    AliRsnDaughter(const AliRsnDaughter &copy);
    AliRsnDaughter(AliESDtrack *track, Bool_t useTPCInnerParam = kFALSE);
    AliRsnDaughter(AliAODTrack *track);
    AliRsnDaughter(AliMCParticle *track);
    virtual ~AliRsnDaughter();
    AliRsnDaughter& operator= (const AliRsnDaughter& copy);

    // 4-momentum
    virtual Double_t E()  const {return TMath::Sqrt(fMass*fMass + P2());}
    virtual Double_t E(Double_t mass) {SetM(mass); return E();}
    virtual Double_t M()  const {return fMass;}
    virtual Double_t P2() const {return Px()*Px() + Py()*Py() + Pz()*Pz();}
    virtual Double_t P()  const {return TMath::Sqrt(P2());}
    virtual Double_t Px() const {return fP[0];}
    virtual Double_t Py() const {return fP[1];}
    virtual Double_t Pz() const {return fP[2];}
    virtual Double_t Pt() const {return TMath::Sqrt(Px()*Px() + Py()*Py());}
    virtual Double_t OneOverPt() const {return 1.0 / Pt();}
    virtual Bool_t   PxPyPz(Double_t p[3]) const {p[0] = Px(); p[1] = Py(); p[2] = Pz(); return kTRUE;}
    void             SetPx(Double_t value) {fP[0] = value;}
    void             SetPy(Double_t value) {fP[1] = value;}
    void             SetPz(Double_t value) {fP[2] = value;}
    void             SetP(Double_t px, Double_t py, Double_t pz) {SetPx(px); SetPy(py); SetPz(pz);}
    void             SetM(Double_t m) {fMass = m;}

    // DCA vertex
    virtual Double_t Xv() const {return fV[0];}
    virtual Double_t Yv() const {return fV[1];}
    virtual Double_t Zv() const {return fV[2];}
    virtual Double_t Dr() const {return TMath::Sqrt(Xv()*Xv() + Yv()*Yv());}
    virtual Bool_t   XvYvZv(Double_t x[3]) const {x[0] = Xv(); x[1] = Yv(); x[2] = Zv(); return kTRUE;}
    void             SetVx(Double_t value) {fV[0] = value;}
    void             SetVy(Double_t value) {fV[1] = value;}
    void             SetVz(Double_t value) {fV[2] = value;}
    void             SetV(Double_t vx, Double_t vy, Double_t vz) {SetVx(vx); SetVy(vy); SetVz(vz);}
    void             ShiftZero(Double_t x, Double_t y, Double_t z){fV[0]-=x;fV[1]-=y;fV[2]-=z;}

    // Orientation
    virtual Double_t Phi() const {return TMath::ATan2(Py(), Px());}
    virtual Double_t Theta() const {return TMath::ATan2(Pt(), Pz());}
    virtual Double_t Eta() const {return -TMath::Log(TMath::Tan(0.5*Theta()));}
    virtual Double_t Y() const {return TMath::Log((E() + Pz()) / (E() - Pz()));}

    // Charge
    virtual Short_t Charge() const {return fCharge;}
    void            SetCharge(Short_t value) {fCharge = value;}

    // PID
    virtual const Double_t* PID() const {return fPIDWeight;}
    const Double_t*         PIDProb() const {return fPIDProb;}
    void                    SetPIDProb(Int_t i, Double_t value);
    void                    SetPIDWeight(Int_t i, Double_t value);
    static void             SetPIDMethod(EPIDMethod method) {fgPIDMethod = method;}
    void                    AssignRealisticPID();
    AliRsnPID::EType        PIDType(Double_t &prob) const;

    // check that contains a given ESD flag
    Bool_t  CheckFlag(ULong_t flag) {return ((fFlags & flag) == flag);}

    // information getters from objects
    Bool_t  Adopt(AliESDtrack *track, EPIDType pidType = kEsd,Double_t divValue = 0.0,Bool_t useTPCInnerParam = kFALSE);
    Bool_t  Adopt(AliAODTrack *track);
    Bool_t  Adopt(AliMCParticle *track);

    // position in stack/array
    Int_t   Index() const {return fIndex;}
    Int_t   Label() const {return fLabel;}
    Int_t   GetLabel() const {return -1;}
    void    SetIndex(Int_t value) {fIndex = value;}
    void    SetLabel(Int_t value) {fLabel = value;}

    // Utilities
    void    Print(Option_t *option = "ALL") const;
    void    InitMCInfo();
    Bool_t  InitMCInfo(TParticle *particle);
    Bool_t  InitMCInfo(AliMCParticle *mcParticle);

    // MC info
    AliRsnMCInfo* GetMCInfo() const { return fMCInfo; }

    // sorting (with respect to Pt)
    virtual Bool_t IsSortable() const {return kTRUE;}
    virtual Int_t  Compare(const TObject* obj) const;

  private:

    Int_t              fIndex;    // index of source object (ESD/AOD/MC) in its collection
    Int_t              fLabel;    // label assigned to the track (act. by GEANT3)

    Short_t            fCharge;   // charge sign
    ULong_t            fFlags;    // status flags

    Double_t           fP[3];     // vector momentum (x, y, z)
    Double_t           fV[3];     // DCA vertex (x, y, z)
    Double_t           fMass;     // mass (assigned externally)

    AliRsnPID::EType   fRealisticPID;                   // PID from Bayesian probs (largest one)
    Double_t           fPIDProb[AliRsnPID::kSpecies];   // PID probabilities (Bayesian comp.)
    Double_t           fPIDWeight[AliRsnPID::kSpecies]; // PID weights

    AliRsnMCInfo      *fMCInfo;     // reference to particle object (if any)

    static EPIDMethod  fgPIDMethod; // flag to define how the PID is computed for this object

    static void        GetESDPID(AliESDtrack *track,Double_t *pid,EPIDType pidType = kEsd,Double_t divValue=-1.0,Double_t val=-1.0);

    ClassDef(AliRsnDaughter, 4);
};

#endif
