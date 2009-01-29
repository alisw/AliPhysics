//
// Class AliRsnMCInfo
//
// Contains informations from the MonteCarlo particle is associated to a track.
// It is used when looking at "perfect" PID and at "true" pairs, but the user
// does not need to access its methods.
//
// authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          M. Vala (martin.vala@cern.ch)
//

#ifndef ALIRSNMCINFO_H
#define ALIRSNMCINFO_H

#include <TParticle.h>
#include <TMath.h>

class AliRsnMCInfo : public TObject
{
  public:

    AliRsnMCInfo();
    AliRsnMCInfo(const AliRsnMCInfo &copy);

    ~AliRsnMCInfo();
    void Adopt(TParticle *part);

    // 4-momentum
    Double_t E()  const {return fEnergy;}
    Double_t E(Double_t mass) {return TMath::Sqrt(mass*mass + P2());}
    Double_t M()  const {return TMath::Sqrt(fEnergy*fEnergy - P2());}
    Double_t P2() const {return Px()*Px() + Py()*Py() + Pz()*Pz();}
    Double_t P()  const {return TMath::Sqrt(P2());}
    Double_t Px() const {return fP[0];}
    Double_t Py() const {return fP[1];}
    Double_t Pz() const {return fP[2];}
    Double_t Pt() const {return TMath::Sqrt(Px() *Px() + Py() *Py());}
    Double_t OneOverPt() const {return 1.0 / Pt();}
    Bool_t   PxPyPz(Double_t p[3]) const {p[0] = Px(); p[1] = Py(); p[2] = Pz(); return kTRUE;}

    Double_t Phi() const {return TMath::ATan2(Py(), Px());}
    Double_t Theta() const {return TMath::ATan2(Pt(), Pz());}
    Double_t Eta() const {return -TMath::Log(TMath::Tan(0.5*Theta()));}
    Double_t Y() const {return 0.5 * TMath::Log((E() + Pz()) / (E() - Pz()));}

    void     SetPx(Double_t value) {fP[0] = value;}
    void     SetPy(Double_t value) {fP[1] = value;}
    void     SetPz(Double_t value) {fP[2] = value;}
    void     SetP(Double_t px, Double_t py, Double_t pz) {SetPx(px); SetPy(py); SetPz(pz);}
    void     SetE(Double_t e) {fEnergy = e;}

    Double_t Vx() const {return fV[0];}
    Double_t Vy() const {return fV[1];}
    Double_t Vz() const {return fV[2];}
    Double_t Dr() const {return TMath::Sqrt(Vx()*Vx() + Vy()*Vy());}
    void     ShiftZero(Double_t x, Double_t y, Double_t z){fV[0]-=x;fV[1]-=y;fV[2]-=z;}

    void     SetVx(Double_t value) {fV[0] = value;}
    void     SetVy(Double_t value) {fV[1] = value;}
    void     SetVz(Double_t value) {fV[2] = value;}
    void     SetV(Double_t x, Double_t y, Double_t z) {SetVx(x); SetVy(y); SetVz(z);}

    Int_t    PDG() const {return fPDG;}
    Int_t    Mother() const {return fMother;}
    Short_t  MotherPDG() const {return fMotherPDG;}
    void     SetPDG(Int_t pdg) {fPDG = pdg;}
    void     SetMother(Int_t mlabel) {fMother = mlabel;}
    void     SetMotherPDG(Int_t pdg) {fMotherPDG = (Short_t) pdg;}

  private:

    Double_t  fP[3];          // MC momentum
    Double_t  fV[3];          // MC position
    Double_t  fEnergy;        // MC energy
    Int_t     fPDG;           // PDG code
    Int_t     fMother;        // GEANT label of mother particle
    Short_t   fMotherPDG;     // PDG code of mother particle

    ClassDef(AliRsnMCInfo,1)
};

#endif
