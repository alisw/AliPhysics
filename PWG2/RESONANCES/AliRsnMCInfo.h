#ifndef AliRsnMCInfo_h
#define AliRsnMCInfo_h

#include <AliVParticle.h>
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
    virtual Double_t E()  const {return fEnergy;}
    virtual Double_t E(Double_t mass) {return TMath::Sqrt(mass*mass + P2());}
    virtual Double_t M()  const {return TMath::Sqrt(fEnergy*fEnergy - P2());}
    virtual Double_t P2() const {return Px()*Px() + Py()*Py() + Pz()*Pz();}
    virtual Double_t P()  const {return TMath::Sqrt(P2());}
    virtual Double_t Px() const {return fP[0];}
    virtual Double_t Py() const {return fP[1];}
    virtual Double_t Pz() const {return fP[2];}
    virtual Double_t Pt() const {return TMath::Sqrt(Px() *Px() + Py() *Py());}
    virtual Double_t OneOverPt() const {return 1.0 / Pt();}
    virtual Bool_t   PxPyPz(Double_t p[3]) const {p[0] = Px(); p[1] = Py(); p[2] = Pz(); return kTRUE;}

    virtual Double_t Phi() const {return TMath::ATan2(Py(), Px());}
    virtual Double_t Theta() const {return TMath::ATan2(Pt(), Pz());}
    virtual Double_t Eta() const {return -TMath::Log(TMath::Tan(0.5*Theta()));}
    virtual Double_t Y() const {return TMath::Log((E() + Pz()) / (E() - Pz()));}

    void             SetPx(Double_t value) {fP[0] = value;}
    void             SetPy(Double_t value) {fP[1] = value;}
    void             SetPz(Double_t value) {fP[2] = value;}
    void             SetP(Double_t px, Double_t py, Double_t pz) {SetPx(px); SetPy(py); SetPz(pz);}
    void             SetE(Double_t e) {fEnergy = e;}

    Int_t     PDG() const {return fPDG;}
    Int_t     Mother() const {return fMother;}
    Short_t   MotherPDG() const {return fMotherPDG;}
    void      SetPDG(Int_t pdg) {fPDG = pdg;}
    void      SetMother(Int_t mlabel) {fMother = mlabel;}
    void      SetMotherPDG(Int_t pdg) {fMotherPDG = (Short_t) pdg;}

  private:

    Double_t  fP[3];          // MC momentum
    Double_t  fEnergy;        // MC energy
    Int_t     fPDG;           // PDG code
    Int_t     fMother;        // GEANT label of mother particle
    Short_t   fMotherPDG;     // PDG code of mother particle

    ClassDef(AliRsnMCInfo,1)
};

#endif
