#ifndef ALIPARTSIMPLEFORCORR_H
#define ALIPARTSIMPLEFORCORR_H

#include "AliVParticle.h"
#include "TString.h"
#include "AliLog.h"
#include "TParticle.h"

class AliPartSimpleForCorr : public AliVParticle
{
  public:
        AliPartSimpleForCorr();

        AliPartSimpleForCorr(Short_t charge, Float_t eta, Float_t phi, Float_t pt, Int_t ID, Int_t ID1, Int_t ID2, Short_t candidate, Double_t multiplicity);

        AliPartSimpleForCorr(Float_t eta, Float_t phi, Double_t multiplicity);

        AliPartSimpleForCorr(Float_t eta, Float_t phi, Float_t pt, Double_t mass);

        virtual ~AliPartSimpleForCorr() {}

        virtual Double_t Px() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Double_t Py() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Double_t Pz() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Double_t Pt() const { return fpT; }
        virtual Double_t P() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Bool_t PxPyPz(Double_t[3]) const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Double_t Xv() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Double_t Yv() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Double_t Zv() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Bool_t XvYvZv(Double_t[3]) const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Double_t OneOverPt() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Double_t Phi() const { return fPhi; }
        virtual Double_t Theta() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Double_t E() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Double_t M() const { return fMass; }
        virtual Double_t Eta() const { return fEta; }
        virtual Double_t Y() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Short_t Charge() const { return fCharge; }
        virtual Int_t GetLabel() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Int_t PdgCode() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual const Double_t *PID() const {
          AliFatal("Not implemented");
          return 0;
        }
        virtual Short_t       WhichCandidate() const { return fCandidate; }
        virtual Int_t         GetID() const { return fID; }
        virtual Int_t         GetIDFirstDaughter() const { return fID1; }
        virtual Int_t         GetIDSecondDaughter() const { return fID2; }
        virtual Double_t      Multiplicity() const { return fMultiplicity; }

  private:

        Short_t     fCharge;
        Float_t     fEta;
        Float_t     fPhi;
        Float_t     fpT;
        Int_t       fID;
        Short_t     fCandidate;
        Double_t    fMultiplicity;
        Int_t       fID1;
        Int_t       fID2;
        Double_t    fMass;

        ClassDef(AliPartSimpleForCorr, 2);
};
#endif
