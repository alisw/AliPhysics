//
// Class AliRsnPairParticle
//
// Implementation of a pair of tracks, for several purposes
// - computing the total 4-momentum & inv. mass for output histos filling
// - evaluating cut checks on the pair of particles
// - evaluating any kind of kinematic value over their sum
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNPAIRPARTICLE_H
#define ALIRSNPAIRPARTICLE_H

#include <TMath.h>

#include "AliRsnDaughter.h"

class AliRsnPairParticle : public TObject
{
  public:

    AliRsnPairParticle();
    AliRsnPairParticle(const AliRsnPairParticle &obj);
    AliRsnPairParticle& operator=(const AliRsnPairParticle &obj);
    virtual ~AliRsnPairParticle();

    Double_t          GetInvMass(Double_t m1, Double_t m2);
    Double_t          GetInvMassMC(Double_t m1, Double_t m2);

    Double_t          GetEtot(Double_t m1, Double_t m2) const;
    Double_t          GetP2() const {return (fPTot[0]*fPTot[0] + fPTot[1]*fPTot[1] + fPTot[2]*fPTot[2]);}
    Double_t          GetPt2() const {return (fPTot[0]*fPTot[0] + fPTot[1]*fPTot[1]);}
    Double_t          GetP() const {return TMath::Sqrt(GetP2());}
    Double_t          GetPx() const {return fPTot[0];}
    Double_t          GetPy() const {return fPTot[1];}
    Double_t          GetPz() const {return fPTot[2];}
    Double_t          GetPt() const {return TMath::Sqrt(GetPt2());}
    Double_t          GetPhi() const {return TMath::Pi() + TMath::ATan2(-fPTot[1], -fPTot[0]);}
    Double_t          GetTheta() const {if (fPTot[2]==0.0){return TMath::PiOver2();}
      else{return TMath::ACos(fPTot[2]/GetP());}}
    Double_t          GetEta() const {return -TMath::Log(TMath::Tan(0.5*GetTheta()));}
    Double_t          GetY(Double_t m1, Double_t m2) const {return 0.5*TMath::Log((GetEtot(m1,m2)+fPTot[2])/(GetEtot(m1,m2)-fPTot[2]));}

    Double_t          GetEtotMC(Double_t m1, Double_t m2) const;
    Double_t          GetP2MC() const {return (fPTotMC[0]*fPTotMC[0] + fPTotMC[1]*fPTotMC[1] + fPTotMC[2]*fPTotMC[2]);}
    Double_t          GetPt2MC() const {return (fPTotMC[0]*fPTotMC[0] + fPTotMC[1]*fPTotMC[1]);}
    Double_t          GetPMC() const {return TMath::Sqrt(GetP2MC());}
    Double_t          GetPxMC() const {return fPTotMC[0];}
    Double_t          GetPyMC() const {return fPTotMC[1];}
    Double_t          GetPzMC() const {return fPTotMC[2];}
    Double_t          GetPtMC() const {return TMath::Sqrt(GetPt2MC());}
    Double_t          GetPhiMC() const {return TMath::Pi() + TMath::ATan2(-fPTotMC[1], -fPTotMC[0]);}
    Double_t          GetThetaMC() const {if (fPTotMC[2]==0.0){return TMath::PiOver2();}
      else{return TMath::ACos(fPTotMC[2]/GetPMC());}}
    Double_t          GetEtaMC() const {return -TMath::Log(TMath::Tan(0.5*GetThetaMC()));}
    Double_t          GetYMC(Double_t m1, Double_t m2) const {return 0.5*TMath::Log((GetEtotMC(m1,m2)+fPTotMC[2])/(GetEtotMC(m1,m2)-fPTotMC[2]));}

    Double_t          GetAngle() const;

    AliRsnDaughter*   GetDaughter(const Int_t &index) const {return fDaughter[index];}

    Bool_t            IsLabelEqual() {return abs(fDaughter[0]->GetLabel()) == abs(fDaughter[1]->GetLabel());}
    Bool_t            IsIndexEqual() {return (fDaughter[0]->GetID() == fDaughter[1]->GetID());}
    Bool_t            IsTruePair(Int_t refPDG = 0);
    Int_t             CommonMother();

    void              SetPair(AliRsnDaughter *daughter1, AliRsnDaughter *daughter2);
    void              ResetPair();
    void              RotateTrack(Int_t i, Double_t angle, Bool_t isDegrees = kTRUE);
    void              PrintInfo(const Option_t *option = "ALL");

  private:

    Double_t         fPTot[3];          // total momentum computed with rec. values
    Double_t         fPTotMC[3];        // total momentum computed with MC values
    Double_t         fPTrack[2][3];     // rec. momentum of single tracks
    Double_t         fPTrackMC[2][3];   // MC momentum of single tracks

    Int_t            fMotherLabel[2];   // GEANT label of tracks
    Int_t            fMotherPDG[2];     // PDG code of mother of tracks

    AliRsnDaughter  *fDaughter[2];      // elements of the pair

    ClassDef(AliRsnPairParticle,1)
};

#endif
