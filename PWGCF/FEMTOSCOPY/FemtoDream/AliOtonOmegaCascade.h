/*
 * AliOtonOmegaCascade.h
 *
 *  Created on: Jan 11, 2018
 *      Author: gu74req
 */

#ifndef ALIOTONOMEGACASCADE_H_
#define ALIOTONOMEGACASCADE_H_
#include "Rtypes.h"
#include "TVector3.h"
#include "TMath.h"

#include "AliFemtoDreamBasePart.h"
#include "AliFemtoDreamTrack.h"
#include "AliAODEvent.h"
#include "AliAODcascade.h"
#include "AliESDEvent.h"
#include "AliESDcascade.h"
#include "TDatabasePDG.h"

class AliOtonOmegaCascade : public AliFemtoDreamBasePart {
 public:
  AliOtonOmegaCascade();
  virtual ~AliOtonOmegaCascade();
  void SetCascade(AliAODEvent *evt, AliAODcascade *casc);
  void SetCascade(AliESDEvent *evt, AliMCEvent *mcEvent, AliESDcascade *casc);
  TString ClassName() {
    return "Cascade";
  }
  ;
  float GetXiMass() {
    return fXiMass;
  }
  ;
  float GetOmegaMass() {
    return fOmegaMass;
  }
  ;
  float GetMass(int PDGValue) const {
    if (TMath::Abs(PDGValue) == 3312) {
      return fXiMass;
    } else if (TMath::Abs(PDGValue) == 3334) {
      return fOmegaMass;
    } else {
      return 0.;
    }
  }
  ;
  float GetDCAXiPrimVtx() {
    return fDCAXiPrimVtx;
  }
  ;
  float GetXiDCADaug() {
    return fDCAXiDaug;
  }
  ;
  float GetXiTransverseRadius() {
    return fTransRadius;
  }
  ;
  float GetXiRapidity() {
    return fRapXi;
  }
  ;
  float GetOmegaRapidity() {
    return fRapOmega;
  }
  ;
  float GetXiAlpha() {
    return fAlphaXi;
  }
  ;
  float GetPtArmXi() {
    return fPtArmXi;
  }
  ;
  float GetDecayLength() {
    return fLength;
  }
  ;
  float BachDCAPrimVtx() {
    return fDCABachPrimVtx;
  }
  ;
  float Getv0Mass() {
    return fMassv0;
  }
  ;
  float Getv0DCADaug() {
    return fv0DCADaug;
  }
  ;
  float Getv0DCAPrimVtx() {
    return fv0DCAPrimVtx;
  }
  ;
  float Getv0TransverseRadius() {
    return fv0TransRadius;
  }
  ;
  float Getv0CPA() {
    return fv0CPA;
  }
  ;
  float Getv0PosToPrimVtx() {
    return fv0PosToPrimVtx;
  }
  ;
  float Getv0NegToPrimVtx() {
    return fv0NegToPrimVtx;
  }
  ;
  float Getv0XiPointingAngle() {
    return fv0ToXiPointAngle;
  }
  ;
  float Getv0DecayLength() {
    return fv0Length;
  }
  ;
  float GetDCAv0Xi() {
    return fDCAv0Xi;
  }
  ;
  void SetPDGDaugPos(int PDGCode) {
    fPosDaug->SetPDGCode(PDGCode);
  }
  ;
  void SetPDGDaugNeg(int PDGCode) {
    fNegDaug->SetPDGCode(PDGCode);
  }
  ;
  void SetPDGDaugBach(int PDGCode) {
    fBach->SetPDGCode(PDGCode);
  }
  ;
  TVector3 Getv0P() {
    return fv0Momentum;
  }
  ;
  float Getv0Pt() {
    return fv0Pt;
  }
  ;
  void Setv0PDGCode(int PDG) {
    fv0PDG = PDG;
  }
  ;
  double E2(int pdgCode, double Ptot2);
  double Ptot2V0();
  double Ptot2Casc();
  double MassCasc(int PDGv0, int PDGBach);
  double Massv0(int PDGPos, int PDGNeg);
  double MassXi();
  double MassOmega();
  double CosPointingAngle(double *primVtx, double *secondVtx);
  float DcaXiToPrimVertex(double *primVtx, double *secondVtx);
  AliFemtoDreamTrack *GetPosDaug() const {
    return fPosDaug;
  }
  ;
  AliFemtoDreamTrack *GetNegDaug() const {
    return fNegDaug;
  }
  ;
  AliFemtoDreamTrack *GetBach() const {
    return fBach;
  }
  ;
 private:
  AliOtonOmegaCascade &operator=(const AliOtonOmegaCascade &obj);
  AliOtonOmegaCascade(const AliOtonOmegaCascade&);
  void Reset();
  void SetMCMotherInfo(AliAODEvent *evt, AliAODcascade *casc);
  AliFemtoDreamTrack *fPosDaug;
  AliFemtoDreamTrack *fNegDaug;
  AliFemtoDreamTrack *fBach;
  float fXiMass;
  float fOmegaMass;
  float fDCAXiDaug;
  float fTransRadius;
  float fRapXi;
  float fRapOmega;
  float fAlphaXi;
  float fPtArmXi;
  float fDCAXiPrimVtx;
  float fLength;
  float fMassv0;
  float fv0DCADaug;
  float fv0DCAPrimVtx;
  TVector3 fv0Momentum;
  float fv0Pt;
  float fDCABachPrimVtx;
  float fv0TransRadius;
  float fv0CPA;
  float fv0PosToPrimVtx;
  float fv0NegToPrimVtx;
  float fv0ToXiPointAngle;
  float fv0Length;
  float fDCAv0Xi;
  int fv0PDG;ClassDef(AliOtonOmegaCascade,2)
};

inline double AliOtonOmegaCascade::E2(int pdgCode, double Ptot2) {
  double mass = TDatabasePDG::Instance()->GetParticle(pdgCode)->Mass();
  return mass * mass + Ptot2;
}

inline double AliOtonOmegaCascade::Ptot2V0() {
  return (fPosDaug->GetMomentum() + fNegDaug->GetMomentum()).Mag2();
}

inline double AliOtonOmegaCascade::Ptot2Casc() {
  return GetMomentum().Mag2();
}

inline double AliOtonOmegaCascade::MassCasc(int PDGv0, int PDGBach) {
  double Ev0 = ::sqrt(E2(PDGv0, Ptot2V0()));
  double EBach = ::sqrt(E2(PDGBach, fBach->GetMomentum().Mag2()));
  double Ptot2 = Ptot2Casc();
  return ::sqrt(::pow(Ev0 + EBach, 2) - Ptot2);
}

inline double AliOtonOmegaCascade::Massv0(int PDGPos, int PDGNeg) {
  double EPos = ::sqrt(E2(PDGPos, fPosDaug->GetMomentum().Mag2()));
  double ENeg = ::sqrt(E2(PDGNeg, fNegDaug->GetMomentum().Mag2()));
  return ::sqrt(
      (EPos + ENeg) * (EPos + ENeg)
          - (fPosDaug->GetMomentum() + fNegDaug->GetMomentum()).Mag2());
}

inline double AliOtonOmegaCascade::MassXi() {
  double Ev0 = ::sqrt(E2(3122, Ptot2V0()));
  double EBach = ::sqrt(E2(211, fBach->GetMomentum().Mag2()));
  double Ptot2 = Ptot2Casc();
  return ::sqrt(::pow(Ev0 + EBach, 2) - Ptot2);
}

inline double AliOtonOmegaCascade::MassOmega() {
  double Ev0 = TMath::Sqrt(E2(3122, Ptot2V0()));
  double EBach = TMath::Sqrt(E2(321, fBach->GetMomentum().Mag2()));
  double Ptot2 = Ptot2Casc();
  return ::sqrt(::pow(Ev0 + EBach, 2) - Ptot2);
}

inline double AliOtonOmegaCascade::CosPointingAngle(double *primVtx,
                                                     double *secondVtx) {
  TVector3 lVectPrimVtxToXi(secondVtx[0] - primVtx[0],
                            secondVtx[1] - primVtx[1],
                            secondVtx[2] - primVtx[2]);
  return GetMomentum().Dot(lVectPrimVtxToXi)
      / TMath::Sqrt(GetMomentum().Mag2() * lVectPrimVtxToXi.Mag2());
}

inline float AliOtonOmegaCascade::DcaXiToPrimVertex(double *primVtx,
                                                     double *secondVtx) {
  Double_t dx = (primVtx[1] - secondVtx[1]) * GetMomentum().Z()
      - (primVtx[2] - secondVtx[2]) * GetMomentum().Y();
  Double_t dy = (primVtx[2] - secondVtx[2]) * GetMomentum().X()
      - (primVtx[0] - secondVtx[0]) * GetMomentum().Z();
  Double_t dz = (primVtx[0] - secondVtx[0]) * GetMomentum().Y()
      - (primVtx[1] - secondVtx[1]) * GetMomentum().X();
  return TMath::Sqrt((dx * dx + dy * dy + dz * dz) / GetMomentum().Mag2());
}

#endif /* ALIOTONOMEGACASCADE_H_ */
