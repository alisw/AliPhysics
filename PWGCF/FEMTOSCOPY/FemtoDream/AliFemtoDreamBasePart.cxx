/*
 * AliFemtoPPbpbLamBaseParticle.cxx
 *
 *  Created on: Aug 29, 2017
 *      Author: gu74req
 */

#include "AliFemtoDreamBasePart.h"
#include "AliMCParticle.h"
#include "AliSigma0ParticleV0.h"

#include <iostream>
ClassImp(AliFemtoDreamBasePart) AliFemtoDreamBasePart::AliFemtoDreamBasePart()
    : fIsReset(false),
      fGTI(0),
      fVGTI(0),
      fTrackBufferSize(0),
      fP(),
      fMCP(),
      fPt(0),
      fMCPt(0),
      fP_TPC(0),
      fEta(0),
      fTheta(0),
      fMCTheta(0),
      fPhi(0),
      fPhiAtRadius(),
      fXYZAtRadius(),
      fMCPhi(0),
      fIDTracks(0),
      fCharge(0),
      fCPA(0),
      fInvMass(0),
      fOrigin(kUnknown),
      fPDGCode(0),
      fMCPDGCode(0),
      fPDGMotherWeak(0),
      fMotherID(0),
      fID(0),
      fMotherPDG(-1),
      fEvtNumber(0),
      fIsMC(false),
      fUse(true),
      fIsSet(true),
      fEvtMultiplicity(-1) {}

AliFemtoDreamBasePart::AliFemtoDreamBasePart(const AliFemtoDreamBasePart &part)
    : fIsReset(part.fIsReset),
      fGTI(part.fGTI),
      fVGTI(part.fVGTI),
      fTrackBufferSize(part.fTrackBufferSize),
      fP(part.fP),
      fMCP(part.fMCP),
      fPt(part.fPt),
      fMCPt(part.fMCPt),
      fP_TPC(part.fP_TPC),
      fEta(part.fEta),
      fTheta(part.fTheta),
      fMCTheta(part.fMCTheta),
      fPhi(part.fPhi),
      fPhiAtRadius(part.fPhiAtRadius),
      fXYZAtRadius(part.fXYZAtRadius),
      fMCPhi(part.fMCPhi),
      fIDTracks(part.fIDTracks),
      fCharge(part.fCharge),
      fCPA(part.fCPA),
      fInvMass(part.fInvMass),
      fOrigin(part.fOrigin),
      fPDGCode(part.fPDGCode),
      fMCPDGCode(part.fMCPDGCode),
      fPDGMotherWeak(part.fPDGMotherWeak),
      fMotherID(part.fMotherID),
      fID(part.fID),
      fMotherPDG(part.fMotherPDG),
      fEvtNumber(part.fEvtNumber),
      fIsMC(part.fIsMC),
      fUse(part.fUse),
      fIsSet(part.fIsSet),
      fEvtMultiplicity(part.fEvtMultiplicity) {}

AliFemtoDreamBasePart &AliFemtoDreamBasePart::operator=(
    const AliFemtoDreamBasePart &obj) {
  if (this == &obj) {
    return *this;
  }
  fIsReset = obj.fIsReset;
  fGTI = 0;
  fVGTI = 0;
  fTrackBufferSize = obj.fTrackBufferSize;
  fP = obj.fP;
  fPt = obj.fPt;
  fMCP = obj.fMCP;
  fMCPt = obj.fMCPt;
  fP_TPC = obj.fP_TPC;
  fEta = obj.fEta;
  fTheta = obj.fTheta;
  fMCTheta = obj.fMCTheta;
  fPhi = obj.fPhi;
  fPhiAtRadius = obj.fPhiAtRadius;
  fXYZAtRadius = obj.fXYZAtRadius;
  fMCPhi = obj.fMCPhi;
  fIDTracks = obj.fIDTracks;
  fCharge = obj.fCharge;
  fCPA = obj.fCPA;
  fInvMass = obj.fInvMass;
  fOrigin = obj.fOrigin;
  fPDGCode = obj.fPDGCode;
  fMCPDGCode = obj.fMCPDGCode;
  fPDGMotherWeak = obj.fPDGMotherWeak;
  fMotherID = obj.fMotherID;
  fID = obj.fID;
  fMotherPDG = obj.fMotherPDG;
  fEvtNumber = obj.fEvtNumber;
  fIsMC = obj.fIsMC;
  fUse = obj.fUse;
  fIsSet = obj.fIsSet;
  fEvtMultiplicity = obj.fEvtMultiplicity;
  return (*this);
}

AliFemtoDreamBasePart::AliFemtoDreamBasePart(
    const AliSigma0ParticlePhotonMother &mother, const AliMCEvent *mcEvent)
    : fIsReset(false),
      fGTI(0),
      fVGTI(0),
      fTrackBufferSize(0),
      fP(TVector3(mother.GetPx(), mother.GetPy(), mother.GetPz())),
      fMCP(TVector3(mother.GetPxMC(), mother.GetPyMC(), mother.GetPzMC())),
      fPt(mother.GetPt()),
      fMCPt(mother.GetPtMC()),
      fP_TPC(0),
      fEta(),
      fTheta(),
      fMCTheta(),
      fPhi(),
      fPhiAtRadius(0),
      fXYZAtRadius(0),
      fMCPhi(),
      fIDTracks(),
      fCharge(0),
      fCPA(0),
      fInvMass(0),
      fOrigin(kUnknown),
      fPDGCode(mother.GetPDGCode()),
      fMCPDGCode(mother.GetPDGCode()),
      fPDGMotherWeak(0),
      fMotherID(),
      fID(0),
      fMotherPDG(0),
      fEvtNumber(0),
      fIsMC((mother.GetMCLabel() > 0)),
      fUse(true),
      fIsSet(true),
      fEvtMultiplicity(-1) {
  fEta.push_back(mother.GetEta());
  fEta.push_back(mother.GetV0().GetPosDaughter().GetEta());
  fEta.push_back(mother.GetPhoton().GetPosDaughter().GetEta());
  fEta.push_back(mother.GetV0().GetNegDaughter().GetEta());
  fEta.push_back(mother.GetPhoton().GetNegDaughter().GetEta());

  fTheta.push_back(mother.GetTheta());
  fTheta.push_back(mother.GetV0().GetPosDaughter().GetTheta());
  fTheta.push_back(mother.GetPhoton().GetPosDaughter().GetTheta());
  fTheta.push_back(mother.GetV0().GetNegDaughter().GetTheta());
  fTheta.push_back(mother.GetPhoton().GetNegDaughter().GetTheta());

  fMCTheta.push_back(mother.GetThetaMC());
  fMCTheta.push_back(mother.GetV0().GetPosDaughter().GetThetaMC());
  fMCTheta.push_back(mother.GetPhoton().GetPosDaughter().GetThetaMC());
  fMCTheta.push_back(mother.GetV0().GetNegDaughter().GetThetaMC());
  fMCTheta.push_back(mother.GetPhoton().GetNegDaughter().GetThetaMC());

  fPhiAtRadius.push_back(mother.GetV0().GetPosDaughter().GetPhiStar());
  fPhiAtRadius.push_back(mother.GetPhoton().GetPosDaughter().GetPhiStar());
  fPhiAtRadius.push_back(mother.GetV0().GetNegDaughter().GetPhiStar());
  fPhiAtRadius.push_back(mother.GetPhoton().GetNegDaughter().GetPhiStar());

  fPhi.push_back(mother.GetPhi());
  fPhi.push_back(mother.GetV0().GetPosDaughter().GetPhi());
  fPhi.push_back(mother.GetPhoton().GetPosDaughter().GetPhi());
  fPhi.push_back(mother.GetV0().GetNegDaughter().GetPhi());
  fPhi.push_back(mother.GetPhoton().GetNegDaughter().GetPhi());

  fMCPhi.push_back(mother.GetPhiMC());
  fMCPhi.push_back(mother.GetV0().GetPosDaughter().GetPhiMC());
  fMCPhi.push_back(mother.GetPhoton().GetPosDaughter().GetPhiMC());
  fMCPhi.push_back(mother.GetV0().GetNegDaughter().GetPhiMC());
  fMCPhi.push_back(mother.GetPhoton().GetNegDaughter().GetPhiMC());

  fIDTracks.push_back(mother.GetV0().GetTrackLabelPos());
  fIDTracks.push_back(mother.GetV0().GetTrackLabelNeg());
  fIDTracks.push_back(mother.GetPhoton().GetTrackLabelPos());
  fIDTracks.push_back(mother.GetPhoton().GetTrackLabelNeg());

  if (mcEvent) {
    auto *partSigma0 =
        static_cast<AliMCParticle *>(mcEvent->GetTrack(mother.GetMCLabel()));
    if (!partSigma0) {
      // if the particle does not exist it is most probably combinatorial
      // background
      SetParticleOrigin(AliFemtoDreamBasePart::kFake);
      fIsSet = false;
      fUse = false;
    } else {
      auto *sigmaMother = static_cast<AliMCParticle *>(
          mcEvent->GetTrack(partSigma0->GetMother()));
      if (sigmaMother) {
        SetPDGMotherWeak(sigmaMother->PdgCode());
        fMotherID = partSigma0->GetMother();
      }
      if (partSigma0->IsPhysicalPrimary() &&
          !(partSigma0->IsSecondaryFromWeakDecay())) {
        SetParticleOrigin(AliFemtoDreamBasePart::kPhysPrimary);
      } else if (partSigma0->IsSecondaryFromWeakDecay() &&
                 !(partSigma0->IsSecondaryFromMaterial())) {
        SetParticleOrigin(AliFemtoDreamBasePart::kWeak);
      } else if (partSigma0->IsSecondaryFromMaterial()) {
        SetParticleOrigin(AliFemtoDreamBasePart::kMaterial);
      } else {
        SetParticleOrigin(AliFemtoDreamBasePart::kUnknown);
      }
    }
  }
}

AliFemtoDreamBasePart::AliFemtoDreamBasePart(
    const AliSigma0ParticleV0 &daughter, const AliMCEvent *mcEvent)
    : fIsReset(false),
      fGTI(0),
      fVGTI(0),
      fTrackBufferSize(0),
      fP(TVector3(daughter.GetPx(), daughter.GetPy(), daughter.GetPz())),
      fMCP(
          TVector3(daughter.GetPxMC(), daughter.GetPyMC(), daughter.GetPzMC())),
      fPt(daughter.GetPt()),
      fMCPt(daughter.GetPtMC()),
      fP_TPC(0),
      fEta(),
      fTheta(),
      fMCTheta(),
      fPhi(),
      fPhiAtRadius(0),
      fXYZAtRadius(0),
      fMCPhi(),
      fIDTracks(),
      fCharge(0),
      fCPA(0),
      fInvMass(0),
      fOrigin(kUnknown),
      fPDGCode(daughter.GetPDGcode()),
      fMCPDGCode(daughter.GetPDGcode()),
      fPDGMotherWeak(0),
      fMotherID(-1),
      fID(0),
      fMotherPDG(0),
      fEvtNumber(0),
      fIsMC(-1),
      fUse(true),
      fIsSet(true),
      fEvtMultiplicity(-1) {
  fEta.push_back(daughter.GetPosDaughter().GetEta());
  fEta.push_back(daughter.GetNegDaughter().GetEta());

  fTheta.push_back(daughter.GetPosDaughter().GetTheta());
  fTheta.push_back(daughter.GetNegDaughter().GetTheta());

  fMCTheta.push_back(daughter.GetPosDaughter().GetThetaMC());
  fMCTheta.push_back(daughter.GetNegDaughter().GetThetaMC());

  fPhiAtRadius.push_back(daughter.GetPosDaughter().GetPhiStar());
  fPhiAtRadius.push_back(daughter.GetNegDaughter().GetPhiStar());

  fPhi.push_back(daughter.GetPosDaughter().GetPhi());
  fPhi.push_back(daughter.GetNegDaughter().GetPhi());

  fCharge.push_back(1);
  fCharge.push_back(-1);

  fMCPhi.push_back(daughter.GetPosDaughter().GetPhiMC());
  fMCPhi.push_back(daughter.GetNegDaughter().GetPhiMC());

  fIDTracks.push_back(daughter.GetTrackLabelPos());
  fIDTracks.push_back(daughter.GetTrackLabelNeg());

  if (mcEvent) {
    auto *partv0 = static_cast<AliMCParticle *>(
        mcEvent->GetTrack(daughter.GetMCLabelV0()));
    if (!partv0) {
      // if the particle does not exist it is most probably combinatorial
      // background
      SetParticleOrigin(AliFemtoDreamBasePart::kFake);
      fIsSet = false;
      fUse = false;
    } else {
      auto *v0Mother =
          static_cast<AliMCParticle *>(mcEvent->GetTrack(partv0->GetMother()));
      if (v0Mother) {
        SetPDGMotherWeak(v0Mother->PdgCode());
        fMotherID = partv0->GetMother();
      }
      if (partv0->IsPhysicalPrimary() &&
          !(partv0->IsSecondaryFromWeakDecay())) {
        SetParticleOrigin(AliFemtoDreamBasePart::kPhysPrimary);
      } else if (partv0->IsSecondaryFromWeakDecay() &&
                 !(partv0->IsSecondaryFromMaterial())) {
        SetParticleOrigin(AliFemtoDreamBasePart::kWeak);
      } else if (partv0->IsSecondaryFromMaterial()) {
        SetParticleOrigin(AliFemtoDreamBasePart::kMaterial);
      } else {
        SetParticleOrigin(AliFemtoDreamBasePart::kUnknown);
      }
    }
  }
}

AliFemtoDreamBasePart::AliFemtoDreamBasePart(
    const AliAODConversionPhoton *gamma, const AliVTrack *posTrack,
    const AliVTrack *negTrack, const AliVEvent *inputEvent)
    : fIsReset(false),
      fGTI(0),
      fVGTI(0),
      fTrackBufferSize(0),
      fP(TVector3(gamma->GetPx(), gamma->GetPy(), gamma->GetPz())),
      fMCP(),
      fPt(gamma->GetPhotonPt()),
      fMCPt(0),
      fP_TPC(0),
      fEta(),
      fTheta(),
      fMCTheta(),
      fPhi(),
      fPhiAtRadius(0),
      fXYZAtRadius(0),
      fMCPhi(),
      fIDTracks(),
      fCharge(0),
      fCPA(0),
      fInvMass(gamma->GetPhotonMass()),
      fOrigin(kUnknown),
      fPDGCode(),
      fMCPDGCode(),
      fPDGMotherWeak(0),
      fMotherID(-1),
      fID(0),
      fMotherPDG(0),
      fEvtNumber(0),
      fIsMC(-1),
      fUse(true),
      fIsSet(true),
      fEvtMultiplicity(-1) {
  double momV0[3] = {0, 0, 0};
  momV0[0] = gamma->Px();
  momV0[1] = gamma->Py();
  momV0[2] = gamma->Pz();

  // Recalculated V0 Position vector
  double PosV0[3] = {
      gamma->GetConversionX() - inputEvent->GetPrimaryVertex()->GetX(),
      gamma->GetConversionY() - inputEvent->GetPrimaryVertex()->GetY(),
      gamma->GetConversionZ() - inputEvent->GetPrimaryVertex()->GetZ()};

  double momV02 = fP[0] * fP[0] + fP[1] * fP[1] + fP[2] * fP[2];
  double PosV02 =
      PosV0[0] * PosV0[0] + PosV0[1] * PosV0[1] + PosV0[2] * PosV0[2];

  double cosinePointingAngle =
      (momV02 * PosV02 > 0.0)
          ? (PosV0[0] * momV0[0] + PosV0[1] * momV0[1] + PosV0[2] * momV0[2]) /
                TMath::Sqrt(momV02 * PosV02)
          : -999.f;
  fCPA = cosinePointingAngle;

  const int posLabel = gamma->GetTrackLabelPositive();
  const int negLabel = gamma->GetTrackLabelNegative();
  fIDTracks.push_back(posLabel);
  fIDTracks.push_back(negLabel);

  fEta.push_back(gamma->GetPhotonEta());
  fEta.push_back(posTrack->Eta());
  fEta.push_back(negTrack->Eta());

  fTheta.push_back(gamma->GetPhotonTheta());
  fTheta.push_back(posTrack->Theta());
  fTheta.push_back(negTrack->Theta());

  fPhi.push_back(gamma->GetPhotonPhi());
  fPhi.push_back(posTrack->Phi());
  fPhi.push_back(negTrack->Phi());

  std::vector<float> phiAtRadiiPos;
  std::vector<float> phiAtRadiiNeg;
  PhiAtRadii(posTrack, inputEvent->GetMagneticField(), phiAtRadiiPos);
  PhiAtRadii(negTrack, inputEvent->GetMagneticField(), phiAtRadiiNeg);
  fPhiAtRadius.push_back(phiAtRadiiPos);
  fPhiAtRadius.push_back(phiAtRadiiNeg);

  fCharge.push_back(posTrack->Charge() + negTrack->Charge());
  fCharge.push_back(posTrack->Charge());
  fCharge.push_back(negTrack->Charge());
}

AliFemtoDreamBasePart::~AliFemtoDreamBasePart() {
  fGTI= nullptr;
  fVGTI = nullptr;
  fEta.clear();
  fTheta.clear();
  fMCTheta.clear();
  fPhi.clear();
  fPhiAtRadius.clear();
  fXYZAtRadius.clear();
  fMCPhi.clear();
  fIDTracks.clear();
  fCharge.clear();
}

void AliFemtoDreamBasePart::SetMCParticle(AliAODMCParticle *mcPart,
                                          AliMCEvent *evt) {
  this->SetMCPt(mcPart->Pt());
  this->SetMCMomentum(mcPart->Px(), mcPart->Py(), mcPart->Pz());
  this->SetMCTheta(mcPart->Theta());
  this->SetMCPhi(mcPart->Phi());
  this->fIsSet = true;
  this->SetUse(true);
  this->fIsReset = false;
  // find the original mother in the mc stack.
  int motherID = mcPart->GetMother();
  int lastMother = motherID;
  while (motherID != -1) {
    lastMother = motherID;
    motherID = evt->GetTrack(motherID)->GetMother();
  }
  this->SetMotherID(lastMother);
}

void AliFemtoDreamBasePart::ResetMCInfo() {
  this->SetPt(0);
  // a change
  this->SetEta(0);
  this->SetTheta(0);
  this->SetPhi(0);
  this->SetCharge(0);
  this->SetPDGCode(0);
  this->SetMotherID(0);
  this->fIsSet = false;
  this->SetUse(false);
  this->fIsReset = true;
}

void AliFemtoDreamBasePart::PhiAtRadii(const AliVTrack *track,
                                       const float bfield,
                                       std::vector<float> &tmpVec) {
  float TPCradii[9] = {85., 105., 125., 145., 165., 185., 205., 225., 245.};
  float phi0 = track->Phi();
  ;
  float pt = track->Pt();
  float chg = track->Charge();
  for (int radius = 0; radius < 9; radius++) {
    tmpVec.push_back(phi0 - TMath::ASin(0.1 * chg * bfield * 0.3 *
                                        TPCradii[radius] * 0.01 / (2. * pt)));
  }
}
