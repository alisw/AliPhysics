#include "AliSigma0ParticleV0.h"
#include "TMath.h"

ClassImp(AliSigma0ParticleV0)

    //____________________________________________________________________________________________________
    AliSigma0ParticleV0::AliSigma0ParticleV0()
    : AliSigma0ParticleBase(),
      fTrackLabelPos(-1),
      fTrackLabelNeg(-1),
      fMCLabelPos(-1),
      fMCLabelNeg(-1),
      fMCLabelV0(-1),
      fTrackPos(),
      fTrackNeg(),
      fCosAlpha(0),
      fRecMass(0),
      fPDGMass(0) {}
//____________________________________________________________________________________________________
AliSigma0ParticleV0::AliSigma0ParticleV0(AliESDv0 *v0, const AliESDtrack *pos,
                                         const AliESDtrack *neg,
                                         const AliESDVertex *vertex,
                                         const int pdg,
                                         const float magneticField,
                                         AliMCEvent *mcEvent)
    : AliSigma0ParticleBase(),
      fTrackLabelPos(-1),
      fTrackLabelNeg(-1),
      fMCLabelPos(-1),
      fMCLabelNeg(-1),
      fMCLabelV0(-1),
      fTrackPos(),
      fTrackNeg(),
      fCosAlpha(0),
      fRecMass(0),
      fPDGMass(0) {
  fP[0] = v0->Px();
  fP[1] = v0->Py();
  fP[2] = v0->Pz();

  fPDGCode = pdg;
  fPt = v0->Pt();
  fTrackLabel = v0->GetLabel();
  fPhi = v0->Phi();
  fEta = v0->Eta();
  if (pdg != 22) {
    v0->ChangeMassHypothesis(pdg);
    fMass = v0->GetEffMass();
    fRecMass = v0->GetEffMass();
  }

  fTrackLabelPos = pos->GetID();
  fTrackLabelNeg = neg->GetID();

  Double_t xPV = vertex->GetX();
  Double_t yPV = vertex->GetY();
  Double_t zPV = vertex->GetZ();
  fCosAlpha = v0->GetV0CosineOfPointingAngle(xPV, yPV, zPV);

  AliSigma0ParticleBase fCandidatePos(pos, 0, magneticField);
  AliSigma0ParticleBase fCandidateNeg(neg, 0, magneticField);
  if (mcEvent) {
    AliMCParticle *mcParticlePos =
        static_cast<AliMCParticle *>(mcEvent->GetTrack(pos->GetLabel()));
    if (mcParticlePos) fCandidatePos.ProcessMCInfo(mcParticlePos, mcEvent);
    fMCLabelPos = pos->GetLabel();
    AliMCParticle *mcParticleNeg =
        static_cast<AliMCParticle *>(mcEvent->GetTrack(neg->GetLabel()));
    if (mcParticleNeg) fCandidateNeg.ProcessMCInfo(mcParticleNeg, mcEvent);
    fMCLabelNeg = neg->GetLabel();
  }

  fTrackPos = fCandidatePos;
  fTrackNeg = fCandidateNeg;

  AliSigma0ParticleBase::SetUse(true);
}

//____________________________________________________________________________________________________
AliSigma0ParticleV0::AliSigma0ParticleV0(const AliAODConversionPhoton *gamma,
                                         const AliESDtrack *pos,
                                         const AliESDtrack *neg,
                                         const AliVEvent *inputEvent)
    : AliSigma0ParticleBase(),
      fTrackLabelPos(-1),
      fTrackLabelNeg(-1),
      fMCLabelPos(-1),
      fMCLabelNeg(-1),
      fMCLabelV0(-1),
      fTrackPos(),
      fTrackNeg(),
      fCosAlpha(0),
      fRecMass(0),
      fPDGMass(0) {
  AliSigma0ParticleBase::SetUse(true);

  fP[0] = gamma->GetPx();
  fP[1] = gamma->GetPy();
  fP[2] = gamma->GetPz();

  fPDGCode = 22;
  fPt = gamma->GetPhotonPt();
  fTrackLabel = gamma->GetV0Index();
  fPhi = gamma->GetPhotonPhi();
  fEta = gamma->GetPhotonEta();
  fMass = gamma->GetPhotonMass();
  fRecMass = gamma->M();

  fTrackLabelPos = gamma->GetTrackLabelPositive();
  fTrackLabelNeg = gamma->GetTrackLabelNegative();

  // commpute CPA
  double momV0[3] = {0, 0, 0};
  momV0[0] = gamma->Px();
  momV0[1] = gamma->Py();
  momV0[2] = gamma->Pz();

  double PosV0[3] = {
      gamma->GetConversionX() - inputEvent->GetPrimaryVertex()->GetX(),
      gamma->GetConversionY() - inputEvent->GetPrimaryVertex()->GetY(),
      gamma->GetConversionZ() -
          inputEvent->GetPrimaryVertex()
              ->GetZ()};  // Recalculated V0 Position vector

  double momV02 =
      momV0[0] * momV0[0] + momV0[1] * momV0[1] + momV0[2] * momV0[2];
  double PosV02 =
      PosV0[0] * PosV0[0] + PosV0[1] * PosV0[1] + PosV0[2] * PosV0[2];

  double cosinePointingAngle =
      (momV02 * PosV02 > 0.0)
          ? (PosV0[0] * momV0[0] + PosV0[1] * momV0[1] + PosV0[2] * momV0[2]) /
                TMath::Sqrt(momV02 * PosV02)
          : -999.f;
  fCosAlpha = cosinePointingAngle;

  fMCLabelPos = gamma->GetMCLabelPositive();
  fMCLabelNeg = gamma->GetMCLabelNegative();

  AliSigma0ParticleBase fCandidatePos(pos, 0, inputEvent->GetMagneticField());
  AliSigma0ParticleBase fCandidateNeg(neg, 0, inputEvent->GetMagneticField());
  fTrackPos = fCandidatePos;
  fTrackNeg = fCandidateNeg;

  AliSigma0ParticleBase::SetUse(true);
}

//____________________________________________________________________________________________________
AliSigma0ParticleV0 &AliSigma0ParticleV0::operator=(
    const AliSigma0ParticleV0 &obj) {
  //  Assignment operator
  if (this == &obj) return *this;

  fP[0] = obj.GetPx();
  fP[1] = obj.GetPy();
  fP[2] = obj.GetPz();
  fPMC[0] = obj.GetPxMC();
  fPMC[1] = obj.GetPyMC();
  fPMC[2] = obj.GetPzMC();

  fPDGCode = obj.GetPDGcode();
  fMass = obj.GetMass();
  fQ = obj.GetQ();
  fPt = obj.GetPt();
  fTrackLabel = obj.GetTrackLabel();
  fPhi = obj.GetPhi();
  fEta = obj.GetEta();

  fUse = obj.GetIsUse();

  fTrackLabelPos = obj.GetTrackLabelPos();
  fTrackLabelNeg = obj.GetTrackLabelNeg();
  fCosAlpha = obj.GetCosineAlpha();

  fMCLabelPos = obj.GetMCLabelPos();
  fMCLabelNeg = obj.GetMCLabelNeg();
  fMCLabelV0 = obj.GetMCLabelV0();
  fTrackPos = obj.GetPosDaughter();
  fTrackNeg = obj.GetNegDaughter();

  fRecMass = obj.GetRecMass();
  fPDGMass = obj.GetPDGMass();

  return (*this);
}

//____________________________________________________________________________________________________
float AliSigma0ParticleV0::GetArmenterosAlpha() const {
  TVector3 posTrack(fTrackPos.GetPx(), fTrackPos.GetPy(), fTrackPos.GetPz());
  TVector3 negTrack(fTrackNeg.GetPx(), fTrackNeg.GetPy(), fTrackNeg.GetPz());
  TVector3 lambda(GetPx(), GetPy(), GetPz());
  return 1. - 2. / (1. + posTrack.Dot(lambda) / negTrack.Dot(lambda));
}

//____________________________________________________________________________________________________
float AliSigma0ParticleV0::GetArmenterosQt() const {
  // Transverse momentum of lambda w.r.t. to total momentum
  TVector3 lambdaP(fTrackPos.GetPx(), fTrackPos.GetPy(), fTrackPos.GetPz());
  TVector3 sigmaP(GetPx(), GetPy(), GetPz());
  return lambdaP.Perp(sigmaP);
}

//____________________________________________________________________________________________________
int AliSigma0ParticleV0::MatchToMC(const AliMCEvent *mcEvent,
                                   const int PIDmother,
                                   const std::vector<int> PIDdaughters) {
  // Adopted from the function in AliAODRecoDecay

  std::vector<int> trackLabels = {{fMCLabelPos, fMCLabelNeg}};
  std::vector<int> labMom = {{0, 0}};
  std::vector<bool> pdgUsed = {{false, false}};
  int lab, labMother, pdgMother, pdgPart;
  AliMCParticle *part = nullptr;
  AliMCParticle *mother = nullptr;
  double pxSumDgs = 0., pySumDgs = 0., pzSumDgs = 0.;

  // loop on daughter labels
  for (int i = 0; i < 2; ++i) {
    labMom[i] = -1;
    lab = TMath::Abs(trackLabels[i]);
    if (lab < 0) {
      return -1;
    }
    part = static_cast<AliMCParticle *>(mcEvent->GetTrack(lab));
    if (!part) {
      return -1;
    }

    // check the PDG of the daughter, if requested
    pdgPart = part->PdgCode();
    for (int j = 0; j < 2; ++j) {
      if (!pdgUsed[j] && pdgPart == PIDdaughters[j]) {
        pdgUsed[j] = true;
        break;
      }
    }

    mother = part;
    while (mother->GetMother() >= 0) {
      labMother = mother->GetMother();
      mother = static_cast<AliMCParticle *>(mcEvent->GetTrack(labMother));
      if (!mother) {
        break;
      }
      pdgMother = mother->PdgCode();
      if (pdgMother == PIDmother) {
        labMom[i] = labMother;
        // keep sum of daughters' momenta, to check for mom conservation
        pxSumDgs += part->Px();
        pySumDgs += part->Py();
        pzSumDgs += part->Pz();
        break;
      } else if (pdgMother > PIDmother || pdgMother < 10) {
        break;
      }
    }
    if (labMom[i] == -1) {
      return -1;
    }
  }

  // check if the candidate is signal
  labMother = labMom[0];
  // all labels have to be the same and !=-1
  for (int i = 0; i < 2; ++i) {
    if (labMom[i] == -1) return -1;
    if (labMom[i] != labMother) return -1;
  }

  // check that all daughter PDGs are matched
  for (int i = 0; i < 2; ++i) {
    if (pdgUsed[i] == false) return -1;
  }

  // the above works only for non-resonant decays,
  // it's better to check for mom conservation
  mother = static_cast<AliMCParticle *>(mcEvent->GetTrack(labMother));
  Double_t pxMother = mother->Px();
  Double_t pyMother = mother->Py();
  Double_t pzMother = mother->Pz();
  // within 0.1%
  if ((TMath::Abs(pxMother - pxSumDgs) / (TMath::Abs(pxMother) + 1.e-13)) >
          0.00001 &&
      (TMath::Abs(pyMother - pySumDgs) / (TMath::Abs(pyMother) + 1.e-13)) >
          0.00001 &&
      (TMath::Abs(pzMother - pzSumDgs) / (TMath::Abs(pzMother) + 1.e-13)) >
          0.00001) {
    return -1;
  }
  fPMC[0] = pxMother;
  fPMC[1] = pyMother;
  fPMC[2] = pzMother;

  fMCLabelV0 = labMother;
  return labMother;
}
