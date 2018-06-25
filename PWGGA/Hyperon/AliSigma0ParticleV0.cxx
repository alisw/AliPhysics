#include "AliSigma0ParticleV0.h"

ClassImp(AliSigma0ParticleV0)

    //____________________________________________________________________________________________________
    AliSigma0ParticleV0::AliSigma0ParticleV0()
    : AliSigma0ParticleBase(),
      fTrackLabelPos(-1),
      fTrackLabelNeg(-1),
      fMCLabelPos(-1),
      fMCLabelNeg(-1),
      fTrackPos(),
      fTrackNeg(),
      fCosAlpha(0),
      fRecMass(0),
      fPDGMass(0),
      fPID() {}
//____________________________________________________________________________________________________
AliSigma0ParticleV0::AliSigma0ParticleV0(AliESDv0 *v0, const AliESDtrack *pos,
                                         const AliESDtrack *neg,
                                         const AliESDVertex *vertex,
                                         const int pdg, const int pid,
                                         const float magneticField,
                                         AliMCEvent *mcEvent)
    : AliSigma0ParticleBase(),
      fTrackLabelPos(-1),
      fTrackLabelNeg(-1),
      fMCLabelPos(-1),
      fMCLabelNeg(-1),
      fTrackPos(),
      fTrackNeg(),
      fCosAlpha(0),
      fRecMass(0),
      fPDGMass(0),
      fPID() {
  fP[0] = v0->Px();
  fP[1] = v0->Py();
  fP[2] = v0->Pz();
  fPMC[0] = -1.;
  fPMC[1] = -1.;
  fPMC[2] = -1.;

  fPDGCode = pdg;
  fPt = v0->Pt();
  fTrackLabel = v0->GetLabel();
  fPhi = v0->Phi();
  fEta = v0->Eta();
  if (pid == 1)
    v0->ChangeMassHypothesis(3122);
  else
    v0->ChangeMassHypothesis(-3122);
  fMass = v0->GetEffMass();
  fRecMass = v0->GetEffMass();
  fPID = pid;

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

  fUse = true;
}

//____________________________________________________________________________________________________
AliSigma0ParticleV0::AliSigma0ParticleV0(const AliAODConversionPhoton &gamma,
                                         const AliVEvent *inputEvent)
    : AliSigma0ParticleBase(),
      fTrackLabelPos(-1),
      fTrackLabelNeg(-1),
      fMCLabelPos(-1),
      fMCLabelNeg(-1),
      fTrackPos(),
      fTrackNeg(),
      fCosAlpha(0),
      fRecMass(0),
      fPDGMass(0),
      fPID(0) {
  fP[0] = gamma.GetPx();
  fP[1] = gamma.GetPy();
  fP[2] = gamma.GetPz();
  fPMC[0] = -1.;
  fPMC[1] = -1.;
  fPMC[2] = -1.;

  fPDGCode = 22;
  fPt = gamma.GetPhotonPt();
  fTrackLabel = gamma.GetV0Index();
  fPhi = gamma.GetPhotonPhi();
  fEta = gamma.GetPhotonEta();
  fMass = gamma.GetPhotonMass();
  fRecMass = gamma.GetMass();

  fTrackLabelPos = gamma.GetTrackLabelPositive();
  fTrackLabelNeg = gamma.GetTrackLabelNegative();

  // commpute CPA
  double momV0[3] = {0, 0, 0};
  momV0[0] = gamma.Px();
  momV0[1] = gamma.Py();
  momV0[2] = gamma.Pz();

  double PosV0[3] = {
      gamma.GetConversionX() - inputEvent->GetPrimaryVertex()->GetX(),
      gamma.GetConversionY() - inputEvent->GetPrimaryVertex()->GetY(),
      gamma.GetConversionZ() -
          inputEvent->GetPrimaryVertex()
              ->GetZ()};  // Recalculated V0 Position vector

  double momV02 =
      momV0[0] * momV0[0] + momV0[1] * momV0[1] + momV0[2] * momV0[2];
  double PosV02 =
      PosV0[0] * PosV0[0] + PosV0[1] * PosV0[1] + PosV0[2] * PosV0[2];

  double cosinePointingAngle =
      (momV02 * PosV02 > 0.0)
          ? (PosV0[0] * momV0[0] + PosV0[1] * momV0[1] + PosV0[2] * momV0[2]) /
                std::sqrt(momV02 * PosV02)
          : -999.f;
  fCosAlpha = cosinePointingAngle;

  fMCLabelPos = gamma.GetMCLabelPositive();
  fMCLabelNeg = gamma.GetMCLabelNegative();

  //  fTrackPos
  //  fTrackNeg

  fUse = true;
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
