#include "AliSigma0ParticlePhotonMother.h"

ClassImp(AliSigma0ParticlePhotonMother)

    //____________________________________________________________________________________________________
    AliSigma0ParticlePhotonMother::AliSigma0ParticlePhotonMother()
    : AliSigma0ParticleBase(),
      fType(-1),
      fRecMass(0),
      fV0(),
      fPhoton(),
      fPhoton2() {}

//____________________________________________________________________________________________________
AliSigma0ParticlePhotonMother::AliSigma0ParticlePhotonMother(
    const AliSigma0ParticleV0 &lambdaCandidate,
    const AliAODConversionPhoton &photonCandidate, const AliVEvent *inputEvent)
    : AliSigma0ParticleBase(),
      fType(1),
      fRecMass(0),
      fV0(),
      fPhoton(),
      fPhoton2() {
  TLorentzVector track1, track2;
  track1.SetXYZM(lambdaCandidate.GetPx(), lambdaCandidate.GetPy(),
                 lambdaCandidate.GetPz(), lambdaCandidate.GetRecMass());
  track2.SetXYZM(photonCandidate.GetPx(), photonCandidate.GetPy(),
                 photonCandidate.GetPz(), photonCandidate.M());
  TLorentzVector trackSum = track1 + track2;

  fP[0] = trackSum.Px();
  fP[1] = trackSum.Py();
  fP[2] = trackSum.Pz();
  fPMC[0] = -1.;
  fPMC[1] = -1.;
  fPMC[2] = -1.;

  //  fPDGCode = pdg;
  fPt = std::sqrt(fP[0] * fP[0] + fP[1] * fP[1]);
  //  fTrackLabel = v0.GetID();
  fPhi = trackSum.Phi();
  fEta = trackSum.Eta();
  fRecMass = trackSum.M();
  // see https://arxiv.org/pdf/1703.04639.pdf
  fMass = trackSum.M() - lambdaCandidate.GetRecMass() +
          lambdaCandidate.GetPDGMass();

  fUse = true;

  fV0 = lambdaCandidate;
  AliSigma0ParticleV0 phot(photonCandidate, inputEvent);
  fPhoton = phot;
}

//____________________________________________________________________________________________________
AliSigma0ParticlePhotonMother &AliSigma0ParticlePhotonMother::operator=(
    const AliSigma0ParticlePhotonMother &obj) {
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

  fType = obj.fType;

  fUse = obj.GetIsUse();

  fRecMass = obj.GetRecMass();
  fV0 = obj.GetV0();
  fPhoton = obj.GetPhoton();
  fPhoton2 = obj.GetPhoton2();

  return (*this);
}

//____________________________________________________________________________________________________
bool AliSigma0ParticlePhotonMother::IsTrueSigma(AliMCEvent *mcEvent) const {
  // get photon mother
  AliMCParticle *photonElePos =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(fPhoton.GetMCLabelPos()));
  AliMCParticle *photonEleNeg =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(fPhoton.GetMCLabelNeg()));
  if (photonElePos == nullptr || photonEleNeg == nullptr) return false;
  if (photonElePos->GetMother() != photonEleNeg->GetMother()) return false;
  AliMCParticle *photonMC = static_cast<AliMCParticle *>(
      mcEvent->GetTrack(photonEleNeg->GetMother()));
  if (photonMC == nullptr) return false;
  AliMCParticle *photonMother =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(photonMC->GetMother()));
  if (photonMother == nullptr) return false;

  // get lambda mother
  AliMCParticle *lambdaPos =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(fV0.GetMCLabelPos()));
  AliMCParticle *lambdaNeg =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(fV0.GetMCLabelNeg()));
  if (lambdaPos == nullptr || lambdaNeg == nullptr) return false;
  if (!(std::abs(lambdaPos->GetMother()) == 211 &&
        std::abs(lambdaNeg->GetMother()) == 2212) ||
      !(std::abs(lambdaPos->GetMother()) == 2211 &&
        std::abs(lambdaNeg->GetMother()) == 211))
    return false;
  AliMCParticle *lambdaMC =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(lambdaPos->GetMother()));
  if (lambdaMC == nullptr) return false;
  AliMCParticle *lambdaMother =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(lambdaMC->GetMother()));
  if (lambdaMother == nullptr) return false;

  if (photonMother->PdgCode() != lambdaMother->PdgCode()) return false;
  if (std::abs(photonMother->PdgCode()) != 3212) return false;
  return true;
}

//____________________________________________________________________________________________________
float AliSigma0ParticlePhotonMother::GetArmenterosAlpha() const {
  TVector3 daughter;
  if (fType == 1)
    daughter.SetXYZ(fPhoton.GetPx(), fPhoton.GetPy(), fPhoton.GetPz());
  if (fType == 2) {
    TLorentzVector track1, track2;
    track1.SetXYZM(fPhoton.GetPx(), fPhoton.GetPy(), fPhoton.GetPz(),
                   fPhoton.GetMass());
    track2.SetXYZM(fPhoton2.GetPx(), fPhoton2.GetPy(), fPhoton2.GetPz(),
                   fPhoton2.GetMass());
    TLorentzVector trackSum = track1 + track2;
    daughter.SetXYZ(trackSum.Px(), trackSum.Py(), trackSum.Pz());
  }

  TVector3 lambdaP(fV0.GetPx(), fV0.GetPy(), fV0.GetPz());
  TVector3 sigmaP(GetPx(), GetPy(), GetPz());
  return 1. - 2. / (1. + daughter.Dot(sigmaP) / lambdaP.Dot(sigmaP));
}

//____________________________________________________________________________________________________
float AliSigma0ParticlePhotonMother::GetArmenterosQt() const {
  // Transverse momentum of lambda w.r.t. to total momentum
  TVector3 lambdaP(fV0.GetPx(), fV0.GetPy(), fV0.GetPz());
  TVector3 sigmaP(GetPx(), GetPy(), GetPz());
  return lambdaP.Perp(sigmaP);
}

//____________________________________________________________________________________________________
double AliSigma0ParticlePhotonMother::GetRapidity() const {
  double energy =
      std::sqrt(GetPt() * GetPt() + GetPz() * GetPz() + GetMass() * GetMass());
  if (energy != std::fabs(GetPz()))
    return 0.5 * std::log((energy + GetPz()) / (energy - GetPz()));
  return (GetPz() >= 0) ? 1.e30 : -1.e30;
}
