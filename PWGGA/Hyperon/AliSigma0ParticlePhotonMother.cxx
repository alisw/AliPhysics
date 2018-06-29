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
    const AliSigma0ParticleV0 &photonCandidate, const AliVEvent *inputEvent)
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
                 photonCandidate.GetPz(), photonCandidate.GetRecMass());
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
  fPhoton = photonCandidate;
}

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
int AliSigma0ParticlePhotonMother::MatchToMC(
    const AliMCEvent *mcEvent, const int PIDmother,
    const std::vector<int> PIDdaughters) const {
  const int labV0 = fV0.GetMCLabelV0();
  const int labPhoton = fPhoton.GetMCLabelV0();
  if (labV0 < 0 || labPhoton < 0) return -1;

  AliMCParticle *partV0 =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(labV0));
  AliMCParticle *partPhoton =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(labPhoton));
  if (!partV0 || !partPhoton) return -1;

  const int pidV0 = partV0->PdgCode();
  const int pidPhoton = partPhoton->PdgCode();
  if (!((pidV0 == PIDdaughters[0] && pidPhoton == PIDdaughters[1]) ||
        (pidV0 == PIDdaughters[1] && pidPhoton == PIDdaughters[0]))) {
    return -1;
  }

  const int labMotherV0 = partV0->GetMother();
  const int labMotherPhoton = partPhoton->GetMother();
  if (labMotherV0 < 0 || labMotherPhoton < 0) return -1;

  AliMCParticle *partMotherV0 =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(labMotherV0));
  AliMCParticle *partMotherPhoton =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(labMotherPhoton));
  if (!partMotherV0 || !partMotherPhoton) return -1;

  const int pdgMotherV0 = partMotherV0->PdgCode();
  const int pdgMotherPhoton = partMotherPhoton->PdgCode();
  if ((pdgMotherV0 != pdgMotherPhoton) || pdgMotherV0 != PIDmother) {
    return -1;
  }

  return labMotherV0;
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
