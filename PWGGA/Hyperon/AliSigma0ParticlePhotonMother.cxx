#include "AliSigma0ParticlePhotonMother.h"

ClassImp(AliSigma0ParticlePhotonMother)

    //____________________________________________________________________________________________________
    AliSigma0ParticlePhotonMother::AliSigma0ParticlePhotonMother()
    : AliSigma0ParticleBase(),
      fType(-1),
      fRecMassPhoton(0),
      fRecMassLambda(0),
      fRecMass(0),
      fMCLabel(-1),
      fPDGCode(-1),
      fV0(),
      fPhoton() {}

//____________________________________________________________________________________________________
AliSigma0ParticlePhotonMother::AliSigma0ParticlePhotonMother(
    const AliSigma0ParticleV0 &lambdaCandidate,
    const AliSigma0ParticleV0 &photonCandidate, const AliVEvent *inputEvent)
    : AliSigma0ParticleBase(),
      fType(1),
      fRecMassPhoton(0),
      fRecMassLambda(0),
      fRecMass(0),
      fMCLabel(-1),
      fPDGCode(-1),
      fV0(),
      fPhoton() {
  TLorentzVector track1, track2;
  track1.SetXYZM(lambdaCandidate.GetPx(), lambdaCandidate.GetPy(),
                 lambdaCandidate.GetPz(), lambdaCandidate.GetRecMass());
  track2.SetXYZM(photonCandidate.GetPx(), photonCandidate.GetPy(),
                 photonCandidate.GetPz(), photonCandidate.GetRecMass());
  TLorentzVector trackSum = track1 + track2;

  TLorentzVector track1MC, track2MC;
  track1MC.SetXYZM(lambdaCandidate.GetPxMC(), lambdaCandidate.GetPyMC(),
                   lambdaCandidate.GetPzMC(), lambdaCandidate.GetRecMass());
  track2MC.SetXYZM(photonCandidate.GetPxMC(), photonCandidate.GetPyMC(),
                   photonCandidate.GetPzMC(), photonCandidate.GetRecMass());
  TLorentzVector trackSumMC = track1MC + track2MC;

  fP[0] = trackSum.Px();
  fP[1] = trackSum.Py();
  fP[2] = trackSum.Pz();
  fPMC[0] = trackSumMC.Px();
  fPMC[1] = trackSumMC.Px();
  fPMC[2] = trackSumMC.Px();

  //  fPDGCode = pdg;
  fPt = std::sqrt(fP[0] * fP[0] + fP[1] * fP[1]);
  //  fTrackLabel = v0.GetID();
  fPhi = trackSum.Phi();
  fEta = trackSum.Eta();
  fRecMass = trackSum.M();
  // see https://arxiv.org/pdf/1703.04639.pdf
  fRecMassLambda = trackSum.M() - photonCandidate.GetRecMass();
  fRecMassPhoton = trackSum.M() - lambdaCandidate.GetRecMass() +
                   lambdaCandidate.GetPDGMass();
  fMass = trackSum.M() - lambdaCandidate.GetRecMass() +
          lambdaCandidate.GetPDGMass() - photonCandidate.GetRecMass();

  fUse = true;

  fV0 = lambdaCandidate;
  fPhoton = photonCandidate;
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

  fMCLabel = obj.GetMCLabel();
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

  return (*this);
}

//____________________________________________________________________________________________________
int AliSigma0ParticlePhotonMother::MatchToMC(
    const AliMCEvent *mcEvent, const int PIDmother,
    const std::vector<int> PIDdaughters, int &pdgLambdaMother,
    int &pdgPhotonMother) {
  const int labV0 = fV0.GetMCLabelV0();
  const int labPhoton = fPhoton.GetMCLabelV0();
  if (labV0 < 0 || labPhoton < 0) return -1;

  AliMCParticle *partV0 =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(labV0));
  AliMCParticle *partPhoton =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(labPhoton));
  if (!partV0 || !partPhoton) return -1;

  const int labMotherV0 = partV0->GetMother();
  const int labMotherPhoton = partPhoton->GetMother();

  AliMCParticle *partMotherV0 =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(labMotherV0));
  AliMCParticle *partMotherPhoton =
      static_cast<AliMCParticle *>(mcEvent->GetTrack(labMotherPhoton));
  if (!partMotherV0 || !partMotherPhoton) return -1;

  pdgLambdaMother = partMotherV0->PdgCode();
  pdgPhotonMother = partMotherPhoton->PdgCode();
  if ((pdgLambdaMother != PIDmother) || pdgPhotonMother != PIDmother) {
    return -1;
  }

  fMCLabel = labMotherV0;
  fPDGCode = pdgLambdaMother;

  fPMC[0] = partMotherV0->Px();
  fPMC[1] = partMotherV0->Py();
  fPMC[2] = partMotherV0->Pz();

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
      TMath::Sqrt(GetPt() * GetPt() + GetPz() * GetPz() + GetMass() * GetMass());
  if (energy != TMath::Abs(GetPz()))
    return 0.5 * TMath::Log((energy + GetPz()) / (energy - GetPz()));
  return (GetPz() >= 0) ? 1.e30 : -1.e30;
}
