#include "AliSigma0ParticleBase.h"
#include <iostream>
#include "TLorentzVector.h"

ClassImp(AliSigma0ParticleBase)

    //____________________________________________________________________________________________________
    AliSigma0ParticleBase::AliSigma0ParticleBase()
    : fP(),
      fPMC(),
      fPDGCode(0),
      fPDGCodeMother(-1),
      fMass(0),
      fQ(0),
      fPt(0),
      fTrackLabel(0),
      fMClabel(0),
      fPhi(0),
      fEta(0),
      fTheta(0),
      fPhiMC(0),
      fThetaMC(0),
      fCharge(0),
      fDCAz(0.f),
      fDCAr(0.f),
      fUse(true),
      fPhistar(),
      fAveragePhistar(0.f) {
  fP[0] = -0.;
  fP[1] = -0.;
  fP[2] = -0.;
  fPMC[0] = -0.;
  fPMC[1] = -0.;
  fPMC[2] = -0.;
}

//____________________________________________________________________________________________________
AliSigma0ParticleBase::AliSigma0ParticleBase(const AliESDtrack *track, int pdg,
                                             const float magneticField,
                                             int filterbit)
    : fP(),
      fPMC(),
      fPDGCode(0),
      fPDGCodeMother(-1),
      fMass(0),
      fQ(0),
      fPt(0),
      fTrackLabel(0),
      fMClabel(0),
      fPhi(0),
      fEta(0),
      fTheta(0),
      fPhiMC(0),
      fThetaMC(0),
      fCharge(0),
      fDCAz(0.f),
      fDCAr(0.f),
      fUse(true),
      fPhistar(),
      fAveragePhistar(0.f) {
  double trackMom[3];
  track->GetPxPyPz(trackMom);
  fP[0] = trackMom[0];
  fP[1] = trackMom[1];
  fP[2] = trackMom[2];
  fPMC[0] = -0.;
  fPMC[1] = -0.;
  fPMC[2] = -0.;

  fCharge = track->Charge();
  fPDGCode = pdg;
  fQ = track->Charge();
  fPt = track->Pt();
  fTrackLabel = (filterbit == 128) ? -track->GetID() - 1 : track->GetID();
  fPhi = track->Phi();
  fEta = track->Eta();
  fTheta = track->Theta();
  fUse = true;

  float TPCradii[9] = {85., 105., 125., 145., 165., 185., 205., 225., 245.};
  for (int radius = 0; radius < 9; radius++) {
    fPhistar.push_back(ComputePhiStar(track, magneticField, TPCradii[radius]));
  }

  float TPCradiiForAverage[33] = {
      85.,  90.,  95.,  100., 105., 110., 115., 120., 125., 130., 135.,
      140., 145., 150., 155., 160., 165., 170., 175., 180., 185., 190.,
      195., 200., 205., 210., 215., 220., 225., 230., 235., 240., 245.};

  for (int radius = 0; radius < 33; radius++) {
    fAveragePhistar +=
        ComputePhiStar(track, magneticField, TPCradiiForAverage[radius]);
  }
  fAveragePhistar /= (33.f + 1e-32);
}

//____________________________________________________________________________________________________
AliSigma0ParticleBase &AliSigma0ParticleBase::operator=(
    const AliSigma0ParticleBase &obj) {
  // Assignment operator
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
  fTrackLabel = -obj.GetTrackLabel() - 1;  // for filterbit 128
  fPhi = obj.GetPhi();
  fEta = obj.GetEta();
  fTheta = obj.GetTheta();
  fPhiMC = obj.GetPhiMC();
  fThetaMC = obj.GetThetaMC();

  fUse = obj.GetIsUse();

  fPhistar = obj.GetPhiStar();
  fAveragePhistar = obj.GetAveragePhiStar();

  return (*this);
}

//____________________________________________________________________________________________________
double AliSigma0ParticleBase::ComputeRelK(const AliSigma0ParticleBase &part2,
                                          bool debug) const {
  TLorentzVector track1, track2;
  track1.SetXYZM(fP[0], fP[1], fP[2], fMass);
  track2.SetXYZM(part2.GetPx(), part2.GetPy(), part2.GetPz(), part2.GetMass());

  TLorentzVector trackSum = track1 + track2;
  TLorentzVector track1cms = track1;
  TLorentzVector track2cms = track2;

  double beta = trackSum.Beta();
  double betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
  double betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
  double betaz = beta * cos(trackSum.Theta());

  track1cms.Boost(-betax, -betay, -betaz);
  track2cms.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK = track1cms - track2cms;
  double relK = 0.5 * trackRelK.P();

  if (debug) {
    std::cout << relK << " " << track1.Px() << " " << track1.Py() << " "
              << track1.Pz() << " " << track1.M() << " " << track2.Px() << " "
              << track2.Py() << " " << track2.Pz() << " " << track2.M() << "\n";
  }
  return relK;
}

//____________________________________________________________________________________________________
double AliSigma0ParticleBase::ComputeRelKMC(
    const AliSigma0ParticleBase &part2) const {
  TLorentzVector track1, track2;
  track1.SetXYZM(fPMC[0], fPMC[1], fPMC[2], fMass);
  track2.SetXYZM(part2.GetPxMC(), part2.GetPyMC(), part2.GetPzMC(),
                 part2.GetMass());

  TLorentzVector trackSum = track1 + track2;
  TLorentzVector track1cms = track1;
  TLorentzVector track2cms = track2;

  double beta = trackSum.Beta();
  double betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
  double betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
  double betaz = beta * cos(trackSum.Theta());

  track1cms.Boost(-betax, -betay, -betaz);
  track2cms.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK = track1cms - track2cms;
  double relK = 0.5 * trackRelK.P();

  return relK;
}

//____________________________________________________________________________________________________
float AliSigma0ParticleBase::ComputePhiStar(const AliESDtrack *track,
                                            const float magneticField,
                                            const float radius) const {
  const float phi0 = track->Phi();  // angle at primary vertex
  const float pt = track->GetTPCInnerParam()->Pt();
  const float charge = track->Charge();

  // To use the following equation:
  // pt must be given in GeV/c
  // bfield in T
  // chg in units of electric charge
  // radius in m

  // 0.3 is a conversion factor for pt and bfield can be plugged in in terms of
  // GeV/c and electric charge, 0.1 converts the magnetic field to Tesla, 0.01
  // transforms the radius from cm to m
  return phi0 + TMath::ASin(0.1 * charge * magneticField * 0.3 * radius * 0.01 /
                            (2. * pt));
}

//____________________________________________________________________________________________________
void AliSigma0ParticleBase::ProcessMCInfo(AliMCParticle *mcParticle,
                                          AliMCEvent *mcEvent) {
  fPMC[0] = mcParticle->Px();
  fPMC[1] = mcParticle->Py();
  fPMC[2] = mcParticle->Pz();
  fPDGCode = mcParticle->PdgCode();
  fMClabel = mcParticle->GetLabel();
  fPhiMC = mcParticle->Phi();
  fThetaMC = mcParticle->Theta();

  if (mcParticle->GetMother() != 0) {
    AliMCParticle *mcMother = static_cast<AliMCParticle *>(
        mcEvent->GetTrack(mcParticle->GetMother()));
    fPDGCodeMother = mcMother->PdgCode();
  }
}
