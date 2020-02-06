/*
 * AliFemtoPPbpbLamBaseParticle.cxx
 *
 *  Created on: Aug 29, 2017
 *      Author: gu74req
 */

#include "AliFemtoDreamBasePart.h"
#include "AliMCParticle.h"

#include <iostream>
ClassImp(AliFemtoDreamBasePart)
AliFemtoDreamBasePart::AliFemtoDreamBasePart(const int part)
    : fIsReset(false),
      fGTI(0),
      fVGTI(0),
      fTrackBufferSize(0),
      fP(part),
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
      fEvtMultiplicity(-1) {
}

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
      fEvtMultiplicity(part.fEvtMultiplicity) {
}

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
    const AliAODConversionPhoton *gamma, const AliVTrack *posTrack,
    const AliVTrack *negTrack, const AliVEvent *inputEvent)
    : fIsReset(false),
      fGTI(0),
      fVGTI(0),
      fTrackBufferSize(0),
      fP(3),
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
  double momV0[3] = { 0, 0, 0 };
  momV0[0] = gamma->Px();
  momV0[1] = gamma->Py();
  momV0[2] = gamma->Pz();
  SetMomentum(0,{gamma->GetPx(), gamma->GetPy(), gamma->GetPz()});
  SetMomentum(1,{posTrack->Px(), posTrack->Py(), posTrack->Pz()});
  SetMomentum(2,{negTrack->Px(), negTrack->Py(), negTrack->Pz()});
  // Recalculated V0 Position vector
  double PosV0[3] = { gamma->GetConversionX()
      - inputEvent->GetPrimaryVertex()->GetX(), gamma->GetConversionY()
      - inputEvent->GetPrimaryVertex()->GetY(), gamma->GetConversionZ()
      - inputEvent->GetPrimaryVertex()->GetZ() };

  double momV02 = GetMomentum().X() * GetMomentum().X()
      + GetMomentum().Y() * GetMomentum().Y()
      + GetMomentum().Z() * GetMomentum().Z();
  double PosV02 = PosV0[0] * PosV0[0] + PosV0[1] * PosV0[1]
      + PosV0[2] * PosV0[2];

  double cosinePointingAngle =
      (momV02 * PosV02 > 0.0) ?
          (PosV0[0] * momV0[0] + PosV0[1] * momV0[1] + PosV0[2] * momV0[2])
              / TMath::Sqrt(momV02 * PosV02) :
          -999.f;
  fCPA = cosinePointingAngle;

  const int posLabel = posTrack->GetID();
  const int negLabel = negTrack->GetID();
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
  fGTI = nullptr;
  fVGTI = nullptr;
  fP.clear();
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
  this->SetMCPDGCode(
      static_cast<AliAODMCParticle*>(evt->GetTrack(lastMother))->GetPdgCode());
}

void AliFemtoDreamBasePart::SetMCParticleRePart(AliAODMCParticle *mcPart) {
  this->SetPt(mcPart->Pt());
  this->SetMomentum(0, mcPart->Px(), mcPart->Py(), mcPart->Pz());
  this->SetTheta(mcPart->Theta());
  this->SetEta(mcPart->Eta());
  this->SetPhi(mcPart->Phi());
  this->fIsSet = true;
  this->SetUse(true);
  this->fIsReset = false;
  this->SetCharge(0);
}

void AliFemtoDreamBasePart::ResetMCInfo() {
  this->SetPt(0);
  this->SetMomentum(0,0,0,0);
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
  float TPCradii[9] = { 85., 105., 125., 145., 165., 185., 205., 225., 245. };
  float phi0 = track->Phi();
  ;
  float pt = track->Pt();
  float chg = track->Charge();
  for (int radius = 0; radius < 9; radius++) {
    tmpVec.push_back(
        phi0
            - TMath::ASin(
                0.1 * chg * bfield * 0.3 * TPCradii[radius] * 0.01
                    / (2. * pt)));
  }
}

void AliFemtoDreamBasePart::DumpParticleInformation() {

  auto dumpfloatVector = [](std::vector<float> &vec) {
    for (auto it :vec) {
      std::cout << it << " ";
    }
    std::cout << "\n";
  };

  auto dumpintVector = [](std::vector<int> &vec) {
    for (auto it :vec) {
      std::cout << it << " ";
    }
    std::cout << "\n";
  };

  std::cout << "Dumping the particle information\n";
  std::cout << "Momentum- x: " << GetMomentum().X() << " y: "
            << GetMomentum().Y() << " z: " << GetMomentum().Z() << "\n";
  std::cout << "Momentum (MC)  - x: " << fMCP.X() << " y: " << fMCP.Y()
            << " z: " << fMCP.Z() << "\n";
  std::cout << "pT: " << fPt << "\n";
  std::cout << "p TPC " << fP_TPC << "\n";
  std::cout << "pT (MC): " << fMCPt << "\n";

  std::cout << "Eta - entries " << fEta.size() << "\n";
  dumpfloatVector(fEta);
  std::cout << "Theta - entries " << fTheta.size() << "\n";
  dumpfloatVector(fTheta);
  std::cout << "Phi - entries " << fPhi.size() << "\n";
  dumpfloatVector(fPhi);
  std::cout << "Theta (MC) - entries " << fMCTheta.size() << "\n";
  dumpfloatVector(fMCTheta);
  std::cout << "Phi (MC) - entries " << fMCPhi.size() << "\n";
  dumpfloatVector(fMCPhi);

  std::cout << "Track IDs - entries " << fIDTracks.size() << "\n";
  dumpintVector(fIDTracks);
  std::cout << "Charge - entries " << fCharge.size() << "\n";
  dumpintVector(fCharge);

  std::cout << "CPA: " << fCPA << "\n";
  std::cout << "Invariant mass " << fInvMass << "\n";
  std::cout << "Origin " << fOrigin << "\n";
  std::cout << "PDG code " << fPDGCode << "\n";
  std::cout << "MC PDG code " << fMCPDGCode << "\n";
  std::cout << "PDG Mother weak " << fPDGMotherWeak << "\n";
  std::cout << "PDG Mother " << fMotherPDG << "\n";
  std::cout << "Mother ID " << fMotherID << "\n";
  std::cout << "ID " << fID << "\n";
  std::cout << "Use particle " << fUse << "\n";
  std::cout << "Is set " << fIsSet << "\n";
}
