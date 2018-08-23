/*
 * AliFemtoPPbpbLamBaseParticle.cxx
 *
 *  Created on: Aug 29, 2017
 *      Author: gu74req
 */

#include "AliFemtoDreamBasePart.h"
#include "AliSigma0ParticleV0.h"

#include <iostream>
ClassImp(AliFemtoDreamBasePart)
AliFemtoDreamBasePart::AliFemtoDreamBasePart()
    : fIsReset(false),
      fGTI(0),
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
      fPhiAtRadius(0),
      fMCPhi(0),
      fIDTracks(0),
      fCharge(0),
      fCPA(0),
      fOrigin(kUnknown),
      fPDGCode(0),
      fMCPDGCode(0),
      fPDGMotherWeak(0),
      fMotherID(0),
      fEvtNumber(0),
      fIsMC(false),
      fUse(true),
      fIsSet(true),
      fEvtMultiplicity(-1) {
}

AliFemtoDreamBasePart::AliFemtoDreamBasePart(const AliFemtoDreamBasePart& part)
    : fIsReset(part.fIsReset),
      fGTI(part.fGTI),
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
      fMCPhi(part.fMCPhi),
      fIDTracks(part.fIDTracks),
      fCharge(part.fCharge),
      fCPA(part.fCPA),
      fOrigin(part.fOrigin),
      fPDGCode(part.fPDGCode),
      fMCPDGCode(part.fMCPDGCode),
      fPDGMotherWeak(part.fPDGMotherWeak),
      fMotherID(part.fMotherID),
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
  fTrackBufferSize = obj.fTrackBufferSize;
  fP = obj.fP;
  fMCPt = obj.fMCPt;
  fP_TPC = obj.fP_TPC;
  fEta = obj.fEta;
  fTheta = obj.fTheta;
  fMCTheta = obj.fMCTheta;
  fPhi = obj.fPhi;
  fPhiAtRadius = obj.fPhiAtRadius;
  fMCPhi = obj.fMCPhi;
  fIDTracks = obj.fIDTracks;
  fCharge = obj.fCharge;
  fCPA = obj.fCPA;
  fOrigin = obj.fOrigin;
  fPDGCode = obj.fPDGCode;
  fMCPDGCode = obj.fMCPDGCode;
  fPDGMotherWeak = obj.fPDGMotherWeak;
  fMotherID = obj.fMotherID;
  fEvtNumber = obj.fEvtNumber;
  fIsMC = obj.fIsMC;
  fUse = obj.fUse;
  fIsSet = obj.fIsSet;
  fEvtMultiplicity = obj.fEvtMultiplicity;
  return (*this);
}

AliFemtoDreamBasePart::AliFemtoDreamBasePart(
    const AliSigma0ParticlePhotonMother &mother)
    : fIsReset(false),
      fGTI(0),
      fTrackBufferSize(0),
      fP(TVector3(mother.GetPx(), mother.GetPy(), mother.GetPz())),
      fMCP(TVector3(mother.GetPxMC(), mother.GetPyMC(), mother.GetPzMC())),
      fPt(mother.GetPt()),
      fMCPt(0),
      fP_TPC(0),
      fEta(),
      fTheta(0),
      fMCTheta(0),
      fPhi(),
      fPhiAtRadius(0),
      fMCPhi(0),
      fIDTracks(),
      fCharge(0),
      fCPA(0),
      fOrigin(kUnknown),
      fPDGCode(0),
      fMCPDGCode(0),
      fPDGMotherWeak(0),
      fMotherID(0),
      fEvtNumber(0),
      fIsMC(false),
      fUse(true),
      fIsSet(true),
      fEvtMultiplicity(-1) {
  fEta.push_back(mother.GetV0().GetPosDaughter().GetEta());
  fEta.push_back(mother.GetPhoton().GetPosDaughter().GetEta());
  fEta.push_back(mother.GetV0().GetNegDaughter().GetEta());
  fEta.push_back(mother.GetPhoton().GetNegDaughter().GetEta());

  fPhi.push_back(mother.GetV0().GetPosDaughter().GetPhi());
  fPhi.push_back(mother.GetPhoton().GetPosDaughter().GetPhi());
  fPhi.push_back(mother.GetV0().GetNegDaughter().GetPhi());
  fPhi.push_back(mother.GetPhoton().GetNegDaughter().GetPhi());

  fIDTracks.push_back(mother.GetV0().GetTrackLabelPos());
  fIDTracks.push_back(mother.GetV0().GetTrackLabelNeg());
  fIDTracks.push_back(mother.GetPhoton().GetTrackLabelPos());
  fIDTracks.push_back(mother.GetPhoton().GetTrackLabelNeg());
}

AliFemtoDreamBasePart::~AliFemtoDreamBasePart() {
}

void AliFemtoDreamBasePart::SetMCParticle(AliAODMCParticle *mcPart,
                                          AliMCEvent *evt) {
  this->SetPt(mcPart->Pt());
  this->SetMomentum(mcPart->Px(), mcPart->Py(), mcPart->Pz());
  this->SetEta(mcPart->Eta());
  this->SetTheta(mcPart->Theta());
  this->SetPhi(mcPart->Phi());
  this->SetCharge(mcPart->Charge());
  this->SetPDGCode(mcPart->GetPdgCode());
  this->fIsSet = true;
  this->SetUse(true);
  this->fIsReset = false;
  //find the original mother in the mc stack.
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
  //a change
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
