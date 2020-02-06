/*
 * AliOtonOmegaCascade.cxx
 *
 *  Created on: Jan 11, 2018
 *      Author: gu74req
 */
#include <TMath.h>

#include "AliOtonOmegaCascade.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODMCParticle.h"

ClassImp(AliOtonOmegaCascade)
AliOtonOmegaCascade::AliOtonOmegaCascade()
    : AliFemtoDreamBasePart(4),
      fPosDaug(new AliFemtoDreamTrack()),
      fNegDaug(new AliFemtoDreamTrack()),
      fBach(new AliFemtoDreamTrack()),
      fXiMass(0),
      fOmegaMass(0),
      fDCAXiDaug(0),
      fTransRadius(0),
      fRapXi(0),
      fRapOmega(0),
      fAlphaXi(0),
      fPtArmXi(0),
      fDCAXiPrimVtx(0),
      fLength(0),
      fMassv0(0),
      fv0DCADaug(0),
      fv0DCAPrimVtx(0),
      fv0Momentum(),
      fv0Pt(0),
      fDCABachPrimVtx(0),
      fv0TransRadius(0),
      fv0CPA(0),
      fv0PosToPrimVtx(0),
      fv0NegToPrimVtx(0),
      fv0ToXiPointAngle(0),
      fv0Length(0),
      fDCAv0Xi(0),
      fv0PDG(0) {
}

AliOtonOmegaCascade::~AliOtonOmegaCascade() {
  if (fPosDaug) {
    delete fPosDaug;
  }
  if (fNegDaug) {
    delete fNegDaug;
  }
  if (fBach) {
    delete fBach;
  }
}

void AliOtonOmegaCascade::SetCascade(AliAODEvent *evt, AliAODcascade *casc) {
  //xi business
  Reset();
  fIsReset = false;
  fIsSet = true;
  this->SetCharge(casc->ChargeXi());
  this->SetMomentum(0, casc->MomXiX(), casc->MomXiY(), casc->MomXiZ());
  this->SetPt(GetMomentum().Pt());
  double PrimVtx[3] = { 99., 99., 99 };
  double decayPosXi[3] = { casc->DecayVertexXiX(), casc->DecayVertexXiY(), casc
      ->DecayVertexXiZ() };
  fXiMass = casc->MassXi();
  fOmegaMass = casc->MassOmega();
  fDCAXiDaug = casc->DcaXiDaughters();
  evt->GetPrimaryVertex()->GetXYZ(PrimVtx);
  this->SetEvtNumber(evt->GetRunNumber());
  fCPA = casc->CosPointingAngleXi(PrimVtx[0], PrimVtx[1], PrimVtx[2]);
  fTransRadius = TMath::Sqrt(
      decayPosXi[0] * decayPosXi[0] + decayPosXi[1] * decayPosXi[1]);
  fRapXi = casc->RapXi();
  fRapOmega = casc->RapOmega();
  this->SetEta(casc->Eta());
  this->SetTheta(GetMomentum().Theta());
  this->SetPhi(GetMomentum().Phi());
  fAlphaXi = casc->AlphaXi();
  fPtArmXi = casc->PtArmXi();
  fDCAXiPrimVtx = casc->DcaXiToPrimVertex(PrimVtx[0], PrimVtx[1], PrimVtx[2]);
  fLength = TMath::Sqrt(
      TMath::Power((decayPosXi[0] - PrimVtx[0]), 2)
          + TMath::Power((decayPosXi[1] - PrimVtx[1]), 2)
          + TMath::Power((decayPosXi[2] - PrimVtx[2]), 2));
  if (this->GetMomentum().Mag() > 0) {
    static double mass =
        TDatabasePDG::Instance()->GetParticle(fPDGCode)->Mass();
    fLength = fLength * mass / this->GetMomentum().Mag();
  } else {
    fLength = -1.;
  }

  //From here daughter business
  AliAODTrack *nTrackXi = dynamic_cast<AliAODTrack*>(casc->GetDaughter(1));
  AliAODTrack *pTrackXi = dynamic_cast<AliAODTrack*>(casc->GetDaughter(0));
  AliAODTrack *bachTrackXi = dynamic_cast<AliAODTrack*>(casc->GetDecayVertexXi()
      ->GetDaughter(0));
  fNegDaug->SetTrack(nTrackXi);
  fPosDaug->SetTrack(pTrackXi);
  fBach->SetTrack(bachTrackXi);
  fNegDaug->SetMomentum(0, casc->MomNegX(), casc->MomNegY(), casc->MomNegZ());
  fPosDaug->SetMomentum(0, casc->MomPosX(), casc->MomPosY(), casc->MomPosZ());
  fBach->SetMomentum(0, casc->MomBachX(), casc->MomBachY(), casc->MomBachZ());
  this->SetMomentum(1, casc->MomNegX(), casc->MomNegY(), casc->MomNegZ());
  this->SetMomentum(2, casc->MomPosX(), casc->MomPosY(), casc->MomPosZ());
  this->SetMomentum(3, casc->MomBachX(), casc->MomBachY(), casc->MomBachZ());

  double posMom[3] = { 0. };
  double negMom[3] = { 0. };
  double bachMom[3] = { 0. };

  fNegDaug->GetMomentum().GetXYZ(negMom);
  fPosDaug->GetMomentum().GetXYZ(posMom);
  fBach->GetMomentum().GetXYZ(bachMom);

//  fMass = MassCasc();
  this->SetIDTracks(nTrackXi->GetID());
  this->SetIDTracks(pTrackXi->GetID());
  this->SetIDTracks(bachTrackXi->GetID());

  this->SetEta(fNegDaug->GetMomentum().Eta());
  this->SetEta(fPosDaug->GetMomentum().Eta());
  this->SetEta(fBach->GetMomentum().Eta());

  this->SetTheta(fNegDaug->GetMomentum().Theta());
  this->SetTheta(fPosDaug->GetMomentum().Theta());
  this->SetTheta(fBach->GetMomentum().Theta());

  this->SetPhi(fNegDaug->GetMomentum().Phi());
  this->SetPhi(fPosDaug->GetMomentum().Phi());
  this->SetPhi(fBach->GetMomentum().Phi());

  this->SetCharge(nTrackXi->Charge());
  this->SetCharge(pTrackXi->Charge());
  this->SetCharge(bachTrackXi->Charge());
  if (fNegDaug->IsSet()) {
    this->SetPhiAtRadius(fNegDaug->GetPhiAtRaidius().at(0));
  }
  if (fPosDaug->IsSet()) {
    this->SetPhiAtRadius(fPosDaug->GetPhiAtRaidius().at(0));
  }
  if (fBach->IsSet()) {
    this->SetPhiAtRadius(fBach->GetPhiAtRaidius().at(0));
  }

  //v0 business
  if (bachTrackXi->Charge() < 0) {  //Xi minus
    fMassv0 = casc->MassLambda();
  } else {  //Xi plus case
    //Workaround to identify also Xi+. We usually only set the PDG code once,
    //assuming a Xi -, meaning we expect a Lambda, not an Antilambda.
    fMassv0 = casc->MassAntiLambda();
  }
  fv0Momentum.SetXYZ(casc->MomV0X(), casc->MomV0Y(), casc->MomV0Z());
  fv0Pt = casc->MomV0X() * casc->MomV0X() + casc->MomV0Y() * casc->MomV0Y();
  fv0DCADaug = casc->DcaV0Daughters();
  fv0DCAPrimVtx = casc->DcaV0ToPrimVertex();
  double decayPosV0[3] = { casc->DecayVertexV0X(), casc->DecayVertexV0Y(), casc
      ->DecayVertexV0Z() };
  fv0TransRadius = TMath::Sqrt(
      decayPosV0[0] * decayPosV0[0] + decayPosV0[1] * decayPosV0[1]);
  fv0CPA = casc->CosPointingAngle(evt->GetPrimaryVertex());

  fDCABachPrimVtx = casc->DcaBachToPrimVertex();
  fv0PosToPrimVtx = casc->DcaPosToPrimVertex();
  fv0NegToPrimVtx = casc->DcaNegToPrimVertex();
  fv0ToXiPointAngle = casc->CosPointingAngle(casc->GetDecayVertexXi());
  float LamPDGMass = 1.115683;
  fv0Length = TMath::Sqrt(
      TMath::Power((decayPosV0[0] - decayPosXi[0]), 2)
          + TMath::Power((decayPosV0[1] - decayPosXi[1]), 2)
          + TMath::Power((decayPosV0[2] - decayPosXi[2]), 2));
  //  Float_t lctauV0 = -1.;
  if (fv0Momentum.Mag() > 0) {
    fv0Length = fv0Length * LamPDGMass / fv0Momentum.Mag();
  } else {
    fv0Length = -1.;
  }
  fDCAv0Xi = TMath::Sqrt(
      TMath::Power((decayPosV0[0] - decayPosXi[0]), 2)
          + TMath::Power((decayPosV0[1] - decayPosXi[1]), 2));
  if (fIsMC) {
    SetMCMotherInfo(evt, casc);
  }
}

void AliOtonOmegaCascade::SetCascade(AliESDEvent *evt, AliMCEvent *mcEvent,
                                      AliESDcascade *casc) {
  Reset();
  fIsReset = false;
  this->fIsSet = true;
  int idxPosFromV0Dghter = casc->GetPindex();
  int idxNegFromV0Dghter = casc->GetNindex();
  int idxBachFromCascade = casc->GetBindex();

  bool IsOmegaTrack = true;
  AliESDtrack *esdCascadePos = evt->GetTrack(idxPosFromV0Dghter);
  fPosDaug->SetTrack(esdCascadePos,mcEvent,-1,false,IsOmegaTrack);
  AliESDtrack *esdCascadeNeg = evt->GetTrack(idxNegFromV0Dghter);
  fNegDaug->SetTrack(esdCascadeNeg,mcEvent,-1,false,IsOmegaTrack);
  AliESDtrack *esdCascadeBach = evt->GetTrack(idxBachFromCascade);
  fBach->SetTrack(esdCascadeBach,mcEvent,-1,false,IsOmegaTrack);

  // Identification of the V0 within the esdCascade (via both daughter track indices)
  AliESDv0 * currentV0 = 0x0;
  int idxV0FromCascade = -1;
  for (int iV0 = 0; iV0 < evt->GetNumberOfV0s(); ++iV0) {
    currentV0 = evt->GetV0(iV0);
    // Make sure not to use online v0 here
    // (the cascade is reconstructed with offline v0 only)
    if (currentV0->GetOnFlyStatus() == kTRUE) continue;
    int posCurrentV0 = currentV0->GetPindex();
    int negCurrentV0 = currentV0->GetNindex();
    if (posCurrentV0 == idxPosFromV0Dghter
        && negCurrentV0 == idxNegFromV0Dghter) {
      idxV0FromCascade = iV0;
      break;
    }
  }
  if (idxV0FromCascade < 0) {
    AliWarning("Cascade - no matching for the V0, skipping ... \n");
    return;
  }
  double posMom[3] = { 0. };
  double negMom[3] = { 0. };
  double bachMom[3] = { 0. };
  //the momenta of the daughters have to be taken at the v0 vertex,
  //the bachelor momenta at the cascade vertex
  currentV0->GetPPxPyPz(posMom[0], posMom[1], posMom[2]);
  fPosDaug->SetMomentum(0, posMom[0], posMom[1], posMom[2]);
  fPosDaug->SetPt(posMom[0]*posMom[0]+posMom[1]*posMom[1]);
  currentV0->GetNPxPyPz(negMom[0], negMom[1], negMom[2]);
  fNegDaug->SetMomentum(0, negMom[0], negMom[1], negMom[2]);
  fNegDaug->SetPt(negMom[0]*negMom[0]+negMom[1]*negMom[1]);
  casc->GetBPxPyPz(bachMom[0], bachMom[1], bachMom[2]);
  fBach->SetMomentum(0, bachMom[0], bachMom[1], bachMom[2]);
  fBach->SetPt(bachMom[0]*bachMom[0]+bachMom[1]*bachMom[1]);

  this->SetMomentum(1, negMom[0], negMom[1], negMom[2]);
  this->SetMomentum(2, posMom[0], posMom[1], posMom[2]);
  this->SetMomentum(3, bachMom[0], bachMom[1], bachMom[2]);

  TVector3 xiMom = fPosDaug->GetMomentum();
  xiMom += fNegDaug->GetMomentum();
  xiMom += fBach->GetMomentum();

  this->SetMomentum(0, xiMom.X(), xiMom.Y(), xiMom.Z());
  this->SetPt(GetMomentum().Pt());

  double PrimVtx[3] = { 99., 99., 99 };
  double decayPosXi[3] = { 0. };
  casc->GetXYZcascade(decayPosXi[0], decayPosXi[1], decayPosXi[2]);
  evt->GetPrimaryVertex()->GetXYZ(PrimVtx);

//  fMass = MassCasc();
  fXiMass = MassXi();
  fOmegaMass = MassOmega();
  fDCAXiDaug = casc->GetDcaXiDaughters();
  this->SetEvtNumber(evt->GetRunNumber());

  fTransRadius = TMath::Sqrt(
      decayPosXi[0] * decayPosXi[0] + decayPosXi[1] * decayPosXi[1]);
  this->SetCharge(casc->Charge());
  fRapXi = casc->RapXi();
  fRapOmega = casc->RapOmega();
  this->SetEta(casc->Eta());
  this->SetTheta(casc->Theta());
  this->SetPhi(casc->Phi());
  fAlphaXi = casc->AlphaXi();
  fPtArmXi = casc->PtArmXi();
  fCPA = casc->GetCascadeCosineOfPointingAngle(PrimVtx[0], PrimVtx[1],
                                               PrimVtx[2]);
  fDCAXiPrimVtx = DcaXiToPrimVertex(PrimVtx, decayPosXi);
  fLength = TMath::Sqrt(
      TMath::Power((decayPosXi[0] - PrimVtx[0]), 2)
          + TMath::Power((decayPosXi[1] - PrimVtx[1]), 2)
          + TMath::Power((decayPosXi[2] - PrimVtx[2]), 2));
  if (this->GetMomentum().Mag() > 0) {
    static double mass =
        TDatabasePDG::Instance()->GetParticle(fPDGCode)->Mass();
    fLength = fLength * mass / this->GetMomentum().Mag();
  } else {
    fLength = -1.;
  }
  this->SetIDTracks(esdCascadeNeg->GetID());
  this->SetIDTracks(esdCascadePos->GetID());
  this->SetIDTracks(esdCascadeBach->GetID());

  this->SetEta(esdCascadeNeg->Eta());
  this->SetEta(esdCascadePos->Eta());
  this->SetEta(esdCascadeBach->Eta());

  this->SetTheta(esdCascadeNeg->Theta());
  this->SetTheta(esdCascadePos->Theta());
  this->SetTheta(esdCascadeBach->Theta());

  this->SetPhi(esdCascadeNeg->Phi());
  this->SetPhi(esdCascadePos->Phi());
  this->SetPhi(esdCascadeBach->Phi());

  this->SetCharge(esdCascadeNeg->Charge());
  this->SetCharge(esdCascadePos->Charge());
  this->SetCharge(esdCascadeBach->Charge());

  if (fNegDaug->IsSet()) {
    this->SetPhiAtRadius(fNegDaug->GetPhiAtRaidius().at(0));
  }
  if (fPosDaug->IsSet()) {
    this->SetPhiAtRadius(fPosDaug->GetPhiAtRaidius().at(0));
  }
  if (fBach->IsSet()) {
    this->SetPhiAtRadius(fBach->GetPhiAtRaidius().at(0));
  }

  TVector3 v0Mom = fPosDaug->GetMomentum() + fNegDaug->GetMomentum();
  fv0Momentum.SetXYZ(v0Mom.X(), v0Mom.Y(), v0Mom.Z());
  fv0Pt = fv0Momentum.Pt();
  if (esdCascadeBach->Charge() < 0) {
    fMassv0 = Massv0(fPosDaug->GetPDGCode(), fNegDaug->GetPDGCode());
  } else {
    //Workaround to identify also Xi+. We usually only set the PDG code once,
    //assuming a Xi -, meaning we expect a Lambda, not an Antilambda.
    fMassv0 = Massv0(fNegDaug->GetPDGCode(), fPosDaug->GetPDGCode());
  }
  fv0DCADaug = casc->GetDcaV0Daughters();
  fv0DCAPrimVtx = currentV0->GetD(evt->GetPrimaryVertex()->GetX(),
                                  evt->GetPrimaryVertex()->GetY(),
                                  evt->GetPrimaryVertex()->GetZ());
  double decayPosV0[3] = { 0. };
  currentV0->GetXYZ(decayPosV0[0], decayPosV0[1], decayPosV0[2]);
  fv0TransRadius = TMath::Sqrt(
      decayPosV0[0] * decayPosV0[0] + decayPosV0[1] * decayPosV0[1]);
  fv0CPA = currentV0->GetV0CosineOfPointingAngle(
      evt->GetPrimaryVertex()->GetX(), evt->GetPrimaryVertex()->GetY(),
      evt->GetPrimaryVertex()->GetZ());

  fDCABachPrimVtx = TMath::Abs(
      esdCascadeBach->GetD(evt->GetPrimaryVertex()->GetX(),
                           evt->GetPrimaryVertex()->GetY(),
                           evt->GetMagneticField()));
  fv0PosToPrimVtx = TMath::Abs(
      esdCascadePos->GetD(evt->GetPrimaryVertex()->GetX(),
                          evt->GetPrimaryVertex()->GetY(),
                          evt->GetMagneticField()));
  fv0NegToPrimVtx = TMath::Abs(
      esdCascadeNeg->GetD(evt->GetPrimaryVertex()->GetX(),
                          evt->GetPrimaryVertex()->GetY(),
                          evt->GetMagneticField()));
  fv0ToXiPointAngle = casc->GetV0CosineOfPointingAngle(decayPosXi[0],
                                                       decayPosXi[1],
                                                       decayPosXi[2]);
  fv0Length = TMath::Sqrt(
      TMath::Power((decayPosV0[0] - decayPosXi[0]), 2)
          + TMath::Power((decayPosV0[1] - decayPosXi[1]), 2)
          + TMath::Power((decayPosV0[2] - decayPosXi[2]), 2));
  if (fv0Momentum.Mag() > 0) {
    float LamPDGMass = TDatabasePDG::Instance()->GetParticle(fv0PDG)->Mass();
    fv0Length = fv0Length * LamPDGMass / fv0Momentum.Mag();
  } else {
    fv0Length = -1.;
  }
  fDCAv0Xi = TMath::Sqrt(
      TMath::Power((decayPosV0[0] - decayPosXi[0]), 2)
          + TMath::Power((decayPosV0[1] - decayPosXi[1]), 2));

  return;
}

void AliOtonOmegaCascade::Reset() {
  if (!fIsReset) {
//    fMass = 0;
    fXiMass = 0;
    fOmegaMass = 0;
    fDCAXiPrimVtx = 0;
    fDCAXiDaug = 0;
    fTransRadius = 0;
    fRapXi = -99;
    fRapOmega = -99;
    fAlphaXi = 99;
    fPtArmXi = 99;
    fLength = 0;
    fMassv0 = 0;
    fv0DCADaug = 0;
    fv0DCAPrimVtx = 0;
    fv0Momentum.SetXYZ(0, 0, 0);
    fv0Pt = 0;
    fDCABachPrimVtx = 0;
    fv0TransRadius = 0;
    fv0CPA = 0;
    fv0PosToPrimVtx = 0;
    fv0NegToPrimVtx = 0;
    fv0ToXiPointAngle = 0;
    fv0Length = 0;
    fDCAv0Xi = 0;
    GetMomentum(0).SetXYZ(0, 0, 0);
    GetMomentum(1).SetXYZ(0, 0, 0);
    GetMomentum(2).SetXYZ(0, 0, 0);
    GetMomentum(3).SetXYZ(0, 0, 0);
    fMCP.SetXYZ(0, 0, 0);
    fPt = 0;
    fMCPt = 0;
    fMCPt = 0;
    fP_TPC = 0;
    fEta.clear();
    fTheta.clear();
    fMCTheta.clear();
    fPhi.clear();
    fPhiAtRadius.clear();
    fMCPhi.clear();
    fIDTracks.clear();
    fCharge.clear();
    fCPA = 0;
    fOrigin = AliFemtoDreamBasePart::kUnknown;
    //we don't want to reset the fPDGCode
    fMCPDGCode = 0;
    fPDGMotherWeak = 0;
    //we don't want to reset isMC
    fUse = false;
    fIsSet = false;
    fIsReset = true;
  }
}

void AliOtonOmegaCascade::SetMCMotherInfo(AliAODEvent *evt,
                                           AliAODcascade *casc) {
  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(evt->FindListObject(
      AliAODMCParticle::StdBranchName()));
  if (!mcarray) {
    AliFatal("No MC Array found");
  }
  if (fBach->IsSet() && fPosDaug->IsSet() && fNegDaug->IsSet()) {
    //look if the bachelor is from a weak decay and find the label of the
    //mother
    int labelBachMother = -1;
    if (fBach->GetParticleOrigin() == AliFemtoDreamBasePart::kWeak) {
      //      std::cout << "Bachelor ID" << fBach->GetIDTracks().at(0) << "\n";
      AliAODTrack* bachTrk = dynamic_cast<AliAODTrack*>(casc->GetDecayVertexXi()
          ->GetDaughter(0));
      int labelBach = bachTrk->GetLabel();
      labelBachMother =
          ((AliAODMCParticle*) mcarray->At(labelBach))->GetMother();
      if (labelBachMother < 0) {
        //This should not happen, just in case somethings very fishy
        this->SetParticleOrigin(AliFemtoDreamBasePart::kFake);
      } else {
//        std::cout << "Found a Track weakling and his mother is: " << labelBachMother <<'\n';
//        std::cout << "PDG Bach: "<<((AliAODMCParticle*)mcarray->At(labelBach))->GetPdgCode() << "\n";
//        std::cout << "PDG Bach-Mother: "<<((AliAODMCParticle*)mcarray->At(labelBachMother))->GetPdgCode() << "\n";
        //look if a v0 exists, and if this v0 stems from a weak decay, and
        //get the label of the mother particle
        int labelv0Mother = -1;
        int PDGDaug[2];
        //order doesn't matter for the MatchToMC
        PDGDaug[0] = TMath::Abs(fPosDaug->GetPDGCode());
        PDGDaug[1] = TMath::Abs(fNegDaug->GetPDGCode());
        int labelv0 = casc->MatchToMC(TMath::Abs(fv0PDG), mcarray, 2, PDGDaug);
        if (labelv0 < 0) {
          //both daughters exist as an MC Particle, therfore this can only
          //be a fake v0, and therefore a fake candidate overall
          this->SetParticleOrigin(AliFemtoDreamBasePart::kFake);
        } else {
//          std::cout <<"Potential v0 weakling with ID: " << labelv0 << " found \n";
          AliAODMCParticle* mcv0 = (AliAODMCParticle*) mcarray->At(labelv0);
          if (!mcv0) {
            this->fIsSet = false;
          } else {
            if (mcv0->IsSecondaryFromWeakDecay()
                && !(mcv0->IsSecondaryFromMaterial())) {
              //the v0 candidate has to be from a weak decay itself
              labelv0Mother = mcv0->GetMother();
//              std::cout << "Indeed a v0 weakling his mothers number is " << labelv0Mother << "\n";
//              std::cout << "PDG v0: "<< mcv0->GetPdgCode() << "\n";
//              std::cout << "PDG v0-Mother: "<<((AliAODMCParticle*)mcarray->At(labelv0Mother))->GetPdgCode() << "\n";
              if (labelv0Mother < 0) {
                //Again, this should not happen, just in case somethings
                //very fishy
                this->fIsSet = false;
              } else {
                //now, for a real candidate bachelor and v0 need to have the
                //same mom, NO ADOPTION OR PATCHWORK :(
                //                std::cout << labelv0Mother << '\t' << labelBachMother << '\n';
                if (labelv0Mother == labelBachMother) {
                  //MOMMY?
                  AliAODMCParticle* mcPart = (AliAODMCParticle*) mcarray->At(
                      labelv0Mother);
//                  std::cout << "MOOOOOM: " << mcPart->GetPdgCode() << "\n";
                  this->SetMCPDGCode(mcPart->GetPdgCode());
                  double mcMom[3] = { 0., 0., 0. };
                  mcPart->PxPyPz(mcMom);
                  this->SetMCMomentum(mcMom[0], mcMom[1], mcMom[2]);
                  this->SetMCPt(mcPart->Pt());
                  this->SetMCPhi(mcPart->Phi());
                  this->SetMCTheta(mcPart->Theta());
                  if (mcPart->IsPhysicalPrimary()
                      && !(mcPart->IsSecondaryFromWeakDecay())) {
                    this->SetParticleOrigin(
                        AliFemtoDreamBasePart::kPhysPrimary);
//                    std::cout << "A primary \n";
                  } else if (mcPart->IsSecondaryFromWeakDecay()
                      && !(mcPart->IsSecondaryFromMaterial())) {
                    this->SetParticleOrigin(AliFemtoDreamBasePart::kWeak);
                    this->SetPDGMotherWeak(
                        ((AliAODMCParticle*) mcarray->At(mcPart->GetMother()))
                            ->PdgCode());
//                    std::cout << "A secondary from" << ((AliAODMCParticle*)mcarray->At(mcPart->GetMother()))->PdgCode() << std::endl;
                  } else if (mcPart->IsSecondaryFromMaterial()) {
                    this->SetParticleOrigin(AliFemtoDreamBasePart::kMaterial);
//                    std::cout << "A Material \n";
                  } else {
                    this->SetParticleOrigin(AliFemtoDreamBasePart::kUnknown);
//                    std::cout << "An Unknown \n";
                  }
                } else {
                  //combinatorial background
                  this->SetParticleOrigin(AliFemtoDreamBasePart::kFake);
                }
              }
            } else {
              //if its not from a secondary decay, this cant be a candidate
              this->SetParticleOrigin(AliFemtoDreamBasePart::kFake);
            }
          }
        }
      }
    } else {
      //if the bachlor is not from a weak decay, then it is a ...
      this->SetParticleOrigin(AliFemtoDreamBasePart::kFake);
    }
  } else {
    //this means the bachelor and daughter tracks have no corresponding MC Particle
    this->fIsSet = false;
  }
}

