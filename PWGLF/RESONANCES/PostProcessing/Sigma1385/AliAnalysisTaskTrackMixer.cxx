/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnalysisTaskTrackMixer
 *
 *  Event Mixer for track
 *  This task will save the pion track in multiple event
 *  for efficient event mixing.
 *
 *  Author: Bong-Hwi Lim
 *
 */

#include <TDatabasePDG.h>
#include <math.h>

#include <iostream>

#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliMultSelectionTask.h"
#include "AliPIDResponse.h"
#include "TChain.h"

// for NanoAOD
#include <AliNanoAODHeader.h>
#include <AliNanoAODTrack.h>

#include "AliAnalysisTaskTrackMixer.h"
#include "THistManager.h"

class AliAnalysisTaskTrackMixer;

ClassImp(AliAnalysisTaskTrackMixer)
    AliAnalysisTaskTrackMixer::AliAnalysisTaskTrackMixer()
    : AliAnalysisTaskSE(),
      fTrackCuts(nullptr),
      fEventCuts(),
      fCheckTPCGeo(kFALSE),
      fTPCActiveLengthCutDeltaY(3.0),
      fTPCActiveLengthCutDeltaZ(220.0),
      fRequireCutGeoNcrNclLength(130),
      fRequireCutGeoNcrNclGeom1Pt(1.5),
      fCutGeoNcrNclFractionNcr(0.85),
      fCutGeoNcrNclFractionNcl(0.7),
      fPIDResponse(0x0),
      fEvt(nullptr),
      fHistos(nullptr),
      fVertex(nullptr),
      fIsAOD(kFALSE),
      fIsNano(kFALSE),
      fFillQAPlot(kTRUE),
      fIsHM(kFALSE),
      fEMpool(0),
      fBinCent(),
      fBinZ(),
      fPosPV(),
      fMagField(0),
      fCent(-1),
      fnMix(10),
      fCentBin(-1),
      fZbin(-1),
      fTrackFilterBit(32.0),
      fTrackTPCNsigCut(3.5),
      fTrackEtaCut(0.8),
      fTrackZVertexCut(2.2),
      fTrackXYVertsigCut(8.0),
      fTPCPIDType(AliPID::kPion),
      fGoodTrackArray(0x0) {
  /// Default constructor
}
//_____________________________________________________________________________
AliAnalysisTaskTrackMixer::AliAnalysisTaskTrackMixer(const char* name)
    : AliAnalysisTaskSE(name),
      fTrackCuts(nullptr),
      fEventCuts(0x0),
      fCheckTPCGeo(kFALSE),
      fTPCActiveLengthCutDeltaY(3.0),
      fTPCActiveLengthCutDeltaZ(220.0),
      fRequireCutGeoNcrNclLength(130),
      fRequireCutGeoNcrNclGeom1Pt(1.5),
      fCutGeoNcrNclFractionNcr(0.85),
      fCutGeoNcrNclFractionNcl(0.7),
      fPIDResponse(0x0),
      fEvt(nullptr),
      fHistos(nullptr),
      fVertex(nullptr),
      fIsAOD(kFALSE),
      fIsNano(kFALSE),
      fFillQAPlot(kTRUE),
      fIsHM(kFALSE),
      fEMpool(0),
      fBinCent(),
      fBinZ(),
      fPosPV(),
      fMagField(0),
      fCent(-1),
      fnMix(10),
      fCentBin(-1),
      fZbin(-1),
      fTrackFilterBit(32),
      fTrackTPCNsigCut(3.5),
      fTrackEtaCut(0.8),
      fTrackZVertexCut(2.2),
      fTrackXYVertsigCut(8.0),
      fTPCPIDType(AliPID::kPion),
      fGoodTrackArray(0x0) {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskTrackMixer::~AliAnalysisTaskTrackMixer() {}
//_____________________________________________________________________________
void AliAnalysisTaskTrackMixer::UserCreateOutputObjects() {
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

  fHistos = new THistManager("TackMixerQAhists");
  std::vector<double> centaxisbin;
  (fIsHM) ? centaxisbin = {0, 0.001, 0.01, 0.05, 0.1}
          : centaxisbin = {
                -1, 0,  1,  5,  10, 15, 20, 30,
                40, 50, 60, 70, 80, 90, 100};  // can be use from pp to PbPb
  fBinCent = AxisVar("Cent", centaxisbin);
  fBinZ = AxisVar("Z", {-10, -5, -3, -1, 1, 3, 5, 10});  // moderate diff

  fEventCuts.AddQAplotsToList(fHistos->GetListOfHistograms());
  if (fIsHM)
    fHistos->CreateTH1("hMultiplicity", "", 100, 0, 0.1, "s");
  else
    fHistos->CreateTH1("hMultiplicity", "", 101, -1, 100, "s");
  if (fFillQAPlot) {
    fHistos->CreateTH2("QA/hTPCPIDPion", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH1("QA/hEtaPion", "", 20, -1.0, 1.0);
    fHistos->CreateTH1("QA/hDCAPVPion", "", 30, 0, 3, "s");
    fHistos->CreateTH1("QA/hDCArPVPion", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QA/hPtPion", "", 200, 0, 20);

    fHistos->CreateTH2("QAcut/hTPCPIDPion", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH1("QAcut/hEtaPion", "", 20, -1.0, 1.0);
    fHistos->CreateTH1("QAcut/hDCAPVPion", "", 30, 0, 3, "s");
    fHistos->CreateTH1("QAcut/hDCArPVPion", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QAcut/hPtPion", "", 200, 0, 20);
  }
  fEMpool.resize(fBinCent.GetNbins() + 1,
                 std::vector<eventpool>(fBinZ.GetNbins() + 1));

  PostData(1, fHistos->GetListOfHistograms());
}
//_____________________________________________________________________________
void AliAnalysisTaskTrackMixer::UserExec(Option_t*) {
  AliVEvent* event = InputEvent();
  if (!event) {
    PostData(1, fHistos->GetListOfHistograms());
    AliInfo("Could not retrieve event");
    return;
  }
  AliNanoAODHeader* nanoHeader =
      dynamic_cast<AliNanoAODHeader*>(fInputEvent->GetHeader());

  event->IsA() == AliESDEvent::Class()
      ? fEvt = dynamic_cast<AliESDEvent*>(event)
      : fEvt = dynamic_cast<AliAODEvent*>(event);
  if (!fIsAOD && (event->IsA() != AliESDEvent::Class()))
    fIsAOD = true;
  if (!fEvt) {
    PostData(1, fHistos->GetListOfHistograms());
    return;
  }

  bool IsEvtSelected{false};
  if (!nanoHeader) {
    IsEvtSelected = fEventCuts.AcceptEvent(event);
    fCent = AliMultSelectionTask::IsINELgtZERO(event)
                ? fEventCuts.GetCentrality()
                : -0.5;
    fPIDResponse = (AliPIDResponse*)(AliAnalysisManager::GetAnalysisManager()
                                         ->GetInputEventHandler())
                       ->GetPIDResponse();
    if (!fPIDResponse)
      AliInfo("No PIDd");
  } else {
    if (!fIsNano)
      fIsNano = kTRUE;
    IsEvtSelected = true;
    fCent = nanoHeader->GetCentr("V0M");
    static int inel_index = -1;
    if (inel_index < 0)
      inel_index = nanoHeader->GetVarIndex("cstINELgt0");
    if ((inel_index > 0) && (nanoHeader->GetVar(inel_index) < 0.5))
      fCent = -0.5;
  }

  if (!IsEvtSelected) {
    PostData(1, fHistos->GetListOfHistograms());
    return;  // event cut
  }

  fHistos->FillTH1("hMultiplicity", (double)fCent);

  if (fIsAOD)
    fVertex = ((AliAODEvent*)fEvt)->GetPrimaryVertex();
  const AliVVertex* pVtx = fEvt->GetPrimaryVertex();
  fPosPV[0] = pVtx->GetX();
  fPosPV[1] = pVtx->GetY();
  fPosPV[2] = pVtx->GetZ();
  fMagField = fEvt->GetMagneticField();

  // Event Mixing pool -----------------------------------------------------
  fZbin = fBinZ.FindBin(fPosPV[2]) - 1;    // Event mixing z-bin
  fCentBin = fBinCent.FindBin(fCent) - 1;  // Event mixing cent bin

  GoodTracksSelection();
  if (fGoodTrackArray.size()) {
    FillTrackToEventPool();  // use only pion track pool.
  }
  PostData(1, fHistos->GetListOfHistograms());
}
//_____________________________________________________________________________
void AliAnalysisTaskTrackMixer::Terminate(Option_t*) {}
//_____________________________________________________________________________
void AliAnalysisTaskTrackMixer::GoodTracksSelection() {
  const UInt_t nTracks = fEvt->GetNumberOfTracks();
  fGoodTrackArray.clear();
  AliVTrack* track;
  Float_t b[2];
  Float_t bCov[3];
  Double_t nTPCNSigPion, pionZ, pionPt, pionSigmaDCA_r, pionDCA_r, lEta,
      lTPCmomentum, lTPCsignal;
  Int_t isTPCGeo;

  for (UInt_t it = 0; it < nTracks; it++) {
    track = (AliVTrack*)fEvt->GetTrack(it);
    if (!track)
      continue;

    // ---------- Track selection begin ----------
    if (!fIsAOD) {
      if (!fTrackCuts->AcceptTrack((AliESDtrack*)track))
        continue;
      if (fCheckTPCGeo)
        isTPCGeo = IsSelectedTPCGeoCut(((AliESDtrack*)track)) ? 1 : 0;
    }  // ESD Case
    else {
      if (!fIsNano) {
        if (!((AliAODTrack*)track)->TestFilterBit(fTrackFilterBit))
          continue;
        if (fCheckTPCGeo)
          isTPCGeo = IsSelectedTPCGeoCut(((AliAODTrack*)track)) ? 1 : 0;
      } else {
        if (!(static_cast<AliNanoAODTrack*>(track)->TestFilterBit(
                fTrackFilterBit)))
          continue;
        if (fCheckTPCGeo) {
          static const Int_t tpcGeo_index =
              AliNanoAODTrackMapping::GetInstance()->GetVarIndex(
                  "cstTPCGeoLength");
          isTPCGeo =
              (static_cast<AliNanoAODTrack*>(track)->GetVar(tpcGeo_index) > 0.5)
                  ? 1
                  : 0;
        }
      }
    }
    // Values
    GetImpactParam(track, b, bCov);
    pionZ = b[1];
    nTPCNSigPion = GetTPCnSigma(track, fTPCPIDType);
    pionPt = track->Pt();
    pionSigmaDCA_r = (0.0026 + 0.0050 / pionPt);
    pionDCA_r = b[0];
    lEta = track->Eta();
    lTPCmomentum = track->GetTPCmomentum();
    lTPCsignal = track->GetTPCsignal();

    if (fFillQAPlot) {
      fHistos->FillTH1("QA/hDCAPVPion", pionZ);
      fHistos->FillTH1("QA/hDCArPVPion", pionDCA_r);
      fHistos->FillTH1("QA/hEtaPion", lEta);
      fHistos->FillTH1("QA/hPtPion", pionPt);
      fHistos->FillTH2("QA/hTPCPIDPion", lTPCmomentum, lTPCsignal);
    }

    // Track Cut
    if (pionZ > fTrackZVertexCut)
      continue;
    if (TMath::Abs(nTPCNSigPion) > fTrackTPCNsigCut)
      continue;
    if (TMath::Abs(lEta) > fTrackEtaCut)
      continue;
    if (pionPt < 0.15)
      continue;
    if (pionDCA_r > pionSigmaDCA_r * fTrackXYVertsigCut)
      continue;
    
    // Save
    fGoodTrackArray.push_back(it);
    if (fFillQAPlot) {
      fHistos->FillTH1("QAcut/hDCAPVPion", pionZ);
      fHistos->FillTH1("QAcut/hDCArPVPion", pionDCA_r);
      fHistos->FillTH1("QAcut/hEtaPion", lEta);
      fHistos->FillTH1("QAcut/hPtPion", pionPt);
      fHistos->FillTH2("QAcut/hTPCPIDPion", lTPCmomentum, lTPCsignal);
    }
  }
}
void AliAnalysisTaskTrackMixer::FillTrackToEventPool() {
  // Fill Selected tracks to event mixing pool
  if ((fCentBin < 0) || (fZbin < 0))
    return;
  AliVTrack* goodtrack;

  tracklist* etl;
  eventpool* ep;
  // Event mixing pool

  ep = &fEMpool[fCentBin][fZbin];
  ep->push_back(tracklist());
  etl = &(ep->back());
  // Fill selected tracks
  for (UInt_t i = 0; i < fGoodTrackArray.size(); i++) {
    goodtrack = (AliVTrack*)fEvt->GetTrack(fGoodTrackArray[i]);
    if (!goodtrack)
      continue;
    etl->push_back((AliVTrack*)goodtrack->Clone());
  }
  if (!fGoodTrackArray.size())
    ep->pop_back();
  if ((int)ep->size() > (int)fnMix) {
    for (auto it : ep->front())
      delete it;
    ep->pop_front();
  }
}
TAxis AliAnalysisTaskTrackMixer::AxisFix(TString name,
                                             int nbin,
                                             Double_t xmin,
                                             Double_t xmax) {
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  TAxis axis(nbin, xmin, xmax);
  axis.SetName(name);
  return axis;
}
TAxis AliAnalysisTaskTrackMixer::AxisVar(TString name,
                                             std::vector<Double_t> bin) {
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  TAxis axis(bin.size() - 1, &bin.front());
  axis.SetName(name);
  return axis;
}
double AliAnalysisTaskTrackMixer::GetTPCnSigma(AliVTrack* track,
                                                   AliPID::EParticleType type) {
  AliNanoAODTrack* nanoT = dynamic_cast<AliNanoAODTrack*>(track);
  if (nanoT) {
    static bool used = false;
    if (!used) {
      AliNanoAODTrack::InitPIDIndex();
      used = true;
    }
    return nanoT->GetVar(
        AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, type));
  } else
    return fPIDResponse->NumberOfSigmasTPC(track, type);
}
void AliAnalysisTaskTrackMixer::GetImpactParam(AliVTrack* track,
                                                   Float_t p[2],
                                                   Float_t cov[3]) {
  AliNanoAODTrack* nanoT = dynamic_cast<AliNanoAODTrack*>(track);
  if (nanoT)
    nanoT->AliNanoAODTrack::GetImpactParameters(p[0], p[1]);
  else
    track->GetImpactParameters(p, cov);
}
Bool_t AliAnalysisTaskTrackMixer::IsSelectedTPCGeoCut(AliAODTrack* track) {
  Bool_t checkResult = kTRUE;

  AliESDtrack esdTrack(track);
  esdTrack.SetTPCClusterMap(track->GetTPCClusterMap());
  esdTrack.SetTPCSharedMap(track->GetTPCSharedMap());
  esdTrack.SetTPCPointsF(track->GetTPCNclsF());

  auto nCrossedRowsTPC = esdTrack.GetTPCCrossedRows();
  auto lengthInActiveZoneTPC = esdTrack.GetLengthInActiveZone(
      0, fTPCActiveLengthCutDeltaY, fTPCActiveLengthCutDeltaZ, fMagField);
  auto cutGeoNcrNclLength = fRequireCutGeoNcrNclLength -
                            TMath::Power(TMath::Abs(esdTrack.GetSigned1Pt()),
                                         fRequireCutGeoNcrNclGeom1Pt);

  if (lengthInActiveZoneTPC < cutGeoNcrNclLength)
    checkResult = false;
  if (nCrossedRowsTPC < fCutGeoNcrNclFractionNcr * cutGeoNcrNclLength)
    checkResult = false;
  if (esdTrack.GetTPCncls() < fCutGeoNcrNclFractionNcl * cutGeoNcrNclLength)
    checkResult = false;

  return checkResult;
}
Bool_t AliAnalysisTaskTrackMixer::IsSelectedTPCGeoCut(AliESDtrack* track) {
  Bool_t checkResult = kTRUE;

  auto nCrossedRowsTPC = track->GetTPCCrossedRows();
  auto lengthInActiveZoneTPC = track->GetLengthInActiveZone(
      0, fTPCActiveLengthCutDeltaY, fTPCActiveLengthCutDeltaZ, fMagField);
  auto cutGeoNcrNclLength = fRequireCutGeoNcrNclLength -
                            TMath::Power(TMath::Abs(track->GetSigned1Pt()),
                                         fRequireCutGeoNcrNclGeom1Pt);

  if (lengthInActiveZoneTPC < cutGeoNcrNclLength)
    checkResult = false;
  if (nCrossedRowsTPC < fCutGeoNcrNclFractionNcr * cutGeoNcrNclLength)
    checkResult = false;
  if (track->GetTPCncls() < fCutGeoNcrNclFractionNcl * cutGeoNcrNclLength)
    checkResult = false;

  return checkResult;
}