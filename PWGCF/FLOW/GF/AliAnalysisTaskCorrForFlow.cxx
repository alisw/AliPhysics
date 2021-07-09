/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliAnalysisTaskCorrForFlow.h"

using namespace std;

ClassImp(AliAnalysisTaskCorrForFlow)

AliAnalysisTaskCorrForFlow::AliAnalysisTaskCorrForFlow() : AliAnalysisTaskSE(),
    fAOD(0),
    fOutputListCharged(0),
    fTracksTrigCharged(0),
    fTracksAss(0),
    fPoolMgr(0),
    fhEventCounter(0),
    fHistPhiEta(0),
    fhTrigTracks(0),
    fhChargedSE(0),
    fhChargedME(0),
    fTrigger(AliVEvent::kINT7),
    fIsHMpp(kFALSE),
    fFilterBit(96),
    fbSign(0),
    fPtMinTrig(0.5),
    fPtMaxTrig(10.0),
    fPtMinAss(0.5),
    fPtMaxAss(1.5),
    fCentMin(0.0),
    fCentMax(10.0),
    fCentrality(-10.0),
    fAbsEtaMax(1.0),
    fPVz(100.0),
    fCentEstimator("V0M"),
    fPoolMaxNEvents(2000),
    fPoolMinNTracks(50000),
    fMinEventsToMix(5),
    fNzVtxBins(10),
    fNCentBins(15),
    fMergingCut(0.0)
{}
//_____________________________________________________________________________
AliAnalysisTaskCorrForFlow::AliAnalysisTaskCorrForFlow(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0),
    fOutputListCharged(0),
    fTracksTrigCharged(0),
    fTracksAss(0),
    fPoolMgr(0),
    fhEventCounter(0),
    fHistPhiEta(0),
    fhTrigTracks(0),
    fhChargedSE(0),
    fhChargedME(0),
    fTrigger(AliVEvent::kINT7),
    fIsHMpp(kFALSE),
    fFilterBit(96),
    fbSign(0),
    fPtMinTrig(0.5),
    fPtMaxTrig(10.0),
    fPtMinAss(0.5),
    fPtMaxAss(1.5),
    fCentMin(0.0),
    fCentMax(10.0),
    fCentrality(-10.0),
    fAbsEtaMax(1.0),
    fPVz(100.0),
    fCentEstimator("V0M"),
    fPoolMaxNEvents(2000),
    fPoolMinNTracks(50000),
    fMinEventsToMix(5),
    fNzVtxBins(10),
    fNCentBins(15),
    fMergingCut(0.0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskCorrForFlow::~AliAnalysisTaskCorrForFlow()
{}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlow::UserCreateOutputObjects()
{
    OpenFile(1);

    // //just for testing
    // fPtBinsTrigCharged = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
    // fPtBinsAss = {0.5, 1.0, 1.5, 2.0, 3.0};
    fzVtxBins = {-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0};
    fCentBins = {0,1,2,3,4,5,10,20,30,40,50,60,70,80,90,100};

    fOutputListCharged = new TList();
    fOutputListCharged->SetOwner(kTRUE);

    fhEventCounter = new TH1D("fhEventCounter","Event Counter",10,0,10);
    fOutputListCharged->Add(fhEventCounter);

    fHistPhiEta = new TH2D("fHistPhiEta", "fHistPhiEta; phi; eta", 100, -0.5, 7, 100, -1.5, 1.5);
    fOutputListCharged->Add(fHistPhiEta);

    fhTrigTracks = new TH2D("fhTrigTracks", "fhTrigTracks; pT (trig); PVz", fPtBinsTrigCharged.size() - 1, fPtBinsTrigCharged.data(), 10, -10, 10);
    fOutputListCharged->Add(fhTrigTracks);

    //mixing
    fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins.data(), fNzVtxBins, fzVtxBins.data());
    if (!fPoolMgr) { AliError("Event Pool manager not created!"); return; }
    fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);

    //sparses

    Int_t nSteps = 1;
    Double_t binning_deta_tpctpc[33] = {-1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,    0.1,  0.2,  0.3,  0.4,  0.5, 0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5, 1.6};
    Double_t binning_dphi[73] = { -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0,       0.087266,  0.174533,  0.261799,  0.349066,  0.436332, 0.523599,  0.610865,  0.698132,  0.785398,  0.872665,  0.959931, 1.047198,  1.134464,  1.221730,  1.308997,  1.396263,  1.483530, 1.570796,  1.658063,  1.745329,  1.832596,  1.919862,  2.007129, 2.094395,  2.181662,  2.268928,  2.356194,  2.443461,  2.530727, 2.617994,  2.705260,  2.792527,  2.879793,  2.967060,  3.054326, 3.141593,  3.228859,  3.316126,  3.403392,  3.490659,  3.577925, 3.665191,  3.752458,  3.839724,  3.926991,  4.014257,  4.101524, 4.188790,  4.276057,  4.363323,  4.450590,  4.537856,  4.625123, 4.712389};
    const Int_t sizePtTrig = fPtBinsTrigCharged.size() - 1;
    const Int_t sizePtAss = fPtBinsAss.size() - 1;
    const Int_t iBinningTPCTPC[5] = {32,72,sizePtTrig,sizePtAss, 10};

    fhChargedSE = new AliTHn("fhChargedSE", "fhChargedSE", nSteps, 5, iBinningTPCTPC);
    fhChargedSE->SetBinLimits(0, binning_deta_tpctpc);
    fhChargedSE->SetBinLimits(1, binning_dphi);
    fhChargedSE->SetBinLimits(2, fPtBinsTrigCharged.data());
    fhChargedSE->SetBinLimits(3, fPtBinsAss.data());
    fhChargedSE->SetBinLimits(4, -10,10);
    fhChargedSE->SetVarTitle(0, "#Delta#eta");
    fhChargedSE->SetVarTitle(1, "#Delta#phi");
    fhChargedSE->SetVarTitle(2, "p_{T} [GeV/c] (trig)");
    fhChargedSE->SetVarTitle(3, "p_{T} [GeV/c] (ass)");
    fhChargedSE->SetVarTitle(4, "PVz");
    fOutputListCharged->Add(fhChargedSE);

    fhChargedME = new AliTHn("fhChargedME", "fhChargedME", nSteps, 5, iBinningTPCTPC);
    fhChargedME->SetBinLimits(0, binning_deta_tpctpc);
    fhChargedME->SetBinLimits(1, binning_dphi);
    fhChargedME->SetBinLimits(2, fPtBinsTrigCharged.data());
    fhChargedME->SetBinLimits(3, fPtBinsAss.data());
    fhChargedME->SetBinLimits(4, -10,10);
    fhChargedME->SetVarTitle(0, "#Delta#eta");
    fhChargedME->SetVarTitle(1, "#Delta#phi");
    fhChargedME->SetVarTitle(2, "p_{T} [GeV/c] (trig)");
    fhChargedME->SetVarTitle(3, "p_{T} [GeV/c] (ass)");
    fhChargedME->SetVarTitle(4, "PVz");
    fOutputListCharged->Add(fhChargedME);

    PostData(1, fOutputListCharged);
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlow::UserExec(Option_t *)
{
    fhEventCounter->Fill("Input",1);

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) { AliError("Event not loaded."); return; }
    if(!IsEventSelected()) { return; }

    fTracksTrigCharged = new TObjArray;
    fTracksAss = new TObjArray;

    Int_t iTracks(fAOD->GetNumberOfTracks());
    if(iTracks < 1 ) {
      fTracksTrigCharged->Clear();
  	  delete fTracksTrigCharged;
      fTracksAss->Clear();
  	  delete fTracksAss;
      AliWarning("No tracks in the event.");
      return;
    }

    for(Int_t i(0); i < iTracks; i++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }

        Double_t trackPt = track->Pt();
        if(trackPt > fPtMinAss && trackPt < fPtMaxAss) fTracksAss->Add((AliAODTrack*)track);
        if(trackPt > fPtMinTrig && trackPt < fPtMaxTrig) {
          fTracksTrigCharged->Add((AliAODTrack*)track);
          fhTrigTracks->Fill(trackPt, fPVz);
        }

        //example histogram
        fHistPhiEta->Fill(track->Phi(), track->Eta());
    }

    if(!fTracksTrigCharged->IsEmpty()){
      FillCorrelations();
      FillCorrelationsMixed();
    }

    fTracksTrigCharged->Clear();
	  delete fTracksTrigCharged;

    fTracksAss->Clear();
	  delete fTracksAss;

    PostData(1, fOutputListCharged);
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlow::Terminate(Option_t *)
{
   if (fPoolMgr) delete fPoolMgr;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlow::IsEventSelected()
{
  fhEventCounter->Fill("EventOK",1);

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();
  if(!(fSelectMask & fTrigger)) { return kFALSE; }
  fhEventCounter->Fill("TriggerOK",1);

  if(fIsHMpp) fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true);
  if(!fEventCuts.AcceptEvent(fAOD)) { return kFALSE; }
  fhEventCounter->Fill("CutsOK",1);

  AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
  if(!multSelection) { return kFALSE; }
  fhEventCounter->Fill("MultOK",1);
  Float_t dPercentile = multSelection->GetMultiplicityPercentile(fCentEstimator);
  if(dPercentile > 100 || dPercentile < 0) { return kFALSE; }
  fhEventCounter->Fill("PercOK",1);

  if(fCentMax > 0.0 && (dPercentile < fCentMin || dPercentile > fCentMax)) { return kFALSE; }
  fhEventCounter->Fill("CentOK",1);
  fCentrality = (Double_t) dPercentile;

  fPVz = fAOD->GetPrimaryVertex()->GetZ();
  if(TMath::Abs(fPVz) >= 10.0) { return kFALSE; }
  fhEventCounter->Fill("PVzOK",1);

  fbSign = (InputEvent()->GetMagneticField() > 0) ? 1 : -1;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlow::IsTrackSelected(const AliAODTrack* track) const
{
  if(!track->TestFilterBit(fFilterBit)) { return kFALSE; }
  if(track->GetTPCNcls() < 70 && fFilterBit != 2) { return kFALSE; }
  if(fAbsEtaMax > 0.0 && TMath::Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
  if(track->Charge() == 0) { return kFALSE; }

  return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlow::RangePhi(Double_t DPhi) {
  if (DPhi < -TMath::Pi() / 2)   DPhi += 2 * TMath::Pi();
  if (DPhi > 3 * TMath::Pi() / 2) DPhi -= 2*TMath::Pi();
  return DPhi;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlow::GetDPhiStar(Double_t phi1, Double_t pt1, Double_t charge1, Double_t phi2, Double_t pt2, Double_t charge2, Double_t radius){
  // calculates delta phi *
  Double_t dPhiStar = phi1 - phi2 - charge1 * fbSign * TMath::ASin(0.075 * radius / pt1) + charge2 * fbSign * TMath::ASin(0.075 * radius / pt2);

  if (dPhiStar > TMath::Pi()) dPhiStar = 2.0*TMath::Pi() - dPhiStar;
  if (dPhiStar < -TMath::Pi()) dPhiStar = -2.0*TMath::Pi() - dPhiStar;

  return dPhiStar;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlow::FillCorrelations()
{
  if(!fhTrigTracks || !fTracksTrigCharged || !fTracksAss) { AliError("Necessary inputs missing, terminating!"); return; }

  Double_t binscont[5];
  binscont[4] = fPVz;

  for(Int_t iTrig(0); iTrig < fTracksTrigCharged->GetEntriesFast(); iTrig++){
    AliAODTrack* track = (AliAODTrack*)fTracksTrigCharged->At(iTrig);
    if(!track) continue;

    Double_t trigPt = track->Pt();
    Double_t trigEta = track->Eta();
    Double_t trigPhi = track->Phi();
    Double_t trigCharge = track->Charge();

    for(Int_t iAss(0); iAss < fTracksAss->GetEntriesFast(); iAss++){
      AliAODTrack* trackAss = (AliAODTrack*)fTracksAss->At(iAss);
      if(!trackAss) continue;

      Double_t assPt = trackAss->Pt();
      Double_t assEta = trackAss->Eta();
      Double_t assPhi = trackAss->Phi();
      Double_t assCharge = trackAss->Charge();

      if(trigPt < assPt) continue;
      if(track->GetID() == trackAss->GetID()) continue;

      binscont[0] = assEta - trigEta;
      binscont[1] = RangePhi(assPhi - trigPhi);

      if(TMath::Abs(binscont[0]) < fMergingCut){
        Double_t dPhiStarLow = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, 0.8);
        Double_t dPhiStarHigh = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, 2.5);

        if(TMath::Abs(dPhiStarLow) < fMergingCut || TMath::Abs(dPhiStarHigh) < fMergingCut) continue;

        const Double_t kLimit = 3.0*fMergingCut;

        if(TMath::Abs(dPhiStarLow) < kLimit || TMath::Abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0 ) {
          Bool_t bIsBelow = kFALSE;
          for(Double_t rad(0.8); rad < 2.51; rad+=0.01){
            Double_t dPhiStar = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, rad);
            if(TMath::Abs(dPhiStar) < fMergingCut) {
              bIsBelow = kTRUE;
              break;
            }
          } // end loop radius
          if(bIsBelow) continue;
        }
      }

      binscont[2] = trigPt;
      binscont[3] = assPt;

      fhChargedSE->Fill(binscont,0);
    }

  }

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlow::FillCorrelationsMixed()
{
  if(!fhTrigTracks || !fTracksTrigCharged || !fTracksAss) { AliError("Necessary inputs missing, terminating!"); return; }

  AliEventPool *pool = fPoolMgr->GetEventPool(fCentrality, fPVz);
  if(!pool) { AliError(Form("No pool found for centrality = %f, zVtx = %f", fCentrality,fPVz)); return; }

  if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() > fMinEventsToMix) {
    Int_t nMix = pool->GetCurrentNEvents();

    Double_t binscont[5];
    binscont[4] = fPVz;

    for(Int_t iTrig(0); iTrig < fTracksTrigCharged->GetEntriesFast(); iTrig++){
      AliAODTrack* track = (AliAODTrack*)fTracksTrigCharged->At(iTrig);
      if(!track) continue;

      Double_t trigPt = track->Pt();
      Double_t trigEta = track->Eta();
      Double_t trigPhi = track->Phi();
      Double_t trigCharge = track->Charge();

      for(Int_t eMix(0); eMix < nMix; eMix++){
        TObjArray *mixEvents = pool->GetEvent(eMix);
        for(Int_t iAss(0); iAss < mixEvents->GetEntriesFast(); iAss++){
          AliAODTrack* trackAss = (AliAODTrack*)mixEvents->At(iAss);
          if(!trackAss) continue;

          Double_t assPt = trackAss->Pt();
          Double_t assEta = trackAss->Eta();
          Double_t assPhi = trackAss->Phi();
          Double_t assCharge = trackAss->Charge();

          if(trigPt < assPt) continue;

          binscont[0] = assEta - trigEta;
          binscont[1] = RangePhi(assPhi - trigPhi);

          if(TMath::Abs(binscont[0]) < fMergingCut){
            Double_t dPhiStarLow = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, 0.8);
            Double_t dPhiStarHigh = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, 2.5);

            if(TMath::Abs(dPhiStarLow) < fMergingCut || TMath::Abs(dPhiStarHigh) < fMergingCut) continue;

            const Double_t kLimit = 3.0*fMergingCut;

            if(TMath::Abs(dPhiStarLow) < kLimit || TMath::Abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0 ) {
              Bool_t bIsBelow = kFALSE;
              for(Double_t rad(0.8); rad < 2.51; rad+=0.01){
                Double_t dPhiStar = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, rad);
                if(TMath::Abs(dPhiStar) < fMergingCut) {
                  bIsBelow = kTRUE;
                  break;
                }
              } // end loop radius
              if(bIsBelow) continue;
            }
          }

          binscont[2] = trigPt;
          binscont[3] = assPt;

          fhChargedME->Fill(binscont,0,1./(Double_t)nMix);
        }
      }
    }
  }

  TObjArray* cloneArray = (TObjArray *)fTracksAss->Clone();
  cloneArray->SetOwner(kTRUE);
  pool->UpdatePool(cloneArray);

  return;
}
