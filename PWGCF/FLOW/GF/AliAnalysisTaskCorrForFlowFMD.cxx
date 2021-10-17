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

#include "AliAnalysisTaskCorrForFlowFMD.h"

using namespace std;

ClassImp(AliAnalysisTaskCorrForFlowFMD);

AliAnalysisTaskCorrForFlowFMD::AliAnalysisTaskCorrForFlowFMD() : AliAnalysisTaskSE(),
    fAOD(0),
    fOutputListCharged(0),
    fInputListEfficiency(0),
    fTracksAss(0),
    fPIDResponse(0),
    fPIDCombined(0),
    fPoolMgr(0),
    fhEventCounter(0),
    fhEventMultiplicity(0),
    fHistPhiEta(0),
    fAnalType(eFMDAFMDC),
    fTrigger(AliVEvent::kINT7),
    fIsHMpp(kFALSE),
    fDoPID(kFALSE),
    fUseNch(kFALSE),
    fUseEfficiency(kFALSE),
    fEfficiencyEtaDependent(kFALSE),
    fUseFMDcut(kTRUE),
    fFilterBit(96),
    fbSign(0),
    fNofTracks(0),
    fNchMin(0),
    fNchMax(100000),
    fPtMinTrig(0.5),
    fPtMaxTrig(10.0),
    fPtMinAss(0.5),
    fPtMaxAss(1.5),
    fFMDcutapar0(1.64755),
    fFMDcutapar1(119.602),
    fFMDcutcpar0(2.73426),
    fFMDcutcpar1(150.31),
    fCentMin(0.0),
    fCentMax(10.0),
    fCentrality(-10.0),
    fAbsEtaMax(0.8),
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
AliAnalysisTaskCorrForFlowFMD::AliAnalysisTaskCorrForFlowFMD(const char* name, Bool_t bUseEff) : AliAnalysisTaskSE(name),
    fAOD(0),
    fOutputListCharged(0),
    fInputListEfficiency(0),
    fTracksAss(0),
    fPIDResponse(0),
    fPIDCombined(0),
    fPoolMgr(0),
    fhEventCounter(0),
    fhEventMultiplicity(0),
    fHistPhiEta(0),
    fAnalType(eFMDAFMDC),
    fTrigger(AliVEvent::kINT7),
    fIsHMpp(kFALSE),
    fDoPID(kFALSE),
    fUseNch(kFALSE),
    fUseEfficiency(bUseEff),
    fEfficiencyEtaDependent(kFALSE),
    fUseFMDcut(kTRUE),
    fFilterBit(96),
    fbSign(0),
    fNofTracks(0),
    fNchMin(0),
    fNchMax(100000),
    fPtMinTrig(0.5),
    fPtMaxTrig(10.0),
    fPtMinAss(0.5),
    fPtMaxAss(1.5),
    fFMDcutapar0(1.64755),
    fFMDcutapar1(119.602),
    fFMDcutcpar0(2.73426),
    fFMDcutcpar1(150.31),
    fCentMin(0.0),
    fCentMax(10.0),
    fCentrality(-10.0),
    fAbsEtaMax(0.8),
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
    if(bUseEff) { DefineInput(1, TList::Class()); }
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskCorrForFlowFMD::~AliAnalysisTaskCorrForFlowFMD()
{}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowFMD::UserCreateOutputObjects()
{
    OpenFile(1);
    PrintSetup();

    if(fAnalType == eFMDAFMDC && fDoPID) { AliWarning("PID on when running FMDA-FMDC. Turning off."); fDoPID = kFALSE; }

    fzVtxBins = {-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0};

    fOutputListCharged = new TList();
    fOutputListCharged->SetOwner(kTRUE);

    fhEventCounter = new TH1D("fhEventCounter","Event Counter",10,0,10);
    fOutputListCharged->Add(fhEventCounter);

    fhEventMultiplicity = new TH1D("fhEventMultiplicity","Event multiplicity; N_{ch}",200,0,200);
    fOutputListCharged->Add(fhEventMultiplicity);

    TString pidName[4] = {"", "_Pion", "_Kaon", "_Proton"};
    for(Int_t i(0); i < 4; i++){
      if(fAnalType != eFMDAFMDC) fhTrigTracks[i] = new TH2D(Form("fhTrigTracks%s",pidName[i].Data()), Form("fhTrigTracks (%s); pT (trig); PVz",pidName[i].Data()), fPtBinsTrigCharged.size() - 1, fPtBinsTrigCharged.data(), 10, -10, 10);
      else fhTrigTracks[i] = new TH2D(Form("fhTrigTracks%s",pidName[i].Data()), Form("fhTrigTracks (%s); #eta; PVz",pidName[i].Data()), 15, 1, 4, 10, -10, 10);
      fOutputListCharged->Add(fhTrigTracks[i]);
    }

    if(fDoPID){
      // PID response
      AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler* inputHandler = (AliInputEventHandler*)mgr->GetInputEventHandler();
      fPIDResponse = inputHandler->GetPIDResponse();
      if(!fPIDResponse) { AliError("AliPIDResponse not found!"); return; }

      fPIDCombined = new AliPIDCombined();
      fPIDCombined->SetDefaultTPCPriors();
      fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF); // setting TPC + TOF mask
    }

    //mixing
    fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins.data(), fNzVtxBins, fzVtxBins.data());
    if (!fPoolMgr) { AliError("Event Pool manager not created!"); return; }
    fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);

    CreateTHnCorrelations();

    if(fUseEfficiency) {
      fInputListEfficiency = (TList*) GetInputData(1);
      if(fEfficiencyEtaDependent && fAbsEtaMax > 0.8) AliWarning("Efficiency loading -- eta can be out of range!");
    }

    PostData(1, fOutputListCharged);
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowFMD::UserExec(Option_t *)
{
    fhEventCounter->Fill("Input",1);

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) { AliError("Event not loaded."); return; }
    if(!IsEventSelected()) { return; }

    Int_t iTracks(fAOD->GetNumberOfTracks());
    if(iTracks < 1 ) {
      AliWarning("No tracks in the event.");
      return;
    }

    fTracksAss = new TObjArray;
    fTracksTrig[0] = new TObjArray;

    if(fUseEfficiency && !AreEfficienciesLoaded()) { return; }

    if(fAnalType != eTPCTPC) {
      if(!PrepareFMDTracks()){
        delete fTracksAss;
        delete fTracksTrig[0];
        PostData(1, fOutputListCharged);
        return;
      }
    }

    if(fDoPID){
      for(Int_t i(1); i < 4; i++){
        fTracksTrig[i] = new TObjArray;
      }
    }

    fNofTracks = 0;
    for(Int_t i(0); i < iTracks; i++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }

        Double_t trackPt = track->Pt();
        if(trackPt > fPtMinAss && trackPt < fPtMaxAss) {
          if(fAnalType == eTPCTPC) fTracksAss->Add((AliAODTrack*)track); // only if associated from TPC
          fNofTracks++;
        }
        if(fAnalType != eFMDAFMDC){
          if(trackPt > fPtMinTrig && trackPt < fPtMaxTrig) {
            fTracksTrig[0]->Add((AliAODTrack*)track);
            fhTrigTracks[0]->Fill(trackPt, fPVz);

            if(fDoPID){
              Int_t trackPid = IdentifyTrack(track);
              if(trackPid > 0 && trackPid < 4){
                fTracksTrig[trackPid]->Add((AliAODTrack*)track);
                fhTrigTracks[trackPid]->Fill(trackPt, fPVz);
              }
            }
          }
        } // POI from TPC
    } // tracks loop end
    fhEventMultiplicity->Fill(fNofTracks);

    if(fUseNch){
      if(fNofTracks < fNchMin || fNofTracks > fNchMax) { return; }
      fhEventCounter->Fill("Nch cut ok ",1);
    }

    if(!fTracksAss->IsEmpty()){
      for(Int_t i(0); i < 4; i++){
        FillCorrelations(i);
        FillCorrelationsMixed(i);

        fTracksTrig[i]->Clear();
        delete fTracksTrig[i];

        if(!fDoPID) break;
      }
    }

    fTracksAss->Clear();
	  delete fTracksAss;

    PostData(1, fOutputListCharged);
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowFMD::Terminate(Option_t *)
{
   if (fPoolMgr) delete fPoolMgr;
   if(fOutputListCharged) delete fOutputListCharged;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowFMD::IsEventSelected()
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
Bool_t AliAnalysisTaskCorrForFlowFMD::IsTrackSelected(const AliAODTrack* track) const
{
  if(!track->TestFilterBit(fFilterBit)) { return kFALSE; }
  if(track->GetTPCNcls() < 70 && fFilterBit != 2) { return kFALSE; }
  if(fAbsEtaMax > 0.0 && TMath::Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
  if(track->Charge() == 0) { return kFALSE; }

  return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlowFMD::RangePhi(Double_t DPhi) {
  if (DPhi < -TMath::Pi() / 2)   DPhi += 2 * TMath::Pi();
  if (DPhi > 3 * TMath::Pi() / 2) DPhi -= 2*TMath::Pi();
  return DPhi;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlowFMD::GetDPhiStar(Double_t phi1, Double_t pt1, Double_t charge1, Double_t phi2, Double_t pt2, Double_t charge2, Double_t radius){
  // calculates delta phi *
  Double_t dPhiStar = phi1 - phi2 - charge1 * fbSign * TMath::ASin(0.075 * radius / pt1) + charge2 * fbSign * TMath::ASin(0.075 * radius / pt2);

  if (dPhiStar > TMath::Pi()) dPhiStar = 2.0*TMath::Pi() - dPhiStar;
  if (dPhiStar < -TMath::Pi()) dPhiStar = -2.0*TMath::Pi() - dPhiStar;

  return dPhiStar;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskCorrForFlowFMD::IdentifyTrack(const AliAODTrack* track) const
{
  // checking detector statuses
  Bool_t bIsTPCok = HasTrackPIDTPC(track);
  Bool_t bIsTOFok = HasTrackPIDTOF(track);

  if(!bIsTPCok) { return -1; }

  Double_t l_Probs[AliPID::kSPECIES];
  Double_t l_MaxProb[] = {0.95,0.85,0.85};
  Bool_t l_TOFUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, l_Probs) & AliPIDResponse::kDetTOF;
  Int_t pidInd = 0;
  for(Int_t i(0); i < AliPID::kSPECIES; i++) pidInd=(l_Probs[i]>l_Probs[pidInd])?i:pidInd;
  Int_t retInd = pidInd-AliPID::kPion+1; //Not interested in e+mu, so realign to 0 -> adding h as 0
  if(retInd<1 || retInd>3) return -1; //Shouldn't be larger than 2, but just to be safe
  if(l_Probs[pidInd] < l_MaxProb[retInd]) return -1;
  //check nsigma cuts
  if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)pidInd))>3) return -1;
  if(bIsTOFok && l_TOFUsed) if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)pidInd))>3) return -1;

  return retInd;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowFMD::HasTrackPIDTPC(const AliAODTrack* track) const
{
  if(!track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
  return (pidStatusTPC == AliPIDResponse::kDetPidOk);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowFMD::HasTrackPIDTOF(const AliAODTrack* track) const
{
  if(!track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
  return ((pidStatusTOF == AliPIDResponse::kDetPidOk) && (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME));
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowFMD::FillCorrelations(const Int_t spec)
{
  if(!fTracksTrig[spec] || !fhTrigTracks[spec] || !fTracksAss) { AliError("Necessary inputs missing, terminating!"); return; }
  if(!fhSE[spec]) { AliError(Form("Output AliTHn missing for %d , terminating!", spec)); return; }

  if(fAnalType == eTPCTPC){
    Double_t binscont[5];
    binscont[2] = fPVz;

    for(Int_t iTrig(0); iTrig < fTracksTrig[spec]->GetEntriesFast(); iTrig++){
      AliAODTrack* track = (AliAODTrack*)fTracksTrig[spec]->At(iTrig);
      if(!track) continue;

      Double_t trigPt = track->Pt();
      Double_t trigEta = track->Eta();
      Double_t trigPhi = track->Phi();
      Double_t trigCharge = track->Charge();
      Double_t trigEff = 1.0;
      if(fUseEfficiency) trigEff = GetEff(trigPt, spec, trigEta);
      binscont[3] = trigPt;

      for(Int_t iAss(0); iAss < fTracksAss->GetEntriesFast(); iAss++){
        AliAODTrack* trackAss = (AliAODTrack*)fTracksAss->At(iAss);
        if(!trackAss) continue;

        Double_t assPt = trackAss->Pt();
        Double_t assEta = trackAss->Eta();
        Double_t assPhi = trackAss->Phi();
        Double_t assCharge = trackAss->Charge();
        Double_t assEff = 1.0;
        if(fUseEfficiency) assEff = GetEff(assPt, 0, assEta);

        if(trigPt < assPt) continue;
        if(track->GetID() == trackAss->GetID()) continue;

        binscont[0] = trigEta - assEta;
        binscont[1] = RangePhi(trigPhi - assPhi);
        binscont[4] = assPt;

        if(TMath::Abs(binscont[0]) < fMergingCut){
          Double_t dPhiStarLow = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, 0.8);
          Double_t dPhiStarHigh = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, 2.5);

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

        fhSE[spec]->Fill(binscont,0,1./(trigEff*assEff));
      }
    }
  } // end TPC - TPC
  else if(fAnalType == eTPCFMDA || fAnalType == eTPCFMDC){
    Double_t binscont[4];
    binscont[2] = fPVz;

    for(Int_t iTrig(0); iTrig < fTracksTrig[spec]->GetEntriesFast(); iTrig++){
      AliAODTrack* track = (AliAODTrack*)fTracksTrig[spec]->At(iTrig);
      if(!track) continue;

      Double_t trigPt = track->Pt();
      Double_t trigEta = track->Eta();
      Double_t trigPhi = track->Phi();
      Double_t trigEff = 1.0;
      if(fUseEfficiency) trigEff = GetEff(trigPt, spec, trigEta);
      binscont[3] = trigPt;

      for(Int_t iAss(0); iAss < fTracksAss->GetEntriesFast(); iAss++){
        AliPartSimpleForCorr* trackAss = (AliPartSimpleForCorr*)fTracksAss->At(iAss);
        if(!trackAss) continue;

        Double_t assEta = trackAss->Eta();
        Double_t assPhi = trackAss->Phi();
        Double_t assMult = trackAss->Multiplicity();

        binscont[0] = trigEta - assEta;
        binscont[1] = RangePhi(trigPhi - assPhi);

        fhSE[spec]->Fill(binscont,0,assMult/(trigEff));
      }
    }
  } // end TPC - FMD
  else{
    Double_t binscont[3];
    binscont[2] = fPVz;

    for(Int_t iTrig(0); iTrig < fTracksTrig[spec]->GetEntriesFast(); iTrig++){
      AliPartSimpleForCorr* track = (AliPartSimpleForCorr*)fTracksTrig[spec]->At(iTrig);
      if(!track) continue;

      Double_t trigEta = track->Eta();
      Double_t trigPhi = track->Phi();
      Double_t trigMult = track->Multiplicity();

      for(Int_t iAss(0); iAss < fTracksAss->GetEntriesFast(); iAss++){
        AliPartSimpleForCorr* trackAss = (AliPartSimpleForCorr*)fTracksAss->At(iAss);
        if(!trackAss) continue;

        Double_t assEta = trackAss->Eta();
        Double_t assPhi = trackAss->Phi();
        Double_t assMult = trackAss->Multiplicity();

        binscont[0] = trigEta - assEta;
        binscont[1] = RangePhi(trigPhi - assPhi);

        fhSE[spec]->Fill(binscont,0,assMult*trigMult);
      }
    }
  } // end FMD - FMD

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowFMD::FillCorrelationsMixed(const Int_t spec)
{
  if(!fTracksTrig[spec] || !fhTrigTracks[spec] || !fTracksAss) { AliError("Necessary inputs missing, terminating!"); return; }

  AliEventPool *pool = fPoolMgr->GetEventPool(fCentrality, fPVz);
  if(!pool) { AliError(Form("No pool found for centrality = %f, zVtx = %f", fCentrality,fPVz)); return; }

  if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() > fMinEventsToMix) {
    Int_t nMix = pool->GetCurrentNEvents();

    if(fAnalType == eTPCTPC){
      Double_t binscont[5];
      binscont[2] = fPVz;

      for(Int_t iTrig(0); iTrig < fTracksTrig[spec]->GetEntriesFast(); iTrig++){
        AliAODTrack* track = (AliAODTrack*)fTracksTrig[spec]->At(iTrig);
        if(!track) continue;

        Double_t trigPt = track->Pt();
        Double_t trigEta = track->Eta();
        Double_t trigPhi = track->Phi();
        Double_t trigCharge = track->Charge();
        binscont[3] = trigPt;

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

            binscont[0] = trigEta - assEta;
            binscont[1] = RangePhi(trigPhi - assPhi);
            binscont[4] = assPt;

            if(TMath::Abs(binscont[0]) < fMergingCut){
              Double_t dPhiStarLow = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, 0.8);
              Double_t dPhiStarHigh = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, 2.5);

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

            fhME[spec]->Fill(binscont,0,1./(Double_t)nMix);
          }
        }
      }
    } // end TPC - TPC
    else if(fAnalType == eTPCFMDA || fAnalType == eTPCFMDC){
      Double_t binscont[4];
      binscont[2] = fPVz;

      for(Int_t iTrig(0); iTrig < fTracksTrig[spec]->GetEntriesFast(); iTrig++){
        AliAODTrack* track = (AliAODTrack*)fTracksTrig[spec]->At(iTrig);
        if(!track) continue;

        Double_t trigPt = track->Pt();
        Double_t trigEta = track->Eta();
        Double_t trigPhi = track->Phi();
        binscont[3] = trigPt;

        for(Int_t eMix(0); eMix < nMix; eMix++){
          TObjArray *mixEvents = pool->GetEvent(eMix);
          for(Int_t iAss(0); iAss < mixEvents->GetEntriesFast(); iAss++){
            AliPartSimpleForCorr* trackAss = (AliPartSimpleForCorr*)mixEvents->At(iAss);
            if(!trackAss) continue;

            Double_t assEta = trackAss->Eta();
            Double_t assPhi = trackAss->Phi();
            Double_t assMult = trackAss->Multiplicity();

            binscont[0] = trigEta - assEta;
            binscont[1] = RangePhi(trigPhi - assPhi);

            fhME[spec]->Fill(binscont,0,assMult/(Double_t)nMix);
          }
        }
      }
    } // end TPC - FMD
    else{
      Double_t binscont[3];
      binscont[2] = fPVz;

      for(Int_t iTrig(0); iTrig < fTracksTrig[spec]->GetEntriesFast(); iTrig++){
        AliPartSimpleForCorr* track = (AliPartSimpleForCorr*)fTracksTrig[spec]->At(iTrig);
        if(!track) continue;

        Double_t trigEta = track->Eta();
        Double_t trigPhi = track->Phi();
        Double_t trigMult = track->Multiplicity();

        for(Int_t eMix(0); eMix < nMix; eMix++){
          TObjArray *mixEvents = pool->GetEvent(eMix);
          for(Int_t iAss(0); iAss < mixEvents->GetEntriesFast(); iAss++){
            AliPartSimpleForCorr* trackAss = (AliPartSimpleForCorr*)mixEvents->At(iAss);
            if(!trackAss) continue;

            Double_t assEta = trackAss->Eta();
            Double_t assPhi = trackAss->Phi();
            Double_t assMult = trackAss->Multiplicity();

            binscont[0] = trigEta - assEta;
            binscont[1] = RangePhi(trigPhi - assPhi);

            fhME[spec]->Fill(binscont,0,(trigMult*assMult)/(Double_t)nMix);
          }
        }
      }
    } // end FMD - FMD

  } // event pool done

  if((!fDoPID && spec == eCharged) || (fDoPID && spec == eProton)){
    TObjArray* cloneArray = (TObjArray *)fTracksAss->Clone();
    cloneArray->SetOwner(kTRUE);
    pool->UpdatePool(cloneArray);
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowFMD::AreEfficienciesLoaded()
{
  if(!fInputListEfficiency) {AliError("Efficiency input list not loaded"); return kFALSE; }
  TString part[4] = {"ch", "pi", "ka", "pr"};
  for(Int_t p(0); p < 4; p++){
    if(!fEfficiencyEtaDependent){
      fhEfficiency[p] = (TH1D*)fInputListEfficiency->FindObject(Form("EffRescaled_%s_eta0",part[p].Data()));
      if(!fhEfficiency[p]) {AliError(Form("Efficiency (%s, not eta dependent) not loaded",part[p].Data())); return kFALSE; }
    }
    else{
      Int_t etaReg[4] = {100, 101, 102, 103};
      for(Int_t eta(0); eta < 4; eta++){
        fhEfficiencyEta[p][eta] = (TH1D*)fInputListEfficiency->FindObject(Form("EffRescaled_%s_eta%d",part[p].Data(), etaReg[eta]));
        if(!fhEfficiencyEta[p][eta]) {AliError(Form("Efficiency (%s, eta region %d) not loaded",part[p].Data(),etaReg[eta])); return kFALSE; }
      }
    }
    if(!fDoPID) break;
  }
  fhEventCounter->Fill("Efficiencies loaded",1);
  return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlowFMD::GetEff(const Double_t dPt, const Int_t spec, const Double_t dEta)
{
  if(!fUseEfficiency) return 1.0;
  if(spec < 0 || spec > 3) { AliError("Efficiency loading -- species out of range! ");}
  if(!fEfficiencyEtaDependent && !fhEfficiency[spec]) { AliFatal("Efficiency not loaded"); return 0.0; }
  if(!fEfficiencyEtaDependent) return fhEfficiency[spec]->GetBinContent(fhEfficiency[spec]->FindFixBin(dPt));
  else{
    Int_t etaReg = -1;
    if(dEta < -0.8) AliError("Efficiency loading -- eta out of range! ");
    else if(dEta < -0.465) etaReg = 0;
    else if(dEta < 0.0) etaReg = 1;
    else if(dEta < 0.465) etaReg = 2;
    else if(dEta < 0.8) etaReg = 3;
    else AliError("Efficiency loading -- eta out of range! ");
    if(!fhEfficiencyEta[spec][etaReg] || etaReg < 0) { AliFatal("Efficiency not loaded"); return 0.0; }
    return fhEfficiencyEta[spec][etaReg]->GetBinContent(fhEfficiencyEta[spec][etaReg]->FindFixBin(dPt));
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowFMD::CreateTHnCorrelations(){
  Int_t nSteps = 1;
  Double_t binning_dphi[73] = { -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0,       0.087266,  0.174533,  0.261799,  0.349066,  0.436332, 0.523599,  0.610865,  0.698132,  0.785398,  0.872665,  0.959931, 1.047198,  1.134464,  1.221730,  1.308997,  1.396263,  1.483530, 1.570796,  1.658063,  1.745329,  1.832596,  1.919862,  2.007129, 2.094395,  2.181662,  2.268928,  2.356194,  2.443461,  2.530727, 2.617994,  2.705260,  2.792527,  2.879793,  2.967060,  3.054326, 3.141593,  3.228859,  3.316126,  3.403392,  3.490659,  3.577925, 3.665191,  3.752458,  3.839724,  3.926991,  4.014257,  4.101524, 4.188790,  4.276057,  4.363323,  4.450590,  4.537856,  4.625123, 4.712389};
  const Int_t sizePtTrig = fPtBinsTrigCharged.size() - 1;
  const Int_t sizePtAss = fPtBinsAss.size() - 1;

  TString nameS[4] = {"fhChargedSE", "fhPidSE_Pion", "fhPidSE_Kaon", "fhPidSE_Proton" };
  TString nameM[4] = {"fhChargedME", "fhPidME_Pion", "fhPidME_Kaon", "fhPidME_Proton" };

  // Double_t binning_dphi_reduce[] = { -1.570796,  -1.221730, -0.872665, -0.523599, -0.174533,  0.174533,  0.523599,
  //      0.872665,  1.221730,   1.570796, 1.919862,  2.268928,   2.617994, 2.967060,  3.316126,  3.665191,
  //      4.014257,  4.363323,   4.712389};
  // Int_t nbinning_dphi_reduce = sizeof(binning_dphi_reduce)/sizeof(Double_t) - 1;

  if(fAnalType == eTPCFMDA || fAnalType == eTPCFMDC){
    Double_t binning_detaFMDTPC[]={-6.,-5.8, -5.6, -5.4, -5.2, -5.0, -4.8, -4.6, -4.4, -4.2, -4., -3.8, -3.6, -3.4, -3.2, -3., -2.8, -2.6, -2.4, -2.2, -2., -1.8, -1.6, -1.4, -1.2, -1., -0.8};
    Double_t binning_detaFMDCTPC[]={ 1., 1.2, 1.4, 1.6, 1.8, 2. , 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4.};
    Int_t ndetatpcfmd;
    if(fAnalType == eTPCFMDA) ndetatpcfmd= sizeof(binning_detaFMDTPC)/sizeof(Double_t) - 1;
    else ndetatpcfmd = sizeof(binning_detaFMDCTPC)/sizeof(Double_t) - 1;

    Int_t iTrackBin_tpcfmdA[] = {26, 72, 10, sizePtTrig};
    Int_t iTrackBin_tpcfmdC[] = {15, 72, 10, sizePtTrig};
    Int_t nTrackBin_tpcfmd = sizeof(iTrackBin_tpcfmdA) / sizeof(Int_t);

    for(Int_t i(0); i < 4; i++){
      if(fAnalType == eTPCFMDA){
        fhSE[i] = new AliTHn(nameS[i], nameS[i], nSteps, nTrackBin_tpcfmd, iTrackBin_tpcfmdA);
        fhME[i] = new AliTHn(nameM[i], nameM[i], nSteps, nTrackBin_tpcfmd, iTrackBin_tpcfmdA);
        fhSE[i]->SetBinLimits(0, -6., -0.8);
        fhME[i]->SetBinLimits(0, -6., -0.8);
      } // TPC - FMDA
      else{
        fhSE[i] = new AliTHn(nameS[i], nameS[i], nSteps, nTrackBin_tpcfmd, iTrackBin_tpcfmdC);
        fhME[i] = new AliTHn(nameM[i], nameM[i], nSteps, nTrackBin_tpcfmd, iTrackBin_tpcfmdC);
        fhSE[i]->SetBinLimits(0, 1., 4.);
        fhME[i]->SetBinLimits(0, 1., 4.);
      } // TPC - FMDC
      fhSE[i]->SetBinLimits(3, fPtBinsTrigCharged.data());
      fhSE[i]->SetVarTitle(3, "p_{T} [GeV/c] (trig)");
      fhME[i]->SetBinLimits(3, fPtBinsTrigCharged.data());
      fhME[i]->SetVarTitle(3, "p_{T} [GeV/c] (trig)");

      if(!fDoPID) break;
    }


  } // end TPC - FMD
  else if(fAnalType == eFMDAFMDC){
    Int_t iTrackBin_fmdAfmdC[] = {48, 72, 10};
    // Int_t iTrackBin_fmdAfmdC[] = {48, 36, 10};
    Int_t nTrackBin_fmdAfmdC = sizeof(iTrackBin_fmdAfmdC) / sizeof(Int_t);

    for(Int_t i(0); i < 4; i++){
      fhSE[i] = new AliTHn(nameS[i], nameS[i], nSteps, nTrackBin_fmdAfmdC, iTrackBin_fmdAfmdC);
      fhSE[i]->SetBinLimits(0, 3.4,8.2);

      fhME[i] = new AliTHn(nameM[i], nameM[i], nSteps, nTrackBin_fmdAfmdC, iTrackBin_fmdAfmdC);
      fhME[i]->SetBinLimits(0, 3.4,8.2);

      if(!fDoPID) break;
    }
  } // end FMD - FMD
  else {
    Double_t binning_deta_tpctpc[33] = {-1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,    0.1,  0.2,  0.3,  0.4,  0.5, 0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5, 1.6};
    const Int_t iBinningTPCTPC[5] = {32,72,sizePtTrig,sizePtAss, 10};

    for(Int_t i(0); i < 4; i++){
      fhSE[i] = new AliTHn(nameS[i], nameS[i], nSteps, 5, iBinningTPCTPC);
      fhSE[i]->SetBinLimits(0, binning_deta_tpctpc);
      fhSE[i]->SetBinLimits(3, fPtBinsTrigCharged.data());
      fhSE[i]->SetBinLimits(4, fPtBinsAss.data());
      fhSE[i]->SetVarTitle(3, "p_{T} [GeV/c] (trig)");
      fhSE[i]->SetVarTitle(4, "p_{T} [GeV/c] (ass)");

      fhME[i] = new AliTHn(nameM[i], nameM[i], nSteps, 5, iBinningTPCTPC);
      fhME[i]->SetBinLimits(0, binning_deta_tpctpc);
      fhME[i]->SetBinLimits(3, fPtBinsTrigCharged.data());
      fhME[i]->SetBinLimits(4, fPtBinsAss.data());
      fhME[i]->SetVarTitle(3, "p_{T} [GeV/c] (trig)");
      fhME[i]->SetVarTitle(4, "p_{T} [GeV/c] (ass)");

      if(!fDoPID) break;
    }

  } // end TPC - TPC

  // all
  for(Int_t i(0); i < 4; i++){
    fhSE[i]->SetBinLimits(1, binning_dphi);
    fhSE[i]->SetBinLimits(2, -10,10);
    fhSE[i]->SetVarTitle(0, "#Delta#eta");
    fhSE[i]->SetVarTitle(1, "#Delta#phi");
    fhSE[i]->SetVarTitle(2, "PVz [cm]");
    fOutputListCharged->Add(fhSE[i]);

    fhME[i]->SetBinLimits(1, binning_dphi);
    fhME[i]->SetBinLimits(2, -10,10);
    fhME[i]->SetVarTitle(0, "#Delta#eta");
    fhME[i]->SetVarTitle(1, "#Delta#phi");
    fhME[i]->SetVarTitle(2, "PVz [cm]");
    fOutputListCharged->Add(fhME[i]);

    if(!fDoPID) break;
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowFMD::PrepareFMDTracks(){
  if(!fTracksAss) { AliError("Problem with fTracksAss, terminating!"); return kFALSE; }
  if(fAnalType == eFMDAFMDC && !fTracksTrig[0]) { AliError("Problem with fTracksTrig (no PID), terminating!"); return kFALSE; }

  AliAODForwardMult* aodForward=static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
  if(!aodForward) { AliError("Problem with aodForward, terminating!"); return kFALSE; }

  const TH2D& d2Ndetadphi = aodForward->GetHistogram();
  Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
  Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();

  Float_t nFMD_fwd_hits=0.;
  Float_t nFMD_bwd_hits=0.;

  for (Int_t iEta = 1; iEta <= nEta; iEta++)
  {
    Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
    if (!valid) continue;
    Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);

    for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++)
    {
      // Bin content is most probable number of particles!
      Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
      Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);

      if(mostProbableN > 0) {
    	   if(eta > 0){
    	     nFMD_fwd_hits+=mostProbableN;
           if(eta > 1.8 && eta < 4.8){
             if(fAnalType == eTPCFMDA) fTracksAss->Add(new AliPartSimpleForCorr(eta,phi,mostProbableN));
             if(fAnalType == eFMDAFMDC) {
               fTracksTrig[0]->Add(new AliPartSimpleForCorr(eta,phi,mostProbableN));
               fhTrigTracks[0]->Fill(eta,fPVz,mostProbableN);
             }
           }
    	   } // eta positive
         else
         {
    	     nFMD_bwd_hits+=mostProbableN;
           if(eta < -1.8 && eta > -3.2){
             if(fAnalType == eTPCFMDC || fAnalType == eFMDAFMDC) fTracksAss->Add(new AliPartSimpleForCorr(eta,phi,mostProbableN));
           }
    	   } // eta negative
    	 } // most probable > 0
    } // end phi
  } // end eta

  if(fUseFMDcut){
    if(nFMD_fwd_hits==0 || nFMD_bwd_hits==0) {
      fTracksAss->Clear();
      if(fAnalType == eFMDAFMDC) fTracksTrig[0]->Clear();
      return kFALSE;
    }
    AliAODVZERO *fvzero = fAOD->GetVZEROData();
    if(!fvzero) { AliError("Problem with VZEROData, terminating!"); return kFALSE; }
    Float_t nV0A_hits = fvzero->GetMTotV0A();
    Float_t nV0C_hits = fvzero->GetMTotV0C();

    if((nV0A_hits<(fFMDcutapar0*nFMD_fwd_hits-fFMDcutapar1)) || (nV0C_hits<(fFMDcutcpar0*nFMD_bwd_hits-fFMDcutcpar1))){
      fTracksAss->Clear();
      if(fAnalType == eFMDAFMDC) fTracksTrig[0]->Clear();
      return kFALSE;
    }
    fhEventCounter->Fill("FMD cuts OK",1);
  }

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowFMD::PrintSetup(){
  printf("\n\n\n ************** Parameters ************** \n");
  printf("\t fAnalType: (Int_t) %d\n", fAnalType);
  printf("\t fDoPID: (Bool_t) %s\n", fDoPID ? "kTRUE" : "kFALSE");
  printf("\t fUseNch: (Bool_t) %s\n", fUseNch ? "kTRUE" : "kFALSE");
  printf("\t fIsHMpp: (Bool_t) %s\n", fIsHMpp ? "kTRUE" : "kFALSE");
  printf("\t fUseEfficiency: (Bool_t) %s\n",  fUseEfficiency ? "kTRUE" : "kFALSE");
  printf("\t fEfficiencyEtaDependent: (Bool_t) %s\n", fEfficiencyEtaDependent ? "kTRUE" : "kFALSE");
  printf(" **************************** \n");
  printf("\t fAbsEtaMax: (Double_t) %f\n", fAbsEtaMax);
  printf("\t fPtMinTrig -- fPtMaxTrig: (Double_t) %f -- %f\n", fPtMinTrig, fPtMaxTrig);
  printf("\t fPtMinAss -- fPtMaxAss: (Double_t) %f -- %f\n", fPtMinAss, fPtMaxAss);
  printf("\t fCentMin -- fCentMax: (Double_t) %f -- %f\n", fCentMin, fCentMax);
  printf(" **************************** \n");
  printf("\t fUseFMDcut: (Bool_t) %s\n", fUseFMDcut ? "kTRUE" : "kFALSE");
  printf("\t fFMDcutapar0 -- fFMDcutapar1: (Double_t) %f -- %f\n", fFMDcutapar0, fFMDcutapar1);
  printf("\t fFMDcutcpar0 -- fFMDcutcpar1: (Double_t) %f -- %f\n", fFMDcutcpar0, fFMDcutcpar1);
}
