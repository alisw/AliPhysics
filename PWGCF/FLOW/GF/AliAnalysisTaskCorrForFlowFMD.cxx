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
    fAnalType(eFMDAFMDC),
    fColSystem(sPPb),
    fTrigger(AliVEvent::kINT7),
    fIsMC(kFALSE),
    fIsTPCgen(kFALSE),
    fIsFMDgen(kFALSE),
    fIsHMpp(kFALSE),
    fDoPID(kFALSE),
    fDoV0(kFALSE),
    fUseNch(kFALSE),
    fUseEfficiency(kFALSE),
    fUseFMDcut(kTRUE),
    fUseOppositeSidesOnly(kFALSE),
    fUseCentralityCalibration(kFALSE),
    fSkipCorr(kFALSE),
    fIsAntiparticleCheck(kFALSE),
    fDoAntiparticleOnly(kFALSE),
    fRejectHighPtEvents(kFALSE),
    fBoostAMPT(kFALSE),
    fFilterBit(96),
    fbSign(0),
    fRunNumber(-1),
    fNofTracks(0),
    fNofMinHighPtTracksForRejection(0),
    fNchMin(0),
    fNchMax(100000),
    fNbinsMinv(60),
    fnTPCcrossedRows(70),
    fNOfSamples(1.0),
    fSampleIndex(0.0),
    fPtMinTrig(0.5),
    fPtMaxTrig(10.0),
    fPtMinAss(0.5),
    fPtMaxAss(1.5),
    fFMDcutapar0(1.64755),
    fFMDcutapar1(119.602),
    fFMDcutcpar0(2.73426),
    fFMDcutcpar1(150.31),
    fFMDAacceptanceCutLower(1.8),
    fFMDAacceptanceCutUpper(4.8),
    fFMDCacceptanceCutLower(1.8),
    fFMDCacceptanceCutUpper(3.2),
    fCentMin(0.0),
    fCentMax(10.0),
    fCentrality(-10.0),
    fAbsEtaMax(0.8),
    fPVz(100.0),
    fPVzCut(10.0),
    fTPCclMin(70.),
    fCutDCAz(0.),
    fCutDCAxySigma(0.),
    fCutTPCchi2pCl(0.),
    fPIDbayesPion(0.95),
    fPIDbayesKaon(0.85),
    fPIDbayesProton(0.85),
    fV0ratioClusters(0.8),
    fV0dcaK0ToPV(0.06),
    fV0dcaNegLambdaToPV(0.25),
    fV0dcaDaugtersK0(1.),
    fV0dcaDaugtersLambda(1.),
    fK0radiusMin(0.5),
    fK0radiusMax(200.),
    fLambdaradiusMin(0.5),
    fLambdaradiusMax(200.),
    fCutCPAK0s(0.97),
    fCutCPALambda(0.995),
    fCutTauK0s(0.),
    fCutTauLambda(0.),
    fSigmaTPC(3.),
    fTrackLength(90),
    fMassRejWindowK0(0.005),
    fMassRejWindowLambda(0.01),
    fCentEstimator("V0M"),
    fSystematicsFlag(""),
    fPoolMaxNEvents(2000),
    fPoolMinNTracks(50000),
    fMinEventsToMix(5),
    fNzVtxBins(10),
    fNCentBins(15),
    fMergingCut(0.0)
{}
//_____________________________________________________________________________
AliAnalysisTaskCorrForFlowFMD::AliAnalysisTaskCorrForFlowFMD(const char* name, Bool_t bUseEff, Bool_t bUseCalib) : AliAnalysisTaskSE(name),
    fAOD(0),
    fOutputListCharged(0),
    fInputListEfficiency(0),
    fTracksAss(0),
    fPIDResponse(0),
    fPIDCombined(0),
    fPoolMgr(0),
    fhEventCounter(0),
    fhEventMultiplicity(0),
    fAnalType(eFMDAFMDC),
    fColSystem(sPPb),
    fTrigger(AliVEvent::kINT7),
    fIsMC(kFALSE),
    fIsTPCgen(kFALSE),
    fIsFMDgen(kFALSE),
    fIsHMpp(kFALSE),
    fDoPID(kFALSE),
    fDoV0(kFALSE),
    fUseNch(kFALSE),
    fUseEfficiency(bUseEff),
    fUseFMDcut(kTRUE),
    fUseOppositeSidesOnly(kFALSE),
    fUseCentralityCalibration(bUseCalib),
    fSkipCorr(kFALSE),
    fIsAntiparticleCheck(kFALSE),
    fDoAntiparticleOnly(kFALSE),
    fRejectHighPtEvents(kFALSE),
    fBoostAMPT(kFALSE),
    fFilterBit(96),
    fbSign(0),
    fRunNumber(-1),
    fNofTracks(0),
    fNofMinHighPtTracksForRejection(0),
    fNchMin(0),
    fNchMax(100000),
    fNbinsMinv(60),
    fnTPCcrossedRows(70),
    fNOfSamples(1.0),
    fSampleIndex(0.0),
    fPtMinTrig(0.5),
    fPtMaxTrig(10.0),
    fPtMinAss(0.5),
    fPtMaxAss(1.5),
    fFMDcutapar0(1.64755),
    fFMDcutapar1(119.602),
    fFMDcutcpar0(2.73426),
    fFMDcutcpar1(150.31),
    fFMDAacceptanceCutLower(1.8),
    fFMDAacceptanceCutUpper(4.8),
    fFMDCacceptanceCutLower(1.8),
    fFMDCacceptanceCutUpper(3.2),
    fCentMin(0.0),
    fCentMax(10.0),
    fCentrality(-10.0),
    fAbsEtaMax(0.8),
    fPVz(100.0),
    fPVzCut(10.0),
    fTPCclMin(70.),
    fCutDCAz(0.),
    fCutDCAxySigma(0.),
    fCutTPCchi2pCl(0.),
    fPIDbayesPion(0.95),
    fPIDbayesKaon(0.85),
    fPIDbayesProton(0.85),
    fV0ratioClusters(0.8),
    fV0dcaK0ToPV(0.06),
    fV0dcaNegLambdaToPV(0.25),
    fV0dcaDaugtersK0(1.),
    fV0dcaDaugtersLambda(1.),
    fK0radiusMin(0.5),
    fK0radiusMax(200.),
    fLambdaradiusMin(0.5),
    fLambdaradiusMax(200.),
    fCutCPAK0s(0.97),
    fCutCPALambda(0.995),
    fCutTauK0s(0.),
    fCutTauLambda(0.),
    fSigmaTPC(3.),
    fTrackLength(90),
    fMassRejWindowK0(0.005),
    fMassRejWindowLambda(0.01),
    fCentEstimator("V0M"),
    fSystematicsFlag(""),
    fPoolMaxNEvents(2000),
    fPoolMinNTracks(50000),
    fMinEventsToMix(5),
    fNzVtxBins(10),
    fNCentBins(15),
    fMergingCut(0.0)
{
    DefineInput(0, TChain::Class());
    if(bUseEff) { DefineInput(1, TList::Class()); }
    if(bUseCalib) {
      if(bUseEff) DefineInput(2, TH1D::Class());
      else  DefineInput(1, TH1D::Class());
    }
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
    if(fAnalType == eFMDAFMDC && fDoV0) { AliWarning("V0 on when running FMDA-FMDC. Turning off."); fDoV0 = kFALSE; }
    if(!fDoV0) { fNbinsMinv = 1; }

    fzVtxBins = {-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0};

    fOutputListCharged = new TList();
    fOutputListCharged->SetOwner(kTRUE);

    fhEventCounter = new TH1D("fhEventCounter","Event Counter",10,0,10);
    fOutputListCharged->Add(fhEventCounter);

    fhEventMultiplicity = new TH1D("fhEventMultiplicity","Event multiplicity; N_{ch}",200,0,200);
    fOutputListCharged->Add(fhEventMultiplicity);

    fHistFMDeta = new TH2D("fHistFMDeta", "FMD eta vs. PVz; eta; PVz [cm]", 90, -4, 5, 20, -10, 10);
    fOutputListCharged->Add(fHistFMDeta);

    TString pidName[6] = {"", "_Pion", "_Kaon", "_Proton", "_K0s", "_Lambda"};
    const Int_t sizePtTrig = fPtBinsTrigCharged.size() - 1;
    const Int_t sizeOfSamples = (Int_t) fNOfSamples;
    const Int_t sizeMbins = fNbinsMinv;
    Int_t binsFMD[] = {10, 10};
    Int_t binsPID[] = {10, 10, sizePtTrig};
    Int_t binsV0[] = {10, 10, sizePtTrig, sizeMbins};
    Double_t min[2] = {0.44, 1.08};
    Double_t max[2] = {0.56, 1.15};

    for(Int_t i(0); i < 6; i++){
      if(!fDoPID && i > 0 && i < 4) continue;
      if(!fDoV0 && i > 3) continue;
      if(fIsAntiparticleCheck && (i == 0 || i == 4)) continue;
      if(fAnalType == eFMDAFMDC) fhTrigTracks[i] = new AliTHn(Form("fhTrigTracks%s",pidName[i].Data()), Form("fhTrigTracks (%s)",pidName[i].Data()), 1, 2, binsFMD);
      else if(i < 4) fhTrigTracks[i] = new AliTHn(Form("fhTrigTracks%s",pidName[i].Data()), Form("fhTrigTracks (%s)",pidName[i].Data()), 1, 3, binsPID);
      else if(i > 3) fhTrigTracks[i] = new AliTHn(Form("fhTrigTracks%s",pidName[i].Data()), Form("fhTrigTracks (%s)",pidName[i].Data()), 1, 4, binsV0);
      else AliError("This should not happen! There might be a problem with trigher AliTHn!");
      fhTrigTracks[i]->SetBinLimits(0,-10,10);
      fhTrigTracks[i]->SetBinLimits(1,0,10);
      fhTrigTracks[i]->SetVarTitle(0, "PVz [cm]");
      fhTrigTracks[i]->SetVarTitle(1, "Sample");
      if(fAnalType != eFMDAFMDC){
        fhTrigTracks[i]->SetBinLimits(2,fPtBinsTrigCharged.data());
        fhTrigTracks[i]->SetVarTitle(2, "p_{T} (trig)");
        if(i > 3){
          fhTrigTracks[i]->SetBinLimits(3,min[i-4],max[i-4]);
          fhTrigTracks[i]->SetVarTitle(3, "M_{inv}");
        }
      }
      fOutputListCharged->Add(fhTrigTracks[i]);
    }


    if(fDoV0){
      for(Int_t i(0); i < 2; i++){
        fhV0Counter[i] = new TH1D(Form("fhCounter_%s",pidName[i-4].Data()),"V0 Counter",10,0,10);
        fOutputListCharged->Add(fhV0Counter[i]);
      }
    }

    if(fDoPID || fDoV0){
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

    if(!fSkipCorr) CreateTHnCorrelations();

    for(Int_t i(0); i < 6; i++){
      if(!fDoPID && i > 0 && i < 4) continue;
      if(!fDoV0 && i > 3) continue;
      fhPT[i] = new TH1D(Form("PT%s",pidName[i].Data()), "PT", 1000, 0, 10);
      fhPT[i]->Sumw2();
      fOutputListCharged->Add(fhPT[i]);
      if(fDoV0 && i > 3){
        fhPTvsMinv[i-4] = new TH2D(Form("PT_minv%s",pidName[i].Data()), "PT vs. minv", 1000, 0, 10, sizeMbins, min[i-4],max[i-4]);
        fhPTvsMinv[i-4]->Sumw2();
        fOutputListCharged->Add(fhPTvsMinv[i-4]);
      }
      if(!fSkipCorr) break;
    }

    if(fUseEfficiency) {
      fInputListEfficiency = (TList*) GetInputData(1);
      if(fAbsEtaMax > 0.8) AliWarning("Efficiency loading -- eta can be out of range!");
      if(fSystematicsFlag.IsNull()) fSystematicsFlag = "Ev0_Tr0";
      if(fColSystem == sPPb && !AreEfficienciesLoaded()) { AliError("Efficiencies not loaded!"); return; }
    }

    if(fUseCentralityCalibration){
      if(fUseEfficiency) fhCentCalib = (TH1D*) GetInputData(2);
      else fhCentCalib = (TH1D*) GetInputData(1);
      if(!fhCentCalib) { AliError("Centrality calibration histogram not loaded!"); return; }
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

    fSampleIndex = gRandom->Uniform(0,fNOfSamples);
    fTracksAss = new TObjArray;
    fTracksTrig[0] = new TObjArray;

    if(fDoPID){
      for(Int_t i(1); i < 4; i++){
        fTracksTrig[i] = new TObjArray;
      }
    }
    if(fDoV0){
      for(Int_t i(4); i < 6; i++){
        fTracksTrig[i] = new TObjArray;
      }
    }

    if(fUseEfficiency && fColSystem == sPP && (fRunNumber != fAOD->GetRunNumber()) && !AreEfficienciesLoaded()) { return; }

    // FMD - V0 correlation event cut
    if(fAnalType != eTPCTPC) {
      if(!PrepareFMDTracks()){
        delete fTracksAss;
        delete fTracksTrig[0];
        PostData(1, fOutputListCharged);
        return;
      }
    }

    if(!fIsTPCgen || fUseNch)  {
      if(!PrepareTPCTracks()){
        for(Int_t i(0); i < 6; i++){
          if(!fDoPID && i > 0 && i < 4) continue;
          if(!fDoV0 && i > 3) continue;
          if(fTracksTrig[i]) delete fTracksTrig[i];
        }
        PostData(1, fOutputListCharged);
        return;
      }
    }

    if(fIsMC){
      if(!PrepareMCTracks()){
        for(Int_t i(0); i < 6; i++){
          if(!fDoPID && i > 0 && i < 4) continue;
          if(!fDoV0 && i > 3) continue;
          if(fTracksTrig[i]) delete fTracksTrig[i];
        }
        PostData(1, fOutputListCharged);
        return;
      }
    } // end MC

    if(!fTracksAss->IsEmpty() && !fSkipCorr){
      for(Int_t i(0); i < 6; i++){
        if(!fDoPID && i > 0 && i < 4) continue;
        if(!fDoV0 && i > 3) continue;
        if(fIsAntiparticleCheck && (i == 0 || i == 4)) continue;

        FillCorrelations(i);
        FillCorrelationsMixed(i);

        fTracksTrig[i]->Clear();
        delete fTracksTrig[i];
      }
    }

    if(fUseEfficiency) fRunNumber = fAOD->GetRunNumber();

    fTracksAss->Clear();
	  delete fTracksAss;

    PostData(1, fOutputListCharged);
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowFMD::Terminate(Option_t *)
{
   if(fPoolMgr) delete fPoolMgr;
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
  if(!fUseCentralityCalibration){
    Float_t dPercentile = multSelection->GetMultiplicityPercentile(fCentEstimator);
    if(dPercentile > 100 || dPercentile < 0) { return kFALSE; }
    fhEventCounter->Fill("PercOK",1);
    fCentrality = (Double_t) dPercentile;
  }
  else if(fIsMC){
    AliMCEvent* mcEvent = dynamic_cast<AliMCEvent*>(MCEvent());
    if(!mcEvent) return kFALSE;
    Int_t ntrackv0aprimary=0;

    for(Int_t i(0); i < mcEvent->GetNumberOfTracks(); i++) {
      AliMCParticle* part = (AliMCParticle*)mcEvent->GetTrack(i);
      if(!part->IsPhysicalPrimary()) continue;
      Double_t mceta = part->Eta();
      if(fBoostAMPT) mceta = TransverseBoost(part);

      if(part->Charge()==0)        continue;
      if(mceta>2.8 && mceta<5.1) ntrackv0aprimary++;
    }
    Int_t nbinmult= fhCentCalib->GetXaxis()->FindBin(ntrackv0aprimary);
    fCentrality = (Double_t) fhCentCalib->GetBinContent(nbinmult);
  }
  else{
    AliAODVZERO* fvzero = fAOD->GetVZEROData();
    Double_t sum = 0.;
    Double_t max = 0.;
     for(Int_t i = 32; i < 64; ++i)
     {
       sum += fvzero->GetMultiplicity(i);
       if (fvzero->GetMultiplicity(i) > max) max = fvzero->GetMultiplicity(i);
     }
     sum -= max;

    Int_t nbinmult= fhCentCalib->GetXaxis()->FindBin(sum);
    fCentrality = (Double_t) fhCentCalib->GetBinContent(nbinmult);
  }
  if(fCentrality < fCentMin || fCentrality > fCentMax) { return kFALSE; }
  fhEventCounter->Fill("CentOK",1);

  fPVz = fAOD->GetPrimaryVertex()->GetZ();
  if(TMath::Abs(fPVz) >= fPVzCut) { return kFALSE; }
  fhEventCounter->Fill("PVzOK",1);

  fbSign = (InputEvent()->GetMagneticField() > 0) ? 1 : -1;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowFMD::IsTrackSelected(const AliAODTrack* track) const
{
  if(!track->TestFilterBit(fFilterBit)) { return kFALSE; }
  if(track->GetTPCNcls() < fTPCclMin && fFilterBit != 2) { return kFALSE; }
  if(fAbsEtaMax > 0.0 && TMath::Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
  if(track->Charge() == 0) { return kFALSE; }

  if(fCutDCAz > 0.){
    Double_t vtxXYZ[3], trXYZ[3];
    track->GetXYZ(trXYZ);
    fAOD->GetPrimaryVertex()->GetXYZ(vtxXYZ);
    trXYZ[2] -= vtxXYZ[2];
    if(TMath::Abs(trXYZ[2]) > fCutDCAz) { return kFALSE; }
  }

  if(fCutDCAxySigma > 0.){
    Double_t vtxXYZ[3], trXYZ[3];
    track->GetXYZ(trXYZ);
    fAOD->GetPrimaryVertex()->GetXYZ(vtxXYZ);
    trXYZ[0] -= vtxXYZ[0];
    trXYZ[1] -= vtxXYZ[1];
    Double_t trDcaxy = TMath::Sqrt(trXYZ[0]*trXYZ[0] + trXYZ[1]*trXYZ[1]);
    Double_t cutDcaxy = 0.0015+0.0050/TMath::Power(track->Pt(),1.1);
    if(trDcaxy > fCutDCAxySigma*cutDcaxy) { return kFALSE; }
  }

  if(fCutTPCchi2pCl > 0. && track->GetTPCchi2perCluster() > fCutTPCchi2pCl)  { return kFALSE; }

  return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlowFMD::RangePhiFMD(Double_t DPhi) {
  DPhi = TMath::ATan2(TMath::Sin(DPhi), TMath::Cos(DPhi));
  if (DPhi < (-0.5*TMath::Pi()-0.0001))    DPhi += 2 * TMath::Pi();
  return DPhi;
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
  Double_t l_MaxProb[] = {fPIDbayesPion,fPIDbayesKaon,fPIDbayesProton};
  Bool_t l_TOFUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, l_Probs) & AliPIDResponse::kDetTOF;
  Int_t pidInd = 0;
  for(Int_t i(0); i < AliPID::kSPECIES; i++) pidInd=(l_Probs[i]>l_Probs[pidInd])?i:pidInd;
  Int_t retInd = pidInd-AliPID::kPion+1; //realigning
  if(retInd<1 || retInd>3) return -1;
  if(l_Probs[pidInd] < l_MaxProb[retInd]) return -1;
  //check nsigma cuts
  if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)pidInd))>3) return -1;
  if(bIsTOFok && l_TOFUsed) if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)pidInd))>3) return -1;

  return retInd;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowFMD::IsV0(const AliAODv0* v0) const
{
  if(!v0) {return kFALSE; }
  if(v0->Pt() < fPtMinTrig || v0->Pt() > fPtMaxTrig) { return kFALSE; }

  AliAODTrack* daughterPos = (AliAODTrack*) v0->GetDaughter(0);
  AliAODTrack* daughterNeg = (AliAODTrack*) v0->GetDaughter(1);
  if(!daughterPos || !daughterNeg) { return kFALSE; }

  if(fAbsEtaMax > 0.0 && ((TMath::Abs(daughterPos->Eta()) > fAbsEtaMax) || (TMath::Abs(daughterNeg->Eta()) > fAbsEtaMax)) ) { return kFALSE; }

  if(v0->Charge() != 0) { return kFALSE; }
  if((daughterPos->Charge() == daughterNeg->Charge()) || (TMath::Abs(daughterPos->Charge()) != 1) || (TMath::Abs(daughterNeg->Charge()) != 1)) { return kFALSE; }

  if(v0->GetOnFlyStatus() != kFALSE) { return kFALSE; }
  if(!daughterPos->IsOn(AliAODTrack::kTPCrefit) || !daughterNeg->IsOn(AliAODTrack::kTPCrefit) ) { return kFALSE; }

  // kink rejection
  const AliAODVertex* prodVtxDaughterPos = (AliAODVertex*) daughterPos->GetProdVertex();
  const AliAODVertex* prodVtxDaughterNeg = (AliAODVertex*) daughterNeg->GetProdVertex();
  if( (prodVtxDaughterPos->GetType() == AliAODVertex::kKink ) || (prodVtxDaughterNeg->GetType() == AliAODVertex::kKink ) ) { return kFALSE; }

  // track quality
  if(daughterPos->GetTPCNcls() < fTPCclMin || daughterNeg->GetTPCNcls() < fTPCclMin) { return kFALSE; }
  if(daughterPos->GetTPCNCrossedRows() < fnTPCcrossedRows || daughterNeg->GetTPCNCrossedRows() < fnTPCcrossedRows) { return kFALSE; }

  if(fColSystem==sPP||fColSystem==sPPb){
    if(daughterPos->GetTPCNclsF() < 1 || daughterNeg->GetTPCNclsF() < 1) { return kFALSE; }

    Double_t dRatioCrossFindPos = (Double_t) daughterPos->GetTPCNCrossedRows() / (Double_t) daughterPos->GetTPCNclsF();
    Double_t dRatioCrossFindNeg = (Double_t) daughterNeg->GetTPCNCrossedRows() / (Double_t) daughterNeg->GetTPCNclsF();
    if( dRatioCrossFindPos < fV0ratioClusters || dRatioCrossFindNeg < fV0ratioClusters) { return kFALSE; }
  }else{
    if(daughterPos->GetIntegratedLength()<fTrackLength || daughterNeg->GetIntegratedLength()<fTrackLength ) { return kFALSE; }
    Double_t dRatioCrossLengthPos = (Double_t) daughterPos->GetTPCNCrossedRows() / daughterPos->GetIntegratedLength();
    Double_t dRatioCrossLengthNeg = (Double_t) daughterNeg->GetTPCNCrossedRows() / daughterNeg->GetIntegratedLength();
    if( dRatioCrossLengthPos < fV0ratioLength || dRatioCrossLengthNeg < fV0ratioLength) { return kFALSE; }
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowFMD::IsK0s(const AliAODv0* v0) const
{
  if(fIsAntiparticleCheck) return kFALSE;

  fhV0Counter[0]->Fill("Input",1);
  Double_t dMass = v0->MassK0Short();
  if(dMass < 0.44 || dMass > 0.56) { return kFALSE; }
  fhV0Counter[0]->Fill("Mass OK",1);

  if(v0->RapK0Short() > 0.5) { return kFALSE; }
  fhV0Counter[0]->Fill("RapK0 OK",1);

  //DCA
  if( TMath::Abs(v0->DcaPosToPrimVertex()) < fV0dcaK0ToPV || TMath::Abs(v0->DcaNegToPrimVertex()) < fV0dcaK0ToPV ) { return kFALSE; }

  if( TMath::Abs(v0->DcaV0Daughters()) > fV0dcaDaugtersK0 ) { return kFALSE; }
  //radius
  Double_t dDecayRadius = v0->RadiusV0();
  if( dDecayRadius < fK0radiusMin || dDecayRadius > fK0radiusMax ) { return kFALSE; }

  // cosine of pointing angle (CPA)
  Double_t dCPA = v0->CosPointingAngle(fAOD->GetPrimaryVertex());
  if(dCPA < fCutCPAK0s) { return kFALSE; }
  fhV0Counter[0]->Fill("CPA OK",1);

  // Armenteros-Podolanski plot
  if(fColSystem==sPbPb){
    Double_t dPtArm = v0->PtArmV0();
    Double_t dAlpha = v0->AlphaV0();
    if(dPtArm < (0.2 * TMath::Abs(dAlpha))) { return kFALSE; }
  }
  fhV0Counter[0]->Fill("AP OK",1);

  if(fCutTauK0s > 0.0 && ProperLifetime(v0, 0.497614) > fCutTauK0s) { return kFALSE; }
  fhV0Counter[0]->Fill("LT OK",1);

  // daughter PID
  const AliAODTrack* daughterPos = (AliAODTrack*) v0->GetDaughter(0);
  const AliAODTrack* daughterNeg = (AliAODTrack*) v0->GetDaughter(1);

  if(!HasTrackPIDTPC(daughterPos) || !HasTrackPIDTPC(daughterNeg)) { return kFALSE; }
  fhV0Counter[0]->Fill("TPC PID OK",1);

  if (daughterPos->GetTPCsignalN() < fTPCclMin || daughterNeg->GetTPCsignalN() < fTPCclMin) { return kFALSE; }
  fhV0Counter[0]->Fill("TPC sig OK",1);
  Float_t nSigmaPiPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(daughterPos, AliPID::kPion));
  Float_t nSigmaPiNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(daughterNeg, AliPID::kPion));
  if(nSigmaPiPos > fSigmaTPC || nSigmaPiNeg > fSigmaTPC) { return kFALSE; }
  fhV0Counter[0]->Fill("3Sigma OK",1);

  // cross mass rejection
  if(fColSystem==sPP||fColSystem==sPPb){
    Double_t dMassLambda = v0->MassLambda();
    Double_t dMassALambda = v0->MassAntiLambda();
    if(TMath::Abs(dMassLambda - 1.11568) < fMassRejWindowK0) { return kFALSE; }
    if(TMath::Abs(dMassALambda - 1.11568) < fMassRejWindowK0) { return kFALSE; }
  }
  fhV0Counter[0]->Fill("Mass crosscheck OK",1);

  Double_t binscont[4] = {fPVz, fSampleIndex, v0->Pt(), dMass};
  fhTrigTracks[4]->Fill(binscont,0,1.);
  fTracksTrig[4]->Add(new AliPartSimpleForCorr(v0->Eta(),v0->Phi(),v0->Pt(),dMass));
  if(fSkipCorr){
    fhPT[4]->Fill(v0->Pt());
    fhPTvsMinv[0]->Fill(v0->Pt(),dMass);
  }
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowFMD::IsLambda(const AliAODv0* v0) const
{
  fhV0Counter[1]->Fill("Input",1);
  Bool_t isL = kFALSE;
  Bool_t isAL = kFALSE;
  // inv. mass window
  Double_t dMassLambda = v0->MassLambda();
  Double_t dMassALambda = v0->MassAntiLambda();
  if( dMassLambda > 1.08 && dMassLambda < 1.15) { isL = kTRUE; }
  if( dMassALambda > 1.08 && dMassALambda < 1.15) { isAL = kTRUE; }
  if(!isL && !isAL)  { return kFALSE; }
  fhV0Counter[1]->Fill("Mass OK",1);

  if(v0->RapLambda() > 0.5) { return kFALSE; }
  fhV0Counter[1]->Fill("RapK0 OK",1);

  if(fIsAntiparticleCheck && fDoAntiparticleOnly && isL) return kFALSE;

  //DCA
  if(isL){
    if( TMath::Abs(v0->DcaPosToPrimVertex()) < fV0dcaPosLambdaToPV || TMath::Abs(v0->DcaNegToPrimVertex()) < fV0dcaNegLambdaToPV ) { return kFALSE; }
  }
  if(isAL){
    if( TMath::Abs(v0->DcaPosToPrimVertex()) < fV0dcaNegLambdaToPV || TMath::Abs(v0->DcaNegToPrimVertex()) < fV0dcaPosLambdaToPV ) { return kFALSE; }
  }
  if( TMath::Abs(v0->DcaV0Daughters()) > fV0dcaDaugtersLambda ) { return kFALSE; }
  //radius
  Double_t dDecayRadius = v0->RadiusV0();
  if( dDecayRadius < fLambdaradiusMin || dDecayRadius > fLambdaradiusMax ) { return kFALSE; }

  // cosine of pointing angle (CPA)
  Double_t dCPA = v0->CosPointingAngle(fAOD->GetPrimaryVertex());
  if(dCPA < fCutCPALambda) { return kFALSE; }
  fhV0Counter[1]->Fill("CPA OK",1);

  if(fCutTauLambda > 0.0 && ProperLifetime(v0, 1.11568) > fCutTauLambda) { return kFALSE; }
  fhV0Counter[1]->Fill("LT OK",1);
  // daughter PID
  const AliAODTrack* daughterPos = (AliAODTrack*) v0->GetDaughter(0);
  const AliAODTrack* daughterNeg = (AliAODTrack*) v0->GetDaughter(1);

  if(!HasTrackPIDTPC(daughterPos) || !HasTrackPIDTPC(daughterNeg)) { return kFALSE; }
  fhV0Counter[1]->Fill("TPC PID OK",1);

  if (daughterPos->GetTPCsignalN() < fTPCclMin || daughterNeg->GetTPCsignalN() < fTPCclMin) { return kFALSE; }
  fhV0Counter[1]->Fill("TPC sig OK",1);
  Float_t nSigmaPiPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(daughterPos, AliPID::kPion));
  Float_t nSigmaPiNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(daughterNeg, AliPID::kPion));
  Float_t nSigmaPPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(daughterPos, AliPID::kProton));
  Float_t nSigmaPNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(daughterNeg, AliPID::kProton));
  if(isL && (nSigmaPPos > fSigmaTPC || nSigmaPiNeg > fSigmaTPC) ) { return kFALSE; }
  if(isAL && (nSigmaPiPos > fSigmaTPC || nSigmaPNeg > fSigmaTPC) ) { return kFALSE; }
  fhV0Counter[1]->Fill("3Sigma OK",1);

  if(fColSystem==sPP||fColSystem==sPPb){
    Double_t dMassK0s = v0->MassK0Short();
    if(TMath::Abs(dMassK0s - 0.497614) < fMassRejWindowLambda) { return kFALSE; }
  }
  fhV0Counter[1]->Fill("Mass crosscheck OK",1);

  Double_t dMass = 0.0;
  if(isL) dMass = dMassLambda;
  if(isAL) dMass = dMassALambda;
  if(isL && isAL) dMass = 0.0;

  Double_t binscont[4] = {fPVz, fSampleIndex, v0->Pt(), dMass};
  fhTrigTracks[5]->Fill(binscont,0,1.);
  fTracksTrig[5]->Add(new AliPartSimpleForCorr(v0->Eta(),v0->Phi(),v0->Pt(),dMass));
  if(fSkipCorr){
    fhPT[5]->Fill(v0->Pt());
    fhPTvsMinv[1]->Fill(v0->Pt(),dMass);
  }
  return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlowFMD::ProperLifetime(const AliAODv0* v0, const Double_t massPDG) const
{
  Double_t dPrimVtxCoor[3] = {0.0,0.0,0.0}; // primary vertex position {x,y,z}
  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
  primVtx->GetXYZ(dPrimVtxCoor);
  Double_t length = v0->DecayLengthV0(dPrimVtxCoor);

  return ( (massPDG / v0->P() ) * length );
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowFMD::PrepareV0()
{
  Int_t nOfV0s = fAOD->GetNumberOfV0s();
  if(nOfV0s < 1) { return; }

  for(Int_t i(0); i < nOfV0s; i++){
    AliAODv0* v0 = static_cast<AliAODv0*>(fAOD->GetV0(i));
    if(!v0) continue;
    if(!IsV0(v0)) continue;

    if(!IsK0s(v0) && !IsLambda(v0)) continue;
  }

  return;
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
    Double_t binscont[6];
    binscont[2] = fPVz;
    binscont[3] = fSampleIndex;
    binscont[4] = 1.0;

    for(Int_t iTrig(0); iTrig < fTracksTrig[spec]->GetEntriesFast(); iTrig++){
      AliVParticle* track = dynamic_cast<AliVParticle*>(fTracksTrig[spec]->At(iTrig));
      if(!track) continue;
      AliAODTrack* trackAOD = nullptr;
      if(!fIsMC && spec < 4) trackAOD = (AliAODTrack*)fTracksTrig[spec]->At(iTrig);

      Double_t trigPt = track->Pt();
      Double_t trigEta = track->Eta();
      Double_t trigPhi = track->Phi();
      Double_t trigCharge = track->Charge();
      Double_t trigEff = 1.0;
      if(fUseEfficiency) {
        trigEff = GetEff(trigPt, spec, trigEta);
        if(trigEff < 0.001) continue;
      }
      binscont[5] = trigPt;
      if(spec > 3) binscont[4] = track->M();

      for(Int_t iAss(0); iAss < fTracksAss->GetEntriesFast(); iAss++){
        AliVParticle* trackAss = dynamic_cast<AliVParticle*>(fTracksAss->At(iAss));
        if(!trackAss) continue;
        AliAODTrack* trackAODAss = nullptr;
        if(!fIsMC && spec < 4) trackAODAss = (AliAODTrack*)fTracksAss->At(iAss);

        Double_t assPt = trackAss->Pt();
        Double_t assEta = trackAss->Eta();
        Double_t assPhi = trackAss->Phi();
        Double_t assCharge = trackAss->Charge();
        Double_t assEff = 1.0;
        if(fUseEfficiency) {
          assEff = GetEff(assPt, 0, assEta);
          if(assEff < 0.001) continue;
        }

        if(!fIsMC && spec < 4 && trackAOD->GetID() == trackAODAss->GetID()) continue;

        binscont[0] = trigEta - assEta;
        binscont[1] = RangePhi(trigPhi - assPhi);

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
    Double_t binscont[6];
    binscont[2] = fPVz;
    binscont[3] = fSampleIndex;
    binscont[4] = 1.0;

    for(Int_t iTrig(0); iTrig < fTracksTrig[spec]->GetEntriesFast(); iTrig++){
      AliVParticle* track = dynamic_cast<AliVParticle*>(fTracksTrig[spec]->At(iTrig));
      if(!track) continue;

      Double_t trigPt = track->Pt();
      Double_t trigEta = track->Eta();
      Double_t trigPhi = track->Phi();
      Double_t trigEff = 1.0;
      if(fUseEfficiency) {
        trigEff = GetEff(trigPt, spec, trigEta);
        if(trigEff < 0.001) continue;
      }
      binscont[5] = trigPt;
      if(spec > 3) binscont[4] = track->M();

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
    Double_t binscont[4];
    binscont[2] = fPVz;
    binscont[3] = fSampleIndex;

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
        binscont[1] = RangePhiFMD(trigPhi - assPhi);

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
      Double_t binscont[6];
      binscont[2] = fPVz;
      binscont[3] = fSampleIndex;
      binscont[4] = 1.0;

      for(Int_t iTrig(0); iTrig < fTracksTrig[spec]->GetEntriesFast(); iTrig++){
        AliVParticle* track = dynamic_cast<AliVParticle*>(fTracksTrig[spec]->At(iTrig));
        if(!track) continue;

        Double_t trigPt = track->Pt();
        Double_t trigEta = track->Eta();
        Double_t trigPhi = track->Phi();
        Double_t trigCharge = track->Charge();
        binscont[5] = trigPt;
        if(spec > 3) binscont[4] = track->M();
        Double_t trigEff = 1.0;
        if(fUseEfficiency) {
          trigEff = GetEff(trigPt, spec, trigEta);
          if(trigEff < 0.001) continue;
        }

        for(Int_t eMix(0); eMix < nMix; eMix++){
          TObjArray *mixEvents = pool->GetEvent(eMix);
          for(Int_t iAss(0); iAss < mixEvents->GetEntriesFast(); iAss++){
            AliVParticle* trackAss = dynamic_cast<AliVParticle*>(mixEvents->At(iAss));
            if(!trackAss) continue;

            Double_t assPt = trackAss->Pt();
            Double_t assEta = trackAss->Eta();
            Double_t assPhi = trackAss->Phi();
            Double_t assCharge = trackAss->Charge();
            Double_t assEff = 1.0;
            if(fUseEfficiency) {
              assEff = GetEff(assPt, 0, assEta);
              if(assEff < 0.001) continue;
            }

            binscont[0] = trigEta - assEta;
            binscont[1] = RangePhi(trigPhi - assPhi);

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

            fhME[spec]->Fill(binscont,0,1./((Double_t)nMix*(trigEff*assEff)));
          }
        }
      }
    } // end TPC - TPC
    else if(fAnalType == eTPCFMDA || fAnalType == eTPCFMDC){
      Double_t binscont[6];
      binscont[2] = fPVz;
      binscont[3] = fSampleIndex;
      binscont[4] = 1.0;

      for(Int_t iTrig(0); iTrig < fTracksTrig[spec]->GetEntriesFast(); iTrig++){
        AliVParticle* track = dynamic_cast<AliVParticle*>(fTracksTrig[spec]->At(iTrig));
        if(!track) continue;

        Double_t trigPt = track->Pt();
        Double_t trigEta = track->Eta();
        Double_t trigPhi = track->Phi();
        binscont[5] = trigPt;
        if(spec > 3) binscont[4] = track->M();
        Double_t trigEff = 1.0;
        if(fUseEfficiency) {
          trigEff = GetEff(trigPt, spec, trigEta);
          if(trigEff < 0.001) continue;
        }

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

            fhME[spec]->Fill(binscont,0,assMult/((Double_t)nMix*trigEff));
          }
        }
      }
    } // end TPC - FMD
    else{
      Double_t binscont[4];
      binscont[2] = fPVz;
      binscont[3] = fSampleIndex;

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
            binscont[1] = RangePhiFMD(trigPhi - assPhi);

            fhME[spec]->Fill(binscont,0,(trigMult*assMult)/(Double_t)nMix);
          }
        }
      }
    } // end FMD - FMD

  } // event pool done

  if((!fDoPID && !fDoV0 && spec == eCharged) || (fDoPID && !fDoV0 && spec == eProton) || (fDoV0 && spec == eLambda) ){
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
  TString part[6] = {"ch", "pi", "ka", "pr", "K0s", "Lambda"};
  if(fColSystem == sPPb){
    TString etaReg[8] = {"0020", "0200", "0204", "0402", "0406", "0604", "0608", "0806"};
    for(Int_t p(0); p < 6; p++){
      for(Int_t eta(0); eta < 8; eta++){
        if(fDoV0 && p < 4) continue;
        if(fDoPID && !fDoV0 && p > 3) continue;
        fhEfficiencyEta[p][eta] = (TH2D*)fInputListEfficiency->FindObject(Form("LHC17f2b_%s_Eta_%s_%s_wFD",part[p].Data(), etaReg[eta].Data(),fSystematicsFlag.Data()));
        if(!fhEfficiencyEta[p][eta]) {AliError(Form("Efficiency (%s, eta region %s, flag %s) not loaded",part[p].Data(),etaReg[eta].Data(),fSystematicsFlag.Data())); return kFALSE; }
      }
      if(!fDoPID && !fDoV0) break;
    }
    fhEventCounter->Fill("Efficiencies loaded",1);
    return kTRUE;
  }
  else if(fColSystem == sPP){
    for(Int_t p(0); p < 6; p++){
      if(fDoV0 && p < 4) continue;
      if(fDoPID && !fDoV0 && p > 3) continue;
      fhEfficiency[p] = (TH2D*)fInputListEfficiency->FindObject(Form("LHC%s_%s_%s_wFD",ReturnPPperiod(fAOD->GetRunNumber()).Data(),part[p].Data(),fSystematicsFlag.Data()));
      if(!fhEfficiency[p]) {AliError(Form("Efficiency (run %d, part %s, flag %s) not loaded",fAOD->GetRunNumber(),part[p].Data(),fSystematicsFlag.Data())); return kFALSE; }
      if(!fDoPID && !fDoV0) break;
    }
    fhEventCounter->Fill("Efficiencies loaded",1);
    return kTRUE;
  }

  return kFALSE;
}
//_____________________________________________________________________________
TString AliAnalysisTaskCorrForFlowFMD::ReturnPPperiod(const Int_t runNumber) const
{
  if(runNumber >= 252235 && runNumber <= 264347){ // LHC16
    if(runNumber >= 252235 && runNumber <= 252375) return "17f6";
    if(runNumber >= 253437 && runNumber <= 253591) return "17f9";
    if(runNumber >= 254128 && runNumber <= 254332) return "17d17";
    if(runNumber >= 254604 && runNumber <= 255467) return "17f5";
    if(runNumber >= 255539 && runNumber <= 255618) return "17d3";
    if(runNumber >= 256219 && runNumber <= 256418) return "17e5";
    if(runNumber >= 256941 && runNumber <= 258537) return "18f1";
    if(runNumber >= 258962 && runNumber <= 259888) return "18d8";
    if(runNumber >= 262424 && runNumber <= 264035) return "17d16";
    if(runNumber >= 264076 && runNumber <= 264347) return "17d18";
  }

  if(runNumber >= 270581 && runNumber <= 282704){ // LHC17
    if(runNumber >= 270581 && runNumber <= 270667) return "18d3";
    if(runNumber >= 270822 && runNumber <= 270830) return "17h1";
    if(runNumber >= 270854 && runNumber <= 270865) return "18d3";
    if(runNumber >= 271870 && runNumber <= 273103) return "18c12";
    if(runNumber >= 273591 && runNumber <= 274442) return "17k4";
    if(runNumber >= 274593 && runNumber <= 274671) return "17h11";
    if(runNumber >= 274690 && runNumber <= 276508) return "18c13";
    if(runNumber >= 276551 && runNumber <= 278216) return "18a8";
    if(runNumber >= 278914 && runNumber <= 280140) return "17l5";
    if(runNumber >= 280282 && runNumber <= 281961) return "18a9";
    if(runNumber >= 282528 && runNumber <= 282704) return "18a1";
  }

  if(runNumber >= 285009 && runNumber <= 294925){ // LHC18
    if(runNumber >= 285009 && runNumber <= 285396) return "18g4";
    if(runNumber >= 285978 && runNumber <= 286350) return "18g5";
    if(runNumber >= 286380 && runNumber <= 286937) return "18g6";
    if(runNumber >= 287000 && runNumber <= 287658) return "18h2";
    if(runNumber >= 288619 && runNumber <= 289201) return "18h4"; //g,h,i,j,k
    if(runNumber >= 289240 && runNumber <= 289971) return "18j1";
    if(runNumber >= 290323 && runNumber <= 292839) return "18j4";
    if(runNumber >= 293357 && runNumber <= 293359) return "18k1";
    if(runNumber >= 293475 && runNumber <= 293898) return "18k2";
    if(runNumber >= 294009 && runNumber <= 294925) return "18k3";
  }

  AliWarning("PP period identifier was called and based on the run number did not pick up the correct efficiency. Setting up efficiencies from LHC18j4.");
  return "18j4";
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlowFMD::GetEff(const Double_t dPt, const Int_t spec, const Double_t dEta)
{
  if(!fUseEfficiency) return 1.0;
  if(fDoV0 && spec == eCharged) return 1.0;
  if(fColSystem == sPPb){
    Int_t region = GetEtaRegion(dEta);
    if(region < 0) { AliWarning("Invalid region, returning efficiency 1.0."); return 1.0; }
    if(!fhEfficiencyEta[spec][region]) { AliError("Efficiency histogram not found, returning efficiency 1.0."); return 1.0; }
    return fhEfficiencyEta[spec][region]->GetBinContent(fhEfficiencyEta[spec][region]->FindFixBin(dPt, fCentrality));
  }else{
    if(!fhEfficiency[spec]) { AliError("Efficiency histogram not found, returning efficiency 1.0."); return 1.0; }
    return fhEfficiency[spec]->GetBinContent(fhEfficiency[spec]->FindFixBin(dPt, fCentrality));
  }

  return 1.0;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskCorrForFlowFMD::GetEtaRegion(const Double_t dEta){
  if(TMath::Abs(dEta) > 0.8) { AliWarning("Eta out of range!"); return -1; }
  if(dEta > 0.0){
    if(dEta > 0.6) return 6;
    if(dEta > 0.4) return 4;
    if(dEta > 0.2) return 2;
    return 0;
  }
  else{
    if(dEta < -0.6) return 7;
    if(dEta < -0.4) return 5;
    if(dEta < -0.2) return 3;
    return 1;
  }

  return -1;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowFMD::CreateTHnCorrelations(){
  Int_t nSteps = 1;
  Double_t binning_dphi[73] = { -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0,       0.087266,  0.174533,  0.261799,  0.349066,  0.436332, 0.523599,  0.610865,  0.698132,  0.785398,  0.872665,  0.959931, 1.047198,  1.134464,  1.221730,  1.308997,  1.396263,  1.483530, 1.570796,  1.658063,  1.745329,  1.832596,  1.919862,  2.007129, 2.094395,  2.181662,  2.268928,  2.356194,  2.443461,  2.530727, 2.617994,  2.705260,  2.792527,  2.879793,  2.967060,  3.054326, 3.141593,  3.228859,  3.316126,  3.403392,  3.490659,  3.577925, 3.665191,  3.752458,  3.839724,  3.926991,  4.014257,  4.101524, 4.188790,  4.276057,  4.363323,  4.450590,  4.537856,  4.625123, 4.712389};
  const Int_t sizePtTrig = fPtBinsTrigCharged.size() - 1;
  const Int_t sizeOfSamples = (Int_t) fNOfSamples;
  const Int_t sizeMbins = fNbinsMinv;

  Double_t min[6] = {0.0, 0.0, 0.0, 0.0, 0.44, 1.08};
  Double_t max[6] = {2.0, 2.0, 2.0, 2.0, 0.56, 1.15};

  TString nameS[6] = {"fhChargedSE", "fhPidSE_Pion", "fhPidSE_Kaon", "fhPidSE_Proton", "fhPidSE_K0s", "fhPidSE_Lambda"};
  TString nameM[6] = {"fhChargedME", "fhPidME_Pion", "fhPidME_Kaon", "fhPidME_Proton", "fhPidME_K0s", "fhPidME_Lambda"};

  if(fAnalType == eTPCFMDA || fAnalType == eTPCFMDC){
    Double_t binning_detaFMDTPC[]={-6.,-5.8, -5.6, -5.4, -5.2, -5.0, -4.8, -4.6, -4.4, -4.2, -4., -3.8, -3.6, -3.4, -3.2, -3., -2.8, -2.6, -2.4, -2.2, -2., -1.8, -1.6, -1.4, -1.2, -1., -0.8};
    Double_t binning_detaFMDCTPC[]={ 1., 1.2, 1.4, 1.6, 1.8, 2. , 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4.};

    Int_t iTrackBin_tpcfmdA[] = {26, 72, 10, sizeOfSamples, sizeMbins, sizePtTrig};
    Int_t iTrackBin_tpcfmdC[] = {15, 72, 10, sizeOfSamples, sizeMbins, sizePtTrig};
    Int_t nTrackBin_tpcfmd = sizeof(iTrackBin_tpcfmdA) / sizeof(Int_t);

    for(Int_t i(0); i < 6; i++){
      if(!fDoPID && i > 0 && i < 4) continue;
      if(!fDoV0 && i > 3) continue;
      if(fIsAntiparticleCheck && (i == 0 || i == 4)) continue;

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
      fhSE[i]->SetBinLimits(1, binning_dphi);
      fhME[i]->SetBinLimits(1, binning_dphi);
    }


  } // end TPC - FMD
  else if(fAnalType == eFMDAFMDC){
    // Int_t iTrackBin_fmdAfmdC[] = {48, 72, 10};
    Int_t iTrackBin_fmdAfmdC[] = {24, 20, 10, sizeOfSamples};
    Int_t nTrackBin_fmdAfmdC = sizeof(iTrackBin_fmdAfmdC) / sizeof(Int_t);

    // FMD only for unidentified
    for(Int_t i(0); i < 1; i++){
      fhSE[i] = new AliTHn(nameS[i], nameS[i], nSteps, nTrackBin_fmdAfmdC, iTrackBin_fmdAfmdC);
      fhSE[i]->SetBinLimits(0,3.4,8.2);
      fhSE[i]->SetBinLimits(1,-0.55*TMath::Pi(), 1.45*TMath::Pi());

      fhME[i] = new AliTHn(nameM[i], nameM[i], nSteps, nTrackBin_fmdAfmdC, iTrackBin_fmdAfmdC);
      fhME[i]->SetBinLimits(0,3.4,8.2);
      fhME[i]->SetBinLimits(1,-0.55*TMath::Pi(), 1.45*TMath::Pi());
    }
  } // end FMD - FMD
  else {
    Double_t binning_deta_tpctpc[33] = {-1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,    0.1,  0.2,  0.3,  0.4,  0.5, 0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5, 1.6};
    Int_t iBinningTPCTPC[] = {32,72,10,sizeOfSamples,sizeMbins,sizePtTrig};
    Int_t nTrackBin_tpctpc = sizeof(iBinningTPCTPC) / sizeof(Int_t);

    for(Int_t i(0); i < 6; i++){
      if(!fDoPID && i > 0 && i < 4) continue;
      if(!fDoV0 && i > 3) continue;
      if(fIsAntiparticleCheck && (i == 0 || i == 4)) continue;

      fhSE[i] = new AliTHn(nameS[i], nameS[i], nSteps, nTrackBin_tpctpc, iBinningTPCTPC);
      fhSE[i]->SetBinLimits(0, binning_deta_tpctpc);
      fhSE[i]->SetBinLimits(1, binning_dphi);

      fhME[i] = new AliTHn(nameM[i], nameM[i], nSteps, nTrackBin_tpctpc, iBinningTPCTPC);
      fhME[i]->SetBinLimits(0, binning_deta_tpctpc);
      fhME[i]->SetBinLimits(1, binning_dphi);
    }

  } // end TPC - TPC

  // all
  for(Int_t i(0); i < 6; i++){
    if(fAnalType == eFMDAFMDC && i > 0) break;
    if(!fDoPID && i > 0 && i < 4) continue;
    if(!fDoV0 && i > 3) continue;
    if(fIsAntiparticleCheck && (i == 0 || i == 4)) continue;

    fhSE[i]->SetBinLimits(2, -10,10);
    fhSE[i]->SetBinLimits(3, 0,10);
    fhSE[i]->SetVarTitle(0, "#Delta#eta");
    fhSE[i]->SetVarTitle(1, "#Delta#phi");
    fhSE[i]->SetVarTitle(2, "PVz [cm]");
    fhSE[i]->SetVarTitle(3, "Sample");
    if(fAnalType != eFMDAFMDC){
      fhSE[i]->SetBinLimits(4, min[i], max[i]);
      fhSE[i]->SetBinLimits(5, fPtBinsTrigCharged.data());
      fhSE[i]->SetVarTitle(4, "Mass [GeV] (trig)");
      fhSE[i]->SetVarTitle(5, "p_{T} [GeV/c] (trig)");
    }
    fOutputListCharged->Add(fhSE[i]);

    fhME[i]->SetBinLimits(2, -10,10);
    fhME[i]->SetBinLimits(3, 0,10);
    fhME[i]->SetVarTitle(0, "#Delta#eta");
    fhME[i]->SetVarTitle(1, "#Delta#phi");
    fhME[i]->SetVarTitle(2, "PVz [cm]");
    fhME[i]->SetVarTitle(3, "Sample");
    if(fAnalType != eFMDAFMDC){
      fhME[i]->SetBinLimits(4, min[i], max[i]);
      fhME[i]->SetBinLimits(5, fPtBinsTrigCharged.data());
      fhME[i]->SetVarTitle(4, "Mass [GeV] (trig)");
      fhME[i]->SetVarTitle(5, "p_{T} [GeV/c] (trig)");
    }
    fOutputListCharged->Add(fhME[i]);
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowFMD::PrepareTPCTracks(){
  if(!fAOD) return kFALSE;
  if(!fTracksAss || !fTracksTrig[0] || !fhTrigTracks[0]) {AliError("Cannot prepare TPC tracks!"); return kFALSE; }

  fNofTracks = 0;
  Int_t iNofHighPtParticles = 0;
  Double_t binscont[3] = {fPVz, fSampleIndex, 0.};

  for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {
      AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
      if(!track || !IsTrackSelected(track)) { continue; }

      Double_t trackPt = track->Pt();
      if(trackPt > 10.) iNofHighPtParticles++;
      if(trackPt > fPtMinAss && trackPt < fPtMaxAss) {
        fNofTracks++;
        if(fAnalType == eFMDAFMDC || fIsTPCgen) continue;
        if(fAnalType == eTPCTPC) fTracksAss->Add((AliAODTrack*)track); // only if associated from TPC
      }
      if(fAnalType != eFMDAFMDC){
        Double_t trackEta = track->Eta();
        if(trackPt > fPtMinTrig && trackPt < fPtMaxTrig) {
          if(fUseOppositeSidesOnly){
            if(fAnalType == eTPCFMDA && trackEta > 0.0) continue;
            if(fAnalType == eTPCFMDC && trackEta < 0.0) continue;
          }

          binscont[2] = trackPt;
          if(!fIsAntiparticleCheck){
            fhPT[0]->Fill(trackPt);
            fTracksTrig[0]->Add((AliAODTrack*)track);
            fhTrigTracks[0]->Fill(binscont,0,1.);
          }

          if(fDoPID){
            Int_t trackPid = IdentifyTrack(track);
            if(trackPid > 0 && trackPid < 4){
              if(fIsAntiparticleCheck && fDoAntiparticleOnly && track->Charge() > 0.) continue;
              fTracksTrig[trackPid]->Add((AliAODTrack*)track);
              fhTrigTracks[trackPid]->Fill(binscont,0,1.);
              if(fSkipCorr) fhPT[trackPid]->Fill(trackPt);
            }
          }
        }
      } // POI from TPC
  } // tracks loop end

  if(fUseNch){
    if(fNofTracks < fNchMin || fNofTracks > fNchMax) { return kFALSE; }
    fhEventCounter->Fill("Nch cut ok ",1);
    fhEventMultiplicity->Fill(fNofTracks);
  }

  if(fRejectHighPtEvents){
    fhEventCounter->Fill("Before HPC",1); //HPC = high pt cut
    if(iNofHighPtParticles > fNofMinHighPtTracksForRejection) return kFALSE;
    fhEventCounter->Fill("After HPC",1);
  }

  if(fDoV0){
    PrepareV0();
  }

  return kTRUE;
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

  Double_t binscontFMD[2] = {fPVz,fSampleIndex};

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
           if(fIsFMDgen) continue;
           if(eta > fFMDAacceptanceCutLower && eta < fFMDAacceptanceCutUpper){
             if(fAnalType == eTPCFMDA) {
               fTracksAss->Add(new AliPartSimpleForCorr(eta,phi,mostProbableN));
               fHistFMDeta->Fill(eta,fPVz,mostProbableN);
             }
             if(fAnalType == eFMDAFMDC) {
               fTracksTrig[0]->Add(new AliPartSimpleForCorr(eta,phi,mostProbableN));
               fhTrigTracks[0]->Fill(binscontFMD,0,mostProbableN);
               fHistFMDeta->Fill(eta,fPVz,mostProbableN);
             }
           }
    	   } // eta positive
         else
         {
    	     nFMD_bwd_hits+=mostProbableN;
           if(fIsFMDgen) continue;
           if(eta < -fFMDCacceptanceCutLower && eta > -fFMDCacceptanceCutUpper){
             if(fAnalType == eTPCFMDC || fAnalType == eFMDAFMDC) {
               fTracksAss->Add(new AliPartSimpleForCorr(eta,phi,mostProbableN));
               fHistFMDeta->Fill(eta,fPVz,mostProbableN);
             }
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
Bool_t AliAnalysisTaskCorrForFlowFMD::PrepareMCTracks(){
  if(!fTracksAss || !fTracksTrig[0] || !fhTrigTracks[0]) {AliError("Cannot prepare MCC tracks!"); return kFALSE; }

  AliMCEvent* mcEvent = dynamic_cast<AliMCEvent*>(MCEvent());
  if(!mcEvent) return kFALSE;

  Double_t binscont[3] = {fPVz, fSampleIndex, 0.};
  Double_t binscontFMD[2] = {fPVz, fSampleIndex};

  for(Int_t i(0); i < mcEvent->GetNumberOfTracks(); i++) {
    AliMCParticle* part = (AliMCParticle*)mcEvent->GetTrack(i);
    if(!part->IsPhysicalPrimary()) continue;
    Double_t partEta = part->Eta();
    Double_t partPt = part->Pt();
    Double_t partPhi = part->Phi();
    Double_t partRapidity = part->Y();
    binscont[2] = partPt;

    if(fBoostAMPT) {
      partEta=TransverseBoost(part);
      partRapidity=partRapidity-0.465;
    }

    // TPC region
    if(TMath::Abs(partEta) < 0.8){
      if(!fIsTPCgen) continue;
      Int_t partPDG = TMath::Abs(part->PdgCode());
      Int_t partIdx = -1;
      if(partPDG == 211) partIdx = 1;
      else if(partPDG == 321) partIdx = 2;
      else if(partPDG == 2212) partIdx = 3;
      else if(partPDG == 310) partIdx = 4;
      else if(partPDG == 3122) partIdx = 5;

      if(partIdx < 4 && part->Charge()==0.) continue;

      if(fAnalType == eTPCTPC){
        if(partPt > fPtMinTrig && partPt < fPtMaxTrig){
          fTracksTrig[0]->Add((AliMCParticle*)part);
          fhTrigTracks[0]->Fill(binscont,0,1.);
          if(fDoPID && partIdx > 0 && partIdx < 4){
            fTracksTrig[partIdx]->Add((AliMCParticle*)part);
            fhTrigTracks[partIdx]->Fill(binscont,0,1.);
          }
          if(fDoV0 && partIdx > 3 && TMath::Abs(partRapidity)<0.5){
            Double_t binscontV0[4] = {fPVz, fSampleIndex, part->Pt(), part->M()};
            fhTrigTracks[partIdx]->Fill(binscontV0,0,1.);
            fTracksTrig[partIdx]->Add(new AliPartSimpleForCorr(part->Eta(),part->Phi(),part->Pt(),part->M()));
          }
        }
        if(partPt > fPtMinAss && partPt < fPtMaxAss) fTracksAss->Add((AliMCParticle*)part);
      } // end TPCTPC
      if(fAnalType == eTPCFMDA || fAnalType == eTPCFMDC){
        if(partPt > fPtMinTrig && partPt < fPtMaxTrig){
          fTracksTrig[0]->Add((AliMCParticle*)part);
          fhTrigTracks[0]->Fill(binscont,0,1.);
          if(fDoPID && partIdx > 0 && partIdx < 4){
            fTracksTrig[partIdx]->Add((AliMCParticle*)part);
            fhTrigTracks[partIdx]->Fill(binscont,0,1.);
          }
          if(fDoV0 && partIdx > 3 && TMath::Abs(partRapidity)<0.5){
            Double_t binscontV0[4] = {fPVz, fSampleIndex, part->Pt(), part->M()};
            fhTrigTracks[partIdx]->Fill(binscontV0,0,1.);
            fTracksTrig[partIdx]->Add(new AliPartSimpleForCorr(part->Eta(),part->Phi(),part->Pt(),part->M()));
          }
        }
      }
    } // end eta within 0.8
    else if(partEta > fFMDAacceptanceCutLower && partEta < fFMDAacceptanceCutUpper){
      if(!fIsFMDgen) continue;
      if(fAnalType == eTPCFMDA) {
        fTracksAss->Add(new AliPartSimpleForCorr(partEta,partPhi,1.));
        fHistFMDeta->Fill(partEta,fPVz,1.);
      }
      if(fAnalType == eFMDAFMDC) {
        fTracksTrig[0]->Add(new AliPartSimpleForCorr(partEta,partPhi,1.));
        fhTrigTracks[0]->Fill(binscontFMD,0,1.);
        fHistFMDeta->Fill(partEta,fPVz,1.);
      }
    } // end eta within FMDA range
    else if(partEta < -fFMDCacceptanceCutLower && partEta > -fFMDCacceptanceCutUpper){
      if(!fIsFMDgen) continue;
      if(fAnalType == eTPCFMDC || fAnalType == eFMDAFMDC) {
        fTracksAss->Add(new AliPartSimpleForCorr(partEta,partPhi,1.));
        fHistFMDeta->Fill(partEta,fPVz,1.);
      }
    } // end eta within FMDC range

  } // end MC track loop


  return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlowFMD::TransverseBoost(const AliMCParticle *track){
  Float_t boost=0.465;
  Float_t beta=TMath::TanH(boost);
  Float_t gamma=1./TMath::Sqrt((1.-TMath::Power(beta,2)));

  Float_t energy=track->E();
  Float_t mass=track->M();
  Float_t px=track->Px();
  Float_t py=track->Py();
  Float_t pz=track->Pz();
  Float_t mT=TMath::Sqrt(TMath::Power(energy,2)-TMath::Power(pz,2));
  Float_t eta=track->Eta();
  Float_t rap=track->Y();

  Float_t energy_boosted=gamma*energy-gamma*beta*pz;
  Float_t pz_boosted=-gamma*beta*energy+gamma*pz;
  Float_t mT_boosted=TMath::Sqrt(TMath::Power(energy_boosted,2)-TMath::Power(pz_boosted,2));
  Float_t rap_boosted=rap-boost;
  Float_t numerator=TMath::Sqrt(TMath::Power(mT_boosted,2)*TMath::Power(TMath::CosH(rap_boosted),2)-TMath::Power(mass,2))+mT_boosted*TMath::SinH(rap_boosted);
  Float_t denumerator=TMath::Sqrt(TMath::Power(mT_boosted,2)*TMath::Power(TMath::CosH(rap_boosted),2)-TMath::Power(mass,2))-mT_boosted*TMath::SinH(rap_boosted);
  Double_t eta_boosted = 0.5*TMath::Log(numerator/denumerator);

  return eta_boosted;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowFMD::PrintSetup(){
  printf("\n\n\n ************** Parameters ************** \n");
  printf("\t fAnalType: (Int_t) %d\n", fAnalType);
  printf("\t fColSystem: (Int_t) %d\n", fColSystem);
  printf("\t fIsMC: (Bool_t) %s\n", fIsMC ? "kTRUE" : "kFALSE");
  printf("\t fIsTPCgen: (Bool_t) %s\n", fIsTPCgen ? "kTRUE" : "kFALSE");
  printf("\t fIsFMDgen: (Bool_t) %s\n", fIsFMDgen ? "kTRUE" : "kFALSE");
  printf("\t fDoPID: (Bool_t) %s\n", fDoPID ? "kTRUE" : "kFALSE");
  printf("\t fDoV0: (Bool_t) %s\n", fDoV0 ? "kTRUE" : "kFALSE");
  printf("\t fUseNch: (Bool_t) %s\n", fUseNch ? "kTRUE" : "kFALSE");
  printf("\t fIsHMpp: (Bool_t) %s\n", fIsHMpp ? "kTRUE" : "kFALSE");
  printf("\t fUseEfficiency: (Bool_t) %s\n",  fUseEfficiency ? "kTRUE" : "kFALSE");
  printf("\t fUseOppositeSidesOnly: (Bool_t) %s\n", fUseOppositeSidesOnly ? "kTRUE" : "kFALSE");
  printf("\t fIsAntiparticleCheck: (Bool_t) %s\n", fIsAntiparticleCheck ? "kTRUE" : "kFALSE");
  printf("\t fDoAntiparticleOnly: (Bool_t) %s\n", fDoAntiparticleOnly ? "kTRUE" : "kFALSE");
  printf("\t fRejectHighPtEvents: (Bool_t) %s\n", fRejectHighPtEvents ? "kTRUE" : "kFALSE");
  printf("\t fNOfSamples: (Int_t) %d\n", (Int_t) fNOfSamples);
  printf(" **************************** \n");
  printf("\t fSystematicsFlag: (TString) %s\n", fSystematicsFlag.Data());
  printf("\t fAbsEtaMax: (Double_t) %f\n", fAbsEtaMax);
  printf("\t fPtMinTrig -- fPtMaxTrig: (Double_t) %f -- %f\n", fPtMinTrig, fPtMaxTrig);
  printf("\t fPtMinAss -- fPtMaxAss: (Double_t) %f -- %f\n", fPtMinAss, fPtMaxAss);
  printf("\t fCentMin -- fCentMax: (Double_t) %f -- %f\n", fCentMin, fCentMax);
  printf("\t fPVzCut: (Double_t) %f\n", fPVzCut);
  printf("\t PID cuts (pion, kaon, proton): (Double_t) %f \t %f \t %f\n", fPIDbayesPion, fPIDbayesKaon, fPIDbayesProton);
  printf(" **************************** \n");
  printf("\t fUseFMDcut: (Bool_t) %s\n", fUseFMDcut ? "kTRUE" : "kFALSE");
  printf("\t fFMDcutapar0 -- fFMDcutapar1: (Double_t) %f -- %f\n", fFMDcutapar0, fFMDcutapar1);
  printf("\t fFMDcutcpar0 -- fFMDcutcpar1: (Double_t) %f -- %f\n", fFMDcutcpar0, fFMDcutcpar1);
  printf("\t fFMDacceptanceCut A - lower, upper: (Double_t) %f, %f\n", fFMDAacceptanceCutLower, fFMDAacceptanceCutUpper);
  printf("\t fFMDacceptanceCut C - lower, upper: (Double_t) %f, %f\n", fFMDCacceptanceCutLower, fFMDCacceptanceCutUpper);
}
