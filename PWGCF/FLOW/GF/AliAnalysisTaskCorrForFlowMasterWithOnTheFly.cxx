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

#include "AliAnalysisTaskCorrForFlowMaster_Anders.h"

using namespace std;

ClassImp(AliAnalysisTaskCorrForFlowMaster_Anders);

AliAnalysisTaskCorrForFlowMaster_Anders::AliAnalysisTaskCorrForFlowMaster_Anders() : AliAnalysisTaskSE(),
    fAOD(0),
    fOutputListCharged(0),
    fInputListEfficiency(0),
    fTracksAss(0),
    fPoolMgr(0),
    fhEventCounter(0),
    fhEventMultiplicity(0),
    fAnalType(eTPCFMDA),
    fColSystem(sPPb),
    fTrigger(AliVEvent::kINT7),
    fIsMC(kFALSE),
    fIsTPCgen(kFALSE),
    fIsHMpp(kFALSE),
    fUseNch(kFALSE),
    fUseEventBias(kFALSE),
    fUseEfficiency(kFALSE),
    fUseOppositeSidesOnly(kFALSE),
    fUseCentralityCalibration(kFALSE),
    fSkipCorr(kFALSE),
    fVetoJetEvents(kFALSE),
    fRejectSecondariesFromMC(kFALSE),
    fBoostAMPT(kFALSE),
    fCreateQAPlots(kFALSE),
    fUseLikeSign(kFALSE),
    fUseUnlikeSign(kFALSE),
    fselectjetsinTPC(kFALSE),
    fJetvetoselectionval(0.5),
    fFilterBit(96),
    fbSign(0),
    fRunNumber(-1),
    fNofTracks(0),
    fNofTrackGlobal(0),
    fNofEventGlobal(0),
    fNofMinHighPtTracksForRejection(0),
    fNchMin(0),
    fNchMax(100000),
    fnTPCcrossedRows(70),
    fNumEventBias(2),
    fNOfSamples(1.0),
    fSampleIndex(0.0),
    fPtMinTrig(0.5),
    fPtMaxTrig(10.0),
    fPtMinAss(0.5),
    fPtMaxAss(1.5),
    fPtRefMin(0.2),
    fPtRefMax(3.0),
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
    fSigmaTPC(3.),
    fJetParticleLowPt(5.),
    fCentEstimator("V0M"),
    fSystematicsFlag(""),
    fPoolMaxNEvents(1000),
    fPoolMinNTracks(50000),
    fMinEventsToMix(5),
    fNzVtxBins(10),
    fNCentBins(15),
    fMergingCut(0.0),
    fUsePhiStar(kFALSE),
    fOnTheFly(kFALSE),
    fGenerator("AMPT"),
    fImpactParameterMC(-1.0)
{}
//_____________________________________________________________________________
AliAnalysisTaskCorrForFlowMaster_Anders::AliAnalysisTaskCorrForFlowMaster_Anders(const char* name, Bool_t bUseEff, Bool_t bUseCalib) : AliAnalysisTaskSE(name),
    fAOD(0),
    fOutputListCharged(0),
    fInputListEfficiency(0),
    fTracksAss(0),
    fPoolMgr(0),
    fhEventCounter(0),
    fhEventMultiplicity(0),
    fAnalType(eTPCFMDA),
    fColSystem(sPPb),
    fTrigger(AliVEvent::kINT7),
    fIsMC(kFALSE),
    fIsTPCgen(kFALSE),
    fIsHMpp(kFALSE),
    fUseNch(kFALSE),
    fUseEventBias(kFALSE),
    fUseEfficiency(kFALSE),
    fUseOppositeSidesOnly(kFALSE),
    fUseCentralityCalibration(kFALSE),
    fSkipCorr(kFALSE),
    fVetoJetEvents(kFALSE),
    fRejectSecondariesFromMC(kFALSE),
    fBoostAMPT(kFALSE),
    fCreateQAPlots(kFALSE),
    fUseLikeSign(kFALSE),
    fUseUnlikeSign(kFALSE),
    fselectjetsinTPC(kFALSE),
    fJetvetoselectionval(0.5),
    fFilterBit(96),
    fbSign(0),
    fRunNumber(-1),
    fNofTracks(0),
    fNofTrackGlobal(0),
    fNofEventGlobal(0),
    fNofMinHighPtTracksForRejection(0),
    fNchMin(0),
    fNchMax(100000),
    fnTPCcrossedRows(70),
    fNumEventBias(2),
    fNOfSamples(1.0),
    fSampleIndex(0.0),
    fPtMinTrig(0.5),
    fPtMaxTrig(10.0),
    fPtMinAss(0.5),
    fPtMaxAss(1.5),
    fCentMin(0.0),
    fCentMax(10.0),
    fPtRefMin(0.2),
    fPtRefMax(3.0),
    fCentrality(-10.0),
    fAbsEtaMax(0.8),
    fPVz(100.0),
    fPVzCut(10.0),
    fTPCclMin(70.),
    fCutDCAz(0.),
    fCutDCAxySigma(0.),
    fCutTPCchi2pCl(0.),
    fSigmaTPC(3.),
    fJetParticleLowPt(5.),
    fCentEstimator("V0M"),
    fSystematicsFlag(""),
    fPoolMaxNEvents(2000),
    fPoolMinNTracks(50000),
    fMinEventsToMix(5),
    fNzVtxBins(10),
    fNCentBins(15),
    fMergingCut(0.0),
    fUsePhiStar(kFALSE),
    fOnTheFly(kFALSE),
    fGenerator("AMPT"),
    fImpactParameterMC(-1.0)
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
AliAnalysisTaskCorrForFlowMaster_Anders::~AliAnalysisTaskCorrForFlowMaster_Anders()
{}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowMaster_Anders::UserCreateOutputObjects()
{
    OpenFile(1);
    PrintSetup();


    fzVtxBins = {-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0};

    fOutputListCharged = new TList();
    fOutputListCharged->SetOwner(kTRUE);

    fhEventCounter = new TH1D("fhEventCounter","Event Counter",10,0,10);
    fOutputListCharged->Add(fhEventCounter);

    fhEventMultiplicity = new TH1D("fhEventMultiplicity","Event multiplicity; N_{ch}",200,0,200);
    fOutputListCharged->Add(fhEventMultiplicity);

    if(fOnTheFly)
    {
      printf("Creating OTF objects\n");
      printf("Generator is %s\n",fGenerator.Data());
      if(centralitymap.empty() && fGenerator.EqualTo("AMPT")) {
        vector<double> b = {0.0,3.72,5.23,7.31,8.88,10.20,11.38,12.47,13.50,14.51,100.0};
        vector<double> cent = {0.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,100.0};
        for(size_t i(0); i<b.size(); ++i) centralitymap[b[i]]=cent[i];
      }
      //fIP = new TH1D("fIP","Impact parameter",1000,0.0,30.0);
      printf("OTF objects created\n");
    }

    const Int_t sizePtTrig = fPtBinsTrigCharged.size() - 1;
    const Int_t sizeOfSamples = (Int_t) fNOfSamples;
    Int_t binsTrig[] = {10, 10, sizePtTrig,100,100};
    fhTrigTracks = new AliTHn("fhTrigTracks","fhTrigTracks", 1, 5, binsTrig);
    fhTrigTracks->SetBinLimits(0,-10,10);
    fhTrigTracks->SetBinLimits(1,0,10);
    fhTrigTracks->SetBinLimits(3, 0,TMath::TwoPi());
    fhTrigTracks->SetBinLimits(4, -fAbsEtaMax, fAbsEtaMax);
    
    fhTrigTracks->SetVarTitle(0, "PVz [cm]");
    fhTrigTracks->SetVarTitle(1, "Sample");
    fhTrigTracks->SetBinLimits(2,fPtBinsTrigCharged.data());
    fhTrigTracks->SetVarTitle(2, "p_{T} (trig)");
    fhTrigTracks->SetVarTitle(3,"#eta");
    fhTrigTracks->SetVarTitle(4,"#phi");
    fOutputListCharged->Add(fhTrigTracks);
    


    //mixing
    fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins.data(), fNzVtxBins, fzVtxBins.data());
    if (!fPoolMgr) { AliError("Event Pool manager not created!"); return; }
    fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);

    if(!fSkipCorr) CreateTHnCorrelations();

    if(fCreateQAPlots){
      fhPT = new TH1D("PT", "PT", 1000, 0, 10);
      fhPT->Sumw2();
      fOutputListCharged->Add(fhPT);


      fhPhi = new TH1D("Phi", "Phi", 100, 0, TMath::TwoPi());
      fhPhi->Sumw2();
      fOutputListCharged->Add(fhPhi);


      fhEta = new TH1D("Eta", "Eta", 100, -fAbsEtaMax, fAbsEtaMax);
      fhEta->Sumw2();
      fOutputListCharged->Add(fhEta);

      fhPVz = new TH1D("PVz", "PVz", 100, -fPVzCut, fPVzCut);
      fhPVz->Sumw2();
      fOutputListCharged->Add(fhPVz);
        //Profile histogram for correlation between pt1,pt2 and Dphi for Mixed events
        //For investigating dip at Dphi~0.
      phDphiPt1Pt2_ME = new TProfile2D("DphiPt1Pt2_ME","DphiPt1Pt2_ME",40,fPtRefMin,fPtRefMax,40,fPtRefMin,fPtRefMax);
      fOutputListCharged->Add(phDphiPt1Pt2_ME);
        
      phDphiPt1Pt2_PtCut_ME = new TProfile2D("DphiPt1Pt2_PtCut_ME","DphiPt1Pt2_PtCut_ME",40,fPtRefMin,fPtRefMax,40,fPtRefMin,fPtRefMax);
      fOutputListCharged->Add(phDphiPt1Pt2_PtCut_ME);
        
        phiCorr_ME = new TH2D("phi_v_eta","phi_v_eta",72, 0, TMath::TwoPi(),18, -0.8, 0.8);
        phiCorr_ME->SetTitle("phi_v_eta;#phi;#eta;N");
        fOutputListCharged->Add(phiCorr_ME);
        
        phiCorr_SE = new TH2D("phiCorr_SE","phiCorr_SE",40, 0, TMath::TwoPi(),40, 0, TMath::TwoPi());
        phiCorr_SE->SetTitle("phiCorr_SE;phi_Trig;phi_Asso;N");
        fOutputListCharged->Add(phiCorr_SE);
    }

    
    
    
      
    

    if(fUseEfficiency) {
      fInputListEfficiency = (TList*) GetInputData(1);
      if(fAbsEtaMax > 0.8) AliWarning("Efficiency loading -- eta can be out of range!");
      if(fSystematicsFlag.IsNull()) fSystematicsFlag = "Ev0_Tr0";
      if(fColSystem == sPPb && fAnalType != eFMDAFMDC && !AreEfficienciesLoaded()) { AliError("Efficiencies not loaded!"); return; }
    }

    if(fUseCentralityCalibration){
      if(fUseEfficiency) fhCentCalib = (TH1D*) GetInputData(2);
      else fhCentCalib = (TH1D*) GetInputData(1);
      if(!fhCentCalib) { AliError("Centrality calibration histogram not loaded!"); return; }
    }

    PostData(1, fOutputListCharged);
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowMaster_Anders::UserExec(Option_t *)
{
    //Fill event counter if event loads successfully.
    fhEventCounter->Fill("Input",1);
    
    //Load event.
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    //Skip if event was not loaded.
    if(!fAOD) { AliError("Event not loaded."); return; }
    //Skip if event fails event cuts.
    if(!IsEventSelected()) { return; }
    //Initiate track counter.
    fNofTrackGlobal = 0;
    //Increment event counter.
    fNofEventGlobal++;
    
    //Fill primary z-vertex if QA-plots are to be created.
    if(fCreateQAPlots){
      fhPVz->Fill(fPVz);
    }
    
    //Get the number of tracks in event.
    Int_t iTracks(fAOD->GetNumberOfTracks());
    //Skip this event if there are no tracks.
    if(iTracks < 1 ) {
      AliWarning("No tracks in the event.");
      return;
    }
    
    //Sample random int for event (Necessary for bootstrapping).
    fSampleIndex = gRandom->Uniform(0,fNOfSamples);
    //Initiate arrays for holding trigger and associate tracks.
    fTracksAss = new TObjArray;
    fTracksTrig = new TObjArray;

    //If efficiencies are to be used, but cannot be loaded, return.
    if(fUseEfficiency && fColSystem == sPP && (fRunNumber != fAOD->GetRunNumber()) && !AreEfficienciesLoaded()) { return; }

    //If is not TPC-tracks are to be used, run prepare function.
    if(!fIsTPCgen || fUseNch)  {
      if(!PrepareTPCTracks()){
            if(fTracksTrig) delete fTracksTrig;
            PostData(1, fOutputListCharged);
            return;
        }
    }
    //If is MC-tracks run prepare funciton.
    if(fIsMC||fOnTheFly){
      if(!PrepareMCTracks()){
            if(fTracksTrig) delete fTracksTrig;
            PostData(1, fOutputListCharged);
            return;
      }
    } // end MC
    
    //Correlate particle tracks (SE and ME).
    if(!fTracksAss->IsEmpty() && !fSkipCorr){
        FillCorrelations();
        FillCorrelationsMixed();
        
        //Clear and delete trigger track array.
        fTracksTrig->Clear();
        delete fTracksTrig;
    }
    
    //If use efficiencies get event runnumber.
    if(fUseEfficiency) fRunNumber = fAOD->GetRunNumber();
    
    //Clear and delete associate tracks.
    fTracksAss->Clear();
    delete fTracksAss;

    PostData(1, fOutputListCharged);
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowMaster_Anders::Terminate(Option_t *)
{
  //if(fPoolMgr) delete fPoolMgr;
  //if(fOutputListCharged) delete fOutputListCharged;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowMaster_Anders::IsEventSelected()
{
  fhEventCounter->Fill("EventOK",1);

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();
    if(!(fSelectMask & fTrigger)) {return kFALSE; }
  fhEventCounter->Fill("TriggerOK",1);

  if(fIsHMpp) fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true);
  if(!fEventCuts.AcceptEvent(fAOD)) { return kFALSE; }
  fhEventCounter->Fill("CutsOK",1);

  //Unfortunately, this has to be done in Prepare MC when running on the fly...
  if(!fOnTheFly){
    AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
    if(!multSelection) { return kFALSE; }
    fhEventCounter->Fill("MultOK",1);
    if(!fUseCentralityCalibration){
      Float_t dPercentile = multSelection->GetMultiplicityPercentile(fCentEstimator);
      if(dPercentile > 100 || dPercentile < 0) { return kFALSE; }
      fhEventCounter->Fill("PercOK",1);
      fCentrality = (Double_t) dPercentile;
    }else{
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
  }

  //AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
  //if(!multSelection) { return kFALSE; }
  //fhEventCounter->Fill("MultOK",1);
  //if(!fUseCentralityCalibration){
  //  Float_t dPercentile = multSelection->GetMultiplicityPercentile(fCentEstimator);
  //  if(dPercentile > 100 || dPercentile < 0) { return kFALSE; }
  //  fhEventCounter->Fill("PercOK",1);
  //  fCentrality = (Double_t) dPercentile;
  //}
  // else if(fIsMC){
  //   AliMCEvent* mcEvent = dynamic_cast<AliMCEvent*>(MCEvent());
  //   if(!mcEvent) return kFALSE;
  //   Int_t ntrackv0aprimary=0;
  //
  //   for(Int_t i(0); i < mcEvent->GetNumberOfTracks(); i++) {
  //     AliMCParticle* part = (AliMCParticle*)mcEvent->GetTrack(i);
  //     if(!part->IsPhysicalPrimary()) continue;
  //     Double_t mceta = part->Eta();
  //     if(fBoostAMPT) mceta = TransverseBoost(part);
  //
  //     if(part->Charge()==0)        continue;
  //     if(mceta>2.8 && mceta<5.1) ntrackv0aprimary++;
  //   }
  //   Int_t nbinmult= fhCentCalib->GetXaxis()->FindBin(ntrackv0aprimary);
  //   fCentrality = (Double_t) fhCentCalib->GetBinContent(nbinmult);
  // }
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
Bool_t AliAnalysisTaskCorrForFlowMaster_Anders::IsTrackSelected(const AliAODTrack* track) const
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

  if(fRejectSecondariesFromMC){
    AliMCEvent* mcEvent = dynamic_cast<AliMCEvent*>(MCEvent());
    if(!mcEvent) return kFALSE;
    if(track->GetLabel() < 0) return kFALSE;
    AliMCParticle* part = (AliMCParticle*)mcEvent->GetTrack(track->GetLabel()); //Get
    if(!part) return kFALSE;
    if(!part->IsPhysicalPrimary()) { return kFALSE; }
  }


  return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlowMaster_Anders::RangePhi(Double_t DPhi) {
    if(fabs(DPhi)>2 * TMath::Pi()){
        cout<<"Dphi out of range"<<"\n";
        cout<<DPhi<<"\n";
    }
  if (DPhi < -TMath::Pi() / 2)   DPhi += 2 * TMath::Pi();
  if (DPhi > 3 * TMath::Pi() / 2) DPhi -= 2*TMath::Pi();
    
  return DPhi;
}
//_____________________________________________________________________________
/*Bool_t AliAnalysisTaskCorrForFlowMaster_Anders::CheckDPhiStar(Double_t dEta, Double_t trigPhi, Double_t trigPt, Double_t trigCharge, Double_t assPhi, Double_t assPt, Double_t assCharge){
  Bool_t check = kFALSE;
  if(TMath::Abs(dEta) < fMergingCut){
    Double_t dPhiStarLow = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, 0.8);
    Double_t dPhiStarHigh = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, 2.5);

    const Double_t kLimit = 3.0*fMergingCut;

    if(TMath::Abs(dPhiStarLow) < kLimit || TMath::Abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0 ) {
      for(Double_t rad(0.8); rad < 2.51; rad+=0.01){
        Double_t dPhiStar = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, rad);
        if(TMath::Abs(dPhiStar) < fMergingCut) {
          check = kTRUE;
          return check;
        }
      } // end loop radius
    }
  }
  return check;
}*/
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowMaster_Anders::CheckDPhiStar(Double_t dEta, Double_t trigPhi, Double_t trigPt, Double_t trigCharge, Double_t assPhi, Double_t assPt, Double_t assCharge){
  //Bool_t check = kFALSE;
  if(TMath::Abs(dEta) < fMergingCut){
    Double_t dPhiStarLow = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, 0.8);
    Double_t dPhiStarHigh = GetDPhiStar(trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge, 2.5);

    //const Double_t kLimit = 3.0*fMergingCut;
      //check=kTRUE means cut failed.
      
      //If either high or low limit fails the merging cut reject.
      if((dPhiStarLow>fMergingCut && dPhiStarHigh>fMergingCut)||(dPhiStarLow<-fMergingCut && dPhiStarHigh<-fMergingCut)){
          return kFALSE;
      }else{
          return kTRUE;
      }
      //If dPhi* changes sign between the inner and outer radius reject.
  }
  return kFALSE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlowMaster_Anders::GetDPhiStar(Double_t phi1, Double_t pt1, Double_t charge1, Double_t phi2, Double_t pt2, Double_t charge2, Double_t radius){
  // calculates delta phi *
  Double_t dPhiStar = phi1 - phi2 - charge1 * fbSign * TMath::ASin(0.075 * radius / pt1) + charge2 * fbSign * TMath::ASin(0.075 * radius / pt2);

  if (dPhiStar > TMath::Pi()) dPhiStar = 2.0*TMath::Pi() - dPhiStar;
  if (dPhiStar < -TMath::Pi()) dPhiStar = -2.0*TMath::Pi() - dPhiStar;

  return dPhiStar;
}

//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowMaster_Anders::FillCorrelations()
{
  if(!fTracksTrig || !fhTrigTracks || !fTracksAss) { AliError("Necessary inputs missing, terminating!"); return; }
  if(!fhSE) { AliError("Output AliTHn missing: terminating!"); return; }
  const Int_t sizePtTrig = fPtBinsTrigCharged.size() - 1;
  TH1D *h = new TH1D();
  h->SetBins(sizePtTrig, fPtBinsTrigCharged.data());
  
  Double_t binscont[6];
  binscont[2] = fPVz;
  binscont[3] = fSampleIndex;

  Double_t binscontref[4];
  binscontref[2] = fPVz;
  binscontref[3] = fSampleIndex;
    
    //Changing this part so auto/double correlations are avoided using the track index instead of Pt.

  for(Int_t iTrig(0); iTrig < fTracksTrig->GetEntriesFast(); iTrig++){
    AliVParticle* track = dynamic_cast<AliVParticle*>(fTracksTrig->At(iTrig));
    if(!track) continue;
    AliAODTrack* trackAOD = nullptr;
    if(!fIsMC) trackAOD = (AliAODTrack*)fTracksTrig->At(iTrig);

    Double_t trigPt = track->Pt();
    Double_t trigEta = track->Eta();
    Double_t trigPhi = track->Phi();
    Double_t trigCharge = track->Charge();
    Double_t trigEff = 1.0;
    if(fUseEfficiency) {
      trigEff = GetEff(trigPt, 0, trigEta);
      if(trigEff < 0.001) continue;
    }
    binscont[4] = trigPt;
      
      phiCorr_ME->Fill(trigPhi,trigEta);

    for(Int_t iAss(0); iAss < fTracksAss->GetEntriesFast(); iAss++){
      AliVParticle* trackAss = dynamic_cast<AliVParticle*>(fTracksAss->At(iAss));
      if(!trackAss) continue;
      AliAODTrack* trackAODAss = nullptr;
      trackAODAss = (AliAODTrack*)fTracksAss->At(iAss); //if(!fIsMC) 

      Double_t assPt = trackAss->Pt();
      Double_t assEta = trackAss->Eta();
      Double_t assPhi = trackAss->Phi();
      Double_t assCharge = trackAss->Charge();
      Double_t assEff = 1.0;

      //if(!fIsMC && trackAOD->GetID() == trackAODAss->GetID()) continue;

      if((Int_t)track->GetUniqueID() == (Int_t)trackAss->GetUniqueID()) {continue;}

      //Ref vs ref - fill seperate histogram for the situation were both tracks are from the large ref. range
    //Checking that the Pt's are in the apropriate ranges, but not requiring trigPt>assPt as this is auto/double correlations are accounted for by indexing.
      if(fPtRefMin<trigPt<fPtRefMax && fPtRefMin<assPt<fPtRefMax && trigPt>assPt){//
        if(fUseEfficiency) {
          assEff = GetEff(assPt, 0, assEta);
          if(assEff < 0.001) continue;
        }
        if(fUseLikeSign){
          if(trigCharge*assCharge < 0) continue;
        }
        if(fUseUnlikeSign){
          if(trigCharge*assCharge > 0) continue;
        }

        binscontref[0] = trigEta - assEta;
        binscontref[1] = RangePhi(trigPhi - assPhi);
        //binscont[5] = assPt;
          
          phiCorr_SE->Fill(trigPhi,assPhi);
        
         
        
        if(fUsePhiStar && CheckDPhiStar(binscontref[0], trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge)) continue;

        fhSEref->Fill(binscontref,0,1./(trigEff*assEff));
      }
      //All other situations
      Int_t TrigBin = h->FindBin(trigPt);
      Int_t AssBin = h->FindBin(assPt);

      if(fUseEfficiency) {
        assEff = GetEff(assPt, 0, assEta);
        if(assEff < 0.001) continue;
      }

      if(fUseLikeSign){
        if(trigCharge*assCharge < 0) continue;
      }
      if(fUseUnlikeSign){
        if(trigCharge*assCharge > 0) continue;
      }


      binscont[0] = trigEta - assEta;
      binscont[1] = RangePhi(trigPhi - assPhi);
      binscont[5] = assPt;

     
      

      if(fUsePhiStar && CheckDPhiStar(binscont[0], trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge)) continue;
      
      fhSE->Fill(binscont,0,1./(trigEff*assEff));
    }
  }

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowMaster_Anders::FillCorrelationsMixed()
{
  if(!fTracksTrig || !fhTrigTracks || !fTracksAss) { AliError("Necessary inputs missing, terminating!"); return; }

  AliEventPool *pool = fPoolMgr->GetEventPool(fCentrality, fPVz);
  if(!pool) { AliError(Form("No pool found for centrality = %f, zVtx = %f", fCentrality,fPVz)); return; }
  const Int_t sizePtTrig = fPtBinsTrigCharged.size() - 1;
  TH1D *h = new TH1D();
  h->SetBins(sizePtTrig, fPtBinsTrigCharged.data());
  if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() > fMinEventsToMix) {
    Int_t nMix = pool->GetCurrentNEvents();

    
    Double_t binscont[6];
    binscont[2] = fPVz;
    binscont[3] = fSampleIndex;

    Double_t binscontref[4];
    binscontref[2] = fPVz;
    binscontref[3] = fSampleIndex;

    for(Int_t iTrig(0); iTrig < fTracksTrig->GetEntriesFast(); iTrig++){
      AliVParticle* track = dynamic_cast<AliVParticle*>(fTracksTrig->At(iTrig));
      if(!track) continue;

      Double_t trigPt = track->Pt();
      Double_t trigEta = track->Eta();
      Double_t trigPhi = track->Phi();
      Double_t trigCharge = track->Charge();
      binscont[4] = trigPt;
      Double_t trigEff = 1.0;
      if(fUseEfficiency) {
        trigEff = GetEff(trigPt, 0, trigEta);
        if(trigEff < 0.001) continue;
      }
        
    //Using indexing to avoid autocorrelations, as in same event loop. Also removes, the potentially problematic trigPt>assPt requirement.
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

          //Ref vs ref - fill seperate histogram for the situation were both tracks are from the large ref. range
          if(fPtRefMin<trigPt<fPtRefMax && fPtRefMin<assPt<fPtRefMax && trigPt>assPt){//
            if(fUseEfficiency) {
              assEff = GetEff(assPt, 0, assEta);
              if(assEff < 0.001) continue;
            }
            if(fUseLikeSign){
              if(trigCharge*assCharge < 0) continue;
            }
            if(fUseUnlikeSign){
              if(trigCharge*assCharge > 0) continue;
            }
            binscontref[0] = trigEta - assEta;
            binscontref[1] = RangePhi(trigPhi - assPhi);
            binscont[5] = assPt;
            
              if(fCreateQAPlots){
                  phDphiPt1Pt2_ME->Fill(trigPt,assPt,trigPhi - assPhi);
                  if(trigPt>assPt){
                      phDphiPt1Pt2_PtCut_ME->Fill(trigPt,assPt,trigPhi - assPhi);
                  }
              }
            

            if(fUsePhiStar && CheckDPhiStar(binscontref[0], trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge)) continue;
        
            fhMEref->Fill(binscontref,0,1./((Double_t)nMix*(trigEff*assEff)));
            
          }
          //All other situations
          Int_t TrigBin = h->FindBin(trigPt);
          Int_t AssBin = h->FindBin(assPt);

          
          if(fUseEfficiency) {
            assEff = GetEff(assPt, 0, assEta);
            if(assEff < 0.001) continue;
          }
          if(fUseLikeSign){
            if(trigCharge*assCharge < 0) continue;
          }
          if(fUseUnlikeSign){
            if(trigCharge*assCharge > 0) continue;
          }

          binscont[0] = trigEta - assEta;
          binscont[1] = RangePhi(trigPhi - assPhi);
          binscont[5] = assPt;

          if(fUsePhiStar && CheckDPhiStar(binscont[0], trigPhi, trigPt, trigCharge, assPhi, assPt, assCharge)) continue;

          fhME->Fill(binscont,0,1./((Double_t)nMix*(trigEff*assEff)));
        } //End of Asso loop
      } //End of mixed event loop
    }
  } // event pool done
  TObjArray* cloneArray = (TObjArray *)fTracksAss->Clone();
  cloneArray->SetOwner(kTRUE);
  pool->UpdatePool(cloneArray);
  
  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowMaster_Anders::AreEfficienciesLoaded()
{
  if(!fInputListEfficiency) {AliError("Efficiency input list not loaded"); return kFALSE; }
  TString part[1] = {"ch"};
  if(fColSystem == sPPb){
    TString etaReg[8] = {"0020", "0200", "0204", "0402", "0406", "0604", "0608", "0806"};
    for(Int_t p(0); p < 1; p++){
      for(Int_t eta(0); eta < 8; eta++){
        fhEfficiencyEta[eta] = (TH2D*)fInputListEfficiency->FindObject(Form("LHC17f2b_%s_Eta_%s_%s_wFD",part[p].Data(), etaReg[eta].Data(),fSystematicsFlag.Data()));
        if(!fhEfficiencyEta[eta]) {AliError(Form("Efficiency (%s, eta region %s, flag %s) not loaded",part[p].Data(),etaReg[eta].Data(),fSystematicsFlag.Data())); return kFALSE; }
      }
    }
    fhEventCounter->Fill("Efficiencies loaded",1);
    return kTRUE;
  }
  else if(fColSystem == sPP){
    for(Int_t p(0); p < 1; p++){
      fhEfficiency[p] = (TH2D*)fInputListEfficiency->FindObject(Form("LHC%s_%s_%s_wFD",ReturnPPperiod(fAOD->GetRunNumber()).Data(),part[p].Data(),fSystematicsFlag.Data()));
      if(!fhEfficiency[p]) {AliError(Form("Efficiency (run %d, part %s, flag %s) not loaded",fAOD->GetRunNumber(),part[p].Data(),fSystematicsFlag.Data())); return kFALSE; }
    }
    fhEventCounter->Fill("Efficiencies loaded",1);
    return kTRUE;
  }

  return kFALSE;
}
//_____________________________________________________________________________
TString AliAnalysisTaskCorrForFlowMaster_Anders::ReturnPPperiod(const Int_t runNumber) const
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
Double_t AliAnalysisTaskCorrForFlowMaster_Anders::GetEff(const Double_t dPt, const Int_t spec, const Double_t dEta)
{
  if(!fUseEfficiency) return 1.0;
  if(fColSystem == sPPb){
    Int_t region = GetEtaRegion(dEta);
    if(region < 0) { AliWarning("Invalid region, returning efficiency 1.0."); return 1.0; }
    if(!fhEfficiencyEta[region]) { AliError("Efficiency histogram not found, returning efficiency 1.0."); return 1.0; }
    return fhEfficiencyEta[region]->GetBinContent(fhEfficiencyEta[region]->FindFixBin(dPt, fCentrality));
  }else{
    if(!fhEfficiency[0]) { AliError("Efficiency histogram not found, returning efficiency 1.0."); return 1.0; }
    return fhEfficiency[0]->GetBinContent(fhEfficiency[0]->FindFixBin(dPt, fCentrality));
  }

  return 1.0;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskCorrForFlowMaster_Anders::GetEtaRegion(const Double_t dEta){
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
void AliAnalysisTaskCorrForFlowMaster_Anders::CreateTHnCorrelations(){
  Int_t nSteps = 1;
  const Int_t sizePtTrig = fPtBinsTrigCharged.size() - 1;
  const Int_t sizePtAss = fPtBinsAssCharged.size() - 1;
  const Int_t sizeOfSamples = (Int_t) fNOfSamples;
  Double_t binning_dphi[73] = { -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0,       0.087266,  0.174533,  0.261799,  0.349066,  0.436332, 0.523599,  0.610865,  0.698132,  0.785398,  0.872665,  0.959931, 1.047198,  1.134464,  1.221730,  1.308997,  1.396263,  1.483530, 1.570796,  1.658063,  1.745329,  1.832596,  1.919862,  2.007129, 2.094395,  2.181662,  2.268928,  2.356194,  2.443461,  2.530727, 2.617994,  2.705260,  2.792527,  2.879793,  2.967060,  3.054326, 3.141593,  3.228859,  3.316126,  3.403392,  3.490659,  3.577925, 3.665191,  3.752458,  3.839724,  3.926991,  4.014257,  4.101524, 4.188790,  4.276057,  4.363323,  4.450590,  4.537856,  4.625123, 4.712389};

  TString nameS[1] = {"fhChargedSE"};
  TString nameM[1] = {"fhChargedME"};
  TString nameSref[1] = {"fhChargedSEref"};
  TString nameMref[1] = {"fhChargedMEref"};

  
  
    Double_t binning_deta_tpctpc[37] = {-1.8,-1.7,-1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,    0.1,  0.2,  0.3,  0.4,  0.5, 0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5, 1.6,1.7,1.8};
    Int_t iBinningTPCTPC[] = {32,72,10,sizeOfSamples,sizePtTrig, sizePtAss};
    Int_t nTrackBin_tpctpc = sizeof(iBinningTPCTPC) / sizeof(Int_t);

    Int_t iBinningTPCTPCref[] = {32,72,10,sizeOfSamples};
    Int_t nTrackBin_tpctpcref = sizeof(iBinningTPCTPCref) / sizeof(Int_t);

   

    fhSE = new AliTHn(nameS[0], nameS[0], nSteps, nTrackBin_tpctpc, iBinningTPCTPC);
    fhSE->SetBinLimits(0, binning_deta_tpctpc);
    fhSE->SetBinLimits(1, binning_dphi);

    fhME = new AliTHn(nameM[0], nameM[0], nSteps, nTrackBin_tpctpc, iBinningTPCTPC);
    fhME->SetBinLimits(0, binning_deta_tpctpc);
    fhME->SetBinLimits(1, binning_dphi);
    fhSE->SetBinLimits(2, -10,10);
    fhSE->SetBinLimits(3, 0,10);
    fhSE->SetVarTitle(0, "#Delta#eta");
    fhSE->SetVarTitle(1, "#Delta#phi");
    fhSE->SetVarTitle(2, "PVz [cm]");
    fhSE->SetVarTitle(3, "Sample");
    fhSE->SetBinLimits(4, fPtBinsTrigCharged.data());
    fhSE->SetBinLimits(5, fPtBinsAssCharged.data());
    fhSE->SetVarTitle(4, "p_{T} [GeV/c] (trig)");
    fhSE->SetVarTitle(5, "p_{T} [GeV/c] (ass)");
    
    fOutputListCharged->Add(fhSE);

    fhME->SetBinLimits(2, -10,10);
    fhME->SetBinLimits(3, 0,10);
    fhME->SetVarTitle(0, "#Delta#eta");
    fhME->SetVarTitle(1, "#Delta#phi");
    fhME->SetVarTitle(2, "PVz [cm]");
    fhME->SetVarTitle(3, "Sample");
    fhME->SetBinLimits(4, fPtBinsTrigCharged.data());
    fhME->SetBinLimits(5, fPtBinsAssCharged.data());
    fhME->SetVarTitle(4, "p_{T} [GeV/c] (trig)");
    fhME->SetVarTitle(5, "p_{T} [GeV/c] (ass)");
    
    fOutputListCharged->Add(fhME);

    //Ref vs ref

    fhSEref = new AliTHn(nameSref[0], nameSref[0], nSteps, nTrackBin_tpctpcref, iBinningTPCTPCref);
    fhSEref->SetBinLimits(0, binning_deta_tpctpc);
    fhSEref->SetBinLimits(1, binning_dphi);

    fhMEref = new AliTHn(nameMref[0], nameMref[0], nSteps, nTrackBin_tpctpcref, iBinningTPCTPCref);
    fhMEref->SetBinLimits(0, binning_deta_tpctpc);
    fhMEref->SetBinLimits(1, binning_dphi);
    fhSEref->SetBinLimits(2, -10,10);
    fhSEref->SetBinLimits(3, 0,10);
    fhSEref->SetVarTitle(0, "#Delta#eta");
    fhSEref->SetVarTitle(1, "#Delta#phi");
    fhSEref->SetVarTitle(2, "PVz [cm]");
    fhSEref->SetVarTitle(3, "Sample");
    // fhSEref->SetBinLimits(4, fPtBinsTrigCharged.data());
    // fhSEref->SetBinLimits(5, fPtBinsAssCharged.data());
    // fhSEref->SetVarTitle(4, "p_{T} [GeV/c] (trig)");
    // fhSEref->SetVarTitle(5, "p_{T} [GeV/c] (ass)");
    
    fOutputListCharged->Add(fhSEref);

    fhMEref->SetBinLimits(2, -10,10);
    fhMEref->SetBinLimits(3, 0,10);
    fhMEref->SetVarTitle(0, "#Delta#eta");
    fhMEref->SetVarTitle(1, "#Delta#phi");
    fhMEref->SetVarTitle(2, "PVz [cm]");
    fhMEref->SetVarTitle(3, "Sample");
    // fhMEref->SetBinLimits(4, fPtBinsTrigCharged.data());
    // fhMEref->SetBinLimits(5, fPtBinsAssCharged.data());
    // fhMEref->SetVarTitle(4, "p_{T} [GeV/c] (trig)");
    // fhMEref->SetVarTitle(5, "p_{T} [GeV/c] (ass)");
    
    fOutputListCharged->Add(fhMEref);
  

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowMaster_Anders::PrepareTPCTracks(){
  if(!fAOD) return kFALSE;
  if(!fTracksAss || !fTracksTrig || !fhTrigTracks) {AliError("Cannot prepare TPC tracks!"); return kFALSE; }

  fNofTracks = 0;
  Double_t binscont[5] = {fPVz, fSampleIndex, 0.,0.,0.};

  TObjArray* fTracksJets = nullptr;
  if(fVetoJetEvents) fTracksJets = new TObjArray;
  vector<AliAODTrack*> tempVec(0);
  vector<vector<AliAODTrack*>> vecTrack(fPtBinsTrigCharged.size(), tempVec);
  for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {
    
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
    if(!track || !IsTrackSelected(track)) { continue; }

    Double_t trackPt = track->Pt();



    if(fVetoJetEvents && trackPt > fJetParticleLowPt) fTracksJets->Add((AliAODTrack*)track);

    for(Int_t j(0); j<fPtBinsTrigCharged.size()-1; j++) {
      if(fPtBinsTrigCharged[j+1]>trackPt){
          vecTrack[j].emplace_back(track);
          break;
      }
    }
      
  } // tracks loop end

  for (Int_t i(0); i<vecTrack.size(); i++){
    if(fUseEventBias && vecTrack[i].size()<fNumEventBias){continue;}
    for (Int_t j(0); j<vecTrack[i].size(); j++){
      fNofTrackGlobal++;
      AliAODTrack* track = vecTrack[i][j];
      if(!track || !IsTrackSelected(track)) { continue; }
      Double_t trackPt = track->Pt();
      binscont[2] = trackPt;
        binscont[3] = track->Eta();
        binscont[4] = track->Phi();
      if(trackPt<fPtMaxTrig && trackPt>fPtMinTrig){
        track->SetUniqueID((Int_t)fNofEventGlobal*10000+(Int_t)fNofTrackGlobal);
        fhTrigTracks->Fill(binscont,0,1.);
        fTracksTrig->Add(track);
        fNofTracks++;



        if(fCreateQAPlots){
          fhPT->Fill(trackPt);
          fhPhi->Fill(track->Phi());
          fhEta->Fill(track->Eta());
        }
        
      }
      
      if(trackPt<fPtMaxAss && trackPt>fPtMinAss){fTracksAss->Add(track);}

    }
  }

  if(fUseNch){
    if(fNofTracks < fNchMin || fNofTracks > fNchMax) { return kFALSE; }
    fhEventCounter->Fill("Nch cut ok ",1);
    fhEventMultiplicity->Fill(fNofTracks);
  }

  if(fVetoJetEvents){
    Bool_t foundjetsinTPC = kFALSE;
    fhEventCounter->Fill("Before Jet Veto",1); //HPC = high pt cut
    
    for(Int_t iTrig=0; iTrig < fTracksJets->GetEntriesFast(); iTrig++){
      AliAODTrack* trackTrig = (AliAODTrack*)fTracksJets->At(iTrig);
      if(!trackTrig) continue;
      Double_t trigPhi = trackTrig->Phi();

      for(Int_t iAss=iTrig+1; iAss < fTracksJets->GetEntriesFast(); iAss++){
        AliAODTrack* trackAss = (AliAODTrack*)fTracksJets->At(iAss);
        if(!trackAss) continue;
        Double_t assPhi = trackAss->Phi();

        Double_t deltaPhi = RangePhi(trigPhi - assPhi);
        if(TMath::Abs(deltaPhi - TMath::Pi()) < fJetvetoselectionval) foundjetsinTPC = kTRUE;//fJetvetoselectionval= 0.5
      }
    }
    
      if(fselectjetsinTPC == kTRUE) {//select events with back to back jets in TPC
      if(foundjetsinTPC == kFALSE) return kFALSE;//reject
    }


    if(fselectjetsinTPC == kFALSE) {//reject events with back to back jets in TPC
      if(foundjetsinTPC == kTRUE) return kFALSE;//reject
    }


    fhEventCounter->Fill("After Jet Veto",1);
  }


  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowMaster_Anders::PrepareMCTracks(){
  if(!fTracksAss || !fTracksTrig || !fhTrigTracks) {AliError("Cannot prepare MCC tracks!"); return kFALSE; }

  AliMCEvent* mcEvent = dynamic_cast<AliMCEvent*>(MCEvent());
  if(!mcEvent) return kFALSE;

  if(fOnTheFly){
    AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(mcEvent->GenEventHeader());
    if(!header){return kFALSE;}
    AliCollisionGeometry *headerH;
    TString genName;
    TList *ltgen = (TList*)mcEvent->GetCocktailList();
    if(ltgen){
      for(auto&& listObject: *ltgen){
        genName = Form("%s",listObject->GetName());
        if(genName.Contains("Hijing")){
          headerH = dynamic_cast<AliCollisionGeometry*>(listObject);
          break;
        }
      }
    }else{
      headerH = dynamic_cast<AliCollisionGeometry*>(mcEvent->GenEventHeader());
    }
    if(headerH){
      fImpactParameterMC = headerH->ImpactParameter();
    }
    fCentrality = getGeneratorCentrality();
    if(fCentrality < fCentMin || fCentrality > fCentMax) { return kFALSE; }
    fhEventCounter->Fill("CentOK",1);

  }

  vector<AliMCParticle*> tempVec(0);
  vector<vector<AliMCParticle*>> vecTrack(fPtBinsTrigCharged.size(), tempVec);

  Double_t binscont[3] = {fPVz, fSampleIndex, 0.};

  for(Int_t i(0); i < mcEvent->GetNumberOfTracks(); i++) {
    
    AliMCParticle* part = (AliMCParticle*)mcEvent->GetTrack(i);
    if(!fIsTPCgen) continue;
    if(part->Charge()==0.) continue;
    if(!part->IsPhysicalPrimary()) continue;
    if(fAbsEtaMax > 0.0 && TMath::Abs(part->Eta()) > fAbsEtaMax) {continue;}
    Double_t trackPt = part->Pt();

    for(Int_t j(0); j<fPtBinsTrigCharged.size()-1; j++) {
      if(fPtBinsTrigCharged[j+1]>trackPt){
          vecTrack[j].emplace_back((AliMCParticle*)part);
          break;
      }
    }
  }


  for (Int_t i(0); i<vecTrack.size(); i++){
    if(fUseEventBias && vecTrack[i].size()<fNumEventBias){continue;}
    for (Int_t j(0); j<vecTrack[i].size(); j++){
      fNofTrackGlobal++;
      AliMCParticle* track = vecTrack[i][j];
      Double_t trackPt = track->Pt();
      binscont[2] = trackPt;
      track->SetUniqueID((Int_t)fNofEventGlobal*10000+(Int_t)fNofTrackGlobal);
      if(trackPt<fPtMaxTrig && trackPt>fPtMinTrig){

        fhTrigTracks->Fill(binscont,0,1.);
        
        fTracksTrig->Add((AliMCParticle*)track);
      }
      
      if(trackPt<fPtMaxAss && trackPt>fPtMinAss){
        
        fTracksAss->Add((AliMCParticle*)track);
      }

    }
  }
  


  return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlowMaster_Anders::TransverseBoost(const AliMCParticle *track){
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
void AliAnalysisTaskCorrForFlowMaster_Anders::PrintSetup(){
  printf("\n\n\n ************** Parameters ************** \n");
  printf("\t fAnalType: (Int_t) %d\n", fAnalType);
  printf("\t fColSystem: (Int_t) %d\n", fColSystem);
  printf("\t fIsMC: (Bool_t) %s\n", fIsMC ? "kTRUE" : "kFALSE");
  printf("\t fIsTPCgen: (Bool_t) %s\n", fIsTPCgen ? "kTRUE" : "kFALSE");
  printf("\t fUseNch: (Bool_t) %s\n", fUseNch ? "kTRUE" : "kFALSE");
  printf("\t fIsHMpp: (Bool_t) %s\n", fIsHMpp ? "kTRUE" : "kFALSE");
  printf("\t fUseEfficiency: (Bool_t) %s\n",  fUseEfficiency ? "kTRUE" : "kFALSE");
  printf("\t fUseOppositeSidesOnly: (Bool_t) %s\n", fUseOppositeSidesOnly ? "kTRUE" : "kFALSE");
  printf("\t fRejectSecondariesFromMC: (Bool_t) %s\n", fRejectSecondariesFromMC ? "kTRUE" : "kFALSE");
  printf("\t fNOfSamples: (Int_t) %d\n", (Int_t) fNOfSamples);
  printf("\t fEventBias: (Bool_t) %s\n", fUseEventBias ? "kTRUE" : "kFALSE");
  printf("\t fNumEventBias: (Int_t) %d\n", (Int_t) fNumEventBias);  
  printf("\t fUsePhiStar: (Bool_t) %s\n", fUsePhiStar ? "kTRUE" : "kFALSE");
  printf(" **************************** \n");
  printf("\t fSystematicsFlag: (TString) %s\n", fSystematicsFlag.Data());
  printf("\t fAbsEtaMax: (Double_t) %f\n", fAbsEtaMax);
  printf("\t fPtMinTrig -- fPtMaxTrig: (Double_t) %f -- %f\n", fPtMinTrig, fPtMaxTrig);
  printf("\t fPtMinAss -- fPtMaxAss: (Double_t) %f -- %f\n", fPtMinAss, fPtMaxAss);
  printf("\t fCentMin -- fCentMax: (Double_t) %f -- %f\n", fCentMin, fCentMax);
  printf("\t fPVzCut: (Double_t) %f\n", fPVzCut);
  printf(" **************************** \n");
}
//___________________________________________________________________________
double AliAnalysisTaskCorrForFlowMaster_Anders::getGeneratorCentrality()
{
  vector<double> b;
  if(centralitymap.empty()) AliFatal("Centralitymap is empty!");
  for (auto const& element : centralitymap) b.push_back(element.first);
  vector<double>::iterator it = upper_bound(b.begin(),b.end(),fImpactParameterMC);
  double l_cent = (fImpactParameterMC<0)?-1.0:(centralitymap[b[it-b.begin()]]+centralitymap[b[it-b.begin()-1]])/2.0;
  return l_cent;
}
