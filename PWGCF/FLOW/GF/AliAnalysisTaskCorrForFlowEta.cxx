/**************************************************************************
 *    Author:       Mikkel Petersen                                       *
 *    Framework for calculating di-hadron correlation                     *
 *    for extraction of v_n{2} coefficients of unidentified particles     *
 *    in TPC-FMD correlations while studying the eta decorrelation.       *
 *                                                                        *
 *    If used, modified, or distributed,                                  *
 *    please aknowledge the author of this code.                          *
 *    Large parts of the code was taken from Zuzana Moravcova, as TPCFMD  *
 *    correlations were aldready implemented by her                       *
 **************************************************************************/

#include "AliAnalysisTaskCorrForFlowEta.h"

using namespace std;

ClassImp(AliAnalysisTaskCorrForFlowEta);

AliAnalysisTaskCorrForFlowEta::AliAnalysisTaskCorrForFlowEta() : AliAnalysisTaskSE(),
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
    fIsFMDgen(kFALSE),
    fIsHMpp(kFALSE),
    fUseNch(kFALSE),
    fUseEfficiency(kFALSE),
    fUseFMDcut(kTRUE),
    fUseOppositeSidesOnly(kFALSE),
    fUseCentralityCalibration(kFALSE),
    fSkipCorr(kFALSE),
    fVetoJetEvents(kFALSE),
    fRejectSecondariesFromMC(kFALSE),
    fBoostAMPT(kFALSE),
    fFilterBit(96),
    fbSign(0),
    fRunNumber(-1),
    fNofTracks(0),
    fNofMinHighPtTracksForRejection(0),
    fNchMin(0),
    fNchMax(100000),
    fnTPCcrossedRows(70),
    fNOfSamples(1.0),
    fSampleIndex(0.0),
    fPtMinTrig(0.2),
    fPtMaxTrig(10.0),
    fFMDcutapar0(1.64755),
    fFMDcutapar1(119.602),
    fFMDcutcpar0(2.73426),
    fFMDcutcpar1(150.31),
    fFMDAacceptanceCutLower(1.8),
    fFMDAacceptanceCutUpper(4.8), //4.8
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
    fSigmaTPC(3.),
    fJetParticleLowPt(5.),
    fCentEstimator("V0M"),
    fSystematicsFlag(""),
    fPoolMaxNEvents(2000),
    fPoolMinNTracks(50000),
    fMinEventsToMix(5),
    fNzVtxBins(10),
    fNCentBins(15),
    fNBinsTPCeta(0),
    fNBinsFMDeta(0),
    fNBinsdEta(25),
    fMergingCut(0.0)
{}
//_____________________________________________________________________________
AliAnalysisTaskCorrForFlowEta::AliAnalysisTaskCorrForFlowEta(const char* name, Bool_t bUseEff, Bool_t bUseCalib) : AliAnalysisTaskSE(name),
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
    fIsFMDgen(kFALSE),
    fIsHMpp(kFALSE),
    fUseNch(kFALSE),
    fUseEfficiency(kFALSE),
    fUseFMDcut(kTRUE),
    fUseOppositeSidesOnly(kFALSE),
    fUseCentralityCalibration(kFALSE),
    fSkipCorr(kFALSE),
    fVetoJetEvents(kFALSE),
    fRejectSecondariesFromMC(kFALSE),
    fBoostAMPT(kFALSE),
    fFilterBit(96),
    fbSign(0),
    fRunNumber(-1),
    fNofTracks(0),
    fNofMinHighPtTracksForRejection(0),
    fNchMin(0),
    fNchMax(100000),
    fnTPCcrossedRows(70),
    fNOfSamples(1.0),
    fSampleIndex(0.0),
    fPtMinTrig(0.2),
    fPtMaxTrig(10.0),
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
    fSigmaTPC(3.),
    fJetParticleLowPt(5.),
    fCentEstimator("V0M"),
    fSystematicsFlag(""),
    fPoolMaxNEvents(2000),
    fPoolMinNTracks(50000),
    fMinEventsToMix(5),
    fNzVtxBins(10),
    fNCentBins(15),
    fNBinsTPCeta(0),
    fNBinsFMDeta(0),
    fNBinsdEta(25),
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
AliAnalysisTaskCorrForFlowEta::~AliAnalysisTaskCorrForFlowEta()
{}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowEta::UserCreateOutputObjects()
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

    fHistFMDeta = new TH2D("fHistFMDeta", "FMD eta vs. PVz; eta; PVz [cm]", 90, -4, 5, 20, -10, 10);
    fOutputListCharged->Add(fHistFMDeta);

    TString fmdv0corrNames[4] = {"A_Before","C_Before", "A_After", "C_After"};
    for(Int_t i(0); i < 4; i++){
      fh2FMDvsV0[i] = new TH2D(Form("fh2FMDvsV0%s",fmdv0corrNames[i].Data()), "FMD vs. V0; FMD; V0", 250, 0, 1000, 250, 0, 1000);
      fOutputListCharged->Add(fh2FMDvsV0[i]);
    }

    TString pidName[1] = {""};
    const Int_t sizeOfSamples = (Int_t) fNOfSamples;
    Int_t binsFMD[] = {10, 10};
    Int_t binsTPCFMD[] = {10, 10, 5};

    
    if(fAnalType == eFMDAFMDC) fhTrigTracks = new AliTHn(Form("fhTrigTracks%s",pidName[0].Data()), Form("fhTrigTracks (%s)",pidName[0].Data()), 1, 2, binsFMD);
    else fhTrigTracks = new AliTHn(Form("fhTrigTracks%s",pidName[0].Data()), Form("fhTrigTracks (%s)",pidName[0].Data()), 1, 3, binsTPCFMD);
    fhTrigTracks->SetBinLimits(0,-10,10);
    fhTrigTracks->SetBinLimits(1,0,10);
    fhTrigTracks->SetVarTitle(0, "PVz [cm]");
    fhTrigTracks->SetVarTitle(1, "Sample");
    if(fAnalType != eFMDAFMDC){
        fhTrigTracks->SetBinLimits(2,0, 5);
        fhTrigTracks->SetVarTitle(2, "p_{T} (trig)");
    }
    fOutputListCharged->Add(fhTrigTracks);
    

    //mixing
    fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins.data(), fNzVtxBins, fzVtxBins.data());
    if (!fPoolMgr) { AliError("Event Pool manager not created!"); return; }
    fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);

    if(!fSkipCorr) CreateTHnCorrelations();

    
    fhPT = new TH1D(Form("PT"), "PT", 1000, 0, 10);
    fhPT->Sumw2();
    fOutputListCharged->Add(fhPT);
    
    //if(!fSkipCorr) break;
    

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

    if(fAnalType == eFMDAFMDC && fUseEfficiency){ AliWarning("Efficeincies inserted when running FMDA-FMDC. Turning off the flag."); fUseEfficiency = kFALSE; }

    PostData(1, fOutputListCharged);
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowEta::UserExec(Option_t *)
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
  fTracksTrig = new TObjArray;

  

  if(fUseEfficiency && fColSystem == sPP && (fRunNumber != fAOD->GetRunNumber()) && !AreEfficienciesLoaded()) { return; }
  
  // FMD - V0 correlation event cut
  if(!PrepareFMDTracks()){
    delete fTracksAss;
    delete fTracksTrig;
    PostData(1, fOutputListCharged);
    return;
    
  }
  

  if(!fIsTPCgen || fUseNch)  {
    if(!PrepareTPCTracks()){
      if(fTracksTrig) delete fTracksTrig;
      PostData(1, fOutputListCharged);
      return;
    }
  }

  if(fIsMC){
    if(!PrepareMCTracks()){
      if(fTracksTrig) delete fTracksTrig;
      PostData(1, fOutputListCharged);
      return;
    }
  } // end MC

  if(!fTracksAss->IsEmpty() && !fSkipCorr){
      FillCorrelations();
      FillCorrelationsMixed();
      fTracksTrig->Clear();
      delete fTracksTrig;
  }

  if(fUseEfficiency) fRunNumber = fAOD->GetRunNumber();

  fTracksAss->Clear();
  delete fTracksAss;

  PostData(1, fOutputListCharged);
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowEta::Terminate(Option_t *)
{
   //if(fPoolMgr) delete fPoolMgr;
   //if(fOutputListCharged) delete fOutputListCharged;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowEta::IsEventSelected()
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
Bool_t AliAnalysisTaskCorrForFlowEta::IsTrackSelected(const AliAODTrack* track) const
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
    AliMCParticle* part = (AliMCParticle*)mcEvent->GetTrack(track->GetLabel());
    if(!part) return kFALSE;
    if(!part->IsPhysicalPrimary()) { return kFALSE; }
  }

  return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlowEta::RangePhiFMD(Double_t DPhi) {
  DPhi = TMath::ATan2(TMath::Sin(DPhi), TMath::Cos(DPhi));
  if (DPhi < (-0.5*TMath::Pi()-0.0001))    DPhi += 2 * TMath::Pi();
  return DPhi;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlowEta::RangePhi(Double_t DPhi) {
  if (DPhi < -TMath::Pi() / 2)   DPhi += 2 * TMath::Pi();
  if (DPhi > 3 * TMath::Pi() / 2) DPhi -= 2*TMath::Pi();
  return DPhi;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowEta::FillCorrelations()
{
  if(!fTracksTrig || !fhTrigTracks || !fTracksAss) { AliError("Necessary inputs missing, terminating!"); return; }
  if(!fhSE) { AliError(Form("Output AliTHn missing , terminating!")); return; }
 
  Double_t binscont[6];
  binscont[4] = fPVz;
  binscont[5] = fSampleIndex;

  for(Int_t iTrig(0); iTrig < fTracksTrig->GetEntriesFast(); iTrig++){
    AliVParticle* track = dynamic_cast<AliVParticle*>(fTracksTrig->At(iTrig));
    if(!track) continue;

    Double_t trigPt = track->Pt();
    Double_t trigEta = track->Eta();
    Double_t trigPhi = track->Phi();
    Double_t trigEff = 1.0;
    if(fUseEfficiency) {
      trigEff = GetEff(trigPt, trigEta);
      if(trigEff < 0.001) continue;
    }


    for(Int_t iAss(0); iAss < fTracksAss->GetEntriesFast(); iAss++){
      AliPartSimpleForCorr* trackAss = (AliPartSimpleForCorr*)fTracksAss->At(iAss);
      if(!trackAss) continue;

      Double_t assEta = trackAss->Eta();
      Double_t assPhi = trackAss->Phi();
      Double_t assMult = trackAss->Multiplicity();

      binscont[0] = trigEta;
      binscont[1] = TMath::Abs(assEta);
      binscont[2] = RangePhi(trigPhi - assPhi);
      binscont[3] = TMath::Abs(trigEta-assEta);
      fhSE->Fill(binscont,0,assMult/(trigEff));
    }
  }
 

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowEta::FillCorrelationsMixed()
{
  if(!fTracksTrig || !fhTrigTracks || !fTracksAss) { AliError("Necessary inputs missing, terminating!"); return; }

  AliEventPool *pool = fPoolMgr->GetEventPool(fCentrality, fPVz);
  if(!pool) { AliError(Form("No pool found for centrality = %f, zVtx = %f", fCentrality,fPVz)); return; }

  if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() > fMinEventsToMix) {
    Int_t nMix = pool->GetCurrentNEvents();

    
    Double_t binscont[6];
    binscont[4] = fPVz;
    binscont[5] = fSampleIndex;

    for(Int_t iTrig(0); iTrig < fTracksTrig->GetEntriesFast(); iTrig++){
      AliVParticle* track = dynamic_cast<AliVParticle*>(fTracksTrig->At(iTrig));
      if(!track) continue;

      Double_t trigPt = track->Pt();
      Double_t trigEta = track->Eta();
      Double_t trigPhi = track->Phi();
      Double_t trigEff = 1.0;
      if(fUseEfficiency) {
        trigEff = GetEff(trigPt, trigEta);
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

          binscont[0] = trigEta;
          binscont[1] = TMath::Abs(assEta);
          binscont[2] = RangePhi(trigPhi - assPhi);
          binscont[3] = TMath::Abs(trigEta-assEta);

          fhME->Fill(binscont,0,assMult/((Double_t)nMix*trigEff));
        }
      }
    }
  } // event pool done

  
    TObjArray* cloneArray = (TObjArray *)fTracksAss->Clone();
    cloneArray->SetOwner(kTRUE);
    pool->UpdatePool(cloneArray);
  

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowEta::AreEfficienciesLoaded()
{
  if(!fInputListEfficiency) {AliError("Efficiency input list not loaded"); return kFALSE; }
  TString part[1] = {"ch"};
  if(fColSystem == sPPb){
    TString etaReg[8] = {"0020", "0200", "0204", "0402", "0406", "0604", "0608", "0806"};
    
    for(Int_t eta(0); eta < 8; eta++){
    
        fhEfficiencyEta[0][eta] = (TH2D*)fInputListEfficiency->FindObject(Form("LHC17f2b_%s_Eta_%s_%s_wFD",part[1].Data(), etaReg[eta].Data(),fSystematicsFlag.Data()));
        if(!fhEfficiencyEta[0][eta]) {AliError(Form("Efficiency (%s, eta region %s, flag %s) not loaded",part[0].Data(),etaReg[eta].Data(),fSystematicsFlag.Data())); return kFALSE; }
    }
    fhEventCounter->Fill("Efficiencies loaded",1);
    return kTRUE;
  }
  else if(fColSystem == sPP){
    fhEfficiency = (TH2D*)fInputListEfficiency->FindObject(Form("LHC%s_%s_%s_wFD",ReturnPPperiod(fAOD->GetRunNumber()).Data(),part[0].Data(),fSystematicsFlag.Data()));
    if(!fhEfficiency) {AliError(Form("Efficiency (run %d, part %s, flag %s) not loaded",fAOD->GetRunNumber(),part[0].Data(),fSystematicsFlag.Data())); return kFALSE; }
    
    fhEventCounter->Fill("Efficiencies loaded",1);
    return kTRUE;
  }

  return kFALSE;
}
//_____________________________________________________________________________
TString AliAnalysisTaskCorrForFlowEta::ReturnPPperiod(const Int_t runNumber) const
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
Double_t AliAnalysisTaskCorrForFlowEta::GetEff(const Double_t dPt, const Double_t dEta)
{
  if(!fUseEfficiency) return 1.0;
  if(fColSystem == sPPb){
    Int_t region = GetEtaRegion(dEta);
    if(region < 0) { AliWarning("Invalid region, returning efficiency 1.0."); return 1.0; }
    if(!fhEfficiencyEta[0][region]) { AliError("Efficiency histogram not found, returning efficiency 1.0."); return 1.0; }
    return fhEfficiencyEta[0][region]->GetBinContent(fhEfficiencyEta[0][region]->FindFixBin(dPt, fCentrality));
  }else{
    if(!fhEfficiency) { AliError("Efficiency histogram not found, returning efficiency 1.0."); return 1.0; }
    return fhEfficiency->GetBinContent(fhEfficiency->FindFixBin(dPt, fCentrality));
  }

  return 1.0;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskCorrForFlowEta::GetEtaRegion(const Double_t dEta){
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
void AliAnalysisTaskCorrForFlowEta::CreateTHnCorrelations(){
  Int_t nSteps = 1;
  Double_t binning_dphi[73] = { -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0,       0.087266,  0.174533,  0.261799,  0.349066,  0.436332, 0.523599,  0.610865,  0.698132,  0.785398,  0.872665,  0.959931, 1.047198,  1.134464,  1.221730,  1.308997,  1.396263,  1.483530, 1.570796,  1.658063,  1.745329,  1.832596,  1.919862,  2.007129, 2.094395,  2.181662,  2.268928,  2.356194,  2.443461,  2.530727, 2.617994,  2.705260,  2.792527,  2.879793,  2.967060,  3.054326, 3.141593,  3.228859,  3.316126,  3.403392,  3.490659,  3.577925, 3.665191,  3.752458,  3.839724,  3.926991,  4.014257,  4.101524, 4.188790,  4.276057,  4.363323,  4.450590,  4.537856,  4.625123, 4.712389};
  const Int_t nBinning_dphi = sizeof(binning_dphi)/sizeof(Int_t) -1;
 
  // Double_t etaFMDA[9] = {-3.4, -3.2, -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8};
  // Double_t etaFMDC[9] = {1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4};
  const Int_t netaFMD = fBinsFMDeta.size() -1;
  //const Int_t netaFMD = sizeof(etaFMDA)/sizeof(Int_t) -1;
  //Double_t etaTPC[9] = {-0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8};
  const Int_t netaTPC = fBinsTPCeta.size() -1;
  const Int_t ndeta = fNBinsdEta - 1;
  const Int_t sizeOfSamples = (Int_t) fNOfSamples;

  //Calculating dEta vector
  Double_t MindEta = TMath::Abs(fBinsTPCeta.back() - fBinsFMDeta.front());
  Double_t MaxdEta = TMath::Abs(fBinsTPCeta.front() - fBinsFMDeta.back());
  Double_t Incriment = (MaxdEta-MindEta)/fNBinsdEta;
  Double_t binning_deta[fNBinsdEta];

  for (Int_t i(0); i<fNBinsdEta; i++){
    binning_deta[i] = MindEta+Incriment*i;
  }

  

  Int_t iTrackBin_tpcfmdA[] = {netaTPC, netaFMD, 72, ndeta, 10, sizeOfSamples};
  Int_t nTrackBin_tpcfmd = sizeof(iTrackBin_tpcfmdA) / sizeof(Int_t);

  fhSE = new AliTHn("fhChargedSE", "fhChargedSE", nSteps, nTrackBin_tpcfmd, iTrackBin_tpcfmdA);
  fhME = new AliTHn("fhChargedME", "fhChargedME", nSteps, nTrackBin_tpcfmd, iTrackBin_tpcfmdA);
  
  fhSE->SetBinLimits(0, fBinsTPCeta.data());
  fhSE->SetBinLimits(1, fBinsFMDeta.data());
  fhSE->SetBinLimits(2, binning_dphi);
  fhSE->SetBinLimits(3, binning_deta);
  fhSE->SetBinLimits(4, -10,10);
  fhSE->SetBinLimits(5, 0,10);
  fhSE->SetVarTitle(0, "#eta (TPC)");
  fhSE->SetVarTitle(1, "#eta (FMD)");
  fhSE->SetVarTitle(2, "#Delta#phi");
  fhSE->SetVarTitle(3, "#Delta#eta");
  fhSE->SetVarTitle(4, "PVz [cm]");
  fhSE->SetVarTitle(5, "Sample");

  fOutputListCharged->Add(fhSE);

  fhME->SetBinLimits(0, fBinsTPCeta.data());
  fhME->SetBinLimits(1, fBinsFMDeta.data()); //CGHANGES THERE BACK TO A!
  fhME->SetBinLimits(2, binning_dphi);
  fhME->SetBinLimits(3, binning_deta);
  fhME->SetBinLimits(4, -10,10);
  fhME->SetBinLimits(5, 0,10);
  fhME->SetVarTitle(0, "#eta (TPC)");
  fhME->SetVarTitle(1, "#eta (FMD)");
  fhME->SetVarTitle(2, "#Delta#phi");
  fhME->SetVarTitle(3, "#Delta#eta");
  fhME->SetVarTitle(4, "PVz [cm]");
  fhME->SetVarTitle(5, "Sample");
  
  fOutputListCharged->Add(fhME);
  

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowEta::PrepareTPCTracks(){
  if(!fAOD) return kFALSE;
  if(!fTracksAss || !fTracksTrig || !fhTrigTracks) {AliError("Cannot prepare TPC tracks!"); return kFALSE; }

  fNofTracks = 0;
  Double_t binscont[3] = {fPVz, fSampleIndex, 0.};

  TObjArray* fTracksJets = nullptr;
  if(fVetoJetEvents) fTracksJets = new TObjArray;

  for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {
      AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
      if(!track || !IsTrackSelected(track)) { continue; }

      Double_t trackPt = track->Pt();
      if(fVetoJetEvents && trackPt > fJetParticleLowPt) fTracksJets->Add((AliAODTrack*)track);
      
      fNofTracks++;
      if(fAnalType == eFMDAFMDC || fIsTPCgen) continue;
    
    
      Double_t trackEta = track->Eta();
      if(trackPt > fPtMinTrig && trackPt < fPtMaxTrig) {
        if(fUseOppositeSidesOnly){
          if(fAnalType == eTPCFMDA && trackEta > 0.0) continue;
          if(fAnalType == eTPCFMDC && trackEta < 0.0) continue;
        }

        binscont[2] = trackPt;
        fhPT->Fill(trackPt);
        fTracksTrig->Add((AliAODTrack*)track);
        fhTrigTracks->Fill(binscont,0,1.);
      }
      
  } // tracks loop end

  if(fUseNch){
    if(fNofTracks < fNchMin || fNofTracks > fNchMax) { return kFALSE; }
    fhEventCounter->Fill("Nch cut ok ",1);
    fhEventMultiplicity->Fill(fNofTracks);
  }

  if(fVetoJetEvents){
    Double_t foundSomething = kFALSE;
    fhEventCounter->Fill("Before Jet Veto",1); //HPC = high pt cut
    for(Int_t iTrig(0); iTrig < fTracksJets->GetEntriesFast(); iTrig++){
      AliAODTrack* trackTrig = (AliAODTrack*)fTracksJets->At(iTrig);
      if(!trackTrig) continue;
      Double_t trigPhi = trackTrig->Phi();

      for(Int_t iAss(iTrig+1); iAss < fTracksJets->GetEntriesFast()+1; iAss++){
        AliAODTrack* trackAss = (AliAODTrack*)fTracksJets->At(iAss);
        if(!trackAss) continue;
        Double_t assPhi = trackAss->Phi();

        Double_t deltaPhi = RangePhi(trigPhi - assPhi);
        if(TMath::Abs(deltaPhi - TMath::Pi()) < 0.5) foundSomething = kTRUE;
      }
    }
    if(!foundSomething){ fhEventCounter->Fill("After Jet Veto",1); }
    else { return kFALSE; }
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowEta::PrepareFMDTracks(){
  if(!fTracksAss) { AliError("Problem with fTracksAss, terminating!"); return kFALSE; }
  if(fAnalType == eFMDAFMDC && !fTracksTrig) { AliError("Problem with fTracksTrig (no PID), terminating!"); return kFALSE; }

  
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
    //std::cout << "Valid: "<< valid << std::endl;

    if (!valid) continue;
    Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
    
    for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++)
    {
      // Bin content is most probable number of particles!
      Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
      Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
      //std::cout << "Phi: " <<  phi << "Eta: " << eta << std::endl;
      
      if(mostProbableN > 0) {
    	   if(eta > 0){
    	     nFMD_fwd_hits+=mostProbableN;
           if(fIsFMDgen) continue;
           if(eta > fFMDAacceptanceCutLower && eta < fFMDAacceptanceCutUpper){
              
             if(fAnalType == eTPCFMDA) {
               fTracksAss->Add(new AliPartSimpleForCorr(eta,phi,mostProbableN));
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
      
      if(fAnalType == eFMDAFMDC) fTracksTrig->Clear();
      return kFALSE;
    }
    
    AliAODVZERO *fvzero = fAOD->GetVZEROData();
    
    if(!fvzero) { AliError("Problem with VZEROData, terminating!"); return kFALSE; }
    
    Float_t nV0A_hits = fvzero->GetMTotV0A();
    Float_t nV0C_hits = fvzero->GetMTotV0C();
    
    fh2FMDvsV0[0]->Fill(nFMD_fwd_hits,nV0A_hits);
    fh2FMDvsV0[1]->Fill(nFMD_bwd_hits,nV0C_hits);

    if((nV0A_hits<(fFMDcutapar0*nFMD_fwd_hits-fFMDcutapar1)) || (nV0C_hits<(fFMDcutcpar0*nFMD_bwd_hits-fFMDcutcpar1))){
      fTracksAss->Clear();
      if(fAnalType == eFMDAFMDC) fTracksTrig->Clear();
      return kFALSE;
    }
    fhEventCounter->Fill("FMD cuts OK",1);
    fh2FMDvsV0[2]->Fill(nFMD_fwd_hits,nV0A_hits);
    fh2FMDvsV0[3]->Fill(nFMD_bwd_hits,nV0C_hits);
  }
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlowEta::PrepareMCTracks(){
  if(!fTracksAss || !fTracksTrig || !fhTrigTracks) {AliError("Cannot prepare MCC tracks!"); return kFALSE; }

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
      
      if(part->Charge()==0.) continue;

      
      if(fAnalType == eTPCFMDA || fAnalType == eTPCFMDC){
        if(partPt > fPtMinTrig && partPt < fPtMaxTrig){
          fTracksTrig->Add((AliMCParticle*)part);
          fhTrigTracks->Fill(binscont,0,1.);
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
        fTracksTrig->Add(new AliPartSimpleForCorr(partEta,partPhi,1.));
        fhTrigTracks->Fill(binscontFMD,0,1.);
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
Double_t AliAnalysisTaskCorrForFlowEta::TransverseBoost(const AliMCParticle *track){
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
void AliAnalysisTaskCorrForFlowEta::PrintSetup(){
  printf("\n\n\n ************** Parameters ************** \n");
  printf("\t fAnalType: (Int_t) %d\n", fAnalType);
  printf("\t fColSystem: (Int_t) %d\n", fColSystem);
  printf("\t fIsMC: (Bool_t) %s\n", fIsMC ? "kTRUE" : "kFALSE");
  printf("\t fIsTPCgen: (Bool_t) %s\n", fIsTPCgen ? "kTRUE" : "kFALSE");
  printf("\t fIsFMDgen: (Bool_t) %s\n", fIsFMDgen ? "kTRUE" : "kFALSE");
  printf("\t fUseNch: (Bool_t) %s\n", fUseNch ? "kTRUE" : "kFALSE");
  printf("\t fIsHMpp: (Bool_t) %s\n", fIsHMpp ? "kTRUE" : "kFALSE");
  printf("\t fUseEfficiency: (Bool_t) %s\n",  fUseEfficiency ? "kTRUE" : "kFALSE");
  printf("\t fUseOppositeSidesOnly: (Bool_t) %s\n", fUseOppositeSidesOnly ? "kTRUE" : "kFALSE");
  printf("\t fVetoJetEvents: (Bool_t) %s\n", fVetoJetEvents ? "kTRUE" : "kFALSE");
  printf("\t fRejectSecondariesFromMC: (Bool_t) %s\n", fRejectSecondariesFromMC ? "kTRUE" : "kFALSE");
  printf("\t fNOfSamples: (Int_t) %d\n", (Int_t) fNOfSamples);
  printf(" **************************** \n");
  printf("\t fSystematicsFlag: (TString) %s\n", fSystematicsFlag.Data());
  printf("\t fAbsEtaMax: (Double_t) %f\n", fAbsEtaMax);
  printf("\t fPtMinTrig -- fPtMaxTrig: (Double_t) %f -- %f\n", fPtMinTrig, fPtMaxTrig);
  printf("\t fCentMin -- fCentMax: (Double_t) %f -- %f\n", fCentMin, fCentMax);
  printf("\t fPVzCut: (Double_t) %f\n", fPVzCut);
  printf(" **************************** \n");
  printf("\t fUseFMDcut: (Bool_t) %s\n", fUseFMDcut ? "kTRUE" : "kFALSE");
  printf("\t fFMDcutapar0 -- fFMDcutapar1: (Double_t) %f -- %f\n", fFMDcutapar0, fFMDcutapar1);
  printf("\t fFMDcutcpar0 -- fFMDcutcpar1: (Double_t) %f -- %f\n", fFMDcutcpar0, fFMDcutcpar1);
  printf("\t fFMDacceptanceCut A - lower, upper: (Double_t) %f, %f\n", fFMDAacceptanceCutLower, fFMDAacceptanceCutUpper);
  printf("\t fFMDacceptanceCut C - lower, upper: (Double_t) %f, %f\n", fFMDCacceptanceCutLower, fFMDCacceptanceCutUpper);
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlowEta::DebugPrint(){
  for(Int_t i(0); i<5; i++){
    printf("\n");
  }
  printf("****** REACHED THIS POINT ******");
  for(Int_t i(0); i<5; i++){
    printf("\n");
  }
}