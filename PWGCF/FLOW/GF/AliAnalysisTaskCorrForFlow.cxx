/**************************************************************************
 *    Author:       Zuzana Moravcova                                      *
 *    Framework for calculating di-hadron correlation                     *
 *    for extraction of v_n{2} and v_n[2] coefficients.                   *
 *                                                                        *
 *    If used, modified, or distributed,                                  *
 *    please aknowledge the author of this code.                          *
 **************************************************************************/

#include "AliAnalysisTaskCorrForFlow.h"

using namespace std;

ClassImp(AliAnalysisTaskCorrForFlow);

AliAnalysisTaskCorrForFlow::AliAnalysisTaskCorrForFlow() : AliAnalysisTaskSE(),
    fAOD(0),
    fOutputListCharged(0),
    fInputListEfficiency(0),
    fTracksTrigCharged(0),
    fTracksAss(0),
    fPoolMgr(0),
    fhEventCounter(0),
    fhEventMultiplicity(0),
    fHistPhiEta(0),
    fhTrigTracks(0),
    fhChargedSE(0),
    fhChargedME(0),
    fNOfSamples(1.0),
    fSampleIndex(0.0),
    fTrigger(AliVEvent::kINT7),
    fIsHMpp(kFALSE),
    fUseNch(kFALSE),
    fUseEfficiency(kFALSE),
    fEfficiencyEtaDependent(kFALSE),
    fFilterBit(96),
    fbSign(0),
    fNofTracks(0),
    fNchMin(0),
    fNchMax(100000),
    fPtMinTrig(0.5),
    fPtMaxTrig(10.0),
    fPtMinAss(0.5),
    fPtMaxAss(10.0),
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
    fMergingCut(0.0),
    fSystematicsFlag(""),
    fPVzCut(10.),
    fCutDCAz(0.),
    fCutDCAxySigma(0.),
    fCutTPCchi2pCl(0.),
    fTPCclMin(70.),
    fEtaPolarity(0)

{}
//_____________________________________________________________________________
AliAnalysisTaskCorrForFlow::AliAnalysisTaskCorrForFlow(const char* name, Bool_t bUseEff) : AliAnalysisTaskSE(name),
    fAOD(0),
    fOutputListCharged(0),
    fInputListEfficiency(0),
    fTracksTrigCharged(0),
    fTracksAss(0),
    fPoolMgr(0),
    fhEventCounter(0),
    fhEventMultiplicity(0),
    fHistPhiEta(0),
    fhTrigTracks(0),
    fhChargedSE(0),
    fhChargedME(0),
    fNOfSamples(1.0),
    fSampleIndex(0.0),
    fTrigger(AliVEvent::kINT7),
    fIsHMpp(kFALSE),
    fUseNch(kFALSE),
    fUseEfficiency(bUseEff),
    fEfficiencyEtaDependent(kFALSE),
    fFilterBit(96),
    fbSign(0),
    fNofTracks(0),
    fNchMin(0),
    fNchMax(100000),
    fPtMinTrig(0.5),
    fPtMaxTrig(10.0),
    fPtMinAss(0.5),
    fPtMaxAss(10.0),
    fCentMin(0.0),
    fCentMax(10.0),
    fCentrality(-10.0),
    fAbsEtaMax(0.8),
    fPVz(100.0),
    fCentEstimator("V0M"),
    fPoolMaxNEvents(2000),
    fPoolMinNTracks(50),
    fMinEventsToMix(5),
    fNzVtxBins(10),
    fNCentBins(15),
    fMergingCut(0.0),
    fSystematicsFlag(""),
    fPVzCut(10.),
    fCutDCAz(0.),
    fCutDCAxySigma(0.),
    fCutTPCchi2pCl(0.),
    fTPCclMin(70.),
    fEtaPolarity(0)
{
    DefineInput(0, TChain::Class());
    if(bUseEff) { DefineInput(1, TList::Class()); }
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskCorrForFlow::~AliAnalysisTaskCorrForFlow()
{}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlow::UserCreateOutputObjects()
{
    OpenFile(1);
    PrintSetup();
    // //
    fzVtxBins = {-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0};
    fsampleBins = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    // //
    fOutputListCharged = new TList();
    fOutputListCharged->SetOwner(kTRUE);
    // //
    fhEventCounter = new TH1D("fhEventCounter","Event Counter",10,0,10);
    fOutputListCharged->Add(fhEventCounter);
    //
    fhEventMultiplicity = new TH1D("fhEventMultiplicity","Event multiplicity; N_{ch}",200,0,200);
    fOutputListCharged->Add(fhEventMultiplicity);
    //
    fHistPhiEta = new TH2D("fHistPhiEta", "fHistPhiEta; phi; eta", 100, 0.0, TMath::TwoPi(), 100, -1.0, 1.0);
    fOutputListCharged->Add(fHistPhiEta);
    //
    fhTrigTracks = new TH3D("fhTrigTracks", "fhTrigTracks; pT (trig); PVz; sample", fPtBinsTrigCharged.size() - 1, fPtBinsTrigCharged.data(), fzVtxBins.size()-1, fzVtxBins.data(), fsampleBins.size()-1, fsampleBins.data());
    fOutputListCharged->Add(fhTrigTracks);

    //
    // //mixing
    fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins, fCentBins.data(), fNzVtxBins, fzVtxBins.data());
    if (!fPoolMgr) { AliError("Event Pool manager not created!"); return; }
    fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);
    //
    // //sparses
    Int_t nSteps = 1;
    Double_t binning_deta_tpctpc[33] = {-1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,    0.1,  0.2,  0.3,  0.4,  0.5, 0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5, 1.6};
    Double_t binning_dphi[73] = { -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0,       0.087266,  0.174533,  0.261799,  0.349066,  0.436332, 0.523599,  0.610865,  0.698132,  0.785398,  0.872665,  0.959931, 1.047198,  1.134464,  1.221730,  1.308997,  1.396263,  1.483530, 1.570796,  1.658063,  1.745329,  1.832596,  1.919862,  2.007129, 2.094395,  2.181662,  2.268928,  2.356194,  2.443461,  2.530727, 2.617994,  2.705260,  2.792527,  2.879793,  2.967060,  3.054326, 3.141593,  3.228859,  3.316126,  3.403392,  3.490659,  3.577925, 3.665191,  3.752458,  3.839724,  3.926991,  4.014257,  4.101524, 4.188790,  4.276057,  4.363323,  4.450590,  4.537856,  4.625123, 4.712389};
    const Int_t sizePtTrig = fPtBinsTrigCharged.size() - 1;
    const Int_t sizePtAss = fPtBinsAss.size() - 1;
    const Int_t sizeOfSamples = (Int_t) fNOfSamples;
    const Int_t iBinningTPCTPC[] = {32,72,10,sizeOfSamples,sizePtTrig,sizePtAss};

    fhChargedSE = new AliTHn("fhChargedSE", "fhChargedSE", nSteps, 6, iBinningTPCTPC);
    fhChargedSE->SetBinLimits(0, binning_deta_tpctpc);
    fhChargedSE->SetBinLimits(1, binning_dphi);
    fhChargedSE->SetBinLimits(2, -10,10);
    fhChargedSE->SetBinLimits(3,  0.,10.);
    fhChargedSE->SetBinLimits(4, fPtBinsTrigCharged.data());
    fhChargedSE->SetBinLimits(5, fPtBinsAss.data());
    fhChargedSE->SetVarTitle(0, "#Delta#eta");
    fhChargedSE->SetVarTitle(1, "#Delta#phi");
    fhChargedSE->SetVarTitle(2, "PVz");
    fhChargedSE->SetVarTitle(3, "sample");
    fhChargedSE->SetVarTitle(4, "p_{T} [GeV/c] (trig)");
    fhChargedSE->SetVarTitle(5, "p_{T} [GeV/c] (ass)");
    fOutputListCharged->Add(fhChargedSE);
    //
    fhChargedME = new AliTHn("fhChargedME", "fhChargedME", nSteps, sizeof(iBinningTPCTPC) / sizeof(Int_t), iBinningTPCTPC);
    fhChargedME->SetBinLimits(0, binning_deta_tpctpc);
    fhChargedME->SetBinLimits(1, binning_dphi);
    fhChargedME->SetBinLimits(2, -10.,10.);
    fhChargedME->SetBinLimits(3,  0.,10.);
    fhChargedME->SetBinLimits(4, fPtBinsTrigCharged.data());
    fhChargedME->SetBinLimits(5, fPtBinsAss.data());
    fhChargedME->SetVarTitle(0, "#Delta#eta");
    fhChargedME->SetVarTitle(1, "#Delta#phi");
    fhChargedME->SetVarTitle(2, "PVz");
    fhChargedME->SetVarTitle(3, "sample");
    fhChargedME->SetVarTitle(4, "p_{T} [GeV/c] (trig)");
    fhChargedME->SetVarTitle(5, "p_{T} [GeV/c] (ass)");
    fOutputListCharged->Add(fhChargedME);
    //

      fInputListEfficiency = (TList*) GetInputData(1);
      if(fEfficiencyEtaDependent && fAbsEtaMax > 0.8) AliWarning("Efficiency loading -- eta can be out of range!");
      if(fSystematicsFlag.IsNull()) fSystematicsFlag = "Ev0_Tr0";
      if(!AreEfficienciesLoaded()) { AliError("Efficiencies not loaded!"); return; }


    //
    PostData(1, fOutputListCharged);
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlow::UserExec(Option_t *)
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
    fTracksTrigCharged = new TObjArray;
    fTracksAss = new TObjArray;

    if(fUseEfficiency && !AreEfficienciesLoaded()) { return; }

    fNofTracks = 0;

    for(Int_t i(0); i < iTracks; i++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }

        Double_t trackPt = track->Pt();
        if(trackPt > fPtMinAss && trackPt < fPtMaxAss) {
             if(fEtaPolarity == 0){
             fTracksAss->Add((AliAODTrack*)track);
             fNofTracks++;
             }
             else if(fEtaPolarity == -1 && track->Eta() > 0){
             fTracksAss->Add((AliAODTrack*)track);
             fNofTracks++;
             }
             else if(fEtaPolarity == 1 && track->Eta() < 0){
             fTracksAss->Add((AliAODTrack*)track);
             fNofTracks++;
             }
     }

        if(trackPt > fPtMinTrig && trackPt < fPtMaxTrig) {
          if(fEtaPolarity == 0){
               fTracksTrigCharged->Add((AliAODTrack*)track);
               fhTrigTracks->Fill(trackPt, fPVz, fSampleIndex);
          }
          else if(fEtaPolarity == -1 && track->Eta() < 0){
               fTracksTrigCharged->Add((AliAODTrack*)track);
               fhTrigTracks->Fill(trackPt, fPVz, fSampleIndex);
          }
          else if(fEtaPolarity == 1 && track->Eta() > 0){
               fTracksTrigCharged->Add((AliAODTrack*)track);
               fhTrigTracks->Fill(trackPt, fPVz, fSampleIndex);
          }
        }

        //example histogram
        fHistPhiEta->Fill(track->Phi(), track->Eta());
    }
    fhEventMultiplicity->Fill(fNofTracks);

    if(fUseNch){
      if(fNofTracks < fNchMin || fNofTracks > fNchMax) { return; }
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
    //
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlow::Terminate(Option_t *)
{
   if (fPoolMgr) delete fPoolMgr;
   // if(fOutputListCharged) delete fOutputListCharged;
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
     if(TMath::Abs(fPVz) >= fPVzCut) { return kFALSE; }
     fhEventCounter->Fill("PVzOK",1);

     fbSign = (InputEvent()->GetMagneticField() > 0) ? 1 : -1;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlow::IsTrackSelected(const AliAODTrack* track) const
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
  if(!fTracksTrigCharged || !fhTrigTracks || !fTracksAss) {AliError("Necessary inputs missing, terminating!"); return;}
  if(!fhChargedSE) { AliError(Form("Output AliTHn missing for ch , terminating!")); return; }

    Double_t binscont[6];
    binscont[2] = fPVz;
    binscont[3] = fSampleIndex;

    for(Int_t iTrig(0); iTrig < fTracksTrigCharged->GetEntriesFast(); iTrig++){
      AliVParticle* track = dynamic_cast<AliVParticle*>(fTracksTrigCharged->At(iTrig));
      if(!track) continue;
      AliAODTrack* trackAOD = nullptr;

      Double_t trigPt = track->Pt();
      Double_t trigEta = track->Eta();
      Double_t trigPhi = track->Phi();
      Double_t trigCharge = track->Charge();
      Double_t trigEff = 1.0;
      if(fUseEfficiency) trigEff = GetEff(trigPt, trigEta);
      binscont[4] = trigPt;

      for(Int_t iAss(0); iAss < fTracksAss->GetEntriesFast(); iAss++){
        AliVParticle* trackAss = dynamic_cast<AliVParticle*>(fTracksAss->At(iAss));
        if(!trackAss) continue;
        AliAODTrack* trackAODAss = nullptr;
        trackAODAss = (AliAODTrack*)fTracksAss->At(iAss);

        Double_t assPt = trackAss->Pt();
        Double_t assEta = trackAss->Eta();
        Double_t assPhi = trackAss->Phi();
        Double_t assCharge = trackAss->Charge();
        Double_t assEff = 1.0;
        if(fUseEfficiency) assEff = GetEff(assPt, assEta);

        binscont[0] = trigEta - assEta;
        binscont[1] = RangePhi(trigPhi - assPhi);
        binscont[5] = assPt;

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

        fhChargedSE->Fill(binscont,0,1./(trigEff*assEff));
      }
    }

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrForFlow::FillCorrelationsMixed()
{
  if(!fTracksTrigCharged || !fhTrigTracks || !fTracksAss) { AliError("Necessary inputs missing, terminating!"); return; }

  AliEventPool *pool = fPoolMgr->GetEventPool(fCentrality, fPVz);
  if(!pool) { AliError(Form("No pool found for centrality = %f, zVtx = %f", fCentrality,fPVz)); return; }

  if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() > fMinEventsToMix) {
    Int_t nMix = pool->GetCurrentNEvents();

      Double_t binscont[6];
      binscont[2] = fPVz;
      binscont[3] = fSampleIndex;

      for(Int_t iTrig(0); iTrig < fTracksTrigCharged->GetEntriesFast(); iTrig++){
        AliVParticle* track = dynamic_cast<AliVParticle*>(fTracksTrigCharged->At(iTrig));
        if(!track) continue;

        Double_t trigPt = track->Pt();
        Double_t trigEta = track->Eta();
        Double_t trigPhi = track->Phi();
        Double_t trigCharge = track->Charge();
        binscont[4] = trigPt;

        for(Int_t eMix(0); eMix < nMix; eMix++){
          TObjArray *mixEvents = pool->GetEvent(eMix);
          for(Int_t iAss(0); iAss < mixEvents->GetEntriesFast(); iAss++){
            AliVParticle* trackAss = dynamic_cast<AliVParticle*>(mixEvents->At(iAss));
            if(!trackAss) continue;

            Double_t assPt = trackAss->Pt();
            Double_t assEta = trackAss->Eta();
            Double_t assPhi = trackAss->Phi();
            Double_t assCharge = trackAss->Charge();

            binscont[0] = trigEta - assEta;
            binscont[1] = RangePhi(trigPhi - assPhi);
            binscont[5] = assPt;

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
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrForFlow::AreEfficienciesLoaded()
{
  if(!fInputListEfficiency) {AliError("Efficiency input list not loaded"); return kFALSE; }
  TString etaReg[8] = {"0020", "0200", "0204", "0402", "0406", "0604", "0608", "0806"};
  for(Int_t eta(0); eta < 8; eta++){
      fhEfficiencyEta[eta] = (TH2D*)fInputListEfficiency->FindObject(Form("LHC17f2b_ch_Eta_%s_%s_wFD", etaReg[eta].Data(),fSystematicsFlag.Data()));
      if(!fhEfficiencyEta[eta]) {AliError(Form("Efficiency (ch, eta region %s, flag %s) not loaded",etaReg[eta].Data(),fSystematicsFlag.Data())); return kFALSE; }
  }
    fhEventCounter->Fill("Efficiencies loaded",1);
    return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForFlow::GetEff(const Double_t dPt, const Double_t dEta)
{
  if(!fUseEfficiency) return 1.0;
    Int_t region = GetEtaRegion(dEta);
    if(region < 0) { AliWarning("Invalid region, returning efficiency 1.0."); return 1.0; }
    if(!fhEfficiencyEta[region]) { AliError("Efficiency histogram not found, returning efficiency 1.0."); return 1.0; }
    return fhEfficiencyEta[region]->GetBinContent(fhEfficiencyEta[region]->FindFixBin(dPt, fCentrality));
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskCorrForFlow::GetEtaRegion(const Double_t dEta){
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
void AliAnalysisTaskCorrForFlow::PrintSetup(){
  printf("\n\n\n ************** Parameters ************** \n");
  printf("\t fUseNch: (Bool_t) %s\n",    fUseNch ? "kTRUE" : "kFALSE");
  printf("\t fIsHMpp: (Bool_t) %s\n",    fIsHMpp ? "kTRUE" : "kFALSE");
  printf("\t fUseEfficiency: (Bool_t) %s\n",    fUseEfficiency ? "kTRUE" : "kFALSE");
  printf("\t fEfficiencyEtaDependent: (Bool_t) %s\n",    fEfficiencyEtaDependent ? "kTRUE" : "kFALSE");
  printf("\n **************************** \n");
  printf("\t fAbsEtaMax: (Double_t) %f\n",    fAbsEtaMax);
  printf("\t fPtMinTrig -- fPtMaxTrig: (Double_t) %f -- %f\n",    fPtMinTrig, fPtMaxTrig);
  printf("\t fPtMinAss -- fPtMaxAss: (Double_t) %f -- %f\n",    fPtMinAss, fPtMaxAss);
  printf("\t fCentMin -- fCentMax: (Double_t) %f -- %f\n",    fCentMin, fCentMax);
}
