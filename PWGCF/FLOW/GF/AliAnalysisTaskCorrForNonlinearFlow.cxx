/* -------------------------------------------
 * Maintainer: Mingrui Zhao
 */
#include "AliAnalysisTaskCorrForNonlinearFlow.h"
#include "AliGFWCuts.h"
#include "AliGFWNFCuts.h"
#include "AliGFWWeights.h"
#include "CorrelationCalculator.h"

#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TComplex.h>
#include <TBits.h>
// AliRoot includes
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"

// ROOT includes
#include <TList.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMatrixDSym.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TComplex.h>

// AliRoot includes
#include "AliEventCuts.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisUtils.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliAODITSsaTrackCuts.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliVVertex.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"

// STL includes
#include <iostream>
using std::cout;
using std::endl;
ClassImp(AliAnalysisTaskCorrForNonlinearFlow)

// ---------------------------------------------------------------------------------
AliAnalysisTaskCorrForNonlinearFlow::~AliAnalysisTaskCorrForNonlinearFlow() {
    // Destructor
}

// ---------------------------------------------------------------------------------
void AliAnalysisTaskCorrForNonlinearFlow::UserCreateOutputObjects() {
    // Create output objects
    fListOfObjects = new TList();
    fListOfObjects->SetOwner();

    // Setting for AliEventCuts:
    fEventCuts.AddQAplotsToList(fListOfObjects);

    if (fPeriod.EqualTo("LHC15o")) {
        // Only for LHC15o pass1
    fGFWSelection15o = new AliGFWNFCuts();
    fGFWSelection15o->PrintSetup();
  } else {
    fGFWSelection = new AliGFWCuts();
    fGFWSelection->PrintSetup();
    }

    // mixing
    fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins, fCentBins.data(), fNzVtxBins, fzVtxBins.data());
    if (!fPoolMgr) {
    	AliError("Event Pool manager not created!");
	return;
    }
    fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);

    Int_t nSteps = 1;
    Double_t binning_deta_tpctpc[33] = {-1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,    0.1,  0.2,  0.3,  0.4,  0.5, 0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5, 1.6};
    Double_t binning_dphi[73] = { -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0,       0.087266,  0.174533,  0.261799,  0.349066,  0.436332, 0.523599,  0.610865,  0.698132,  0.785398,  0.872665,  0.959931, 1.047198,  1.134464,  1.221730,  1.308997,  1.396263,  1.483530, 1.570796,  1.658063,  1.745329,  1.832596,  1.919862,  2.007129, 2.094395,  2.181662,  2.268928,  2.356194,  2.443461,  2.530727, 2.617994,  2.705260,  2.792527,  2.879793,  2.967060,  3.054326, 3.141593,  3.228859,  3.316126,  3.403392,  3.490659,  3.577925, 3.665191,  3.752458,  3.839724,  3.926991,  4.014257,  4.101524, 4.188790,  4.276057,  4.363323,  4.450590,  4.537856,  4.625123, 4.712389};
    std::vector<Double_t>   fPtBinsTrigCharged;    // I don't want to set the things outside
    std::vector<Double_t>   fPtBinsAss;            // I don't want to set the things outside
    const Int_t sizePtTrig = fPtBinsTrigCharged.size() - 1;
    const Int_t sizePtAss = fPtBinsAss.size() - 1;
    const Int_t sizeOfSamples = 1; // (Int_t) fNOfSamples; Subsample, don't use it so far
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
    fListOfObjects->Add(fhChargedSE);
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
    fListOfObjects->Add(fhChargedME);
}

// ---------------------------------------------------------------------------------
void AliAnalysisTaskCorrForNonlinearFlow::NotifyRun() {
    if (fAddTPCPileupCuts) {
      Bool_t dummy = fEventCuts.AcceptEvent(InputEvent());
      fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
      fEventCuts.fESDvsTPConlyLinearCut[0] = fESDvsTPConlyLinearCut;
    }
}

// ---------------------------------------------------------------------------------
void AliAnalysisTaskCorrForNonlinearFlow::UserExec(Option_t *) {
  // Mingrui: apply the bootstrap later
  // bootstrap_value = rand.Integer(30);

  // Check if it can pass the trigger
  //..apply physics selection
  UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isTrigselected = false;
  if (fTrigger == 0) {
    isTrigselected = fSelectMask&AliVEvent::kINT7;
    fAliTrigger = AliVEvent::kINT7;
  } else if (fTrigger == 1) {
    isTrigselected = fSelectMask&AliVEvent::kHighMultV0;
    fAliTrigger = AliVEvent::kHighMultV0;
  }
  if(isTrigselected == false) return;

  //..check if I have AOD
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) {
    Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    return;
  }

  if (fLowMultiplicityMode) {
     // Number of AOD tracks before track cuts
     const int nAODTracks = fAOD->GetNumberOfTracks();
     if (nAODTracks > 200) {
       PostData(1,fListOfObjects);
       int outputslot = 2;
       PostData(2, fListOfProfile);
       for (int i = 0; i < 30; i++) {
         outputslot++;
         PostData(outputslot, fListOfProfiles[i]);
       }
       return;
     }
  }

  // Check if it passed the standard AOD selection
  if (!AcceptAOD(fAOD) ) {
    PostData(1,fListOfObjects);
    int outputslot = 2;
    PostData(2, fListOfProfile);
    for (int i = 0; i < 30; i++) {
      outputslot++;
      PostData(outputslot, fListOfProfiles[i]);
    }
    return;
  }
  hEventCount->Fill("after fEventCuts", 1.);

  if (fPeriod.EqualTo("LHC15o")) { // Only for LHC15o pass1
    fGFWSelection15o->ResetCuts();
  } else {
    fGFWSelection->ResetCuts();
  }
  //..filling Vz distribution
  AliVVertex *vtx = fAOD->GetPrimaryVertex();
  float fVtxZ = vtx->GetZ();

  if (fPeriod.EqualTo("LHC15o")) { // Only for LHC15o pass1
	   if (!fGFWSelection15o->AcceptVertex(fAOD)) {
	    PostData(1,fListOfObjects);
	    int outputslot = 2;
	    PostData(2, fListOfProfile);
	    for (int i = 0; i < 30; i++) {
	      outputslot++;
	      PostData(outputslot, fListOfProfiles[i]);
	    }
	    return;
	  }
  } else {
	  if (!fGFWSelection->AcceptVertex(fAOD)) {
	    PostData(1,fListOfObjects);
	    int outputslot = 2;
	    PostData(2, fListOfProfile);
	    for (int i = 0; i < 30; i++) {
	      outputslot++;
	      PostData(outputslot, fListOfProfiles[i]);
	    }
	    return;
	  }
  }

  // checking the run number for aplying weights & loading TList with weights
  //
  if (lastRunNumber != fAOD->GetRunNumber()) {
    lastRunNumber = fAOD->GetRunNumber();
    if (fPeriod.EqualTo("LHC15oKatarina")) {
      if (fNUA && !LoadWeightsKatarina()) {
        AliFatal("Trying to Load Systematics but weights not loaded!");
        return;
      }
      if (fNUE && !LoadPtWeightsKatarina()) {
        AliFatal("PtWeights not loaded!");
        return;
      }

    } else {
      if (fNUA && !LoadWeightsSystematics()) {
        AliFatal("Trying to Load Systematics but weights not loaded!");
        return;
      }
      if (fNUE && !LoadPtWeights()) {
        AliFatal("PtWeights not loaded!");
        return;
      }
    }

  }

  // Here we calcuate the multiplicity distribution
  NTracksCalculation(fInputEvent);

  Int_t nTracks = fInputEvent->GetNumberOfTracks();
  fTracksTrigCharged = new TObjArray;
  fTracksAss = new TObjArray;

  //..for DCA
  double pos[3], vz, vx, vy;
  vz = fInputEvent->GetPrimaryVertex()->GetZ();
  vx = fInputEvent->GetPrimaryVertex()->GetX();
  vy = fInputEvent->GetPrimaryVertex()->GetY();
  double vtxp[3] = {vx, vy, vz};
  

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(iTrack));

    // Require track to be existing and pass the track selection
    if (!track) continue;
    track->GetXYZ(pos);
    if (!AcceptAODTrack(track, pos,vtxp)) continue;
    
    Double_t pt = track->Pt();
    if (pt > fPtMinAss && pt < fPtMaxAss) {
        // Mingrui Polarity ??
        fTracksAss->Add(track);
        // fNofTracksAss++; // number of associate tracks in the event
    }

    if (pt > fPtMinTrig && pt < fPtMaxTrig) {
        fTracksTrigCharged->Add(track);
        // fNofTracksTrig++; // number of trigger tracks in the event
        // fhTracksTrigPt->Fill(pt, fPVz);
    }
  }

  if (!fTracksTrigCharged) {
      FillCorrelations();
      FillCorrelationsMixed();
  }

  fTracksTrigCharged->Clear();
  delete fTracksTrigCharged;

  fTracksAss->Clear();
  delete fTracksAss;

  PostData(1, fListOfObjects);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskCorrForNonlinearFlow::Terminate(Option_t *)
{
    if (fPoolMgr) {
        delete fPoolMgr;
    }
}

//________________________________________________________________________
void AliAnalysisTaskCorrForNonlinearFlow::FillCorrelations() {
    if (!fTracksTrigCharged || !fTracksAss) {
        AliError("Necessary inputs missing, terminating!"); return;
    }
    /* don't check this Mingrui
    if(!fhChargedSE) { 
	AliError(Form("Output AliTHn missing for ch , terminating!")); return; 
        return;
    }
    */

    Double_t binscont[6];
    binscont[2] = fPVz;
    binscont[3] = 1;

    for (Int_t iTrig = 0; iTrig < fTracksTrigCharged->GetEntriesFast(); iTrig++) {
        AliAODTrack* trackTrig = dynamic_cast<AliAODTrack*>(fTracksTrigCharged->At(iTrig));

        Double_t ptTrig = trackTrig->Pt();
        Double_t phiTrig = trackTrig->Phi();
        Double_t etaTrig = trackTrig->Eta();
        Double_t chargeTrig = trackTrig->Charge();
        Double_t trigEff = 1.0; // Efficiency

        for (Int_t iAss = 0; iAss < fTracksAss->GetEntriesFast(); iAss++) {
            AliAODTrack* trackAss = dynamic_cast<AliAODTrack*>(fTracksTrigCharged->At(iAss));

            Double_t ptAss = trackAss->Pt();
            Double_t phiAss = trackAss->Phi();
            Double_t etaAss = trackAss->Eta();
            Double_t chargeAss = trackAss->Charge();
            Double_t assEff = 1.0; // Efficiency

            //..check if the tracks are the same
            // Mingrui: I don't see Zuzana uses this
            // if (trackTrig == trackAss) continue;

            // Here these's a complicated way to check whether the pair pass or not

            fhChargedSE->Fill(binscont, 0, 1./(trigEff*assEff));
        }
    }
}

//________________________________________________________________________
void AliAnalysisTaskCorrForNonlinearFlow::FillCorrelationsMixed() {
    if (!fTracksTrigCharged  || !fTracksAss) {
        AliError("Necessary inputs missing, terminating!"); return;
    }

    Double_t binscont[6];
    binscont[2] = fPVz;
    binscont[3] = 1;

    AliEventPool* pool = fPoolMgr->GetEventPool(fCentrality, fPVz);
    if (!pool) {
    }

    if (pool->IsReady()) {
        int nMix = pool->GetCurrentNEvents();
        for (Int_t iTrig = 0; iTrig < fTracksTrigCharged->GetEntriesFast(); iTrig++) {
            AliAODTrack* trackTrig = dynamic_cast<AliAODTrack*>(fTracksTrigCharged->At(iTrig));

            Double_t ptTrig = trackTrig->Pt();
            Double_t phiTrig = trackTrig->Phi();
            Double_t etaTrig = trackTrig->Eta();
            Double_t chargeTrig = trackTrig->Charge();

            for (Int_t iMix = 0; iMix < nMix; iMix++) {
                TObjArray* mixTracks = pool->GetEvent(iMix);
                for (Int_t iAss = 0; iAss < mixTracks->GetEntriesFast(); iAss++) {
                    AliAODTrack* trackAss = dynamic_cast<AliAODTrack*>(mixTracks->At(iAss));

                    Double_t ptAss = trackAss->Pt();
                    Double_t phiAss = trackAss->Phi();
                    Double_t etaAss = trackAss->Eta();
                    Double_t chargeAss = trackAss->Charge();

                    //..check if the tracks are the same
                    //
                    // if (trackTrig == trackAss) continue;
                    // Here these's a complicated way to check whether the pair pass or not

                    fhChargedME->Fill(binscont, 0, 1./nMix);
                }
            }
        }
    }
    TObjArray* cloneArray = (TObjArray*)fTracksTrigCharged->Clone();
    cloneArray->SetOwner(kTRUE);
    pool->UpdatePool(cloneArray);

    return;
}

Bool_t AliAnalysisTaskCorrForNonlinearFlow::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp) {
  // Pt cut
  if(mtr->Pt() < fMinPt) return kFALSE;
  if(mtr->Pt() > fMaxPt) return kFALSE;

  // DCA cut
  if(ltrackXYZ && vtxp) {
    mtr->GetXYZ(ltrackXYZ);
    ltrackXYZ[0] = ltrackXYZ[0]-vtxp[0];
    ltrackXYZ[1] = ltrackXYZ[1]-vtxp[1];
    ltrackXYZ[2] = ltrackXYZ[2]-vtxp[2];
  } else return kFALSE; //DCA cut is a must for now

  // Additional cut for TPCchi2perCluster
  if (mtr->GetTPCchi2perCluster()>fTPCchi2perCluster) return kFALSE;

  if (fPeriod.EqualTo("LHC15o")) { // Only for LHC15o pass1
    return fGFWSelection15o->AcceptTrack(mtr,ltrackXYZ,0,kFALSE);
  } else {
    return fGFWSelection->AcceptTrack(mtr,ltrackXYZ,0,kFALSE);
  }
}

Bool_t AliAnalysisTaskCorrForNonlinearFlow::AcceptAOD(AliAODEvent *inEv) {
  // LHC15i, LHC15l, LHC16, LHC17, LHC18: means: pp sample
  if (fPeriod.EqualTo("LHC15i") ||
      fPeriod.EqualTo("LHC15l") ||
      fPeriod.EqualTo("LHC16Preview") ||
      fPeriod.EqualTo("LHC17Preview") ||
      fPeriod.EqualTo("LHC18Preview") || 
      fPeriod.EqualTo("LHC16") ||
      fPeriod.EqualTo("LHC17") ||
      fPeriod.EqualTo("LHC18") || 
      fPeriod.EqualTo("LHC16ZM") ||
      fPeriod.EqualTo("LHC17ZM") ||
      fPeriod.EqualTo("LHC18ZM") ) {
    fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true);
  }


  if(!fEventCuts.AcceptEvent(inEv)) return false;

  // Primary vertex
  const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(inEv->GetPrimaryVertex());
  if(!vtx || vtx->GetNContributors() < 1)
    return kFALSE;

  // SPD Vertex
  const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(inEv->GetPrimaryVertexSPD());
  Double_t dMaxResol = 0.25; // suggested from DPG
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return kFALSE;

  if (fPeriod.EqualTo("LHC15o") ||
      fPeriod.EqualTo("LHC15o_pass2") ||
      fPeriod.EqualTo("LHC18qr_pass3") ||
      fPeriod.EqualTo("LHC16qt") ||
      fPeriod.EqualTo("LHC17n") ||
      fPeriod.EqualTo("LHC15oKatarina")) {
    // return false;
  } else {
    // if(fAOD->IsPileupFromSPDInMultBins() ) { return false; }

    // AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
    // if (!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return false; }

    // if(!multSelection->GetThisEventIsNotPileup() || !multSelection->GetThisEventIsNotPileupInMultBins() || !multSelection->GetThisEventHasNoInconsistentVertices() || !multSelection->GetThisEventPassesTrackletVsCluster()) { return false; }

    // Int_t nTracksPrim = fAOD->GetPrimaryVertex()->GetNContributors();
    // if(nTracksPrim < 0.5) { return false; }
  }

  // Vertex Z
  const Double_t aodVtxZ = vtx->GetZ();
  if(TMath::Abs(aodVtxZ) > 10) return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskCorrForNonlinearFlow::LoadWeightsSystematics() {

  // If the period is not pPb LHC16qt
  if (!fPeriod.EqualTo("LHC16qt")) {
    // Only if it is the new LHC16,17,18, We need the period NUA
    if (fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18")) {
      std::string ppperiod = ReturnPPperiod(fAOD->GetRunNumber());
      if(fCurrSystFlag == 0) fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%s", ppperiod.c_str()));
      if(!fWeightsSystematics)
      {
        printf("Weights could not be found in list!\n");
        return kFALSE;
      }
      fWeightsSystematics->CreateNUA();
    } else {
      if(fCurrSystFlag == 0 || fUseDefaultWeight) fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i",fAOD->GetRunNumber()));
      else fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i_SystFlag%i_",fAOD->GetRunNumber(), fCurrSystFlag));
      if(!fWeightsSystematics)
      {
        printf("Weights could not be found in list!\n");
        return kFALSE;
      }
      fWeightsSystematics->CreateNUA();
    }

  } 
  // If it is the pPb LHC16qt
  else {
    int EvFlag = 0, TrFlag = 0;
    if (fCurrSystFlag == 0) EvFlag = 0, TrFlag = 0;
    if (fCurrSystFlag == 1) EvFlag = 0, TrFlag = 1;
    if (fCurrSystFlag == 2) EvFlag = 0, TrFlag = 3;
    if (fCurrSystFlag == 3) EvFlag = 0, TrFlag = 0; // Abandoned
    if (fCurrSystFlag == 4) EvFlag = 0, TrFlag = 2;
    if (fCurrSystFlag == 5) EvFlag = 0, TrFlag = 3;
    if (fCurrSystFlag == 6) EvFlag = 0, TrFlag = 8;
    if (fCurrSystFlag == 7) EvFlag = 0, TrFlag = 9;
    if (fCurrSystFlag == 8) EvFlag = 0, TrFlag = 10;

    if (fCurrSystFlag == 17) EvFlag = 1, TrFlag = 0;
    if (fCurrSystFlag == 18) EvFlag = 2, TrFlag = 0;
    if (fCurrSystFlag == 19) EvFlag = 3, TrFlag = 0;


  fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i_Ev%d_Tr%d",fAOD->GetRunNumber(),EvFlag,TrFlag));
  if(!fWeightsSystematics)
  {
    printf("Weights could not be found in list!\n");
    return kFALSE;
  }
  fWeightsSystematics->CreateNUA();
  }
  return kTRUE;
}

Bool_t AliAnalysisTaskCorrForNonlinearFlow::LoadPtWeights() {
  int EvFlag = 0, TrFlag = 0;

  // If the period is **NOT** pPb LHC16qt
  if ( !fPeriod.EqualTo("LHC16qt") &&
       !fPeriod.EqualTo("LHC16") && !fPeriod.EqualTo("LHC17") && !fPeriod.EqualTo("LHC18") &&
       !fPeriod.EqualTo("LHC16Preview") && !fPeriod.EqualTo("LHC17Preview") && !fPeriod.EqualTo("LHC18Preview")
		  ) {
    if(fCurrSystFlag == 0) fPtWeightsSystematics = (TH1D*)fFlowPtWeightsList->FindObject(Form("EffRescaled_Cent0"));
    else fPtWeightsSystematics = (TH1D*)fFlowPtWeightsList->FindObject(Form("EffRescaled_Cent0_SystFlag%i_", fCurrSystFlag));
    if(!fPtWeightsSystematics)
    {
      printf("PtWeights could not be found in list!\n");
      return kFALSE;
    }
  } 
  // If it is the pPb LHC16qt 
  else {
    if (fCurrSystFlag == 0) EvFlag = 0, TrFlag = 0;
    if (fCurrSystFlag == 1) EvFlag = 0, TrFlag = 1;
    if (fCurrSystFlag == 2) EvFlag = 0, TrFlag = 3;
    if (fCurrSystFlag == 3) EvFlag = 0, TrFlag = 0; // Abandoned
    if (fCurrSystFlag == 4) EvFlag = 0, TrFlag = 2;
    if (fCurrSystFlag == 5) EvFlag = 0, TrFlag = 3;
    if (fCurrSystFlag == 6) EvFlag = 0, TrFlag = 8;
    if (fCurrSystFlag == 7) EvFlag = 0, TrFlag = 9;
    if (fCurrSystFlag == 8) EvFlag = 0, TrFlag = 10;

    if (fCurrSystFlag == 17) EvFlag = 1, TrFlag = 0;
    if (fCurrSystFlag == 18) EvFlag = 2, TrFlag = 0;
    if (fCurrSystFlag == 19) EvFlag = 3, TrFlag = 0;

    if (fPeriod.EqualTo("LHC16qt")) {
        fPtWeightsSystematics = (TH1D*)fFlowPtWeightsList->FindObject(Form("LHC17f2b_ch_Eta_0020_Ev%d_Tr%d", EvFlag, TrFlag));

        cout << "Trying to load" << Form("LHC17f2b_ch_Eta_0020_Ev%d_Tr%d", EvFlag, TrFlag) << endl;
        fPtWeightsFeeddown    = (TH1D*)fFlowFeeddownList->FindObject(Form("LHC17f2b_ch_Eta_0020_Ev%d_Tr%d", EvFlag, TrFlag));
        if(!fPtWeightsSystematics)
        {
            printf("pPb: PtWeights could not be found in list!\n");
            return kFALSE;
        }
    } else {
	std::string period = ReturnPPperiodMC(fAOD->GetRunNumber());
        fPtWeightsSystematics = (TH1D*)fFlowPtWeightsList->FindObject(Form("LHC%s_ch_Ev%d_Tr%d", period.c_str(), EvFlag, TrFlag));
        fPtWeightsFeeddown    = (TH1D*)fFlowFeeddownList->FindObject(Form("LHC%s_ch_Ev%d_Tr%d", period.c_str(), EvFlag, TrFlag));
        if(!fPtWeightsSystematics)
        {
            printf("pp: PtWeights could not be found in list!\n");
            return kFALSE;
        }
    }
  }
  return kTRUE;

}

double AliAnalysisTaskCorrForNonlinearFlow::GetWeightKatarina(double phi, double eta, double vz) {
  double weight = hPhiWeightRun->GetBinContent(hPhiWeightRun->GetXaxis()->FindBin(phi),
      hPhiWeightRun->GetYaxis()->FindBin(eta),
      hPhiWeightRun->GetZaxis()->FindBin(vz));
  return weight;
}

// Load Katarina's weights
Bool_t AliAnalysisTaskCorrForNonlinearFlow::LoadWeightsKatarina() {
  hPhiWeightRun = (TH3F*)fFlowWeightsList->FindObject(Form("fPhiWeight_%0.lf", (double)(fAOD->GetRunNumber())));
  if (!hPhiWeightRun) {
    printf("Weights could not be found in list!\n");
    return kFALSE;
  }
  return kTRUE;
}

double AliAnalysisTaskCorrForNonlinearFlow::GetPtWeightKatarina(double pt, double eta, double vz)
{
  double weight = 1;
  double binPt = hTrackEfficiencyRun->GetXaxis()->FindBin(pt);
  double binEta = hTrackEfficiencyRun->GetYaxis()->FindBin(eta);
  double binVz = hTrackEfficiencyRun->GetZaxis()->FindBin(vz);
  //..take into account error on efficiency: randomly get number from gaussian distribution of eff. where width = error
  double eff = hTrackEfficiencyRun->GetBinContent(binPt, binEta, binVz);
  double error = hTrackEfficiencyRun->GetBinError(binPt, binEta, binVz);

  if((eff < 0.03) || ((error/eff) > 0.1)) weight = 1;
  else {
    TRandom3 r(0);
    double efficiency = 0;
    efficiency = r.Gaus(eff, error);
    weight = 1./efficiency; //..taking into account errors
    //weight = 1./eff;
  }
  return weight;
}

// Load Katarina's pt weights
Bool_t AliAnalysisTaskCorrForNonlinearFlow::LoadPtWeightsKatarina() {
  hTrackEfficiencyRun = (TH3F*)fFlowPtWeightsList->FindObject(Form("eff_LHC15o_HIJING_%.0lf", (double)(fAOD->GetRunNumber())));
  if (!hTrackEfficiencyRun) {
    printf("Pt Weights could not be found in list!\n");
    return kFALSE;
  }
  return kTRUE;
}

Double_t AliAnalysisTaskCorrForNonlinearFlow::GetFlowWeightSystematics(const AliVParticle* track, double fVtxZ, const PartSpecies species) {

  double dPhi = track->Phi();
  double dEta = track->Eta();
  double dVz = fVtxZ;
  double dWeight = 1.0;
  dWeight = fWeightsSystematics->GetNUA(dPhi, dEta, dVz);
  return dWeight;
}

Bool_t AliAnalysisTaskCorrForNonlinearFlow::LoadWeights() {
  // (Re-) Loading of flow vector weights
  // ***************************************************************************
  if(!fFlowWeightsList) { AliError("Flow weights list not found! Terminating!"); return kFALSE; }

  TList* listFlowWeights = nullptr;

  TString fFlowWeightsTag = "";
  if(!fFlowWeightsTag.IsNull()) {
    // using weights Tag if provided (systematics)
    listFlowWeights = (TList*) fFlowWeightsList->FindObject(fFlowWeightsTag.Data());
    if(!listFlowWeights) { AliError(Form("TList with tag '%s' not found!",fFlowWeightsTag.Data())); fFlowWeightsList->ls(); return kFALSE; }
  } else {
    if(!fFlowRunByRunWeights && !fFlowPeriodWeights) {
      // loading run-averaged weights
      listFlowWeights = (TList*) fFlowWeightsList->FindObject("averaged");
      if(!listFlowWeights) { AliError("TList with flow run-averaged weights not found."); fFlowWeightsList->ls(); return kFALSE; }
    } else if(fFlowPeriodWeights){
      // loading period-specific weights
      listFlowWeights = (TList*) fFlowWeightsList->FindObject(ReturnPPperiod(fAOD->GetRunNumber()));
      if(!listFlowWeights) { AliError("Loading period weights failed!"); fFlowWeightsList->ls(); return kFALSE; }
    }
    else {
      // loading run-specific weights
      listFlowWeights = (TList*) fFlowWeightsList->FindObject(Form("%d",fAOD->GetRunNumber()));

      if(!listFlowWeights) {
        // run-specific weights not found for this run; loading run-averaged instead
        AliWarning(Form("TList with flow weights (run %d) not found. Using run-averaged weights instead (as a back-up)", fAOD->GetRunNumber()));
        listFlowWeights = (TList*) fFlowWeightsList->FindObject("averaged");
        if(!listFlowWeights) { AliError("Loading run-averaged weights failed!"); fFlowWeightsList->ls(); return kFALSE; }
      }
    }
  }


  for(Int_t iSpec(0); iSpec <= kRefs; ++iSpec) {
    if(fFlowUse3Dweights) {
      fh3Weights[iSpec] = (TH3D*) listFlowWeights->FindObject(Form("%s3D",GetSpeciesName(PartSpecies(iSpec))));
      if(!fh3Weights[iSpec]) { AliError(Form("Weight 3D (%s) not found",GetSpeciesName(PartSpecies(iSpec)))); return kFALSE; }
    } else {
      fh2Weights[iSpec] = (TH2D*) listFlowWeights->FindObject(GetSpeciesName(PartSpecies(iSpec)));
      if(!fh2Weights[iSpec]) { AliError(Form("Weight 2D (%s) not found",GetSpeciesName(PartSpecies(iSpec)))); return kFALSE; }
    }
  }

  return kTRUE;
}

Double_t AliAnalysisTaskCorrForNonlinearFlow::GetFlowWeight(const AliVParticle* track, double fVtxZ, const PartSpecies species) {
  // if not applying for reconstructed
  // if(!fFlowWeightsApplyForReco && HasMass(species)) { return 1.0; }

  Double_t dWeight = 1.0;
  if(fFlowUse3Dweights) {
    Int_t iBin = fh3Weights[species]->FindFixBin(track->Phi(),track->Eta(),fVtxZ);
    dWeight = fh3Weights[species]->GetBinContent(iBin);
  } else {
    Int_t iBin = fh2Weights[species]->FindFixBin(track->Phi(),track->Eta());
    dWeight = fh2Weights[species]->GetBinContent(iBin);
  }

  if(dWeight <= 0.0) { dWeight = 1.0; }
  return dWeight;
}

const char* AliAnalysisTaskCorrForNonlinearFlow::GetSpeciesName(const PartSpecies species) const {
  const char* name;

  switch(species) {
    case kRefs: name = "Refs"; break;
    case kCharged: name = "Charged"; break;
    case kPion: name = "Pion"; break;
    case kKaon: name = "Kaon"; break;
    case kProton: name = "Proton"; break;
    case kCharUnidentified: name = "UnidentifiedCharged"; break;
    case kK0s: name = "K0s"; break;
    case kLambda: name = "Lambda"; break;
    case kPhi: name = "Phi"; break;
    default: name = "Unknown";
  }

  return name;
}

void AliAnalysisTaskCorrForNonlinearFlow::NTracksCalculation(AliVEvent* aod) {
  const int nAODTracks = aod->GetNumberOfTracks();
  NtrksCounter = 0;
  NTracksCorrected = 0;
  NTracksUncorrected = 0;

  //..for DCA
  double pos[3], vz, vx, vy;
  vz = aod->GetPrimaryVertex()->GetZ();
  vx = aod->GetPrimaryVertex()->GetX();
  vy = aod->GetPrimaryVertex()->GetY();
  double vtxp[3] = {vx, vy, vz};
  float fVtxZ = vz;
  double runNumber = fInputEvent->GetRunNumber();

  //..LOOP OVER TRACKS........
  //........................................
  for(Int_t nt = 0; nt < nAODTracks; nt++)
  {
    AliAODTrack *aodTrk = (AliAODTrack*) fInputEvent->GetTrack(nt);

    if (!aodTrk) {
      continue;
    }

    aodTrk->GetXYZ(pos);
    if (!AcceptAODTrack(aodTrk, pos, vtxp)) continue;
    if(TMath::Abs(aodTrk->Eta()) > fEtaCut) continue;

    /*
    //..get phi-weight for NUA correction
    double weight = 1;
    if (fPeriod.EqualTo("LHC15oKatarina") ) {
      if(fNUA == 1) weight = GetWeightKatarina(aodTrk->Phi(), aodTrk->Eta(), fVtxZ);
    } else {
      if(fNUA == 1) weight = GetFlowWeightSystematics(aodTrk, fVtxZ, kRefs);
    }
    */
    double weightPt = 1;
    if (fPeriod.EqualTo("LHC15oKatarina") ) {
      if(fNUE == 1) weightPt = GetPtWeightKatarina(aodTrk->Pt(), aodTrk->Eta(), fVtxZ);
    } else {
      // This is to make sure the extreme efficiency is only applied when calculating the Qs.
      double extremeEfficiency = fExtremeEfficiency;
      fExtremeEfficiency = 0;
      if(fNUE == 1) weightPt = GetPtWeight(aodTrk->Pt(), aodTrk->Eta(), fVtxZ, runNumber);
      fExtremeEfficiency = extremeEfficiency;
    }

    NTracksUncorrected += 1;
    NTracksCorrected += weightPt;
  } // end loop of all track
  if (!fUseCorrectedNTracks) {
    NtrksCounter = NTracksUncorrected;
  } else {
    NtrksCounter = NTracksCorrected; 
  }
  // hTracksCorrection2d->Fill(NTracksUncorrected, NTracksCorrected);
  // hnCorrectedTracks->Fill(NtrksCounter, NTracksCorrected);
}

