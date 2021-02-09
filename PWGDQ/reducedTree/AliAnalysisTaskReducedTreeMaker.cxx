/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//    Analysis task for creating a reduced data tree                     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include <TChain.h>
#include <TH1D.h>
#include <TH2I.h>
#include <TFile.h>
#include <TBits.h>
#include <TRandom.h>
#include <TTimeStamp.h>

#include <AliAnalysisTaskSE.h>
#include <AliCFContainer.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliESDHeader.h>
#include <AliAODHeader.h>
#include <AliAODTrack.h>
#include <AliAODTrdTrack.h>
#include <AliAODForwardMult.h>
#include <AliForwardUtil.h>
#include <AliTriggerAnalysis.h>
#include <AliESDtrack.h>
#include <AliESDTrdTrack.h>
#include <AliESDtrackCuts.h>
#include <AliVZDC.h>
#include <AliESDv0.h>
#include <AliAODv0.h>
#include <AliESDv0Cuts.h>
#include <AliAODv0KineCuts.h>
#include <AliESDv0KineCuts.h>
#include <AliESDFMD.h>
#include <AliVCluster.h>
#include <AliAODTracklets.h>
#include <AliMultiplicity.h>
#include <AliAODTracklets.h>
#include <AliPIDResponse.h>
#include <AliTPCdEdxInfo.h>
#include <AliMCParticle.h>
#include <AliAODMCParticle.h>
#include "AliAnalysisUtils.h"
#include <AliMultSelection.h>
#include <AliMultEstimator.h>
#include <AliCentrality.h>
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTimeRangeCut.h"
#include "AliDielectronVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliReducedCaloClusterInfo.h"
#include "AliReducedFMDInfo.h"
#include "AliReducedEventPlaneInfo.h"
#include "AliSignalMC.h"
#include "AliDielectronCutGroup.h"
#include "AliDielectronVarCuts.h"
#include "AliAnalysisTaskReducedTreeMaker.h"

#include <iostream>
#include <vector>
#include <algorithm>
using std::cout;
using std::endl;
using std::flush;

ClassImp(AliAnalysisTaskReducedTreeMaker)

//_________________________________________________________________________________
AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker() :
  AliAnalysisTaskSE(),
  fAnalysisUtils(0x0),
  fUseAnalysisUtils(kFALSE),
  fMinVtxContributors(0),
  fMaxVtxZ(100.),
  fCutOnSPDVtxZ(kFALSE),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kAny),
  fRejectPileup(kFALSE),
  fTreeWritingOption(kBaseEventsWithBaseTracks),
  fWriteTree(kTRUE),
  fScaleDownEvents(0.0),
  fWriteSecondTrackArray(kFALSE),
  fWriteBaseTrack(),
  fMinSelectedTracks(),
  fMaxSelectedTracks(),
  fNSelectedFullTracks(),
  fNSelectedBaseTracks(),  
  fEventsList(0x0),
  fEventsHistogram(0x0),
  fTRDEventsHistogram(0x0),
  fEMCalEventsHistogram(0x0),
  fCentEventsList(0x0),
  fTracksHistogram(0x0),
  fMCSignalsHistogram(0x0),
  fFillTrackInfo(kTRUE),
  fFillV0Info(kTRUE),
  fFillGammaConversions(kTRUE),
  fFillK0s(kTRUE),
  fFillLambda(kTRUE),
  fFillALambda(kTRUE),
  fFillCaloClusterInfo(kTRUE),
  fFillFMDInfo(kFALSE),
  fFillEventPlaneInfo(kFALSE),
  fEventPlaneTPCetaGap(1.0),
  fFillMCInfo(kFALSE),
  fFillHFInfo(kFALSE),
  fMCsignals(),
  fMCsignalsWritingOptions(),
  fFillTRDMatchedTracks(kFALSE),
  fFillAllTRDMatchedTracks(kFALSE),
  fTRDtrglayerMaskEl(0x1),
  fEventFilter(0x0),
  fTrackFilter(),
  fFlowTrackFilter(0x0),
  fClusterFilter(0x0),
  fK0sCuts(0x0),
  fLambdaCuts(0x0),
  fGammaConvCuts(0x0),
  fK0sPionCuts(0x0),
  fLambdaProtonCuts(0x0),
  fLambdaPionCuts(0x0),
  fGammaElectronCuts(0x0),
  fV0OpenCuts(0x0),
  fV0StrongCuts(0x0),
  fV0CutsAOD(0x0),
  fFMDhist(0x0),
  fK0sMassRange(),
  fLambdaMassRange(),
  fGammaMassRange(),
  fActiveBranches(""),
  fInactiveBranches(""),
  fTreeFile(0x0),
  fTree(0x0),
  fNevents(0),
  fReducedEvent(0x0),
  fUsedVars(0x0),
  fTimeRangeCut(),
  fTimeRangeReject(kFALSE)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker(const char *name, Bool_t writeTree /*=kTRUE*/) :
  AliAnalysisTaskSE(name),
  fAnalysisUtils(0x0),
  fUseAnalysisUtils(kFALSE),
  fMinVtxContributors(0),
  fMaxVtxZ(100.),
  fCutOnSPDVtxZ(kFALSE),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kAny),
  fRejectPileup(kFALSE),
  fTreeWritingOption(kBaseEventsWithBaseTracks),
  fWriteTree(kTRUE),
  fScaleDownEvents(0.0),
  fWriteSecondTrackArray(kFALSE),
  fWriteBaseTrack(),
  fMinSelectedTracks(),
  fMaxSelectedTracks(),
  fNSelectedFullTracks(),
  fNSelectedBaseTracks(),  
  fEventsList(0x0),
  fEventsHistogram(0x0),
  fTRDEventsHistogram(0x0),
  fEMCalEventsHistogram(0x0),
  fCentEventsList(0x0),
  fTracksHistogram(0x0),
  fMCSignalsHistogram(0x0),
  fFillTrackInfo(kTRUE),
  fFillV0Info(kTRUE),
  fFillGammaConversions(kTRUE),
  fFillK0s(kTRUE),
  fFillLambda(kTRUE),
  fFillALambda(kTRUE),
  fFillCaloClusterInfo(kTRUE),
  fFillFMDInfo(kFALSE),
  fFillEventPlaneInfo(kFALSE),
  fEventPlaneTPCetaGap(1.0),
  fFillMCInfo(kFALSE),
  fFillHFInfo(kFALSE),
  fMCsignals(),
  fMCsignalsWritingOptions(),
  fFillTRDMatchedTracks(kFALSE),
  fFillAllTRDMatchedTracks(kFALSE),
  fTRDtrglayerMaskEl(0x1),
  fEventFilter(0x0),
  fTrackFilter(),
  fFlowTrackFilter(0x0),
  fClusterFilter(0x0),
  fK0sCuts(0x0),
  fLambdaCuts(0x0),
  fGammaConvCuts(0x0),
  fK0sPionCuts(0x0),
  fLambdaProtonCuts(0x0),
  fLambdaPionCuts(0x0),
  fGammaElectronCuts(0x0),
  fV0OpenCuts(0x0),
  fV0StrongCuts(0x0),
  fV0CutsAOD(0x0),
  fFMDhist(0x0),
  fK0sMassRange(),
  fLambdaMassRange(),
  fGammaMassRange(),
  fActiveBranches(""),
  fInactiveBranches(""),
  fTreeFile(0x0),
  fTree(0x0),
  fNevents(0),
  fReducedEvent(0x0),
  fUsedVars(0x0),
  fTimeRangeCut(),
  fTimeRangeReject(kFALSE)
{
  //
  // Constructor
  //
  fK0sMassRange[0] = 0.4; fK0sMassRange[1] = 0.6;
  fLambdaMassRange[0] = 1.08; fLambdaMassRange[1] = 1.15;
  fGammaMassRange[0] = 0.0; fGammaMassRange[1] = 0.1;
  for(Int_t i=0; i<kMaxMCsignals; ++i) fMCsignalsWritingOptions[i] = kBaseTrack;
 
  DefineInput(0, TChain::Class());
  //DefineInput(2,AliAODForwardMult::Class());
  DefineOutput(1, AliReducedBaseEvent::Class());   // reduced information tree
  if(writeTree) {
    DefineOutput(2, TTree::Class());  // reduced information tree
    DefineOutput(3, TList::Class());  // Tlist of event statistics information
  }
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //
  if(fUseAnalysisUtils) fAnalysisUtils = new AliAnalysisUtils();
  if(fTree) return; //already initialised
  
  if(fWriteTree) {
    OpenFile(2);
    fTree = new TTree("DstTree","Reduced ESD information");
  }
  
  // check for tension between fTreeWritingOption and individual choices from AddTrackFilter
  if(fTreeWritingOption==kBaseEventsWithBaseTracks || fTreeWritingOption==kFullEventsWithBaseTracks) {
    for(UInt_t i=0; i<fWriteBaseTrack.size(); i++) {
      if(!fWriteBaseTrack.at(i)) {
        printf("AliAnalysisTaskReducedTreeMaker::UserCreateOutputObjects() WARNING: Full tracks requested for filter %d, but interferes with fTreeWritingOption choice! Only base tracks will be written. \n", i);
        fWriteBaseTrack.at(i) = kTRUE;
      }
    }
  }

  // check for tension between fTreeWritingOption and individual choices for MC signals
  if(fTreeWritingOption==kBaseEventsWithBaseTracks || fTreeWritingOption==kFullEventsWithBaseTracks) {
    for(Int_t i=0; i<kMaxMCsignals; i++) {
      if(fMCsignalsWritingOptions[i]==kFullTrack) {
        printf("AliAnalysisTaskReducedTreeMaker::UserCreateOutputObjects() WARNING: Full tracks requested for MC signal %d, but interferes with fTreeWritingOption choice! Only base tracks will be written. \n", i);
        fMCsignalsWritingOptions[i] = kBaseTrack;
      }
    }
  }

  // print active filters
  for(Int_t i=0; i<fTrackFilter.GetEntries(); i++) {
    cout << "AliAnalysisTaskReducedTreeMaker::UserCreateOutputObjects() filter " << i << ", base track = " << fWriteBaseTrack.at(i);
    cout << ", min tracks: " << fMinSelectedTracks[i] << ", max tracks: " << fMaxSelectedTracks[i] << endl;
  }

  // check if second track array is needed, i.e. fTracks contains full tracks, fTracks2 contains base tracks
  if(fTreeWritingOption==kBaseEventsWithFullTracks || fTreeWritingOption==kFullEventsWithFullTracks) {
    // data
    if(std::find(fWriteBaseTrack.begin(), fWriteBaseTrack.end(), kTRUE) != fWriteBaseTrack.end()) {
      printf("AliAnalysisTaskReducedTreeMaker::UserCreateOutputObjects(): Second track array will be used.\n");
      fWriteSecondTrackArray = kTRUE;
    }

    // MC
    if(fWriteSecondTrackArray==kFALSE && fMCsignals.GetEntries()) {
      for(Int_t i=0; i<fMCsignals.GetEntries(); i++) {
        if(fMCsignalsWritingOptions[i]==kBaseTrack) {
          printf("AliAnalysisTaskReducedTreeMaker::UserCreateOutputObjects(): Second track array will be used (MC signal choice).\n");
          fWriteSecondTrackArray = kTRUE;
          break;
        }
      }
    }
  }

  Int_t track2Option = AliReducedBaseEvent::kNoInit;
  if(fWriteSecondTrackArray) track2Option = AliReducedBaseEvent::kUseBaseTracks;

  switch(fTreeWritingOption) {
     case kBaseEventsWithBaseTracks:
        fReducedEvent = new AliReducedBaseEvent("DstEvent", AliReducedBaseEvent::kUseBaseTracks, track2Option);
        break;
     case kBaseEventsWithFullTracks:
        fReducedEvent = new AliReducedBaseEvent("DstEvent", AliReducedBaseEvent::kUseReducedTracks, track2Option);
        break;
     case kFullEventsWithBaseTracks:
        fReducedEvent = new AliReducedEventInfo("DstEvent", AliReducedBaseEvent::kUseBaseTracks, track2Option);
        break;
     case kFullEventsWithFullTracks:
        fReducedEvent = new AliReducedEventInfo("DstEvent", AliReducedBaseEvent::kUseReducedTracks, track2Option);
        break;
     default:
        break;
  };

  if(fWriteTree) {
    fTree->Branch("Event",&fReducedEvent,16000,99);

    // if user set active branches
    TObjArray* aractive=fActiveBranches.Tokenize(";");
    if(aractive->GetEntries()>0) {fTree->SetBranchStatus("*", 0);}
    for(Int_t i=0; i<aractive->GetEntries(); i++){
	  fTree->SetBranchStatus(aractive->At(i)->GetName(), 1);
    }

    // if user set inactive branches
    TObjArray* arinactive=fInactiveBranches.Tokenize(";");
    for(Int_t i=0; i<arinactive->GetEntries(); i++){
	  fTree->SetBranchStatus(arinactive->At(i)->GetName(), 0);
    }

    // if MC info is not requested, then set the respective branches off
    if(!fFillMCInfo) fTree->SetBranchStatus("fTracks.fMC*", 0);
    
    // switch event plane information branch
    if(!fFillEventPlaneInfo) fTree->SetBranchStatus("fEventPlane.*", 0);
    
    // if calorimeter cluster is filled, switch on cluster ID branch
    if(fFillCaloClusterInfo) fTree->SetBranchStatus("fCaloClusters.fClusterID", 1);
  }
  
  // enable all variables in the VarManager
  fUsedVars = new TBits(AliDielectronVarManager::kNMaxValues);
  for(Int_t i=0;i<AliDielectronVarManager::kParticleMax;++i) fUsedVars->SetBitNumber(i,kTRUE);
  
  AliDielectronVarManager::SetFillMap(fUsedVars);

  // TList for Event statistic histograms
  fEventsList = new TList();
  fEventsList->SetOwner();

  // event statistics histogram
  fEventsHistogram = new TH2I("EventStatistics", "Event statistics", 13, -0.5,12.5,34,-2.5,31.5);
  const Char_t* offlineTriggerNames[34] = {"Total", "No Phys Sel", "MB/INT1", "INT7", "MUON", "HighMult/HighMultSPD", "EMC1", "CINT5/INT5", "CMUS5/MUSPB/INT7inMUON",
     "MuonSingleHighPt7/MUSH7/MUSHPB", "MuonLikeLowPt7/MUL7/MuonLikePB", "MuonUnlikeLowPt7/MUU7/MuonUnlikePB", "EMC7/EMC8", 
     "MUS7/MuonSingleLowPt7", "PHI1", "PHI7/PHI8/PHOSPb", "EMCEJE", "EMCEGA", "Central/HighMultV0", "SemiCentral", "DG/DG5", "ZED", 
     "SPI7/SPI", "INT8", "MuonSingleLowPt8", "MuonSingleHighPt8", "MuonLikeLowPt8", "MuonUnlikeLowPt8", "MuonUnlikeLowPt0/INT6", "UserDefined", 
     "TRD", "MuonCalo/CaloOnly", "FastOnly", "N/A"
  };  
  const Char_t* selectionNames[13] = {"All events", 
     "TR and Physics Selection events (PS)", "Rejected due to PS",
     "PS and Trigger Selected (TS)", "Rejected due to TS",
     "TS and Pileup Checked (PC)", "Rejected due to PC",
     "PC and Event cuts Checked (EC)", "Rejected due to EC",
     "EC and Time Range accepted events (TR)", "Rejected due to TR",
     "Written ev. (track filters passed)", "Written ev. (unbiased)"};
  for(Int_t i=1;i<=34;++i)
     fEventsHistogram->GetYaxis()->SetBinLabel(i, offlineTriggerNames[i-1]);
  for(Int_t i=1;i<=13;++i)
     fEventsHistogram->GetXaxis()->SetBinLabel(i, selectionNames[i-1]);
  fEventsList->Add(fEventsHistogram);

  // TRD event statistics histogram
  const Char_t* offlineTRDTriggerNames[7] = {"Total", "No Phys Sel", "HQU | HSE", "HQU", "HSE", "Nuclei", "Jet"};
  fTRDEventsHistogram = new TH2I("TRDEventStatistics", "TRD Event statistics", 13, -0.5,12.5,12,-2.5,9.5);
  for(Int_t i=1;i<=7;++i)  fTRDEventsHistogram->GetYaxis()->SetBinLabel(i, offlineTRDTriggerNames[i-1]);
  for(Int_t i=1;i<=13;++i) fTRDEventsHistogram->GetXaxis()->SetBinLabel(i, selectionNames[i-1]);
	fEventsList->Add(fTRDEventsHistogram);

  // EMCal event statistics histogram
  const Char_t* offlineEMCalTriggerNames[22] = {"Total", "No Phys Sel", 
    "EMC7 | DMC7", "EMC7", "DMC7", "EMC1 | DMC1", "EMC1", "DMC1", "EMCEGA", "EG1 | DG1", "EG2 | DG2", "EG1", 
    "EG2", "DG1", "DG2", "EMCEJE", "EJ1 | DJ1", "EJ2 | DJ2", "EJ1", "EJ2", "DJ1", "DJ2"};
  fEMCalEventsHistogram = new TH2I("EMCalEventStatistics", "EMCal Event statistics", 13,-0.5,12.5, 22,-2.5,19.5);
  for(Int_t i=1;i<=22;++i) fEMCalEventsHistogram->GetYaxis()->SetBinLabel(i, offlineEMCalTriggerNames[i-1]);
  for(Int_t i=1;i<=13;++i) fEMCalEventsHistogram->GetXaxis()->SetBinLabel(i, selectionNames[i-1]);
  fEventsList->Add(fEMCalEventsHistogram);

  // Event Centrality statistics list of histograms (for the provided trigger mask)
  fCentEventsList = new TList();
  fCentEventsList->SetName("EventStatVsCent");
  fCentEventsList->SetOwner();
  const Int_t nCentEstimators = 3;
  const Char_t* estimatorNames[nCentEstimators] = {"V0M", "ZNA", "CL1"};
  TH2I *centEventsHistogram[nCentEstimators];
  for(Int_t i=0; i<nCentEstimators; ++i) {
    centEventsHistogram[i] = new TH2I(estimatorNames[i], Form("%s estimator",estimatorNames[i]), 13, -0.5,12.5,120,-5.,115.);
    for(Int_t xi=1;xi<=13;++xi) 
      centEventsHistogram[i]->GetXaxis()->SetBinLabel(xi, selectionNames[xi-1]);
    fCentEventsList->Add(centEventsHistogram[i]);
  }
  fEventsList->Add(fCentEventsList);

  // track statistics histogram
  fTracksHistogram = new TH2I("TrackStatistics", "Track statistics", fTrackFilter.GetEntries()+15, -1.5, fTrackFilter.GetEntries()+13.5, 3, -0.5, 2.5);
  fTracksHistogram->GetYaxis()->SetBinLabel(3, "base tracks");
  fTracksHistogram->GetYaxis()->SetBinLabel(2, "full tracks");
  fTracksHistogram->GetYaxis()->SetBinLabel(1, "total");
  fTracksHistogram->GetXaxis()->SetBinLabel(1, "total tracks");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+2, "Pure MC tracks");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+3, "Conversion electrons");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+4, "K^{0}_{s} pions");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+5, "#Lambda protons");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+6, "#Lambda pions");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+7, "total V0 pairs");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+8, "on the fly #gamma");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+9, "on the fly K^{0}_{s}");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+10, "on the fly #Lambda");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+11, "on the fly #bar{#Lambda}");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+12, "offline #gamma");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+13, "offline K^{0}_{s}");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+14, "offline #Lambda");
  fTracksHistogram->GetXaxis()->SetBinLabel(fTrackFilter.GetEntries()+15, "offline #bar{#Lambda}");
  for(Int_t i=0; i<fTrackFilter.GetEntries(); i++)
    fTracksHistogram->GetXaxis()->SetBinLabel(i+2, Form("%s", ((AliAnalysisCuts*)fTrackFilter.At(i))->GetName()));
  fEventsList->Add(fTracksHistogram);
  // Add counters for the V0 prongs (conversion electrons, K0s pions, Lambda protons, Lambda pions)
  //        and V0 pairs (separately on the fly and offline)
  // TODO: indexing has to be done in a less error prone way (define an enum)
  for(Int_t i=0; i<4; i++) {
    fNSelectedFullTracks.push_back(0);
    fNSelectedBaseTracks.push_back(0);
  }
  for(Int_t i=0; i<8; i++) fNSelectedFullTracks.push_back(0);
  // add counters for the MC tracks (full and base tracks)
  fNSelectedFullTracks.push_back(0);
  fNSelectedBaseTracks.push_back(0);
  
  // MC statistics histogram
  fMCSignalsHistogram = new TH2I("MCSignalsStatistics", "Monte-Carlo signals statistics", 
                                 fMCsignals.GetEntries(), -0.5, Double_t(fMCsignals.GetEntries())-0.5, 32, -0.5, 31.5);
  for(Int_t i=1;i<=32;++i) fMCSignalsHistogram->GetYaxis()->SetBinLabel(i, offlineTriggerNames[i-1]);
  for(Int_t i=1;i<=fMCsignals.GetEntries();++i) {
     TString trackTypeStr = "base track";
     if(fMCsignalsWritingOptions[i-1]==kFullTrack) trackTypeStr = "full track";
     fMCSignalsHistogram->GetXaxis()->SetBinLabel(i, Form("%s (%s)", ((AliSignalMC*)fMCsignals.At(i-1))->GetName(), trackTypeStr.Data()));
  }
  if(fFillMCInfo) fEventsList->Add(fMCSignalsHistogram);

  // set a seed for the random number generator
  TTimeStamp ts;
  gRandom->SetSeed(ts.GetNanoSec());
  
  PostData(1, fReducedEvent);
  if(fWriteTree) {
    PostData(2, fTree);
    PostData(3, fEventsList);
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::UserExec(Option_t *option)
{
  //
  // Main loop. Called for every event
  //
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  fNevents++;

  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(!inputHandler) return;

  if(inputHandler->GetPIDResponse()) 
    AliDielectronVarManager::SetPIDResponse(inputHandler->GetPIDResponse());
  else
    AliFatal("This task needs the PID response attached to the input event handler!");
  
  // Was event selected ?
  UInt_t isPhysSel = AliVEvent::kAny;
  UInt_t isPhysAndTrigSel = AliVEvent::kAny;
  UChar_t trdtrgtype=0;
  UInt_t emcaltrgtype=0;
  
  if((isESD && inputHandler->GetEventSelection()) || isAOD){
    isPhysSel = inputHandler->IsEventSelected();
    isPhysAndTrigSel = isPhysSel & fTriggerMask;

    // get type of TRD triggered event
    TString trgClasses = inputHandler->GetEvent()->GetFiredTriggerClasses();
    if(trgClasses.Contains("HQU") || trgClasses.Contains("HSE")) trdtrgtype |= (UChar_t(1)<<0);
    if(trgClasses.Contains("HQU")) trdtrgtype |= (UChar_t(1)<<1);
    if(trgClasses.Contains("HSE")) trdtrgtype |= (UChar_t(1)<<2);
    if(trgClasses.Contains("HNU")) trdtrgtype |= (UChar_t(1)<<3);
    if(trgClasses.Contains("HJT")) trdtrgtype |= (UChar_t(1)<<4);

    // get type of EMCal triggered event
    if(trgClasses.Contains("EMC7") || trgClasses.Contains("DMC7")) emcaltrgtype |= (UInt_t(1)<<0);
    if(trgClasses.Contains("EMC7")) emcaltrgtype |= (UInt_t(1)<<1);
    if(trgClasses.Contains("DMC7")) emcaltrgtype |= (UInt_t(1)<<2);
    if(trgClasses.Contains("EMC1") || trgClasses.Contains("DMC1")) emcaltrgtype |= (UInt_t(1)<<3);
    if(trgClasses.Contains("EMC1")) emcaltrgtype |= (UInt_t(1)<<4);
    if(trgClasses.Contains("DMC1")) emcaltrgtype |= (UInt_t(1)<<5);
    if(trgClasses.Contains("EG1") || trgClasses.Contains("EG2") || 
       trgClasses.Contains("DG1") || trgClasses.Contains("DG2")) emcaltrgtype |= (UInt_t(1)<<6);
    if(trgClasses.Contains("EG1") || trgClasses.Contains("DG1")) emcaltrgtype |= (UInt_t(1)<<7);
    if(trgClasses.Contains("EG2") || trgClasses.Contains("DG2")) emcaltrgtype |= (UInt_t(1)<<8);
    if(trgClasses.Contains("EG1")) emcaltrgtype |= (UInt_t(1)<<9);
    if(trgClasses.Contains("EG2")) emcaltrgtype |= (UInt_t(1)<<10);
    if(trgClasses.Contains("DG1")) emcaltrgtype |= (UInt_t(1)<<11);
    if(trgClasses.Contains("DG2")) emcaltrgtype |= (UInt_t(1)<<12);
    if(trgClasses.Contains("EJ1") || trgClasses.Contains("EJ2") || 
       trgClasses.Contains("DJ1") || trgClasses.Contains("DJ2")) emcaltrgtype |= (UInt_t(1)<<13);
    if(trgClasses.Contains("EJ1") || trgClasses.Contains("DJ1")) emcaltrgtype |= (UInt_t(1)<<14);
    if(trgClasses.Contains("EJ2") || trgClasses.Contains("DJ2")) emcaltrgtype |= (UInt_t(1)<<15);
    if(trgClasses.Contains("EJ1")) emcaltrgtype |= (UInt_t(1)<<16);
    if(trgClasses.Contains("EJ2")) emcaltrgtype |= (UInt_t(1)<<17);
    if(trgClasses.Contains("DJ1")) emcaltrgtype |= (UInt_t(1)<<18);
    if(trgClasses.Contains("DJ2")) emcaltrgtype |= (UInt_t(1)<<19);
  }
  
  // Get centrality object
  const Int_t nCentEstimators = 3;
  const Char_t* estimatorNames[nCentEstimators] = {"V0M", "ZNA", "CL1"};
  Double_t percentileEstimators[nCentEstimators] = {999., 999., 999.};
  AliVEvent* event = InputEvent();
  AliCentrality *centrality = 0x0;
  AliMultSelection* multSelection = 0x0;
  Bool_t isOldCent = kFALSE;
  if(event->GetRunNumber()<200000) isOldCent = kTRUE;
  if(isOldCent) centrality = event->GetCentrality(); // old centrality framework
  else multSelection = (AliMultSelection*)event->FindListObject("MultSelection"); // new centrality framework
  if(!centrality && !multSelection) AliInfo("No centrality object found");
  for(Int_t i=0; i<nCentEstimators; ++i) 
    percentileEstimators[i] = (isOldCent ? centrality->GetCentralityPercentile(Form("%s",estimatorNames[i])) : 
                                           multSelection->GetMultiplicityPercentile(Form("%s",estimatorNames[i])));
        
  // event statistics before any selection
  FillStatisticsHistograms(Bool_t(isPhysSel), isPhysSel, trdtrgtype,  emcaltrgtype, 0., percentileEstimators,  nCentEstimators);
  
  // rejected due to physics selection
  if(fSelectPhysics && !isPhysSel) {
    fEventsHistogram->Fill(2., -1.);
    fEventsHistogram->Fill(2., -2.);
    fTRDEventsHistogram->Fill(2., -1.);
    fTRDEventsHistogram->Fill(2., -2.);
    fEMCalEventsHistogram->Fill(2., -1.);
    fEMCalEventsHistogram->Fill(2., -2.);
	for(Int_t i=0; i<nCentEstimators; ++i)
      ((TH2I*)fCentEventsList->At(i))->Fill(2., percentileEstimators[i]);
    PostData(3, fEventsList);
    return;
  }
  
  // event statistics after physics selection
  // NOTE: if physics selection was not applied (as requested by user) then we can still have events with PS not fulfilled
  FillStatisticsHistograms(Bool_t(isPhysSel), isPhysSel, trdtrgtype,  emcaltrgtype, 1., percentileEstimators,  nCentEstimators);
  
  // event statistics after physics selection and trigger selection
  if(isPhysAndTrigSel)
    FillStatisticsHistograms(kTRUE, isPhysSel, trdtrgtype,  emcaltrgtype, 3., percentileEstimators,  nCentEstimators);
  else {
    // reject events which do not fulfill the requested trigger mask
    FillStatisticsHistograms(Bool_t(isPhysSel), isPhysSel, trdtrgtype,  emcaltrgtype, 4., percentileEstimators,  nCentEstimators);
    PostData(3, fEventsList);
    return;
  }

  // rejected by pileup selection
  if(fRejectPileup && InputEvent()->IsPileupFromSPD(3,0.8,3.,2.,5.)) {
    FillStatisticsHistograms(Bool_t(isPhysSel), isPhysSel, trdtrgtype,  emcaltrgtype, 6., percentileEstimators,  nCentEstimators);
    PostData(3, fEventsList);
    return;
  }
  else {
    // accepted pileup selection
    FillStatisticsHistograms(Bool_t(isPhysSel), isPhysSel, trdtrgtype,  emcaltrgtype, 5., percentileEstimators,  nCentEstimators);
  }  
  
  // user defined event filter
  if(fEventFilter && !fEventFilter->IsSelected(InputEvent())) {
    // event statistics for events failing selection cuts
    FillStatisticsHistograms(Bool_t(isPhysSel), isPhysSel, trdtrgtype,  emcaltrgtype, 8., percentileEstimators,  nCentEstimators);
    PostData(3, fEventsList);
    return;
  }
  else {
    FillStatisticsHistograms(Bool_t(isPhysSel), isPhysSel, trdtrgtype,  emcaltrgtype, 7., percentileEstimators,  nCentEstimators);
  }
  
  // Apply the time range cut needed to reject events with bad TPC pid in some chambers
  // For details see: https://indico.cern.ch/event/830757/contributions/3479738/attachments/1873483/3083952/IArsene_DPG_2019July3.pdf
  fTimeRangeCut.InitFromEvent(InputEvent());
  if(fTimeRangeCut.CutEvent(InputEvent())) {
     if(fTimeRangeReject) {
        // if rejecting the bad events from writing, fill the appropriate event stats bins and return 
        FillStatisticsHistograms(Bool_t(isPhysSel), isPhysSel, trdtrgtype,  emcaltrgtype, 10., percentileEstimators,  nCentEstimators);
        PostData(3, fEventsList);
        return;
     }
     fReducedEvent->fEventTag |= (ULong64_t(1)<<15);
  }
  // fill the event statistics for events passing the Time Range cut
  FillStatisticsHistograms(Bool_t(isPhysSel), isPhysSel, trdtrgtype,  emcaltrgtype, 9., percentileEstimators,  nCentEstimators);
  
  if(fFillMCInfo) {
     Bool_t hasMC=AliDielectronMC::Instance()->HasMC();
     if(hasMC) {
       AliDielectronMC::Instance()->SetCheckHF(fFillHFInfo);
       AliDielectronMC::Instance()->ConnectMCEvent();
       AliDielectronVarManager::SetEvent(AliDielectronMC::Instance()->GetMCEvent());
     }
  }
  AliDielectronVarManager::SetEvent(InputEvent());

  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField(bz);
  
  //Fill event wise information
  fReducedEvent->ClearEvent();  
  FillEventInfo();
  if(fFillCaloClusterInfo) FillCaloClusters();
  if(fFillFMDInfo) FillFMDInfo(isAOD);

  // reset track counters
  // TODO: indexing needs to be in a less error prone way (define an enum)
  for(Int_t i=0; i<fTrackFilter.GetEntries()+4; ++i) {
    fNSelectedFullTracks[i] = 0;
    fNSelectedBaseTracks[i] = 0;
  }
  // additional counters for the V0 pairs
  for(Int_t i=0; i<8; ++i) fNSelectedFullTracks[fTrackFilter.GetEntries()+4+i] = 0;
  // counters for the MC tracks
  fNSelectedFullTracks[fTrackFilter.GetEntries()+4+8] = 0;
  fNSelectedBaseTracks[fTrackFilter.GetEntries()+4] = 0;
  
  // NOTE: It is important that FillV0PairInfo() is called before FillTrackInfo()
  if(fFillMCInfo) FillMCTruthInfo();
  if(fFillV0Info && isESD) FillV0PairInfo();
  if(fFillV0Info && isAOD) FillV0PairInfoAOD();
  if(fFillTrackInfo) FillTrackInfo();
  
  if(fWriteTree) {
    // write a random sample of events with a downscale of fScaleDownEvents
    Bool_t unbiasedEvent = kFALSE;
    if(gRandom->Rndm()<fScaleDownEvents) {
      unbiasedEvent = kTRUE;
      fReducedEvent->fEventTag |= (ULong64_t(1)<<14);                    // mark unbiased events
      FillStatisticsHistograms(Bool_t(isPhysSel), isPhysSel, trdtrgtype, emcaltrgtype, 12.0, percentileEstimators, nCentEstimators);
    }
    
    // if the event was not already selected to be written, check that it fullfills the conditions
    Bool_t trackFilterAccepted = kTRUE;
    for(Int_t i=0; i<fTrackFilter.GetEntries(); ++i) {
      if(fMinSelectedTracks[i]>=0 && (fNSelectedFullTracks[i]+fNSelectedBaseTracks[i])<fMinSelectedTracks[i]) {
        trackFilterAccepted = kFALSE; break;
      }
      if(fMaxSelectedTracks[i]>=0 && (fNSelectedFullTracks[i]+fNSelectedBaseTracks[i])>fMaxSelectedTracks[i]) {
        trackFilterAccepted = kFALSE; break;
      }
    }
    if(trackFilterAccepted) 
      FillStatisticsHistograms(Bool_t(isPhysSel), isPhysSel, trdtrgtype, emcaltrgtype, 11., percentileEstimators, nCentEstimators);
    
    if(trackFilterAccepted || unbiasedEvent) {
      fTree->Fill();
      FillTrackStatisticsHistogram();
    }
  }  // end if(writeTree)

  PostData(1, fReducedEvent);
  if(fWriteTree) {
    PostData(2, fTree);
    PostData(3, fEventsList);
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillStatisticsHistograms(Bool_t isSelected, UInt_t physSel, 
                                    UChar_t trdTrigMap, UInt_t emcalTrigMap, Double_t xbin, Double_t* percentiles, Int_t nEstimators) {
   //
   // Fill event statistics histograms
   //
   if(isSelected) {             //  count the events that passed the selection 
      for(Int_t i=0; i<32; ++i) 
         if(physSel & (UInt_t(1)<<i)) fEventsHistogram->Fill(xbin, Double_t(i));
      for(Int_t i=0; i<5; ++i)
         if(trdTrigMap & (UChar_t(1) << i)) fTRDEventsHistogram->Fill(xbin, Double_t(i));
      for(Int_t i=0; i<20; ++i)
         if(emcalTrigMap & (UInt_t(1) << i)) fEMCalEventsHistogram->Fill(xbin, Double_t(i));
   }
   else {                                                   //  count the events which did not pass
      fEventsHistogram->Fill(xbin, -1.);
      fTRDEventsHistogram->Fill(xbin, -1.);
      fEMCalEventsHistogram->Fill(xbin, -1.);
   }
   // count all events,  passing or not the condition
   fEventsHistogram->Fill(xbin, -2.);
   fTRDEventsHistogram->Fill(xbin, -2.);
   fEMCalEventsHistogram->Fill(xbin, -2.);
  
   // count events as a function of centrality
   for(Int_t i=0; i<nEstimators; ++i) ((TH2I*)fCentEventsList->At(i))->Fill(xbin, percentiles[i]);
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::AddTrackFilter(AliAnalysisCuts* const filter, Bool_t option/*=kFALSE*/, Int_t minSel/*=-1*/, Int_t maxSel/*=-1*/)
{
  //
  // add track filter to track filter list
  //
  if(fTrackFilter.GetEntries()<kNMaxTrackFilters) {
    fTrackFilter.Add(filter);
    fWriteBaseTrack.push_back(option);
    fMinSelectedTracks.push_back(minSel);
    fMaxSelectedTracks.push_back(maxSel);
    fNSelectedFullTracks.push_back(0);
    fNSelectedBaseTracks.push_back(0);
  } else {
    printf("AliAnalysisTaskReducedTreeMaker::AddTrackFilter() WARNING: Track filter list full (%d entries), will not add another filter!\n", fTrackFilter.GetEntries());
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::AddCaloClusterFilter(AliAnalysisCuts * const filter)
{
  //
  // add cluster filter to cluster filter list
  //
  if(fClusterFilter.GetEntries()<32) fClusterFilter.Add(filter);
  else printf("AliAnalysisTaskReducedTreeMaker::AddCaloClusterFilter() WARNING: Cluster filter list full (%d entries), will not add another filter!\n", fClusterFilter.GetEntries());
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeMaker::IsTrackSelected(AliVParticle* track, Double_t* values, std::vector<Bool_t>& filterDecision)
{
  //
  // check if track is selected and write filter decision to vector
  //
  Bool_t trackIsSelected = kFALSE;
  for(Int_t i=0; i<fTrackFilter.GetEntries(); i++) {
    AliAnalysisCuts* filter = (AliAnalysisCuts*)fTrackFilter.At(i);
    Bool_t                                                  filterIsSelected = kFALSE;
    if(filter->IsA()==AliDielectronCutGroup::Class())       filterIsSelected = (dynamic_cast<AliDielectronCutGroup*>(filter))->IsSelected(track, values);
    else if (filter->IsA()==AliDielectronVarCuts::Class())  filterIsSelected = (dynamic_cast<AliDielectronVarCuts*>(filter))->IsSelected(values);
    else                                                    filterIsSelected = filter->IsSelected(track);
    if(filterIsSelected) {
      filterDecision.push_back(kTRUE);
      trackIsSelected = kTRUE;
    } else {
      filterDecision.push_back(kFALSE);
    }
  }
  return trackIsSelected;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeMaker::IsSelectedTrackRequestedBaseTrack(std::vector<Bool_t> filterDecision, Bool_t usedForV0Or)
{
  //
  // compare passed track filter and corresponding choice of base or full track
  // full track wins if there is some overlap
  //
  Bool_t isBaseTrack = kTRUE;
  for(Int_t i=0; i<fTrackFilter.GetEntries(); i++) {
    if(filterDecision[i] && !fWriteBaseTrack[i]) {
      isBaseTrack = kFALSE;
      break;
    }
  }
  if(isBaseTrack && usedForV0Or) {
    if(fTreeWritingOption==kBaseEventsWithFullTracks || fTreeWritingOption==kFullEventsWithFullTracks)
      isBaseTrack = kFALSE;
  }
  return isBaseTrack;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeMaker::IsClusterSelected(AliVCluster* cluster, std::vector<Bool_t>& filterDecision) {
  //
  // check if cluster is selected and write filter decision to vector
  //
  Bool_t clusterIsSelected = kFALSE;
  for(Int_t i=0; i<fClusterFilter.GetEntries(); i++) {
    AliAnalysisCuts* filter = (AliAnalysisCuts*)fClusterFilter.At(i);
    if(filter->IsSelected(cluster)) {
      filterDecision.push_back(kTRUE);
      clusterIsSelected = kTRUE;
    } else {
      filterDecision.push_back(kFALSE);
    }
  }
  return clusterIsSelected;
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::SetTrackFilterQualityFlags(AliReducedBaseTrack* track, std::vector<Bool_t> filterDecision)
{
  //
  // set track quality flags for passed track filters
  //
  for(Int_t i=0; i<fTrackFilter.GetEntries(); i++) {
    if(filterDecision[i])
      track->SetQualityFlag(32+i); // AliReduceBaseTrack::fQualityFlags BIT 32+i (0<=i<fTrackFilter.GetEntries())
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillTrackStatisticsHistogram()
{
  //
  // fill track statistics histogram
  //
  // Total number of tracks (1st bin)
  fTracksHistogram->Fill(-1.0, 0.0, fReducedEvent->fTracks->GetEntries()+fReducedEvent->fTracks2->GetEntries());
  if(fTreeWritingOption==kBaseEventsWithBaseTracks || fTreeWritingOption==kFullEventsWithBaseTracks)
    fTracksHistogram->Fill(-1.0, 2.0, fReducedEvent->fTracks->GetEntries());
  else {
    fTracksHistogram->Fill(-1.0, 1.0, fReducedEvent->fTracks->GetEntries());
    fTracksHistogram->Fill(-1.0, 2.0, fReducedEvent->fTracks2->GetEntries());
  }
  // counters based on track filters
  for(Int_t i=0; i<fTrackFilter.GetEntries(); i++) { 
    fTracksHistogram->Fill(i, 0.0, fNSelectedBaseTracks[i]+fNSelectedFullTracks[i]);
    fTracksHistogram->Fill(i, 1.0, fNSelectedFullTracks[i]);
    fTracksHistogram->Fill(i, 2.0, fNSelectedBaseTracks[i]);
  }
  // counters for tracks belonging to V0s (first 4 elements), on the fly V0 pairs (next 4 elements) and offline V0 pairs (next 4 elements)
  for(Int_t i=0; i<4; i++) {
    fTracksHistogram->Fill(fTrackFilter.GetEntries()+1+i, 0.0, fNSelectedFullTracks[fTrackFilter.GetEntries()+i]+fNSelectedBaseTracks[fTrackFilter.GetEntries()+i]);
    fTracksHistogram->Fill(fTrackFilter.GetEntries()+1+i, 1.0, fNSelectedFullTracks[fTrackFilter.GetEntries()+i]);
    fTracksHistogram->Fill(fTrackFilter.GetEntries()+1+i, 2.0, fNSelectedBaseTracks[fTrackFilter.GetEntries()+i]);
  }
  // counter for the total number of V0 pairs selected (also written if the fCandidates branch is not switched off)
  fTracksHistogram->Fill(fTrackFilter.GetEntries()+5, 0.0, fReducedEvent->fCandidates->GetEntries());
  for(Int_t i=0; i<8; i++)
    fTracksHistogram->Fill(fTrackFilter.GetEntries()+6+i, 0.0, fNSelectedFullTracks[fTrackFilter.GetEntries()+4+i]);
  // MC tracks counters
  fTracksHistogram->Fill(fTrackFilter.GetEntries(), 0.0, fNSelectedFullTracks[fTrackFilter.GetEntries()+12]+fNSelectedBaseTracks[fTrackFilter.GetEntries()+4]);
  fTracksHistogram->Fill(fTrackFilter.GetEntries(), 1.0, fNSelectedFullTracks[fTrackFilter.GetEntries()+12]);
  fTracksHistogram->Fill(fTrackFilter.GetEntries(), 2.0, fNSelectedBaseTracks[fTrackFilter.GetEntries()+4]);
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillEventInfo() 
{
  //
  // fill reduced event information
  //
  AliVEvent* event = InputEvent();
  // Was event selected ?
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD = (event->IsA()==AliESDEvent::Class());
  Bool_t isAOD = (event->IsA()==AliAODEvent::Class());
  
  AliESDEvent* esdEvent = 0x0;
  if(isESD) esdEvent = static_cast<AliESDEvent*>(event);
  AliAODEvent* aodEvent = 0x0;
  if(isAOD) aodEvent = static_cast<AliAODEvent*>(event);
  
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  UInt_t isSelected = AliVEvent::kAny;
  if(inputHandler){
    if((isESD && inputHandler->GetEventSelection()) || isAOD){
      isSelected = inputHandler->IsEventSelected();
      isSelected&=fTriggerMask;
    }
  }
  
  if(fUseAnalysisUtils) {
    if(fAnalysisUtils->IsVertexSelected2013pA(event))  // 2013 p-Pb event selection    
      fReducedEvent->fEventTag |= (ULong64_t(1)<<0);
    fAnalysisUtils->SetMinPlpContribMV(5); fAnalysisUtils->SetMaxPlpChi2MV(5.);
    fAnalysisUtils->SetCheckPlpFromDifferentBCMV(kTRUE);
    fAnalysisUtils->SetMinWDistMV(15.);
    if(fAnalysisUtils->IsPileUpMV(event))              // multi-vertexer pileup
      fReducedEvent->fEventTag |= (ULong64_t(1)<<1);
    fAnalysisUtils->SetCheckPlpFromDifferentBCMV(kFALSE);
    if(fAnalysisUtils->IsPileUpMV(event))              // multi-vertexer pileup, with no BC check
      fReducedEvent->fEventTag |= (ULong64_t(1)<<2);
    fAnalysisUtils->SetCheckPlpFromDifferentBCMV(kTRUE);
    fAnalysisUtils->SetMinWDistMV(10.);
    if(fAnalysisUtils->IsPileUpMV(event))              // multi-vertexer pileup, with min weighted distance 10
      fReducedEvent->fEventTag |= (ULong64_t(1)<<3);
    fAnalysisUtils->SetMinWDistMV(5.);
    if(fAnalysisUtils->IsPileUpMV(event))              // multi-vertexer pileup, with min weighted distance 5
      fReducedEvent->fEventTag |= (ULong64_t(1)<<4);
  }
    
  if(event->IsPileupFromSPD(3,0.6,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<5);
  if(event->IsPileupFromSPD(4,0.6,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<6);
  if(event->IsPileupFromSPD(5,0.6,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<7);
  if(event->IsPileupFromSPD(6,0.6,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<8);
  if(event->IsPileupFromSPD(3,0.8,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<9);
  if(event->IsPileupFromSPD(4,0.8,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<10);
  if(event->IsPileupFromSPD(5,0.8,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<11);
  if(event->IsPileupFromSPD(6,0.8,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<12);
    
  fReducedEvent->fRunNo       = event->GetRunNumber();
  AliVVertex* eventVtx = 0x0;
  if(isESD) eventVtx = const_cast<AliESDVertex*>(esdEvent->GetPrimaryVertexTracks());
  if(isAOD) eventVtx = const_cast<AliAODVertex*>(aodEvent->GetPrimaryVertex());
  if(eventVtx) {
     fReducedEvent->fVtx[0] = (isESD ? ((AliESDVertex*)eventVtx)->GetX() : ((AliAODVertex*)eventVtx)->GetX());
     fReducedEvent->fVtx[1] = (isESD ? ((AliESDVertex*)eventVtx)->GetY() : ((AliAODVertex*)eventVtx)->GetY());
     fReducedEvent->fVtx[2] = (isESD ? ((AliESDVertex*)eventVtx)->GetZ() : ((AliAODVertex*)eventVtx)->GetZ());
     fReducedEvent->fNVtxContributors = eventVtx->GetNContributors();
  }
  
  AliCentrality *centrality = 0x0;
  AliMultSelection* multSelection = 0x0;
  if(event->GetRunNumber()<200000) {
     centrality = event->GetCentrality();
     if(centrality) {
       fReducedEvent->fCentrality[0] = centrality->GetCentralityPercentile("V0M");
       fReducedEvent->fCentrality[1] = centrality->GetCentralityPercentile("CL1");
       fReducedEvent->fCentrality[2] = centrality->GetCentralityPercentile("TRK");
       fReducedEvent->fCentrality[3] = centrality->GetCentralityPercentile("ZEMvsZDC");
       fReducedEvent->fCentrality[4] = centrality->GetCentralityPercentile("V0A");
       fReducedEvent->fCentrality[5] = centrality->GetCentralityPercentile("V0C");
       fReducedEvent->fCentrality[6] = centrality->GetCentralityPercentile("ZNA");
       fReducedEvent->fCentQuality   = centrality->GetQuality();
     }
  }
  else {
     multSelection = (AliMultSelection*)event->FindListObject("MultSelection");
     if(multSelection) {
        fReducedEvent->fCentrality[0] = multSelection->GetMultiplicityPercentile("V0M");
        fReducedEvent->fCentrality[1] = multSelection->GetMultiplicityPercentile("CL1");
        fReducedEvent->fCentrality[2] = multSelection->GetMultiplicityPercentile("TRK");
        fReducedEvent->fCentrality[3] = multSelection->GetMultiplicityPercentile("ZEMvsZDC");
        fReducedEvent->fCentrality[4] = multSelection->GetMultiplicityPercentile("V0A");
        fReducedEvent->fCentrality[5] = multSelection->GetMultiplicityPercentile("V0C");
        fReducedEvent->fCentrality[6] = multSelection->GetMultiplicityPercentile("ZNA");
        fReducedEvent->fCentQuality   = multSelection->GetEvSelCode();
    }
  }
  fReducedEvent->fNtracks[0] = event->GetNumberOfTracks();
  
  // In case we want to write just basic event information, we stop here
  if(fTreeWritingOption==kBaseEventsWithBaseTracks || fTreeWritingOption==kBaseEventsWithFullTracks) 
     return;
  
  AliReducedEventInfo* eventInfo = dynamic_cast<AliReducedEventInfo*>(fReducedEvent);
  if(!eventInfo) return;
  
  if(multSelection) {
     eventInfo->fMultiplicityEstimatorPercentiles[0] = multSelection->GetMultiplicityPercentile("OnlineV0M");
     eventInfo->fMultiplicityEstimatorPercentiles[1] = multSelection->GetMultiplicityPercentile("OnlineV0A");
     eventInfo->fMultiplicityEstimatorPercentiles[2] = multSelection->GetMultiplicityPercentile("OnlineV0C");
     eventInfo->fMultiplicityEstimatorPercentiles[3] = multSelection->GetMultiplicityPercentile("ADM");
     eventInfo->fMultiplicityEstimatorPercentiles[4] = multSelection->GetMultiplicityPercentile("ADA");
     eventInfo->fMultiplicityEstimatorPercentiles[5] = multSelection->GetMultiplicityPercentile("ADC");
     eventInfo->fMultiplicityEstimatorPercentiles[6] = multSelection->GetMultiplicityPercentile("SPDClusters");
     eventInfo->fMultiplicityEstimatorPercentiles[7] = multSelection->GetMultiplicityPercentile("SPDTracklets");
     eventInfo->fMultiplicityEstimatorPercentiles[8] = multSelection->GetMultiplicityPercentile("RefMult05");
     eventInfo->fMultiplicityEstimatorPercentiles[9] = multSelection->GetMultiplicityPercentile("RefMult08");
     eventInfo->fMultiplicityEstimatorPercentiles[10] = multSelection->GetMultiplicityPercentile("V0M");
     eventInfo->fMultiplicityEstimatorPercentiles[11] = multSelection->GetMultiplicityPercentile("V0A");
     eventInfo->fMultiplicityEstimatorPercentiles[12] = multSelection->GetMultiplicityPercentile("V0C");
     AliMultEstimator* estimator = 0x0;
     estimator = multSelection->GetEstimator("OnlineV0M"); if(estimator) eventInfo->fMultiplicityEstimators[0] = estimator->GetValue();
     estimator = multSelection->GetEstimator("OnlineV0A"); if(estimator) eventInfo->fMultiplicityEstimators[1] = estimator->GetValue();
     estimator = multSelection->GetEstimator("OnlineV0C"); if(estimator) eventInfo->fMultiplicityEstimators[2] = estimator->GetValue();
     estimator = multSelection->GetEstimator("ADM"); if(estimator) eventInfo->fMultiplicityEstimators[3] = estimator->GetValue();
     estimator = multSelection->GetEstimator("ADA"); if(estimator) eventInfo->fMultiplicityEstimators[4] = estimator->GetValue();
     estimator = multSelection->GetEstimator("ADC"); if(estimator) eventInfo->fMultiplicityEstimators[5] = estimator->GetValue();
     estimator = multSelection->GetEstimator("SPDClusters"); if(estimator) eventInfo->fMultiplicityEstimators[6] = estimator->GetValue();
     estimator = multSelection->GetEstimator("SPDTracklets"); if(estimator) eventInfo->fMultiplicityEstimators[7] = estimator->GetValue();
     estimator = multSelection->GetEstimator("RefMult05"); if(estimator) eventInfo->fMultiplicityEstimators[8] = estimator->GetValue();
     estimator = multSelection->GetEstimator("RefMult08"); if(estimator) eventInfo->fMultiplicityEstimators[9] = estimator->GetValue();   
     estimator = multSelection->GetEstimator("V0M"); if(estimator) eventInfo->fMultiplicityEstimators[10] = estimator->GetValue();
     estimator = multSelection->GetEstimator("V0A"); if(estimator) eventInfo->fMultiplicityEstimators[11] = estimator->GetValue();
     estimator = multSelection->GetEstimator("V0C"); if(estimator) eventInfo->fMultiplicityEstimators[12] = estimator->GetValue();  
  }
  
  eventInfo->fSPDntracklets = GetSPDTrackletMultiplicity(event, -1.0, 1.0);
  for(Int_t ieta=0; ieta<32; ++ieta)
    eventInfo->fSPDntrackletsEta[ieta] = GetSPDTrackletMultiplicity(event, -1.6+0.1*ieta, -1.6+0.1*(ieta+1));
  
  if(eventVtx){
    Double_t covTracks[6];
    eventVtx->GetCovarianceMatrix(covTracks);
    for(Int_t i=0;i<6;++i) {
      eventInfo->fVtxCovMatrix[i] = covTracks[i];
    }
  }
  
  AliVVertex* eventVtxSPD = 0x0;
  if(isESD) eventVtxSPD = const_cast<AliESDVertex*>(esdEvent->GetPrimaryVertexSPD());
  if(isAOD) eventVtxSPD = const_cast<AliAODVertex*>(aodEvent->GetPrimaryVertexSPD());
  if(eventVtxSPD) {
     eventInfo->fVtxSPD[0] = (isESD ? ((AliESDVertex*)eventVtxSPD)->GetX() : ((AliAODVertex*)eventVtxSPD)->GetX());
     eventInfo->fVtxSPD[1] = (isESD ? ((AliESDVertex*)eventVtxSPD)->GetY() : ((AliAODVertex*)eventVtxSPD)->GetY());
     eventInfo->fVtxSPD[2] = (isESD ? ((AliESDVertex*)eventVtxSPD)->GetZ() : ((AliAODVertex*)eventVtxSPD)->GetZ());
     eventInfo->fNVtxSPDContributors = eventVtxSPD->GetNContributors();
  }
  
  // ------------------------------------------------------------------------------------------------------------------
  // Improved cut on the distance between SPD and track vertices 
  // See Francesco Prino's slides during Physics Forum from 5 october 2016, slide 33
  //  based on input from Ruben Shahoyan and Alex Dobrin
  //------------------------------------------------------------------------------------------------------------------
  Bool_t vertexDistanceSelected = kTRUE;
  if(!eventVtx) vertexDistanceSelected = kFALSE;
  if(!eventVtxSPD) vertexDistanceSelected = kFALSE;
  if(vertexDistanceSelected) {
     if(eventVtx->GetNContributors()<2 || eventVtxSPD->GetNContributors()<1) vertexDistanceSelected = kFALSE;   
  }
  if(vertexDistanceSelected) {
     Double_t covTracks[6], covSPD[6];
     eventVtx->GetCovarianceMatrix(covTracks);
     eventVtxSPD->GetCovarianceMatrix(covSPD);
     Double_t dz = eventVtx->GetZ() - eventVtxSPD->GetZ();
     Double_t errTot = TMath::Sqrt(covTracks[5]+covSPD[5]);
     Double_t errTrk = TMath::Sqrt(covTracks[5]);
     Double_t nsigTot = TMath::Abs(dz)/errTot;
     Double_t nsigTrk = TMath::Abs(dz)/errTrk;
     if(TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrk>20) vertexDistanceSelected = kFALSE;
  }
  if(vertexDistanceSelected) fReducedEvent->fEventTag |= (ULong64_t(1)<<13);
  //-------------------------------------------------------------------------------------------------------------------------
  
  eventInfo->fBC          = event->GetBunchCrossNumber();
  eventInfo->fEventType   = event->GetEventType();
  eventInfo->fOnlineTriggerMask = event->GetTriggerMask();
  eventInfo->fOnlineTriggerMaskNext50 = event->GetTriggerMaskNext50();
  eventInfo->fTriggerMask = inputHandler->IsEventSelected();
  eventInfo->fTriggerClass = event->GetFiredTriggerClasses();
  eventInfo->fIsPhysicsSelection = (isSelected!=0 ? kTRUE : kFALSE);
  eventInfo->fIsSPDPileup = event->IsPileupFromSPD(3,0.8,3.,2.,5.);
  eventInfo->fIsSPDPileupMultBins = event->IsPileupFromSPDInMultBins();
  
  if(isESD) {
    eventInfo->fEventNumberInFile = esdEvent->GetEventNumberInFile();
    eventInfo->fL0TriggerInputs = esdEvent->GetHeader()->GetL0TriggerInputs();
    eventInfo->fL1TriggerInputs = esdEvent->GetHeader()->GetL1TriggerInputs();
    eventInfo->fL2TriggerInputs = esdEvent->GetHeader()->GetL2TriggerInputs();

    TString trgClasses = esdEvent->GetFiredTriggerClasses();
    if((trgClasses.Contains("HQU")) && (trgClasses.Contains("HSE"))) eventInfo->fTRDfired = 3;
    else {
	if(trgClasses.Contains("HQU")) eventInfo->fTRDfired = 1;
	if(trgClasses.Contains("HSE")) eventInfo->fTRDfired = 2;
    }

    eventInfo->fIRIntClosestIntMap[0] = esdEvent->GetHeader()->GetIRInt1ClosestInteractionMap();
    eventInfo->fIRIntClosestIntMap[1] = esdEvent->GetHeader()->GetIRInt2ClosestInteractionMap();
    eventVtx = const_cast<AliESDVertex*>(esdEvent->GetPrimaryVertexTPC());
    if(eventVtx) {
      eventInfo->fVtxTPC[0] = ((AliESDVertex*)eventVtx)->GetX();
      eventInfo->fVtxTPC[1] = ((AliESDVertex*)eventVtx)->GetY();
      eventInfo->fVtxTPC[2] = ((AliESDVertex*)eventVtx)->GetZ();
      eventInfo->fNVtxTPCContributors = eventVtx->GetNContributors();
    }
    eventInfo->fTimeStamp     = esdEvent->GetTimeStamp();
    eventInfo->fNpileupSPD    = esdEvent->GetNumberOfPileupVerticesSPD();
    eventInfo->fNpileupTracks = esdEvent->GetNumberOfPileupVerticesTracks();
    eventInfo->fNPMDtracks    = esdEvent->GetNumberOfPmdTracks();
    eventInfo->fNTRDtracks    = esdEvent->GetNumberOfTrdTracks();
    eventInfo->fNTRDtracklets = esdEvent->GetNumberOfTrdTracklets();
    eventInfo->fNTPCclusters  = esdEvent->GetNumberOfTPCClusters();
    eventInfo->fNtracksTPCout = esdEvent->GetNTPCTrackBeforeClean();
    
    for(Int_t ilayer=0; ilayer<2; ++ilayer)
      eventInfo->fSPDFiredChips[ilayer] = esdEvent->GetMultiplicity()->GetNumberOfFiredChips(ilayer);
    for(Int_t ilayer=0; ilayer<6; ++ilayer)
      eventInfo->fITSClusters[ilayer] = esdEvent->GetMultiplicity()->GetNumberOfITSClusters(ilayer);
    eventInfo->fSPDnSingle = esdEvent->GetMultiplicity()->GetNumberOfSingleClusters();
    
    // ZDC information
    AliESDZDC* zdc = esdEvent->GetESDZDC();
    if(zdc) {
      eventInfo->fZDCnTotalEnergy[0] = zdc->GetZN2TowerEnergy()[0];
      eventInfo->fZDCnTotalEnergy[1] = zdc->GetZN1TowerEnergy()[0];
      eventInfo->fZDCpTotalEnergy[0] = zdc->GetZP2TowerEnergy()[0];
      eventInfo->fZDCpTotalEnergy[1] = zdc->GetZP1TowerEnergy()[0];
      for(Int_t i=0; i<5; ++i)  eventInfo->fZDCnEnergy[i]   = zdc->GetZN1TowerEnergy()[i];
      for(Int_t i=5; i<10; ++i)  eventInfo->fZDCnEnergy[i]   = zdc->GetZN2TowerEnergy()[i-5];
      for(Int_t i=0; i<5; ++i)  eventInfo->fZDCpEnergy[i]   = zdc->GetZP1TowerEnergy()[i];
      for(Int_t i=5; i<10; ++i)  eventInfo->fZDCpEnergy[i]   = zdc->GetZP2TowerEnergy()[i-5];      
    }
    
    // T0 information
    const AliESDTZERO* tzero = esdEvent->GetESDTZERO();
    if(tzero) {
      eventInfo->fT0start = tzero->GetT0();
      eventInfo->fT0zVertex = tzero->GetT0zVertex();
      for(Int_t i = 0;i<24;i++)
        eventInfo->fT0amplitude[i] = tzero->GetT0amplitude()[i];
      for(Int_t i = 0;i<3;i++)
        eventInfo->fT0TOF[i] = tzero->GetT0TOF()[i];
      for(Int_t i = 0;i<3;i++)
        eventInfo->fT0TOFbest[i] = tzero->GetT0TOFbest()[i];
      eventInfo->fT0pileup = tzero->GetPileupFlag();
      eventInfo->fT0sattelite = tzero->GetSatellite();
    }
    eventInfo->fDiamondDim[0] = esdEvent->GetDiamondX(); 
    eventInfo->fDiamondDim[1] = esdEvent->GetDiamondY();
    eventInfo->fDiamondDim[2] = esdEvent->GetDiamondZ();
    Float_t cov[3]; esdEvent->GetDiamondCovXY(cov);
    for(Int_t icomp=0; icomp<3; icomp++) 
      eventInfo->fDiamondCov[icomp] = cov[icomp];
  }
  if(isAOD) {
    eventInfo->fIRIntClosestIntMap[0] = aodEvent->GetHeader()->GetIRInt1ClosestInteractionMap();
    eventInfo->fIRIntClosestIntMap[1] = aodEvent->GetHeader()->GetIRInt2ClosestInteractionMap();
    eventInfo->fEventNumberInFile = aodEvent->GetEventNumberInFile();
    eventInfo->fL0TriggerInputs = aodEvent->GetHeader()->GetL0TriggerInputs();
    eventInfo->fL1TriggerInputs = aodEvent->GetHeader()->GetL1TriggerInputs();
    eventInfo->fL2TriggerInputs = aodEvent->GetHeader()->GetL2TriggerInputs();

    TString trgClasses = aodEvent->GetFiredTriggerClasses();
    eventInfo->fTRDfired = 0;
    if((trgClasses.Contains("HQU")) && (trgClasses.Contains("HSE"))) eventInfo->fTRDfired = 3;
    else {
	if(trgClasses.Contains("HQU")) eventInfo->fTRDfired = 1;
	if(trgClasses.Contains("HSE")) eventInfo->fTRDfired = 2;
    }

    eventInfo->fTimeStamp     = aodEvent->GetTimeStamp();
    eventInfo->fNpileupSPD    = aodEvent->GetNumberOfPileupVerticesSPD();
    eventInfo->fNpileupTracks = aodEvent->GetNumberOfPileupVerticesTracks();
    eventInfo->fNPMDtracks    = aodEvent->GetNPmdClusters();
    eventInfo->fNTRDtracks    = aodEvent->GetNumberOfTrdTracks();
    eventInfo->fNTRDtracklets = 0;
    eventInfo->fNTPCclusters  = aodEvent->GetNumberOfTPCClusters();
    
    eventVtx = const_cast<AliAODVertex*>(aodEvent->GetPrimaryVertexTPC());
    if(eventVtx) {
       eventInfo->fVtxTPC[0] = ((AliAODVertex*)eventVtx)->GetX();
       eventInfo->fVtxTPC[1] = ((AliAODVertex*)eventVtx)->GetY();
       eventInfo->fVtxTPC[2] = ((AliAODVertex*)eventVtx)->GetZ();
       eventInfo->fNVtxTPCContributors = eventVtx->GetNContributors();
    }
    
    for(Int_t ilayer=0; ilayer<2; ++ilayer)
      eventInfo->fSPDFiredChips[ilayer] = aodEvent->GetMultiplicity()->GetNumberOfFiredChips(ilayer);
    for(Int_t ilayer=0; ilayer<6; ++ilayer)
       eventInfo->fITSClusters[ilayer] = aodEvent->GetMultiplicity()->GetNumberOfITSClusters(ilayer);
  
    // ZDC information
    AliAODZDC* zdc = aodEvent->GetZDCData();
    if(zdc) {
       eventInfo->fZDCnTotalEnergy[0] = zdc->GetZNATowerEnergy()[0];
       eventInfo->fZDCnTotalEnergy[1] = zdc->GetZNCTowerEnergy()[0];
       eventInfo->fZDCpTotalEnergy[0] = zdc->GetZPATowerEnergy()[0];
       eventInfo->fZDCpTotalEnergy[1] = zdc->GetZPCTowerEnergy()[0];
      for(Int_t i=0; i<5; ++i)  eventInfo->fZDCnEnergy[i]   = zdc->GetZNATowerEnergy()[i];
      for(Int_t i=5; i<10; ++i)  eventInfo->fZDCnEnergy[i]   = zdc->GetZNCTowerEnergy()[i-5];
      for(Int_t i=0; i<5; ++i)  eventInfo->fZDCpEnergy[i]   = zdc->GetZPATowerEnergy()[i];
      for(Int_t i=5; i<10; ++i)  eventInfo->fZDCpEnergy[i]   = zdc->GetZPCTowerEnergy()[i-5];
    }
    
    // T0 information
    AliAODTZERO* tzero = aodEvent->GetTZEROData();
    if(tzero) {
      eventInfo->fT0start = -999.;   // not available
      eventInfo->fT0zVertex = tzero->GetT0zVertex();
      for(Int_t i = 0;i<26;i++)
        eventInfo->fT0amplitude[i] = tzero->GetAmp(i);
      for(Int_t i = 0;i<3;i++)
        eventInfo->fT0TOF[i] = tzero->GetT0TOF()[i];
      for(Int_t i = 0;i<3;i++)
        eventInfo->fT0TOFbest[i] = tzero->GetT0TOFbest()[i];
      eventInfo->fT0pileup = tzero->GetPileupFlag();
      eventInfo->fT0sattelite = tzero->GetSatellite();
    }
    eventInfo->fDiamondDim[0] = aodEvent->GetDiamondX(); 
    eventInfo->fDiamondDim[1] = aodEvent->GetDiamondY();
    eventInfo->fDiamondDim[2] = aodEvent->GetDiamondZ();
    Float_t cov[3]; aodEvent->GetDiamondCovXY(cov);
    for(Int_t icomp=0; icomp<3; icomp++) 
      eventInfo->fDiamondCov[icomp] = cov[icomp];
  }
  
  // V0 information
  AliVVZERO* vzero = event->GetVZEROData();
  for(Int_t i=0;i<64;++i) 
    eventInfo->fVZEROMult[i] = vzero->GetMultiplicity(i);  
  Float_t multVZERO = 0.0;
  for(Int_t i=0;i<32;++i) multVZERO +=  vzero->GetMultiplicity(i);
  eventInfo->fVZEROTotalMult[1] = multVZERO;
  multVZERO = 0.0;
  for(Int_t i=32;i<64;++i) multVZERO +=  vzero->GetMultiplicity(i);
  eventInfo->fVZEROTotalMult[0] = multVZERO;
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillCaloClusters() {
  //
  // Fill info about the calorimeter clusters
  //
  AliVEvent* event = InputEvent();
  Int_t nclusters = event->GetNumberOfCaloClusters();

  AliReducedEventInfo* eventInfo = dynamic_cast<AliReducedEventInfo*>(fReducedEvent);
  if(!eventInfo) return;
  
  eventInfo->fNCaloClusters = 0;
  for(Int_t iclus=0; iclus<nclusters; ++iclus) {
    AliVCluster* cluster = event->GetCaloCluster(iclus);

    Bool_t clusterFilterDecision = kTRUE;
    std::vector<Bool_t> individualFilterDecisions;
    if (fClusterFilter.GetEntries()>0)  clusterFilterDecision = IsClusterSelected(cluster, individualFilterDecisions);
    if (!clusterFilterDecision) continue;
    
    TClonesArray& clusters = *(eventInfo->fCaloClusters);
    AliReducedCaloClusterInfo *reducedCluster=new(clusters[eventInfo->fNCaloClusters]) AliReducedCaloClusterInfo();
    
    reducedCluster->fClusterID = iclus;
    reducedCluster->fType    = (cluster->IsEMCAL() ? AliReducedCaloClusterInfo::kEMCAL : AliReducedCaloClusterInfo::kPHOS);
    reducedCluster->fEnergy  = cluster->E();
    reducedCluster->fTrackDx = cluster->GetTrackDx();
    reducedCluster->fTrackDz = cluster->GetTrackDz();
    reducedCluster->fM20     = cluster->GetM20();
    reducedCluster->fM02     = cluster->GetM02();
    reducedCluster->fDispersion = cluster->GetDispersion();
    reducedCluster->fNMatchedTracks = cluster->GetNTracksMatched();
    cluster->GetPosition(reducedCluster->fPosition);
    reducedCluster->fTOF = cluster->GetTOF();
    reducedCluster->fNCells = cluster->GetNCells();
    eventInfo->fNCaloClusters += 1;
  }  // end loop over clusters
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillFMDInfo(Bool_t isAOD) {
  //
  // Fill FMD information  
  //
  Float_t m;
  AliReducedEventInfo *eventInfo = dynamic_cast<AliReducedEventInfo*>(fReducedEvent);
  if(!eventInfo) return;
  TClonesArray &fmd = *(eventInfo->GetFMD());

  if(isAOD) {
    AliAODEvent *aodEvent = static_cast<AliAODEvent*>(InputEvent());
    TObject *obj = aodEvent->FindListObject("Forward");
    if (!obj) return;
    AliAODForwardMult *aodForward = static_cast<AliAODForwardMult*>(obj);
    TH2D &d2Ndetadphi = aodForward->GetHistogram();
    Int_t nFMD = -1;
    // Loop over Eta
    for (Int_t iEta = 1; iEta <= d2Ndetadphi.GetNbinsX(); iEta++) {
      // Loop over phi
      for (Int_t iPhi = 1; iPhi <= d2Ndetadphi.GetNbinsY(); iPhi++) {
        m = d2Ndetadphi.GetBinContent(iEta, iPhi);
        if(m<1E-6) continue;
        nFMD++;
        AliReducedFMDInfo *reducedFMD = (AliReducedFMDInfo*) fmd.ConstructedAt(nFMD);
        reducedFMD->fMultiplicity = m;
        reducedFMD->fId = iEta * d2Ndetadphi.GetNbinsY() + iPhi;
      }
    }
  } 
  else {
    AliAODEvent* aodEvent = AliForwardUtil::GetAODEvent(this);
    if (!aodEvent) {cout<<"didn't get AOD"<<endl; return;}
    TH2D* histos[5];
    histos[0] = static_cast<TH2D*>(aodEvent->FindListObject("FMD1I_cache"));
    histos[1] = static_cast<TH2D*>(aodEvent->FindListObject("FMD2I_cache"));
    histos[2] = static_cast<TH2D*>(aodEvent->FindListObject("FMD2O_cache"));
    histos[3] = static_cast<TH2D*>(aodEvent->FindListObject("FMD3I_cache"));
    histos[4] = static_cast<TH2D*>(aodEvent->FindListObject("FMD3O_cache"));
    // Loop over eta
    Int_t nFMD = -1;
    for (Int_t ih = 0; ih < 5; ih++) {
      if(!histos[ih]) continue;
      for (Int_t iEta = 1; iEta <= histos[ih]->GetNbinsX(); iEta++) {
        // Loop over phi
        for (Int_t iPhi = 1; iPhi <= histos[ih]->GetNbinsY(); iPhi++) {
        m = histos[ih]->GetBinContent(iEta, iPhi);
        if(m<1E-6) continue;
        nFMD++;
        AliReducedFMDInfo *reducedFMD = new(fmd[nFMD]) AliReducedFMDInfo();
        reducedFMD->fMultiplicity = m;
        reducedFMD->fId = iEta*histos[ih]->GetNbinsY()+iPhi;
        if(ih == 2 || ih == 4) reducedFMD->fId *= -1;
        }
      }
    }
  }
}

//________________________________________________________________________________________
Double_t AliAnalysisTaskReducedTreeMaker::Rapidity(Double_t r, Double_t z){
  //
  // calculate eta based on radius from beampipe r and distance from interaction point z
  //
  Double_t x = r/z;
  if(z<0) x = x*-1;

  Double_t eta = -1.*TMath::Log((TMath::Sqrt(x*x+1)-1)/x);

  if(z<0) eta = eta*-1;
  return eta;
}

//________________________________________________________________________________________
Double_t AliAnalysisTaskReducedTreeMaker::Radius(Double_t eta, Double_t z){
  //
  // calculate radius from beampipe based on distance from interaction point z and eta
  //
  Double_t r = 2*TMath::Power(TMath::E(), eta)*z/(TMath::Power(TMath::E(), 2.0*eta)-1);
  return r;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeMaker::CheckPDGcode(AliMCEvent* event, Int_t ipart, AliSignalMC* mcSignal) {
   //
   // Check that the particle satisfies the PDG code criteria specified in the mcSignal
   // Work on just 1 pronged MC signals here
   // Method: All of the defined generations of the prong must fulfill the defined PDG criteria
   // 
   if(mcSignal->GetNProngs()>1) return kFALSE;
   
   // loop over all generations
   AliVParticle* currentGenerationParticle = event->GetTrack(ipart);
   Int_t currentGenerationLabel = ipart;
   for(UInt_t ig=0; ig<mcSignal->GetNGenerations(); ++ig) {      
      // test the PDG code of this particle
      // In case the MC history finished (no current particle), test the MC signal using the not assigned PDG.
      // If there is no PDG requested in this generation, the MC test can still pass
      if(!mcSignal->TestPDG(0, ig, currentGenerationParticle ? currentGenerationParticle->PdgCode() : AliSignalMC::kPDGnotAssigned)) 
         return kFALSE;
      
      // get the next generation
      currentGenerationLabel = (currentGenerationParticle ? currentGenerationParticle->GetMother() : 0);
      currentGenerationParticle = (currentGenerationParticle ? event->GetTrack(currentGenerationLabel) : 0x0);
   }
   return kTRUE;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeMaker::CheckParticleSource(AliMCEvent* event, Int_t ipart, AliSignalMC* mcSignal) {
   //
   // Check that the particle satisfies the source criteria specified in the mcSignal
   // Work on just 1 pronged MC signals here
   // Method: All of the defined generations of the prong must fulfill the defined source bit map
   //               For a given generation, all the sources for which corresponding bits are enabled, must be fulfilled 
   // 
   if(mcSignal->GetNProngs()>1) return kFALSE;
   
   // loop over all generations
   AliVParticle* currentGenerationParticle = event->GetTrack(ipart);
   Int_t currentGenerationLabel = ipart;
   for(UInt_t ig=0; ig<mcSignal->GetNGenerations(); ++ig) {
      if(!mcSignal->GetSources(0,ig)) continue;       // no sources requested
      if(!currentGenerationParticle) return kFALSE;   // if there are sources requested, but MC history finished, evaluate to FALSE
      
      // check all implemented sources
      UInt_t decision = 0;
      // use logical XOR between the presence of a given source and the exclude flag
      if(mcSignal->CheckSourceBit(0,ig, AliSignalMC::kPhysicalPrimary)) { 
         if(mcSignal->GetSourceExclude(0,ig,AliSignalMC::kPhysicalPrimary) != event->IsPhysicalPrimary(currentGenerationLabel)) 
            decision |= (UInt_t(1) << AliSignalMC::kPhysicalPrimary);
      }
      if(mcSignal->CheckSourceBit(0,ig, AliSignalMC::kFromBGEvent)) { 
         if(mcSignal->GetSourceExclude(0,ig,AliSignalMC::kFromBGEvent) != event->IsFromBGEvent(currentGenerationLabel)) 
            decision |= (UInt_t(1) << AliSignalMC::kFromBGEvent);
      }
      if(mcSignal->CheckSourceBit(0,ig, AliSignalMC::kSecondaryFromWeakDecay)) { 
         if(mcSignal->GetSourceExclude(0,ig,AliSignalMC::kSecondaryFromWeakDecay) != event->IsSecondaryFromWeakDecay(currentGenerationLabel)) 
            decision |= (UInt_t(1) << AliSignalMC::kSecondaryFromWeakDecay);
      }
      if(mcSignal->CheckSourceBit(0,ig, AliSignalMC::kSecondaryFromMaterial)) { 
         if(mcSignal->GetSourceExclude(0,ig,AliSignalMC::kSecondaryFromMaterial) != event->IsSecondaryFromMaterial(currentGenerationLabel)) 
            decision |= (UInt_t(1) << AliSignalMC::kSecondaryFromMaterial);
      }
      if(mcSignal->CheckSourceBit(0,ig, AliSignalMC::kFromSubsidiaryEvent)) { 
         if(mcSignal->GetSourceExclude(0,ig,AliSignalMC::kFromSubsidiaryEvent) != event->IsFromSubsidiaryEvent(currentGenerationLabel)) 
            decision |= (UInt_t(1) << AliSignalMC::kFromSubsidiaryEvent);
      }
      if(mcSignal->CheckSourceBit(0,ig, AliSignalMC::kRadiativeDecay)) { 
         if(mcSignal->GetSourceExclude(0,ig,AliSignalMC::kRadiativeDecay) != (currentGenerationParticle->GetNDaughters()>2)) 
            decision |= (UInt_t(1) << AliSignalMC::kRadiativeDecay);
      }
      if(mcSignal->CheckSourceBit(0,ig, AliSignalMC::kFirstInStack)) { 
         if(mcSignal->GetSourceExclude(0,ig,AliSignalMC::kFirstInStack) != (ipart==0)) 
            decision |= (UInt_t(1) << AliSignalMC::kFirstInStack);
      }
      if(mcSignal->CheckSourceBit(0,ig, AliSignalMC::kSecondInStack)) { 
         if(mcSignal->GetSourceExclude(0,ig,AliSignalMC::kSecondInStack) != (ipart==1)) 
            decision |= (UInt_t(1) << AliSignalMC::kSecondInStack);
      }
      if(mcSignal->CheckSourceBit(0,ig, AliSignalMC::kFirstTenInStack)) { 
         if(mcSignal->GetSourceExclude(0,ig,AliSignalMC::kFirstTenInStack) != (ipart<10)) 
            decision |= (UInt_t(1) << AliSignalMC::kFirstTenInStack);
      }
      
      if(!decision) return kFALSE;
      decision &= mcSignal->GetSources(0,ig);
      if(mcSignal->GetUseANDonSourceBits(0,ig) && (decision != mcSignal->GetSources(0,ig))) return kFALSE;  // not all req sources are fullfilled
      
      // get the next generation
      currentGenerationLabel = (currentGenerationParticle ? currentGenerationParticle->GetMother() : 0);
      currentGenerationParticle = (currentGenerationParticle ? event->GetTrack(currentGenerationLabel) : 0x0);
   }
   return kTRUE;
}

//_________________________________________________________________________________
UInt_t AliAnalysisTaskReducedTreeMaker::MatchMCsignals(Int_t iparticle) {
   //
   // check whether the defined MC signals match this particle
   //
   if(!AliDielectronMC::Instance()->HasMC()) return 0;
   
   Int_t nMCsignals = fMCsignals.GetEntries();
   if(!nMCsignals) return 0;
   
   AliMCEvent* event = AliDielectronMC::Instance()->GetMCEvent();
   
   UInt_t mcSignalsMap = 0;
   for(Int_t isig=0; isig<nMCsignals; ++isig) {
      Bool_t mcMatch = CheckPDGcode(event, iparticle, (AliSignalMC*)fMCsignals.At(isig)) && 
                                    CheckParticleSource(event, iparticle, (AliSignalMC*)fMCsignals.At(isig));
      
      if(mcMatch)
         mcSignalsMap |= (UInt_t(1)<<isig);
   }
   return mcSignalsMap;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeMaker::CheckMCtruthWriteFormat(UInt_t bitMap) {
   //
   // For the bits which are on, check which writing options were requested
   // If both base and full track formats are requested, the track will be written as full track
   // Return TRUE if base track format is chosen, and FALSE otherwise 
   //
   Bool_t writeBaseTrack = kTRUE;
   for(Int_t iflag=0;iflag<32;++iflag) {
      if(!(bitMap & (UInt_t(1)<<iflag))) continue;
      if(fMCsignalsWritingOptions[iflag]==kFullTrack) writeBaseTrack = kFALSE;
   }
   return writeBaseTrack;
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillMCTruthInfo() 
{
   //
   // fill MC truth info
   //
   Bool_t hasMC = AliDielectronMC::Instance()->HasMC();
   if(!hasMC) return;
   Int_t nMCsignals = fMCsignals.GetEntries();
   if(!nMCsignals) return;
   AliInputEventHandler* inputHandler = (AliInputEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
   
   AliMCEvent* event = AliDielectronMC::Instance()->GetMCEvent();
   if(!event) return;
   
   AliReducedEventInfo* eventInfo = NULL; 
   if(fTreeWritingOption==kFullEventsWithBaseTracks || fTreeWritingOption==kFullEventsWithFullTracks) {
     eventInfo = dynamic_cast<AliReducedEventInfo*>(fReducedEvent);
     const AliVVertex* mcVtx = event->GetPrimaryVertex();
     if(mcVtx){
       eventInfo->fVtxMC[0] = mcVtx->GetX();
       eventInfo->fVtxMC[1] = mcVtx->GetY();
       eventInfo->fVtxMC[2] = mcVtx->GetZ();
     } 
   }
   
   // We loop over all particles in the MC event
   for(Int_t i=0; i<event->GetNumberOfTracks(); ++i) {
      AliVParticle* particle = event->GetTrack(i);
      if(!particle) continue;
      // fill MC truth number of charged particles
      if(eventInfo){
        if(particle->IsPhysicalPrimary() && particle->Charge()){
          Float_t eta = particle->Eta();
          Float_t etaAbs = TMath::Abs(eta);
          if(etaAbs < 1.6) eventInfo->fNch[0]++;
          if(etaAbs < 1.0) eventInfo->fNch[2]++;
          else if(eta > 2.8 && eta < 5.1) eventInfo->fNch[4]++;   // V0A
          else if(eta > -3.7 && eta < -1.7) eventInfo->fNch[6]++; // V0C
          
          // check if particle is J/psi daughter
          AliVParticle* mother = event->GetTrack(particle->GetMother());
          Bool_t jpsiDaughter = (mother && mother->PdgCode() == 443);
          if(!jpsiDaughter){
            if(etaAbs < 1.6) eventInfo->fNch[1]++;
            if(etaAbs < 1.0) eventInfo->fNch[3]++;
            else if(eta > 2.8 && eta < 5.1 ) eventInfo->fNch[5]++;   // V0A
            else if(eta > -3.7 && eta < -1.7 ) eventInfo->fNch[7]++; // V0C
          }
        }
      }
      
      UInt_t mcSignalsMap = MatchMCsignals(i);    // check which MC signals match this particle and fill the bit map
      if(!mcSignalsMap) continue;
      
      // fill MC statistics summary
      for(Int_t iTrig=0;iTrig<32;++iTrig) {
         if(inputHandler->IsEventSelected() & (UInt_t(1)<<iTrig)) {
            for(Int_t iSig=0;iSig<fMCsignals.GetEntries();++iSig) {
               if(mcSignalsMap & (UInt_t(1)<<iSig)) fMCSignalsHistogram->Fill(Double_t(iSig), Double_t(iTrig));
            }
         }
      }
      
      Bool_t writeBaseTrack = kFALSE;      // if false write full track format
      writeBaseTrack = CheckMCtruthWriteFormat(mcSignalsMap);  // check which track format (base/full) should be used
      // write the track in the first track array if the format is full track
      // if the track format is base track then write it on either the first or the second array, depending on the tree writing options
      Bool_t useFirstTrackArray = kTRUE;
      if(writeBaseTrack) {
         if(fTreeWritingOption==kBaseEventsWithFullTracks || fTreeWritingOption==kFullEventsWithFullTracks)
            useFirstTrackArray = kFALSE;
      }
      
      if(writeBaseTrack) fNSelectedBaseTracks[fTrackFilter.GetEntries()+4] += 1;
      else fNSelectedFullTracks[fTrackFilter.GetEntries()+12] += 1;
      
      TClonesArray* trackArrPointer = fReducedEvent->fTracks;
      if(!useFirstTrackArray) trackArrPointer = fReducedEvent->fTracks2;
      TClonesArray& tracks = *(trackArrPointer);
      Int_t currentTrackIdx = tracks.GetEntries();
      
      AliReducedBaseTrack* reducedParticle=NULL;
      if(writeBaseTrack) 
         reducedParticle=new(tracks[currentTrackIdx]) AliReducedBaseTrack();
      else
         reducedParticle=new(tracks[currentTrackIdx]) AliReducedTrackInfo();
      
      reducedParticle->fMCFlags = mcSignalsMap;
      reducedParticle->fIsMCTruth = kTRUE;
      reducedParticle->PxPyPz(particle->Px(), particle->Py(), particle->Pz());
      reducedParticle->Charge(particle->Charge());
   
      if(writeBaseTrack) continue;
      
      AliReducedTrackInfo* trackInfo = dynamic_cast<AliReducedTrackInfo*>(reducedParticle);
      if(!trackInfo) continue;
      
      trackInfo->fMCLabels[0] = particle->GetLabel();
      trackInfo->fMCPdg[0] = particle->PdgCode();
      trackInfo->fMCMom[0] = particle->Px();
      trackInfo->fMCMom[1] = particle->Py();
      trackInfo->fMCMom[2] = particle->Pz();
      trackInfo->fMCFreezeout[0] = particle->Xv();
      trackInfo->fMCFreezeout[1] = particle->Yv();
      trackInfo->fMCFreezeout[2] = particle->Zv();
      
      AliVParticle* mother = event->GetTrack(particle->GetMother());
      if(mother) {
         trackInfo->fMCLabels[1] = mother->GetLabel();
         trackInfo->fMCPdg[1] = mother->PdgCode();
         
         AliVParticle* grandmother = event->GetTrack(mother->GetMother());
         if(grandmother) {
            trackInfo->fMCLabels[2] = grandmother->GetLabel();
            trackInfo->fMCPdg[2] = grandmother->PdgCode();
            
            AliVParticle* grandgrandmother = event->GetTrack(grandmother->GetMother());
            if(grandgrandmother) {
               trackInfo->fMCLabels[3] = grandgrandmother->GetLabel();
               trackInfo->fMCPdg[3] = grandgrandmother->PdgCode();
            }
         }
      }
      if(fFillHFInfo) trackInfo->fHFProc = AliDielectronMC::Instance()->GetHFProcess(particle->GetLabel());
        
      fReducedEvent->fNtracks[1] += 1;  
   }
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillTrackInfo() 
{
  //
  // fill reduced track information
  //
  AliVEvent* event = InputEvent();
  Bool_t isESD = (event->IsA()==AliESDEvent::Class());
  Bool_t isAOD = (event->IsA()==AliAODEvent::Class());
  AliESDEvent* esdEvent = 0x0;
  AliAODEvent* aodEvent = 0x0;
  if(isESD) esdEvent = static_cast<AliESDEvent*>(InputEvent());  
  if(isAOD) aodEvent = static_cast<AliAODEvent*>(InputEvent());  
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  AliPIDResponse* pidResponse = inputHandler->GetPIDResponse();
  
  Bool_t hasMC = AliDielectronMC::Instance()->HasMC();
  
  // find all the tracks which belong to a V0 stored in the reduced event
  UShort_t trackIdsV0[4][20000]={{0}};
  UShort_t trackIdsPureV0[4][20000]={{0}};
  Int_t nV0LegsTagged[4] = {0}; Int_t nPureV0LegsTagged[4] = {0};
  Bool_t leg1Found[4]; Bool_t leg2Found[4];
  for(Int_t iv0=0;iv0<fReducedEvent->fNV0candidates[1];++iv0) {
    AliReducedPairInfo* pair = fReducedEvent->GetV0Pair(iv0);
    if(!pair) continue;
    Int_t pairId = 0; Bool_t isPureV0 = kFALSE;
    if(pair->fCandidateId==AliReducedPairInfo::kGammaConv) {
      pairId=0;
      if(pair->IsPureV0Gamma()) isPureV0 = kTRUE;
    }
    if(pair->fCandidateId==AliReducedPairInfo::kK0sToPiPi) {
      pairId=1;
      if(pair->IsPureV0K0s()) isPureV0 = kTRUE;
    }
    if(pair->fCandidateId==AliReducedPairInfo::kLambda0ToPPi) {
      pairId=2;
      if(pair->IsPureV0Lambda()) isPureV0 = kTRUE;
    }
    if(pair->fCandidateId==AliReducedPairInfo::kALambda0ToPPi) {
      pairId=3;
      if(pair->IsPureV0ALambda()) isPureV0 = kTRUE;
    }
    
    leg1Found[pairId] = kFALSE; leg2Found[pairId] = kFALSE;
    for(Int_t it=0;it<nV0LegsTagged[pairId];++it) {
      if(trackIdsV0[pairId][it]==pair->fLegIds[0]) leg1Found[pairId]=kTRUE;
      if(trackIdsV0[pairId][it]==pair->fLegIds[1]) leg2Found[pairId]=kTRUE;
    }
    // if the legs of this V0 were not already stored then add them now to the list
    if(!leg1Found[pairId]) {trackIdsV0[pairId][nV0LegsTagged[pairId]] = pair->fLegIds[0]; ++nV0LegsTagged[pairId];}
    if(!leg2Found[pairId]) {trackIdsV0[pairId][nV0LegsTagged[pairId]] = pair->fLegIds[1]; ++nV0LegsTagged[pairId];}
    
    if(isPureV0) {
      leg1Found[pairId] = kFALSE; leg2Found[pairId] = kFALSE;
      for(Int_t it=0;it<nPureV0LegsTagged[pairId];++it) {
        if(trackIdsPureV0[pairId][it]==pair->fLegIds[0]) leg1Found[pairId]=kTRUE;
        if(trackIdsPureV0[pairId][it]==pair->fLegIds[1]) leg2Found[pairId]=kTRUE;
      }
      // if the legs of this pure V0 were not already stored then add them now to the list
      if(!leg1Found[pairId]) {trackIdsPureV0[pairId][nPureV0LegsTagged[pairId]] = pair->fLegIds[0]; ++nPureV0LegsTagged[pairId];}
      if(!leg2Found[pairId]) {trackIdsPureV0[pairId][nPureV0LegsTagged[pairId]] = pair->fLegIds[1]; ++nPureV0LegsTagged[pairId];}
    }
  }
        
  // check for tracks matched in TRD 
  Int_t trackIdsTRD[20000]={0};
  Int_t trackTRDGTUtracklets[20000]={0};
  Int_t trackTRDGTUlayermask[20000]={0};
  Double_t trackTRDGTUpt[20000]={0};
  Float_t trackTRDGTUsagitta[20000]={2};
  Int_t trackTRDGTUPID[20000]={0};
  Int_t nTracksTRD = 0;
  if(fFillTRDMatchedTracks) {
    for(Int_t itrackTRD=0; itrackTRD<event->GetNumberOfTrdTracks(); ++itrackTRD) {
       if(isESD) {
         AliESDTrdTrack* trdTrack = (AliESDTrdTrack*)esdEvent->GetTrdTrack(itrackTRD);
         if(!trdTrack) {
            cout << "############## Bad TRD track found" << endl;
            continue;
         }
         AliESDtrack* tempESDtrack = dynamic_cast<AliESDtrack*>(trdTrack->GetTrackMatch());
         if(!tempESDtrack) continue;
         Int_t trackID = tempESDtrack->GetID();
         Bool_t found = kFALSE;
         for(Int_t k=0; k<nTracksTRD; ++k) {
           if(trackID==trackIdsTRD[k]) {
              found = kTRUE;   
              break;
           }
         }
         if(!found) {
	       trackIdsTRD[nTracksTRD] = trackID;
	       trackTRDGTUtracklets[nTracksTRD] = trdTrack->GetNTracklets();
	       if((trdTrack->GetLayerMask() & fTRDtrglayerMaskEl) != fTRDtrglayerMaskEl) trackTRDGTUlayermask[nTracksTRD] = 0;
           else trackTRDGTUlayermask[nTracksTRD] = 1;
	       trackTRDGTUpt[nTracksTRD] = trdTrack->Pt();
	       Int_t b = trdTrack->GetB();
	       Int_t c = trdTrack->GetC();
	       trackTRDGTUsagitta[nTracksTRD] = GetInvPtDevFromBC(b,c);
	       trackTRDGTUPID[nTracksTRD] = trdTrack->GetPID();
	       nTracksTRD++;
         }
       }
       if(isAOD) {
          AliAODTrdTrack* trdTrack = (AliAODTrdTrack*)aodEvent->GetTrdTrack(itrackTRD);
          if(!trdTrack) {
             cout << "############## Bad TRD track found" << endl;
             continue;
          }
          AliAODTrack* tempAODtrack = dynamic_cast<AliAODTrack*>(trdTrack->GetTrackMatch());
          if(!tempAODtrack) continue;
          Int_t trackID = tempAODtrack->GetID();
          Bool_t found = kFALSE;
          for(Int_t k=0; k<nTracksTRD; ++k) {
             if(trackID==trackIdsTRD[k]) {
                found = kTRUE;   
                break;
             }
          }
          if(!found) {
	        trackIdsTRD[nTracksTRD] = trackID;
	        trackTRDGTUtracklets[nTracksTRD] = trdTrack->GetNTracklets();
	        if ((trdTrack->GetLayerMask() & fTRDtrglayerMaskEl) != fTRDtrglayerMaskEl) trackTRDGTUlayermask[nTracksTRD] = 0;
	        else trackTRDGTUlayermask[nTracksTRD] = 1;
	        trackTRDGTUpt[nTracksTRD] =  trdTrack->Pt();
	        trackTRDGTUPID[nTracksTRD] = trdTrack->GetPID();
	        nTracksTRD++;
          }
       }
    }  // end loop over TRD tracks
  }  // end if(fFillTRDMatchedTracks)

  AliReducedEventPlaneInfo* evPlane = 0x0;  
  if(fFillEventPlaneInfo) {
     evPlane = new AliReducedEventPlaneInfo();     
     for(Int_t i=0; i<3; i++) {
       evPlane->fEventPlaneStatus[AliReducedEventPlaneInfo::kTPC][i] = AliReducedEventPlaneInfo::kRaw;
       evPlane->fEventPlaneStatus[AliReducedEventPlaneInfo::kTPCptWeights][i] = AliReducedEventPlaneInfo::kRaw;
       evPlane->fEventPlaneStatus[AliReducedEventPlaneInfo::kTPCpos][i] = AliReducedEventPlaneInfo::kRaw;
       evPlane->fEventPlaneStatus[AliReducedEventPlaneInfo::kTPCneg][i] = AliReducedEventPlaneInfo::kRaw;
       evPlane->fEventPlaneStatus[AliReducedEventPlaneInfo::kTPCsideA][i] = AliReducedEventPlaneInfo::kRaw;
       evPlane->fEventPlaneStatus[AliReducedEventPlaneInfo::kTPCsideC][i] = AliReducedEventPlaneInfo::kRaw;
     }
  }
  
  // prepare to run over the track loop
  AliESDtrack* esdTrack=0;
  AliAODTrack* aodTrack=0;
  Int_t ntracks=event->GetNumberOfTracks();
  Int_t trackId = 0; 
  Bool_t usedForV0[4] = {kFALSE}; 
  Bool_t usedForPureV0[4] = {kFALSE};
  Bool_t usedForV0Or = kFALSE;
  Float_t pileupTrackArrayP[20000];
  Float_t pileupTrackArrayM[20000];
  Int_t pileupCounterP = 0, pileupCounterM = 0;
  Float_t pileupTrackArrayP2[20000];
  Float_t pileupTrackArrayM2[20000];
  Int_t pileupCounterP2 = 0, pileupCounterM2 = 0;
  AliTPCdEdxInfo tpcdEdxInfo;
  
  AliVVertex* eventVtx = 0x0;
  if(isESD) eventVtx = const_cast<AliESDVertex*>(esdEvent->GetPrimaryVertexTracks());
  if(isAOD) eventVtx = const_cast<AliAODVertex*>(aodEvent->GetPrimaryVertex());
  
  for(Int_t itrack=0; itrack<ntracks; ++itrack){
     
    AliVParticle *particle=event->GetTrack(itrack);
    if(!particle) continue;
    if(isESD) {
      esdTrack=static_cast<AliESDtrack*>(particle);
      trackId = esdTrack->GetID();
    }
    if(isAOD) {
      aodTrack=static_cast<AliAODTrack*>(particle);
      trackId = aodTrack->GetID();
    }
        
    // check whether this track belongs to a V0 stored in the reduced event
    usedForV0Or = kFALSE;
    for(Int_t i=0; i<4; ++i) {
      usedForV0[i] = kFALSE;
      for(Int_t ii=0; ii<nV0LegsTagged[i]; ++ii) {
        if(UShort_t(trackId)==trackIdsV0[i][ii]) {
          usedForV0[i] = kTRUE;
          break;
        }
      }
      usedForV0Or = usedForV0Or || usedForV0[i];
      usedForPureV0[i] = kFALSE;
      for(Int_t ii=0; ii<nPureV0LegsTagged[i]; ++ii) {
        if(UShort_t(trackId)==trackIdsPureV0[i][ii]) {
          usedForPureV0[i] = kTRUE;
          break;
        }
      }
      usedForV0Or = usedForV0Or || usedForPureV0[i];
    }
    
    // check whether this track is matched in TRD
    Bool_t matchedInTRD = kFALSE;
    Int_t indexmatchedtrackinTRD=-1;
    if(fFillTRDMatchedTracks) {
      for(Int_t kk=0; kk<nTracksTRD; ++kk) {
        if(trackId==trackIdsTRD[kk]) {
	      matchedInTRD = kTRUE;
	      indexmatchedtrackinTRD = kk;
	      break;
        }   
      }
    }
    
    ULong_t status = (isESD ? esdTrack->GetStatus() : aodTrack->GetStatus());
    
    AliReducedEventInfo* eventInfo=0x0;
    if(fTreeWritingOption==kFullEventsWithBaseTracks || fTreeWritingOption==kFullEventsWithFullTracks) {
      eventInfo = dynamic_cast<AliReducedEventInfo*>(fReducedEvent);
      for(Int_t ibit=0; ibit<32; ++ibit) {
         if(status & (ULong_t(1)<<ibit)) {
            eventInfo->fNtracksPerTrackingFlag[ibit] += 1;
         }
      }
      if(isAOD && (status & AliVTrack::kTPCout)) eventInfo->fNtracksTPCout += 1;
      
      if(fFillEventPlaneInfo) {
         if(!fFlowTrackFilter ||
            (fFlowTrackFilter && fFlowTrackFilter->IsSelected(particle))) {
            
            Double_t x = TMath::Cos(particle->Phi());
            Double_t y = TMath::Sin(particle->Phi());
            Double_t localQvec[3][2];
            localQvec[0][0] = x; localQvec[0][1] = y;
            localQvec[1][0] = (2.0*TMath::Power(x,2.0)-1); localQvec[1][1] = 2.0*x*y;
            localQvec[2][0] = (4.0*TMath::Power(x,3.0)-3.0*x); localQvec[2][1] = (3.0*y-4.0*TMath::Power(y,3.0));
            for(Int_t iharmonic=0; iharmonic<=2; iharmonic++) {
              for(Int_t icomp=0; icomp<=1; icomp++) {
                evPlane->fQvector[AliReducedEventPlaneInfo::kTPC][iharmonic][icomp] += localQvec[iharmonic][icomp];  
                evPlane->fQvector[AliReducedEventPlaneInfo::kTPCptWeights][iharmonic][icomp] += (particle->Pt()<2.0 ? particle->Pt() : 2.0)*localQvec[iharmonic][icomp];  
                if(particle->Charge()<0.0) 
                  evPlane->fQvector[AliReducedEventPlaneInfo::kTPCneg][iharmonic][icomp] += localQvec[iharmonic][icomp];  
                if(particle->Charge()>0.0) 
                  evPlane->fQvector[AliReducedEventPlaneInfo::kTPCpos][iharmonic][icomp] += localQvec[iharmonic][icomp];  
                if(particle->Eta() < -0.5*fEventPlaneTPCetaGap)
                  evPlane->fQvector[AliReducedEventPlaneInfo::kTPCsideC][iharmonic][icomp] += localQvec[iharmonic][icomp];    
                if(particle->Eta() > 0.5*fEventPlaneTPCetaGap)
                  evPlane->fQvector[AliReducedEventPlaneInfo::kTPCsideA][iharmonic][icomp] += localQvec[iharmonic][icomp];    
              }  
            }
         }
      }
    }
    
    if((isESD && !esdTrack->IsOn(0x1)) || (!isESD && !aodTrack->IsOn(0x1))) {
      Float_t dcaXY, dcaZ;
      if(isESD) esdTrack->GetImpactParameters(dcaXY, dcaZ); 
      if(!isESD) aodTrack->GetImpactParameters(dcaXY, dcaZ);
      if(TMath::Abs(dcaXY)<3.0 && TMath::Abs(dcaZ)>4.0) {
         Double_t tgl = particle->Pz() / particle->Pt();
         if(tgl > 0.1) pileupTrackArrayP[++pileupCounterP] = (isESD ? esdTrack->GetZ() : aodTrack->GetZ());
         if(tgl < -0.1) pileupTrackArrayM[++pileupCounterM] = (isESD ? esdTrack->GetZ() : aodTrack->GetZ());
         if(TMath::Abs(dcaZ)>10.0) {
            if(tgl > 0.1) pileupTrackArrayP2[++pileupCounterP2] = (isESD ? esdTrack->GetZ() : aodTrack->GetZ());
            if(tgl < -0.1) pileupTrackArrayM2[++pileupCounterM2] = (isESD ? esdTrack->GetZ() : aodTrack->GetZ());
         }
      }
    }
    
    // fill values
    Double_t values[AliDielectronVarManager::kNMaxValues];
    // set the fill map (all 1's) for the AliDielectronVarManager
    AliDielectronVarManager::SetFillMap(fUsedVars);
    AliDielectronVarManager::Fill(particle, values);

    // decide whether to write the track in the tree
    Bool_t writeTrack = kFALSE;
    Bool_t trackFilterDecision = kFALSE;
    std::vector<Bool_t> individualFilterDecisions;
    if(fTrackFilter.GetEntries()==0) trackFilterDecision = kTRUE;
    if(fTrackFilter.GetEntries()>0) trackFilterDecision = IsTrackSelected(particle, values, individualFilterDecisions);
    if(trackFilterDecision) writeTrack = kTRUE;
    if(matchedInTRD) {
       if(fFillAllTRDMatchedTracks) writeTrack = kTRUE;
       else 
          if(trackFilterDecision) writeTrack = kTRUE;     // not needed since the track will be written anyway
    }
    if(usedForV0Or) writeTrack = kTRUE;
    if(!writeTrack) continue;
    
    Bool_t fSelectedTrackIsBaseTrack = IsSelectedTrackRequestedBaseTrack(individualFilterDecisions, usedForV0Or);
    TClonesArray& tracks = (fWriteSecondTrackArray && fSelectedTrackIsBaseTrack) ? *(fReducedEvent->fTracks2) : *(fReducedEvent->fTracks);
    AliReducedBaseTrack* reducedParticle = NULL;
    if (fSelectedTrackIsBaseTrack && fWriteSecondTrackArray)
      reducedParticle=new(tracks[fReducedEvent->NTracks2()]) AliReducedBaseTrack();
    else if (fSelectedTrackIsBaseTrack && !fWriteSecondTrackArray)
      reducedParticle=new(tracks[fReducedEvent->NTracks1()]) AliReducedBaseTrack();
    else
      reducedParticle=new(tracks[fReducedEvent->NTracks1()]) AliReducedTrackInfo();

    // increment track counters
    for(Int_t ifilter=0;ifilter<fTrackFilter.GetEntries(); ++ifilter) {
      if(individualFilterDecisions[ifilter]) {
        if(fSelectedTrackIsBaseTrack) fNSelectedBaseTracks[ifilter] += 1;
        else fNSelectedFullTracks[ifilter] += 1;
      }
    }
    for(Int_t itype=0; itype<4; ++itype) {
      if(usedForV0[itype] || usedForPureV0[itype]) {
        if(fSelectedTrackIsBaseTrack) fNSelectedBaseTracks[fTrackFilter.GetEntries()+itype] += 1;
        else fNSelectedFullTracks[fTrackFilter.GetEntries()+itype] += 1;    
      }
    }
    
    // set track quality flags
    SetTrackFilterQualityFlags(reducedParticle, individualFilterDecisions);

    reducedParticle->PtPhiEta(values[AliDielectronVarManager::kPt],values[AliDielectronVarManager::kPhi],values[AliDielectronVarManager::kEta]);
    reducedParticle->fCharge        = values[AliDielectronVarManager::kCharge];
    
    if(fFlowTrackFilter) {
      // switch on the first bit if this particle was used for the event plane
      if(fFlowTrackFilter->IsSelected(particle)) reducedParticle->fQualityFlags |= (ULong_t(1)<<0);
    }
    for(Int_t iV0type=0;iV0type<4;++iV0type) {
      if(usedForV0[iV0type]) reducedParticle->fQualityFlags |= (ULong_t(1)<<(iV0type+1));
      if(usedForPureV0[iV0type]) reducedParticle->fQualityFlags |= (ULong_t(1)<<(iV0type+8));
    }
    if(matchedInTRD) reducedParticle->fQualityFlags |= (ULong_t(1)<<26);

    if(isESD) {
      reducedParticle->fTrackId          = (UShort_t)esdTrack->GetID();
      for(Int_t idx=0; idx<3; ++idx) if(esdTrack->GetKinkIndex(idx)>0) reducedParticle->fQualityFlags |= (ULong_t(1)<<(5+idx));
      for(Int_t idx=0; idx<3; ++idx) if(esdTrack->GetKinkIndex(idx)<0) reducedParticle->fQualityFlags |= (ULong_t(1)<<(12+idx));
      if(((AliESDVertex*)eventVtx)->UsesTrack(esdTrack->GetID())) reducedParticle->fQualityFlags |= (ULong_t(1)<<27);
   }
   if(isAOD) {
      reducedParticle->fTrackId = aodTrack->GetID();
      for(Int_t idx=0; idx<3; ++idx) if(aodTrack->GetKinkIndex(idx)>0) reducedParticle->fQualityFlags |= (ULong_t(1)<<(5+idx));
      for(Int_t idx=0; idx<3; ++idx) if(aodTrack->GetKinkIndex(idx)<0) reducedParticle->fQualityFlags |= (ULong_t(1)<<(12+idx));
      for(Int_t idx=0; idx<11; ++idx) if(aodTrack->TestFilterBit(BIT(idx))) reducedParticle->SetQualityFlag(15+idx);
      if(((AliAODVertex*)eventVtx)->HasDaughter(aodTrack)) reducedParticle->fQualityFlags |= (ULong_t(1)<<27);
   }
   
    // If we want to write only AliReducedBaseTrack objects, then we stop here
    if(fSelectedTrackIsBaseTrack) {
       if(fFillMCInfo && hasMC) {
          AliVParticle* mcTruth = AliDielectronMC::Instance()->GetMCTrack(particle);
          if(mcTruth)
             reducedParticle->fMCFlags = MatchMCsignals(mcTruth->GetLabel());    // check which MC signals match this particle
       }
       
      fReducedEvent->fNtracks[1] += 1;
      continue;
    }
    
    AliReducedTrackInfo* trackInfo = dynamic_cast<AliReducedTrackInfo*>(reducedParticle);
    if(!trackInfo) continue;
    
    trackInfo->fStatus        = status;
    trackInfo->fMomentumInner = values[AliDielectronVarManager::kPIn];
    trackInfo->fDCA[0]        = values[AliDielectronVarManager::kImpactParXY];
    trackInfo->fDCA[1]        = values[AliDielectronVarManager::kImpactParZ];
    trackInfo->fTrackLength   = values[AliDielectronVarManager::kTrackLength];
    
    trackInfo->fITSclusterMap = (UChar_t)values[AliDielectronVarManager::kITSclusterMap];
    trackInfo->fITSsignal     = values[AliDielectronVarManager::kITSsignal];
    trackInfo->fITSnSig[0]    = values[AliDielectronVarManager::kITSnSigmaEle];
    trackInfo->fITSnSig[1]    = values[AliDielectronVarManager::kITSnSigmaPio];
    trackInfo->fITSnSig[2]    = values[AliDielectronVarManager::kITSnSigmaKao];
    trackInfo->fITSnSig[3]    = values[AliDielectronVarManager::kITSnSigmaPro];
    trackInfo->fITSchi2       = values[AliDielectronVarManager::kITSchi2Cl];
    
    trackInfo->fTPCNcls      = (UChar_t)values[AliDielectronVarManager::kNclsTPC];
    trackInfo->fTPCNclsF     = (UChar_t)values[AliDielectronVarManager::kNFclsTPC];
    trackInfo->fTPCNclsShared = (UChar_t)values[AliDielectronVarManager::kNclsSTPC];
    trackInfo->fTPCCrossedRows = values[AliDielectronVarManager::kNFclsTPCr];
    trackInfo->fTPCsignal    = values[AliDielectronVarManager::kTPCsignal];
    trackInfo->fTPCsignalN   = values[AliDielectronVarManager::kTPCsignalN];
    trackInfo->fTPCnSig[0]   = values[AliDielectronVarManager::kTPCnSigmaEle];
    trackInfo->fTPCnSig[1]   = values[AliDielectronVarManager::kTPCnSigmaPio];
    trackInfo->fTPCnSig[2]   = values[AliDielectronVarManager::kTPCnSigmaKao];
    trackInfo->fTPCnSig[3]   = values[AliDielectronVarManager::kTPCnSigmaPro];
    trackInfo->fTPCClusterMap = EncodeTPCClusterMap(particle, isAOD);
    trackInfo->fTPCchi2       = values[AliDielectronVarManager::kTPCchi2Cl];
    trackInfo->fTPCActiveLength = values[AliDielectronVarManager::kTPCActiveLength];
    trackInfo->fTPCGeomLength = values[AliDielectronVarManager::kTPCGeomLength];
        
    trackInfo->fTOFbeta      = values[AliDielectronVarManager::kTOFbeta];
    trackInfo->fTOFtime      = values[AliDielectronVarManager::kTOFsignal]-pidResponse->GetTOFResponse().GetTimeZero();
    trackInfo->fTOFmismatchProbab = values[AliDielectronVarManager::kTOFmismProb];
    trackInfo->fTOFnSig[0]   = values[AliDielectronVarManager::kTOFnSigmaEle];
    trackInfo->fTOFnSig[1]   = values[AliDielectronVarManager::kTOFnSigmaPio];
    trackInfo->fTOFnSig[2]   = values[AliDielectronVarManager::kTOFnSigmaKao];
    trackInfo->fTOFnSig[3]   = values[AliDielectronVarManager::kTOFnSigmaPro];
    
    trackInfo->fEMCALnSigEle = values[AliDielectronVarManager::kEMCALnSigmaEle];

    Double_t trdProbab[AliPID::kSPECIES]={0.0};
    if(isESD) {
       trackInfo->fMassForTracking = esdTrack->GetMassForTracking();
       
       AliESDEvent* esdEvent = static_cast<AliESDEvent*>(InputEvent());
       AliESDVertex* eventVtx = const_cast<AliESDVertex*>(esdEvent->GetPrimaryVertexTracks());

       AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
       TClass* esdClass = esdTrack->Class();
       if(esdClass->GetMethodAny("GetChi2TPCConstrainedVsGlobal") && fld)
         trackInfo->fChi2TPCConstrainedVsGlobal = esdTrack->GetChi2TPCConstrainedVsGlobal(eventVtx);
       
      const AliExternalTrackParam* tpcInner = esdTrack->GetTPCInnerParam();

      for(Int_t i=0; i<6; ++i) {
         if(esdTrack->HasSharedPointOnITSLayer(i)) trackInfo->fITSSharedClusterMap |= (1<<i);
      }
      
      Float_t xyDCA,zDCA;
      Double_t helixinfo[6];
      if(tpcInner){
        trackInfo->fTPCPhi        = (tpcInner ? tpcInner->Phi() : 0.0);
        trackInfo->fTPCPt         = (tpcInner ? tpcInner->Pt() : 0.0);
        trackInfo->fTPCEta        = (tpcInner ? tpcInner->Eta() : 0.0);
        esdTrack->GetImpactParametersTPC(xyDCA,zDCA);
        trackInfo->fTPCDCA[0]     = xyDCA;
        trackInfo->fTPCDCA[1]     = zDCA;
        
        // helix information (Alex Chauvin)
        tpcInner->GetHelixParameters(helixinfo,InputEvent()->GetMagneticField());
        if(helixinfo[2] < 0) helixinfo[2] = helixinfo[2] + 2*TMath::Pi();
        helixinfo[2] -= TMath::Pi()/2.;
        trackInfo->fHelixCenter[0]= helixinfo[5]+(TMath::Cos(helixinfo[2])*TMath::Abs(1./helixinfo[4])*copysignf(1.0, InputEvent()->GetMagneticField()*values[AliDielectronVarManager::kCharge]));
        trackInfo->fHelixCenter[1]= helixinfo[0]+(TMath::Sin(helixinfo[2])*TMath::Abs(1./helixinfo[4])*copysignf(1.0, InputEvent()->GetMagneticField()*values[AliDielectronVarManager::kCharge]));
        trackInfo->fHelixRadius   = TMath::Abs(1./helixinfo[4]);
      }
      
      if(esdTrack->GetTPCdEdxInfo(tpcdEdxInfo)) {
         for(Int_t i=0;i<4;++i) {
            trackInfo->fTPCdEdxInfoQmax[i] = tpcdEdxInfo.GetSignalMax(i);
            trackInfo->fTPCdEdxInfoQtot[i] = tpcdEdxInfo.GetSignalTot(i);
         }
      }
      
      trackInfo->fTOFdeltaBC    = esdTrack->GetTOFDeltaBC();
      trackInfo->fTOFdx         = esdTrack->GetTOFsignalDx();
      trackInfo->fTOFdz         = esdTrack->GetTOFsignalDz();
      trackInfo->fTOFchi2       = esdTrack->GetTOFchi2();
      
      trackInfo->fTRDntracklets[0] = esdTrack->GetTRDntracklets();
      trackInfo->fTRDntracklets[1] = esdTrack->GetTRDntrackletsPID();
      pidResponse->ComputeTRDProbability(esdTrack,AliPID::kSPECIES,trdProbab,AliTRDPIDResponse::kLQ1D);
      trackInfo->fTRDpid[0]    = trdProbab[AliPID::kElectron];
      trackInfo->fTRDpid[1]    = trdProbab[AliPID::kPion];
      pidResponse->ComputeTRDProbability(esdTrack,AliPID::kSPECIES,trdProbab,AliTRDPIDResponse::kLQ2D);
      trackInfo->fTRDpidLQ2D[0]    = trdProbab[AliPID::kElectron];
      trackInfo->fTRDpidLQ2D[1]    = trdProbab[AliPID::kPion];

      if(fFillTRDMatchedTracks && (indexmatchedtrackinTRD!=-1)) {
	    const Int_t indexTRD         = indexmatchedtrackinTRD;
	    trackInfo->fTRDGTUtracklets  = trackTRDGTUtracklets[indexTRD];
	    trackInfo->fTRDGTUlayermask  = trackTRDGTUlayermask[indexTRD];
	    trackInfo->fTRDGTUpt         = trackTRDGTUpt[indexTRD];
	    trackInfo->fTRDGTUsagitta    = trackTRDGTUsagitta[indexTRD];
	    trackInfo->fTRDGTUPID        = trackTRDGTUPID[indexTRD];
      }

      if(esdTrack->IsEMCAL()) trackInfo->fCaloClusterId = esdTrack->GetEMCALcluster();
      if(esdTrack->IsPHOS()) trackInfo->fCaloClusterId = esdTrack->GetPHOScluster();
      // NOTE: extrapolation depends on radius and PHOS radius differs slightly from EMCal radius
      if (esdTrack->IsExtrapolatedToEMCAL()) trackInfo->fMomentumOnCalo = esdTrack->GetTrackPOnEMCal();
      if (esdTrack->IsExtrapolatedToEMCAL()) trackInfo->fPhiOnCalo = esdTrack->GetTrackPhiOnEMCal();
      if (esdTrack->IsExtrapolatedToEMCAL()) trackInfo->fEtaOnCalo = esdTrack->GetTrackEtaOnEMCal();

      Double_t xyz[3], pxpypz[3];
      Double_t covMat[21];
      esdTrack->GetXYZ(xyz);
      esdTrack->GetPxPyPz(pxpypz);
      esdTrack->GetCovarianceXYZPxPyPz(covMat);
      for(Int_t i=0;i<3;++i) {
        trackInfo->fTrackParam[i] = xyz[i];
        trackInfo->fTrackParam[i+3] = pxpypz[i];
      }
      for(Int_t i=0;i<21;++i) {
        trackInfo->fCovMatrix[i] = covMat[i];
      }
            
      if(fFillMCInfo && hasMC) {
         AliMCParticle* truthParticle = AliDielectronMC::Instance()->GetMCTrack(esdTrack);
         if(truthParticle) {
           trackInfo->fMCFlags = MatchMCsignals(truthParticle->GetLabel());    // check which MC signals match this particle and fill the bit map
                      
           trackInfo->fMCMom[0] = truthParticle->Px();
           trackInfo->fMCMom[1] = truthParticle->Py();
           trackInfo->fMCMom[2] = truthParticle->Pz();
           trackInfo->fMCFreezeout[0] = truthParticle->Xv();
           trackInfo->fMCFreezeout[1] = truthParticle->Yv();
           trackInfo->fMCFreezeout[2] = truthParticle->Zv();
           trackInfo->fMCLabels[0] = esdTrack->GetLabel();
           trackInfo->fMCPdg[0] = truthParticle->PdgCode();
           trackInfo->fMCGeneratorIndex = truthParticle->GetGeneratorIndex();
           
           AliMCParticle* motherTruth = AliDielectronMC::Instance()->GetMCTrackMother(truthParticle);
           if(motherTruth) {
             trackInfo->fMCLabels[1] = truthParticle->GetMother();
             trackInfo->fMCPdg[1] = motherTruth->PdgCode();
          }
          
          AliMCParticle* grandmotherTruth = NULL;
          if(motherTruth) grandmotherTruth = AliDielectronMC::Instance()->GetMCTrackMother(motherTruth);
          if(grandmotherTruth) {
             trackInfo->fMCLabels[2] = motherTruth->GetMother();
             trackInfo->fMCPdg[2] = grandmotherTruth->PdgCode();
          }
           
          AliMCParticle* grandgrandmotherTruth = NULL;
          if(grandmotherTruth) grandgrandmotherTruth = AliDielectronMC::Instance()->GetMCTrackMother(grandmotherTruth);
          if(grandgrandmotherTruth) {
             trackInfo->fMCLabels[3] = grandmotherTruth->GetMother();
             trackInfo->fMCPdg[3] = grandgrandmotherTruth->PdgCode();
          }

	      if(fFillHFInfo)      trackInfo->fHFProc = AliDielectronMC::Instance()->GetHFProcess(truthParticle->GetLabel());
	    }
      }
    }  // end if(isESD)
    if(isAOD) {
      trackInfo->fMassForTracking = aodTrack->GetMassForTracking();
      trackInfo->fChi2TPCConstrainedVsGlobal = aodTrack->GetChi2TPCConstrainedVsGlobal(); 
      
      //trackInfo->fITSSharedClusterMap = aodTrack->GetITSSharedClusterMap();
      for(Int_t i=0; i<6; ++i) {
         if(aodTrack->HasSharedPointOnITSLayer(i)) trackInfo->fITSSharedClusterMap |= (1<<i);
      }
      
      const AliExternalTrackParam* tpcInner = aodTrack->GetInnerParam();
      Float_t xyDCA,zDCA;
      Double_t helixinfo[6];
      if(tpcInner){
        trackInfo->fTPCPhi        = (tpcInner ? tpcInner->Phi() : 0.0);
        trackInfo->fTPCPt         = (tpcInner ? tpcInner->Pt() : 0.0);
        trackInfo->fTPCEta        = (tpcInner ? tpcInner->Eta() : 0.0);
      
        aodTrack->GetImpactParametersTPC(xyDCA,zDCA);
        trackInfo->fTPCDCA[0]     = xyDCA;
        trackInfo->fTPCDCA[1]     = zDCA;
      
        // helix information (Alex Chauvin)
        tpcInner->GetHelixParameters(helixinfo,InputEvent()->GetMagneticField());
        if(helixinfo[2] < 0) helixinfo[2] = helixinfo[2] + 2*TMath::Pi();
        helixinfo[2] -= TMath::Pi()/2.;
        trackInfo->fHelixCenter[0]= helixinfo[5]+(TMath::Cos(helixinfo[2])*TMath::Abs(1./helixinfo[4])*copysignf(1.0, InputEvent()->GetMagneticField()*values[AliDielectronVarManager::kCharge]));
        trackInfo->fHelixCenter[1]= helixinfo[0]+(TMath::Sin(helixinfo[2])*TMath::Abs(1./helixinfo[4])*copysignf(1.0, InputEvent()->GetMagneticField()*values[AliDielectronVarManager::kCharge]));
        trackInfo->fHelixRadius   = TMath::Abs(1./helixinfo[4]);
      }
      
      if(aodTrack->GetTPCdEdxInfo(tpcdEdxInfo)) {
         for(Int_t i=0;i<4;++i) {
            trackInfo->fTPCdEdxInfoQmax[i] = tpcdEdxInfo.GetSignalMax(i);
            trackInfo->fTPCdEdxInfoQtot[i] = tpcdEdxInfo.GetSignalTot(i);
         }
      }
      
      trackInfo->fTOFdz         = aodTrack->GetTOFsignalDz();
      trackInfo->fTOFdeltaBC = eventInfo->fBC - aodTrack->GetTOFBunchCrossing();
      
      trackInfo->fTRDntracklets[0] = aodTrack->GetTRDntrackletsPID();
      trackInfo->fTRDntracklets[1] = aodTrack->GetTRDntrackletsPID();
      pidResponse->ComputeTRDProbability(aodTrack,AliPID::kSPECIES,trdProbab,AliTRDPIDResponse::kLQ1D);
      trackInfo->fTRDpid[0]    = trdProbab[AliPID::kElectron];
      trackInfo->fTRDpid[1]    = trdProbab[AliPID::kPion];
      pidResponse->ComputeTRDProbability(aodTrack,AliPID::kSPECIES,trdProbab,AliTRDPIDResponse::kLQ2D);
      trackInfo->fTRDpidLQ2D[0]    = trdProbab[AliPID::kElectron];
      trackInfo->fTRDpidLQ2D[1]    = trdProbab[AliPID::kPion];
      
      if(aodTrack->IsEMCAL()) trackInfo->fCaloClusterId = aodTrack->GetEMCALcluster();
      if(aodTrack->IsPHOS()) trackInfo->fCaloClusterId = aodTrack->GetPHOScluster();
      // NOTE: extrapolation depends on radius and PHOS radius differs slightly from EMCal radius
      if (aodTrack->IsExtrapolatedToEMCAL()) trackInfo->fMomentumOnCalo = aodTrack->GetTrackPOnEMCal();
      if (aodTrack->IsExtrapolatedToEMCAL()) trackInfo->fPhiOnCalo = aodTrack->GetTrackPhiOnEMCal();
      if (aodTrack->IsExtrapolatedToEMCAL()) trackInfo->fEtaOnCalo = aodTrack->GetTrackEtaOnEMCal();

      Double_t xyz[3], pxpypz[3];
      Double_t covMat[21];
      aodTrack->GetXYZ(xyz);
      aodTrack->GetPxPyPz(pxpypz);
      aodTrack->GetCovarianceXYZPxPyPz(covMat);
      for(Int_t i=0;i<3;++i) {
        trackInfo->fTrackParam[i] = xyz[i];
        trackInfo->fTrackParam[i+3] = pxpypz[i];
      }
      for(Int_t i=0;i<21;++i) {
        trackInfo->fCovMatrix[i] = covMat[i];
      }
        
      if(fFillMCInfo && hasMC) {
         AliAODMCParticle* truthParticle = AliDielectronMC::Instance()->GetMCTrack(aodTrack);
         if(truthParticle) {
            trackInfo->fMCFlags = MatchMCsignals(aodTrack->GetLabel());    // check which MC signals match this particle and fill the bit map
            
            trackInfo->fMCMom[0] = truthParticle->Px();
            trackInfo->fMCMom[1] = truthParticle->Py();
            trackInfo->fMCMom[2] = truthParticle->Pz();
            trackInfo->fMCFreezeout[0] = truthParticle->Xv();
            trackInfo->fMCFreezeout[1] = truthParticle->Yv();
            trackInfo->fMCFreezeout[2] = truthParticle->Zv();
            trackInfo->fMCLabels[0] = aodTrack->GetLabel();
            trackInfo->fMCPdg[0] = truthParticle->PdgCode();
            trackInfo->fMCGeneratorIndex = truthParticle->GetGeneratorIndex();
            
            AliAODMCParticle* motherTruth = AliDielectronMC::Instance()->GetMCTrackMother(truthParticle);
            if(motherTruth) {
               trackInfo->fMCLabels[1] = truthParticle->GetMother();
               trackInfo->fMCPdg[1] = motherTruth->PdgCode();
            }
            
            AliAODMCParticle* grandmotherTruth = NULL;
            if(motherTruth) grandmotherTruth = AliDielectronMC::Instance()->GetMCTrackMother(motherTruth);
            if(grandmotherTruth) {
               trackInfo->fMCLabels[2] = motherTruth->GetMother();
               trackInfo->fMCPdg[2] = grandmotherTruth->PdgCode();
            }
            
            AliAODMCParticle* grandgrandmotherTruth = NULL;
            if(grandmotherTruth) grandgrandmotherTruth = AliDielectronMC::Instance()->GetMCTrackMother(grandmotherTruth);
            if(grandgrandmotherTruth) {
               trackInfo->fMCLabels[3] = grandmotherTruth->GetMother();
               trackInfo->fMCPdg[3] = grandgrandmotherTruth->PdgCode();
            }
         }
      }
    }  // end if(isAOD)

    fReducedEvent->fNtracks[1] += 1;
  }  // end loop over tracks
    
  AliReducedEventInfo* eventInfo = NULL; 
  if(fTreeWritingOption==kFullEventsWithBaseTracks || fTreeWritingOption==kFullEventsWithFullTracks) {
     eventInfo = dynamic_cast<AliReducedEventInfo*>(fReducedEvent);
     eventInfo->SetEventPlane(evPlane);
     eventInfo->fTPCpileupZ[0] = (pileupCounterP>0 ? TMath::Median(pileupCounterP, pileupTrackArrayP) : 0.0);
     eventInfo->fTPCpileupZ[1] = (pileupCounterM>0 ? -1.0*TMath::Median(pileupCounterM, pileupTrackArrayM) : 0.0);
     eventInfo->fTPCpileupContributors[0] = pileupCounterP;
     eventInfo->fTPCpileupContributors[1] = pileupCounterM;
     eventInfo->fTPCpileupZ2[0] = (pileupCounterP2>0 ? TMath::Median(pileupCounterP2, pileupTrackArrayP2) : 0.0);
     eventInfo->fTPCpileupZ2[1] = (pileupCounterM2>0 ? -1.0*TMath::Median(pileupCounterM2, pileupTrackArrayM2) : 0.0);
     eventInfo->fTPCpileupContributors2[0] = pileupCounterP2;
     eventInfo->fTPCpileupContributors2[1] = pileupCounterM2;
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillV0PairInfo() 
{
  //
  // fill reduced pair information
  //
  AliESDEvent* esd = (AliESDEvent*)InputEvent();
  const AliESDVertex *primaryVertex = esd->GetPrimaryVertex();
  AliKFVertex primaryVertexKF(*primaryVertex);
  
  fReducedEvent->fNV0candidates[0] = InputEvent()->GetNumberOfV0s();
  
  if(!(fFillK0s || fFillLambda || fFillALambda || fFillGammaConversions)) return;
    
  if(fV0OpenCuts) {
    fV0OpenCuts->SetEvent(esd);
    fV0OpenCuts->SetPrimaryVertex(&primaryVertexKF);
  }
  if(fV0StrongCuts) {
    fV0StrongCuts->SetEvent(esd);
    fV0StrongCuts->SetPrimaryVertex(&primaryVertexKF);
  }
  
  Int_t pdgV0=0; Int_t pdgP=0; Int_t pdgN=0;
  for(Int_t iV0=0; iV0<InputEvent()->GetNumberOfV0s(); ++iV0) {   // loop over V0s
    AliESDv0 *v0 = esd->GetV0(iV0);
       
    AliESDtrack* legPos = esd->GetTrack(v0->GetPindex());
    AliESDtrack* legNeg = esd->GetTrack(v0->GetNindex());
 
    if(legPos->GetSign() == legNeg->GetSign()) continue;

    Bool_t v0ChargesAreCorrect = (legPos->GetSign()==+1 ? kTRUE : kFALSE);
    legPos = (!v0ChargesAreCorrect ? esd->GetTrack(v0->GetNindex()) : legPos);
    legNeg = (!v0ChargesAreCorrect ? esd->GetTrack(v0->GetPindex()) : legNeg); 
    
    pdgV0=0; pdgP=0; pdgN=0;
    Bool_t goodK0s = kFALSE; Bool_t goodLambda = kFALSE; Bool_t goodALambda = kFALSE; Bool_t goodGamma = kFALSE;
    if(fV0OpenCuts) {
      Bool_t processV0 = fV0OpenCuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
      if(processV0 && TMath::Abs(pdgV0)==310 && TMath::Abs(pdgP)==211 && TMath::Abs(pdgN)==211) {
        goodK0s = kTRUE;
        if(fK0sPionCuts && (!fK0sPionCuts->IsSelected(legPos) || !fK0sPionCuts->IsSelected(legNeg))) goodK0s = kFALSE;
      }
      if(processV0 && pdgV0==3122 && (TMath::Abs(pdgP)==211 || TMath::Abs(pdgP)==2212) && (TMath::Abs(pdgN)==211 || TMath::Abs(pdgN)==2212)) {
        goodLambda = kTRUE;
        if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(legPos)) goodLambda = kFALSE;
        if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(legNeg)) goodLambda = kFALSE;
      }
      if(processV0 && pdgV0==-3122 && (TMath::Abs(pdgP)==211 || TMath::Abs(pdgP)==2212) && (TMath::Abs(pdgN)==211 || TMath::Abs(pdgN)==2212)) {
        goodALambda = kTRUE;
        if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(legNeg)) goodALambda = kFALSE;
        if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(legPos)) goodALambda = kFALSE;
      }
      if(processV0 && TMath::Abs(pdgV0)==22 && TMath::Abs(pdgP)==11 && TMath::Abs(pdgN)==11) {
        goodGamma = kTRUE;
        if(fGammaElectronCuts && (!fGammaElectronCuts->IsSelected(legPos) || !fGammaElectronCuts->IsSelected(legNeg))) goodGamma = kFALSE;
      }
    }
    
    Bool_t veryGoodK0s = kFALSE; Bool_t veryGoodLambda = kFALSE; Bool_t veryGoodALambda = kFALSE; Bool_t veryGoodGamma = kFALSE;
    if(fV0StrongCuts && (goodK0s || goodLambda || goodALambda || goodGamma)) {
      pdgV0=0; pdgP=0; pdgN=0;
      Bool_t processV0 = fV0StrongCuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
      if(processV0 && goodK0s && TMath::Abs(pdgV0)==310 && TMath::Abs(pdgP)==211 && TMath::Abs(pdgN)==211)
        veryGoodK0s = kTRUE;
      if(processV0 && goodLambda && pdgV0==3122 && (TMath::Abs(pdgP)==211 || TMath::Abs(pdgP)==2212) && (TMath::Abs(pdgN)==211 || TMath::Abs(pdgN)==2212))
        veryGoodLambda = kTRUE;
      if(processV0 && goodALambda && pdgV0==-3122 && (TMath::Abs(pdgP)==211 || TMath::Abs(pdgP)==2212) && (TMath::Abs(pdgN)==211 || TMath::Abs(pdgN)==2212))
        veryGoodALambda = kTRUE;
      if(processV0 && goodGamma && TMath::Abs(pdgV0)==22 && TMath::Abs(pdgP)==11 && TMath::Abs(pdgN)==11)
        veryGoodGamma = kTRUE;
    }
              
    if(!((goodK0s && fFillK0s) || 
         (goodLambda && fFillLambda) || 
         (goodALambda && fFillALambda) || 
         (goodGamma && fFillGammaConversions))) continue;
    
    // Fill the V0 information into the tree for 4 hypothesis: K0s, Lambda, Anti-Lambda and gamma conversion
    AliReducedPairInfo* k0sReducedPair     = FillV0PairInfo(v0, AliReducedPairInfo::kK0sToPiPi,     legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    AliReducedPairInfo* lambdaReducedPair  = FillV0PairInfo(v0, AliReducedPairInfo::kLambda0ToPPi,  legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    AliReducedPairInfo* alambdaReducedPair = FillV0PairInfo(v0, AliReducedPairInfo::kALambda0ToPPi, legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    AliReducedPairInfo* gammaReducedPair   = FillV0PairInfo(v0, AliReducedPairInfo::kGammaConv,     legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    
    if(fFillK0s && goodK0s && k0sReducedPair->fMass[0]>fK0sMassRange[0] && k0sReducedPair->fMass[0]<fK0sMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPairInfo *goodK0sPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPairInfo(*k0sReducedPair);
      goodK0sPair->fMass[0] = k0sReducedPair->fMass[0];
      goodK0sPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodK0sPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodK0sPair->fMass[3] = gammaReducedPair->fMass[0];
      if(veryGoodK0s) goodK0sPair->fQualityFlags |= (ULong_t(1)<<1);
      fReducedEvent->fNV0candidates[1] += 1;
      if(v0->GetOnFlyStatus()) fNSelectedFullTracks[fTrackFilter.GetEntries()+4+1] += 1;
      else fNSelectedFullTracks[fTrackFilter.GetEntries()+4+5] += 1;
    } else {goodK0s=kFALSE;}
    if(fFillLambda && goodLambda && lambdaReducedPair->fMass[0]>fLambdaMassRange[0] && lambdaReducedPair->fMass[0]<fLambdaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPairInfo *goodLambdaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPairInfo(*lambdaReducedPair);
      goodLambdaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodLambdaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodLambdaPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodLambdaPair->fMass[3] = gammaReducedPair->fMass[0];
      if(veryGoodLambda) goodLambdaPair->fQualityFlags |= (ULong_t(1)<<2);
      fReducedEvent->fNV0candidates[1] += 1;
      if(v0->GetOnFlyStatus()) fNSelectedFullTracks[fTrackFilter.GetEntries()+4+2] += 1;
      else fNSelectedFullTracks[fTrackFilter.GetEntries()+4+6] += 1;
    } else {goodLambda=kFALSE;}
    if(fFillALambda && goodALambda && alambdaReducedPair->fMass[0]>fLambdaMassRange[0] && alambdaReducedPair->fMass[0]<fLambdaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPairInfo *goodALambdaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPairInfo(*alambdaReducedPair);
      goodALambdaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodALambdaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodALambdaPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodALambdaPair->fMass[3] = gammaReducedPair->fMass[0];
      if(veryGoodALambda) goodALambdaPair->fQualityFlags |= (ULong_t(1)<<3);
      fReducedEvent->fNV0candidates[1] += 1;
      if(v0->GetOnFlyStatus()) fNSelectedFullTracks[fTrackFilter.GetEntries()+4+3] += 1;
      else fNSelectedFullTracks[fTrackFilter.GetEntries()+4+7] += 1;
    } else {goodALambda = kFALSE;}
    if(fFillGammaConversions && goodGamma && gammaReducedPair->fMass[0]>fGammaMassRange[0] && gammaReducedPair->fMass[0]<fGammaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPairInfo *goodGammaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPairInfo(*gammaReducedPair);
      goodGammaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodGammaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodGammaPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodGammaPair->fMass[3] = gammaReducedPair->fMass[0];
      if(veryGoodGamma) goodGammaPair->fQualityFlags |= (ULong_t(1)<<4);
      fReducedEvent->fNV0candidates[1] += 1;
      if(v0->GetOnFlyStatus()) fNSelectedFullTracks[fTrackFilter.GetEntries()+4+4] += 1;
      else fNSelectedFullTracks[fTrackFilter.GetEntries()+4+8] += 1;
    } else {goodGamma=kFALSE;}
    delete k0sReducedPair;
    delete lambdaReducedPair;
    delete alambdaReducedPair;
    delete gammaReducedPair;
  }   // end loop over V0s
}

//_________________________________________________________________________________
AliReducedPairInfo* AliAnalysisTaskReducedTreeMaker::FillV0PairInfo(AliESDv0* v0, Int_t id, 
						    AliESDtrack* legPos, AliESDtrack* legNeg,
						    AliKFVertex* vtxKF, Bool_t chargesAreCorrect) {
  //
  // Create a reduced V0 object and fill it
  //
  AliReducedPairInfo* reducedPair=new AliReducedPairInfo();  
  reducedPair->fCandidateId = id;
  reducedPair->fPairType    = v0->GetOnFlyStatus();    // on the fly status
  reducedPair->fLegIds[0]   = legPos->GetID();
  reducedPair->fLegIds[1]   = legNeg->GetID();
  if(!reducedPair->fPairType) {    // offline
    UInt_t pidPos = AliPID::kPion;
    if(id==AliReducedPairInfo::kLambda0ToPPi) pidPos = AliPID::kProton;
    if(id==AliReducedPairInfo::kGammaConv) pidPos = AliPID::kElectron;
    UInt_t pidNeg = AliPID::kPion;
    if(id==AliReducedPairInfo::kALambda0ToPPi) pidNeg = AliPID::kProton;
    if(id==AliReducedPairInfo::kGammaConv) pidNeg = AliPID::kElectron;
    reducedPair->fMass[0]      = v0->GetEffMass(pidPos, pidNeg);
    reducedPair->fIsCartesian  = kFALSE;
    reducedPair->fP[1]         = v0->Phi();
    if(reducedPair->fP[1]<0.0) reducedPair->fP[1] = 2.0*TMath::Pi() + reducedPair->fP[1];  // converted to [0,2pi]
    reducedPair->fP[0]         = v0->Pt();
    reducedPair->fP[2]         = v0->Eta();
    reducedPair->fLxy          = v0->GetRr();
    reducedPair->fPointingAngle = v0->GetV0CosineOfPointingAngle(vtxKF->GetX(), vtxKF->GetY(), vtxKF->GetZ());
    reducedPair->fChisquare    = v0->GetChi2V0();
  }
  else {
    const AliExternalTrackParam *negHelix=v0->GetParamN();
    const AliExternalTrackParam *posHelix=v0->GetParamP();
    if(!chargesAreCorrect) {
      negHelix = v0->GetParamP();
      posHelix = v0->GetParamN();
    }
    Int_t pdgPos = 211;
    if(id==AliReducedPairInfo::kLambda0ToPPi) pdgPos = 2212;
    if(id==AliReducedPairInfo::kGammaConv) pdgPos = -11;
    Int_t pdgNeg = -211;
    if(id==AliReducedPairInfo::kALambda0ToPPi) pdgNeg = -2212;
    if(id==AliReducedPairInfo::kGammaConv) pdgNeg = 11;
    AliKFParticle negKF(*(negHelix), pdgPos);
    AliKFParticle posKF(*(posHelix), pdgNeg);
    AliKFParticle v0Refit;
    v0Refit += negKF;
    v0Refit += posKF;
    Double_t massFit=0.0, massErrFit=0.0;
    v0Refit.GetMass(massFit,massErrFit);
    reducedPair->fMass[0] = massFit;
    reducedPair->fIsCartesian  = kFALSE;
    reducedPair->fP[1]         = v0Refit.GetPhi();
    if(reducedPair->fP[1]<0.0) reducedPair->fP[1] = 2.0*TMath::Pi() + reducedPair->fP[1];  // converted to [0,2pi]
    reducedPair->fP[0]         = v0Refit.GetPt();
    reducedPair->fP[2]         = v0Refit.GetEta();
    reducedPair->fLxy          = v0Refit.GetPseudoProperDecayTime(*vtxKF, massFit);
    Double_t deltaPos[3];
    deltaPos[0] = v0Refit.GetX() - vtxKF->GetX(); deltaPos[1] = v0Refit.GetY() - vtxKF->GetY(); deltaPos[2] = v0Refit.GetZ() - vtxKF->GetZ();
    Double_t momV02 = v0Refit.GetPx()*v0Refit.GetPx() + v0Refit.GetPy()*v0Refit.GetPy() + v0Refit.GetPz()*v0Refit.GetPz();
    Double_t deltaPos2 = deltaPos[0]*deltaPos[0] + deltaPos[1]*deltaPos[1] + deltaPos[2]*deltaPos[2];
    reducedPair->fPointingAngle = (deltaPos[0]*v0Refit.GetPx() + deltaPos[1]*v0Refit.GetPy() + deltaPos[2]*v0Refit.GetPz()) / 
                                  TMath::Sqrt(momV02*deltaPos2);
    reducedPair->fChisquare = v0Refit.GetChi2();                              
  }
  return reducedPair;
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillV0PairInfoAOD() 
{
   //
   // fill reduced pair information
   //
   AliAODEvent* aod = (AliAODEvent*)InputEvent();
   const AliAODVertex *primaryVertex = aod->GetPrimaryVertex();
   AliKFVertex primaryVertexKF(*primaryVertex);
   
   fReducedEvent->fNV0candidates[0] = InputEvent()->GetNumberOfV0s();
   
   if(!(fFillK0s || fFillLambda || fFillALambda || fFillGammaConversions)) return;
   
   if(fV0CutsAOD) {
      fV0CutsAOD->SetEvent(aod);
      fV0CutsAOD->SetPrimaryVertex(&primaryVertexKF);
   }
   
   Int_t pdgV0=0; Int_t pdgP=0; Int_t pdgN=0;
   for(Int_t iV0=0; iV0<InputEvent()->GetNumberOfV0s(); ++iV0) {   // loop over V0s
      AliAODv0 *v0 = aod->GetV0(iV0);
      if(!v0) continue;
      
      if(!aod->GetTrack(v0->GetPosID())) continue;
      AliAODTrack* legPos = dynamic_cast<AliAODTrack*>(aod->GetTrack(v0->GetPosID()));
      if(!legPos) continue;
      if(!aod->GetTrack(v0->GetNegID())) continue;
      AliAODTrack* legNeg = dynamic_cast<AliAODTrack*>(aod->GetTrack(v0->GetNegID()));
      if(!legNeg) continue;
      
      pdgV0=0; pdgP=0; pdgN=0;
      Bool_t goodK0s = kTRUE; Bool_t goodLambda = kTRUE; Bool_t goodALambda = kTRUE; Bool_t goodGamma = kTRUE;
      if(fV0CutsAOD) {
         goodK0s = kFALSE; goodLambda = kFALSE; goodALambda = kFALSE; goodGamma = kFALSE;
         Bool_t processV0 = fV0CutsAOD->ProcessV0(v0, pdgV0, pdgP, pdgN);
         if(processV0 && TMath::Abs(pdgV0)==310 && TMath::Abs(pdgP)==211 && TMath::Abs(pdgN)==211) {
            goodK0s = kTRUE;
            if(fK0sPionCuts && (!fK0sPionCuts->IsSelected(legPos) || !fK0sPionCuts->IsSelected(legNeg))) goodK0s = kFALSE;
         }
         if(processV0 && pdgV0==3122 && (TMath::Abs(pdgP)==211 || TMath::Abs(pdgP)==2212) && (TMath::Abs(pdgN)==211 || TMath::Abs(pdgN)==2212)) {
            goodLambda = kTRUE;
            if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(legPos)) goodLambda = kFALSE;
            if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(legNeg)) goodLambda = kFALSE;
         }
         if(processV0 && pdgV0==-3122 && (TMath::Abs(pdgP)==211 || TMath::Abs(pdgP)==2212) && (TMath::Abs(pdgN)==211 || TMath::Abs(pdgN)==2212)) {
            goodALambda = kTRUE;
            if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(legNeg)) goodALambda = kFALSE;
            if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(legPos)) goodALambda = kFALSE;
         }
         if(processV0 && TMath::Abs(pdgV0)==22 && TMath::Abs(pdgP)==11 && TMath::Abs(pdgN)==11) {
            goodGamma = kTRUE;
            if(fGammaElectronCuts && (!fGammaElectronCuts->IsSelected(legPos) || !fGammaElectronCuts->IsSelected(legNeg))) goodGamma = kFALSE;
         }
      }
      
      if(!((goodK0s && fFillK0s) || 
         (goodLambda && fFillLambda) || 
         (goodALambda && fFillALambda) || 
         (goodGamma && fFillGammaConversions))) {
         continue;
      }
      
      // Fill the V0 information into the tree for 4 hypothesis: K0s, Lambda, Anti-Lambda and gamma conversion
      AliReducedPairInfo* k0sReducedPair     = FillV0PairInfoAOD(v0, AliReducedPairInfo::kK0sToPiPi,     legPos, legNeg, &primaryVertexKF);
      AliReducedPairInfo* lambdaReducedPair  = FillV0PairInfoAOD(v0, AliReducedPairInfo::kLambda0ToPPi,  legPos, legNeg, &primaryVertexKF);
      AliReducedPairInfo* alambdaReducedPair = FillV0PairInfoAOD(v0, AliReducedPairInfo::kALambda0ToPPi, legPos, legNeg, &primaryVertexKF);
      AliReducedPairInfo* gammaReducedPair   = FillV0PairInfoAOD(v0, AliReducedPairInfo::kGammaConv,     legPos, legNeg, &primaryVertexKF);
      
      if(fFillK0s && goodK0s && k0sReducedPair->fMass[0]>fK0sMassRange[0] && k0sReducedPair->fMass[0]<fK0sMassRange[1]) {
         TClonesArray& tracks = *(fReducedEvent->fCandidates);
         AliReducedPairInfo *goodK0sPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPairInfo(*k0sReducedPair);
         goodK0sPair->fMass[0] = k0sReducedPair->fMass[0];
         goodK0sPair->fMass[1] = lambdaReducedPair->fMass[0];
         goodK0sPair->fMass[2] = alambdaReducedPair->fMass[0];
         goodK0sPair->fMass[3] = gammaReducedPair->fMass[0];
         fReducedEvent->fNV0candidates[1] += 1;
      } else {goodK0s=kFALSE;}
      if(fFillLambda && goodLambda && lambdaReducedPair->fMass[0]>fLambdaMassRange[0] && lambdaReducedPair->fMass[0]<fLambdaMassRange[1]) {
         TClonesArray& tracks = *(fReducedEvent->fCandidates);
         AliReducedPairInfo *goodLambdaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPairInfo(*lambdaReducedPair);
         goodLambdaPair->fMass[0] = k0sReducedPair->fMass[0];
         goodLambdaPair->fMass[1] = lambdaReducedPair->fMass[0];
         goodLambdaPair->fMass[2] = alambdaReducedPair->fMass[0];
         goodLambdaPair->fMass[3] = gammaReducedPair->fMass[0];
         fReducedEvent->fNV0candidates[1] += 1;
      } else {goodLambda=kFALSE;}
      if(fFillALambda && goodALambda && alambdaReducedPair->fMass[0]>fLambdaMassRange[0] && alambdaReducedPair->fMass[0]<fLambdaMassRange[1]) {
         TClonesArray& tracks = *(fReducedEvent->fCandidates);
         AliReducedPairInfo *goodALambdaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPairInfo(*alambdaReducedPair);
         goodALambdaPair->fMass[0] = k0sReducedPair->fMass[0];
         goodALambdaPair->fMass[1] = lambdaReducedPair->fMass[0];
         goodALambdaPair->fMass[2] = alambdaReducedPair->fMass[0];
         goodALambdaPair->fMass[3] = gammaReducedPair->fMass[0];
         fReducedEvent->fNV0candidates[1] += 1;
      } else {goodALambda = kFALSE;}
      if(fFillGammaConversions && goodGamma && gammaReducedPair->fMass[0]>fGammaMassRange[0] && gammaReducedPair->fMass[0]<fGammaMassRange[1]) {
         TClonesArray& tracks = *(fReducedEvent->fCandidates);
         AliReducedPairInfo *goodGammaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPairInfo(*gammaReducedPair);
         goodGammaPair->fMass[0] = k0sReducedPair->fMass[0];
         goodGammaPair->fMass[1] = lambdaReducedPair->fMass[0];
         goodGammaPair->fMass[2] = alambdaReducedPair->fMass[0];
         goodGammaPair->fMass[3] = gammaReducedPair->fMass[0];
         fReducedEvent->fNV0candidates[1] += 1;
      } else {goodGamma=kFALSE;}
      delete k0sReducedPair;
      delete lambdaReducedPair;
      delete alambdaReducedPair;
      delete gammaReducedPair;
   }   // end loop over V0s
}

//_________________________________________________________________________________
AliReducedPairInfo* AliAnalysisTaskReducedTreeMaker::FillV0PairInfoAOD(AliAODv0* v0, Int_t id, 
                                                                    AliAODTrack* legPos, AliAODTrack* legNeg,
                                                                    AliKFVertex* vtxKF) {
   //
   // Create a reduced V0 object and fill it
   //
   AliReducedPairInfo* reducedPair=new AliReducedPairInfo();  
   reducedPair->fCandidateId = id;
   reducedPair->fPairType    = v0->GetOnFlyStatus();    // on the fly status
   reducedPair->fLegIds[0]   = legPos->GetID();
   reducedPair->fLegIds[1]   = legNeg->GetID();
   
   // NOTE same treatment for both offline and on the fly V0s 
   if(id==AliReducedPairInfo::kLambda0ToPPi) reducedPair->fMass[0] = v0->MassLambda();
   if(id==AliReducedPairInfo::kALambda0ToPPi) reducedPair->fMass[0] = v0->MassAntiLambda();
   if(id==AliReducedPairInfo::kK0sToPiPi) reducedPair->fMass[0] = v0->MassK0Short();
   if(id==AliReducedPairInfo::kGammaConv) reducedPair->fMass[0] = 0.0;
   reducedPair->fIsCartesian  = kFALSE;
   reducedPair->fP[1]         = v0->Phi();
   if(reducedPair->fP[1]<0.0) reducedPair->fP[1] = 2.0*TMath::Pi() + reducedPair->fP[1];  // converted to [0,2pi]
   reducedPair->fP[0]         = v0->Pt();
   reducedPair->fP[2]         = v0->Eta();
   reducedPair->fLxy          = v0->RadiusV0();
   Double_t secVtx[3] = {vtxKF->GetX(), vtxKF->GetY(), vtxKF->GetZ()};
   reducedPair->fPointingAngle = v0->CosPointingAngle(secVtx);
   reducedPair->fChisquare    = v0->Chi2V0();
   
   return reducedPair;
}

//_________________________________________________________________________________
UChar_t AliAnalysisTaskReducedTreeMaker::EncodeTPCClusterMap(AliVParticle* track, Bool_t isAOD) {
  //
  // Encode the TPC cluster map into an UChar_t
  // Divide the 159 bits from the bit map into 8 groups of adiacent clusters
  // For each group enable its corresponding bit if in that group there are more clusters compared to
  // a threshold.
  //
  AliESDtrack* esdTrack=0x0;
  AliAODTrack* aodTrack=0x0;
  if(isAOD)
    aodTrack=static_cast<AliAODTrack*>(track);
  else
    esdTrack=static_cast<AliESDtrack*>(track);
  
  const UChar_t threshold=5;
  TBits tpcClusterMap = (isAOD ? aodTrack->GetTPCClusterMap() : esdTrack->GetTPCClusterMap());
  UChar_t map=0;
  UChar_t n=0;
  UChar_t j=0;
  for(UChar_t i=0; i<8; ++i) {
    n=0;
    for(j=i*20; j<(i+1)*20 && j<159; ++j) n+=tpcClusterMap.TestBitNumber(j);
    if(n>=threshold) map |= (1<<i);
  }
  return map;
}

//_________________________________________________________________________________
Int_t AliAnalysisTaskReducedTreeMaker::GetSPDTrackletMultiplicity(AliVEvent* event, Float_t lowEta, Float_t highEta) {
  //
  // Count the number of SPD tracklets in a given eta range
  //
  if (!event) return -1;
  
  Int_t nTracklets = 0;
  Int_t nAcc = 0;
  
  if(event->IsA() == AliAODEvent::Class()) {
    AliAODTracklets *tracklets = ((AliAODEvent*)event)->GetTracklets();
    nTracklets = tracklets->GetNumberOfTracklets();
    for(Int_t nn=0; nn<nTracklets; ++nn) {
      Double_t theta = tracklets->GetTheta(nn);
      Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
      if(eta < lowEta) continue;
      if(eta > highEta) continue;
      ++nAcc;
    }
  } else if(event->IsA() == AliESDEvent::Class()) {
    nTracklets = ((AliESDEvent*)event)->GetMultiplicity()->GetNumberOfTracklets();
    for(Int_t nn=0; nn<nTracklets; ++nn) {
      Double_t eta = ((AliESDEvent*)event)->GetMultiplicity()->GetEta(nn);
      if(eta < lowEta) continue;
      if(eta > highEta) continue; 
      ++nAcc;
    }
  } else return -1;
  
  return nAcc;
}

//______________________________________________________
Float_t AliAnalysisTaskReducedTreeMaker::GetInvPtDevFromBC(Int_t b, Int_t c)
{
  //
  //returns d(1/Pt) in c/GeV
  //in case of no gtu simulation -> return maximum 0.5
  //
  if(b==0 && c==0) return 0.5;
  Int_t tmp = (((b & 0xfff) << 12) ^ 0x800000) - 0x800000;
  tmp += (c & 0xfff);
  Float_t invPtDev = tmp * 0.000001;
  return invPtDev;
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FinishTaskOutput()
{
  //
  // Finish Task 
  //
  PostData(2, fTree);
}
