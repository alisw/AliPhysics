/*
  Example of AliEffFDContainer usage.
  The PCC framework (AliMCSpectraWeights) used in AliEffFDContainer is written by Patrick Huhn.
  For postprocessing of output to get efficiencies and feed-down corrections, use the
  postprocessing macro at: https://github.com/vvislavi/Feeddown
  If used, please acknowledge the authors of AliMCSpectraWeights as well as myself
  Author: Vytautas Vislavicius
*/

#include "AliAnalysisTaskLWTree.h"
#include "AliEventCuts.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliEventCuts.h"
#include "AliEffFDContainer.h"
#include "AliMultSelection.h"
#include "TClonesArray.h"


ClassImp(AliAnalysisTaskLWTree);

AliAnalysisTaskLWTree::AliAnalysisTaskLWTree():
  AliAnalysisTaskSE(),
  fOutTree(0),
  fLWEvent(0),
  fTPCTracks(0),
  fFMDTracks(0),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
  fVtxZCut(10),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fCentEst("V0M"),
  fTCtoAdd(0),
  fFBtoAdd(0)
{
};
AliAnalysisTaskLWTree::AliAnalysisTaskLWTree(const char *name):
  AliAnalysisTaskSE(name),
  fOutTree(0),
  fLWEvent(0),
  fTPCTracks(0),
  fFMDTracks(0),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
  fVtxZCut(10),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fCentEst("V0M"),
  fTCtoAdd(0),
  fFBtoAdd(0)
{
  DefineOutput(1,TTree::Class());//AliEffFDContainer::Class());
};
AliAnalysisTaskLWTree::~AliAnalysisTaskLWTree() {
  if(fTCtoAdd) { fTCtoAdd->Clear(); delete fTCtoAdd; };
};
void AliAnalysisTaskLWTree::NotifyRun() {
  fEventCuts.OverrideAutomaticTriggerSelection(fTriggerType,true);
};
void AliAnalysisTaskLWTree::UserCreateOutputObjects(){
  OpenFile(1);
  //Initialize LW event and TClonesArrays for tracks
  fTPCTracks = new TClonesArray("AliLWTPCTrack");
  fFMDTracks = new TClonesArray("AliLWFMDTrack");
  fLWEvent   = new AliLWEvent();
  fOutTree   = new TTree("LWTree","LWTree");
  fOutTree->Branch("Event","AliLWEvent",&fLWEvent);
  fOutTree->Branch("TPCTracks","TClonesArray",&fTPCTracks);
  fOutTree->Branch("FMDTracks","TClonesArray",&fFMDTracks);
  PostData(1,fOutTree);
};
void AliAnalysisTaskLWTree::UserExec(Option_t*) {
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;
  //Checking trigger
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(!(fTriggerType&fSelMask)) return;
  //Processing/filling the tree
  //First, check if the event passes event cuts
  AliGFWFlags *lFlags = (AliGFWFlags*)fInputEvent->FindListObject("GFWFlags");
  if(!lFlags) { printf("GFWFlags not found!\n"); return; };
  UInt_t gEventFlag = lFlags->GetEventFlags();
  if(!(gEventFlag&fEvNomFlag)) return; //otherwise, if not the selected event flag, then move on
  AliMultSelection *lMultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  Double_t cent = lMultSel->GetMultiplicityPercentile(fCentEst.Data());
  fLWEvent->Setup(fInputEvent->GetRunNumber(),fAOD->GetPrimaryVertex()->GetZ(),cent,gEventFlag);
  //Looping over TPC tracks
  fTPCTracks->Clear();
  fFMDTracks->Clear();
  UInt_t gTrackFlags;
  AliAODTrack *lTrack;
  Int_t l_tcaInd;
  //Processing TPC tracks
  for(Int_t lTr=0;lTr<lFlags->GetNFiltered();lTr++) {
    gTrackFlags = lFlags->GetTrackFlag(lTr);
    if(!(gTrackFlags&fTrNomFlag)) continue; //Check if we want to accept the track
    Int_t trInd = lFlags->GetTrackIndex(lTr);
    lTrack = (AliAODTrack*)fAOD->GetTrack(trInd);
    l_tcaInd = fTPCTracks->GetEntries();
    new ((*fTPCTracks)[l_tcaInd]) AliLWTPCTrack(lTrack->Pt(),lTrack->Phi(),lTrack->Eta(),gTrackFlags);
  };
  if(fTPCTracks->GetEntries()>1) fTPCTracks->Sort(); //Sort tpc tracks following pT for more efficient postprocessing
  //Processing FMD info
  AliAODForwardMult* aodForward=static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
  if(!aodForward) { return; }
  const TH2D& d2Ndetadphi = aodForward->GetHistogram();
  for(Int_t iEta=1; iEta<=d2Ndetadphi.GetNbinsX(); iEta++) {
    Double_t etaCent = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
    //No need to loop over uncovered regions
    if(etaCent<-3.5) continue;
    else if(etaCent>-1.6 && etaCent<1.7) continue;
    else if(etaCent>5) continue;
    for(Int_t iPhi=1; iPhi<=d2Ndetadphi.GetNbinsY();iPhi++) {
      if(!d2Ndetadphi.GetBinContent(iEta,iPhi)) continue;
      Double_t phiCent = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
      new ((*fFMDTracks)[fFMDTracks->GetEntries()]) AliLWFMDTrack(phiCent,etaCent,d2Ndetadphi.GetBinContent(iEta,iPhi));
    }
  }
  fOutTree->Fill();
  PostData(1,fOutTree);
};
void AliAnalysisTaskLWTree::Terminate(Option_t*) {
};
void AliAnalysisTaskLWTree::SetupFlags() {
  fEvNomFlag=1<<kNominal;
  //Here let's store only FB96 and 768 tracks without further info, at least for now
  fTrNomFlag=1<<kFB96;
  fTrNomFlag+= (1<<kFB768);
/*  switch(ind) {
    default: // also 0
      break;
    //Event flags:
    case 1:
      fEvNomFlag = 1<<kVtx9;
      break;
    case 2:
      fEvNomFlag = 1<<kVtx7;
      break;
    case 3:
      fEvNomFlag = 1<<kVtx5;
      break;
    //Track flags:
    case 4:
      fTrNomFlag = 1<<kFB768;
      break;
    case 5:
      fTrNomFlag = 1<<kDCAz10;
      break;
    case 6:
      fTrNomFlag = 1<<kDCAz05;
      break;
    case 7:
      fTrNomFlag = 1<<kDCA4Sigma;
      break;
    case 8:
      fTrNomFlag = 1<<kDCA10Sigma;
      break;
    case 9:
      fTrNomFlag = 1<<kChiSq2;
      break;
    case 10:
      fTrNomFlag = 1<<kChiSq3;
      break;
    case 11:
      fTrNomFlag = 1<<kNTPC80;
      break;
    case 12:
      fTrNomFlag = 1<<kNTPC90;
      break;
    case 13:
      fTrNomFlag = 1<<kNTPC100;
      break;
    case 14:
      fTrNomFlag = 1<<kFB768Tuned;
      break;
    case 15:
      fTrNomFlag = 1<<kFB96Tuned;
      break;
    case 16:
      fTrNomFlag = 1<<kFB768DCAz;
      break;
    case 17:
      fTrNomFlag = 1<<kFB768DCAxyLow;
      break;
    case 18:
      fTrNomFlag = 1<<kFB768DCAxyHigh;
      break;
    case 19:
      fTrNomFlag = 1<<kFB768ChiSq2;
      break;
    case 20:
      fTrNomFlag = 1<<kFB768ChiSq3;
      break;
    case 21:
      fTrNomFlag = 1<<kFB768nTPC;
      break;
    case 22:
      fTrNomFlag = 1<<kFB96MergedDCA;
      break;
  }*/
}
