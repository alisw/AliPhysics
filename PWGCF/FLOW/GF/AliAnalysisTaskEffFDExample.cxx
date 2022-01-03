/*
  Example of AliEffFDContainer usage.
  The PCC framework (AliMCSpectraWeights) used in AliEffFDContainer is written by Patrick Huhn.
  For postprocessing of output to get efficiencies and feed-down corrections, use the
  postprocessing macro at: https://github.com/vvislavi/Feeddown
  If used, please acknowledge the authors of AliMCSpectraWeights as well as myself
  Author: Vytautas Vislavicius
*/

#include "AliAnalysisTaskEffFDExample.h"
#include "AliEventCuts.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliVVertex.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliEventCuts.h"
#include "AliEffFDContainer.h"

ClassImp(AliAnalysisTaskEffFDExample);

AliAnalysisTaskEffFDExample::AliAnalysisTaskEffFDExample():
  AliAnalysisTaskSE(),
  fIsMC(kFALSE),
  fUseGenPt(kTRUE),
  fMCEvent(0),
  fAddPID(kFALSE),
  fPIDResponse(0),
  fBayesPID(0),
  fBayesPIDProbs({0.95,0.85,0.85}),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
  fVtxZCut(10),
  fEfFd(0),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fContPF(""),
  fPtAxis(0),
  fMultiAxis(0),
  fPtBins(0),
  fNPtBins(0),
  fMultiBins(0),
  fNMultiBins(0),
  fCentEst("V0M"),
  fTCtoAdd(0),
  fFBtoAdd(0)
{
};
AliAnalysisTaskEffFDExample::AliAnalysisTaskEffFDExample(const char *name, Bool_t IsMC, TString pf):
  AliAnalysisTaskSE(name),
  fIsMC(IsMC),
  fUseGenPt(kTRUE),
  fMCEvent(0),
  fAddPID(kFALSE),
  fPIDResponse(0),
  fBayesPID(0),
  fBayesPIDProbs({0.95,0.85,0.85}),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
  fVtxZCut(10),
  fEfFd(0),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fContPF(pf),
  fPtAxis(0),
  fMultiAxis(0),
  fPtBins(0),
  fNPtBins(0),
  fMultiBins(0),
  fNMultiBins(0),
  fCentEst("V0M"),
  fTCtoAdd(0),
  fFBtoAdd(0)
{
  DefineInput(1,TH3D::Class());
  DefineOutput(1,AliEffFDContainer::Class());//AliEffFDContainer::Class());
};
AliAnalysisTaskEffFDExample::~AliAnalysisTaskEffFDExample() {
  if(fTCtoAdd) { fTCtoAdd->Clear(); delete fTCtoAdd; };
};
void AliAnalysisTaskEffFDExample::NotifyRun() {
  fEventCuts.OverrideAutomaticTriggerSelection(fTriggerType,true);
};
void AliAnalysisTaskEffFDExample::UserCreateOutputObjects(){
  OpenFile(1);
  //Setting up default pT and centrality bins
  const Int_t default_NCentBins=11;
  Double_t default_CentBins[default_NCentBins+1] = {0,5,10,20,30,40,50,60,70,80,90,100}; //Last bin to include V0M beyond anchor point
  if(!fMultiAxis) SetMultiBins(default_NCentBins,default_CentBins);

  const Int_t default_NPtBins = 21;
  Double_t default_PtBins[default_NPtBins+1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25,
                                              1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0,
                                              8.0, 10.0};
  if(!fPtAxis) SetPtBins(default_NPtBins,default_PtBins);
  //Creating the EffFD Object
  fEfFd = new AliEffFDContainer(Form("%s_EffAndFD%s",this->GetName(), fContPF.Data()),Form("%s_EffAndFD%s",this->GetName(), fContPF.Data()),fIsMC);
  fEfFd->SetEta(fEtaMin,fEtaMax);
  fEfFd->SetUseGenPt(fUseGenPt);
  if(!fContPF.IsNull()) fEfFd->SetIdentifier(fContPF);
  //Setting up centrality & pt bins
  Double_t *ptbins = GetBinsFromAxis(fPtAxis);
  fEfFd->SetPtBins(fPtAxis->GetNbins(),ptbins);
  Double_t *multibins = GetBinsFromAxis(fMultiAxis);
  fEfFd->SetCentralityBins(fMultiAxis->GetNbins(),multibins);
  //Centrality estimator, default is V0M
  fEfFd->SetCentralityEstimator(fCentEst);
  //Setting weights in case those are specified
  //In principle here we don't know if we'll be working with AODs, so just setting them up and if we end up working with ESDs, they won't be used
  TH3D *hWeights = (TH3D*)GetInputData(1);
  fEfFd->SetMCSWeights(hWeights);
  fEfFd->SetAODSelectionFlags(fEvNomFlag,fTrNomFlag);
  /* //Example how to add custom cuts:
  AliESDtrackCuts *tc = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  fEfFd->AddCut(tc); */
  //Instead, filter bits 32, 64, 256, and 512 are implemented:
  if(fTCtoAdd) for(Int_t i=0;i<fTCtoAdd->GetEntries();i++) fEfFd->AddCut((AliESDtrackCuts*)fTCtoAdd->At(i));
  if(fFBtoAdd) fEfFd->AddCut(fFBtoAdd);
  //Checking & setting up PID, if needed
  if(fAddPID) {
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    fBayesPID = new AliPIDCombined();
    fBayesPID->SetDefaultTPCPriors();
    fBayesPID->SetSelectedSpecies(AliPID::kSPECIES);
    fBayesPID->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
    fEfFd->SetPIDObjects(fPIDResponse,fBayesPID);
    fEfFd->SetBayesianProbs(fBayesPIDProbs);
  };
  PostData(1,fEfFd);
};
void AliAnalysisTaskEffFDExample::UserExec(Option_t*) {
  AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fESD && !fAOD) return;
  //Checking trigger
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(!(fTriggerType&fSelMask)) return ;
  if(fESD) { //If working on ESDs
    //Checking event cuts
    if(!fEventCuts.AcceptEvent(fESD)) return;
    if(TMath::Abs(fESD->GetPrimaryVertex()->GetZ())>fVtxZCut) return;
    //Fetching MC event
    if(fIsMC) {
      fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
      if (!fMCEvent) return;
    }
    //Filling: either ESD + MC (for MC), or ESD (for data)
    if(fIsMC) fEfFd->Fill(*fESD, *fMCEvent);
              else fEfFd->Fill(*fESD);
    PostData(1,fEfFd);
    return;
  } else if(fAOD) { //else if working on AODs
    if(!fIsMC) return; //If running on AODs, only makes sense running on MC ( = efficiencies), since we can't get the DCAxy anyways
    //Event and vertex selection embedded in AliGFWFlags
    // if(!fEventCuts.AcceptEvent(fAOD)) return;
    // if(TMath::Abs(fAOD->GetPrimaryVertex()->GetZ())>fVtxZCut) return;
    //Fetching MC event
    if(fIsMC) {
      fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
      if (!fMCEvent) return;
    }
    //Filling: for AODs, only makes sense to fill for MC
    fEfFd->Fill(*fAOD, *fMCEvent);
    PostData(1,fEfFd);
    return;

  }
};
void AliAnalysisTaskEffFDExample::Terminate(Option_t*) {
};
Double_t *AliAnalysisTaskEffFDExample::GetBinsFromAxis(TAxis *inax) {
  Int_t lBins = inax->GetNbins();
  Double_t *retBins = new Double_t[lBins+1];
  for(Int_t i=0;i<lBins;i++)
    retBins[i] = inax->GetBinLowEdge(i+1);
  retBins[lBins] = inax->GetBinUpEdge(lBins);
  return retBins;
}
void AliAnalysisTaskEffFDExample::SetupFlagsByIndex(Int_t ind) {
  fEvNomFlag=1<<kNominal;
  fTrNomFlag=1<<kFB96;
  switch(ind) {
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
  }
}
