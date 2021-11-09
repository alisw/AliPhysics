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
  fMCEvent(0),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
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
  fMCEvent(0),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
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
  if(!fContPF.IsNull()) fEfFd->SetIdentifier(fContPF);
  //Setting up centrality & pt bins
  Double_t *ptbins = GetBinsFromAxis(fPtAxis);
  fEfFd->SetPtBins(fPtAxis->GetNbins(),ptbins);
  Double_t *multibins = GetBinsFromAxis(fMultiAxis);
  fEfFd->SetCentralityBins(fMultiAxis->GetNbins(),multibins);
  //Centrality estimator, default is V0M
  fEfFd->SetCentralityEstimator(fCentEst);
  /* //Example how to add custom cuts:
  AliESDtrackCuts *tc = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  fEfFd->AddCut(tc); */
  //Instead, filter bits 32, 64, 256, and 512 are implemented:
  if(fTCtoAdd) for(Int_t i=0;i<fTCtoAdd->GetEntries();i++) fEfFd->AddCut((AliESDtrackCuts*)fTCtoAdd->At(i));
  if(fFBtoAdd) fEfFd->AddCut(fFBtoAdd);
  PostData(1,fEfFd);
};
void AliAnalysisTaskEffFDExample::UserExec(Option_t*) {
  AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fESD) return;
  //Checking trigger
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(!(fTriggerType&fSelMask)) return ;
  //Checking event cuts
  if(!fEventCuts.AcceptEvent(fESD)) return;
  //Fetching MC event
  if(fIsMC) {
    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMCEvent) return;
  }
  //Filling: either ESD + MC (for MC), or ESD (for data)
  if(fIsMC) fEfFd->Fill(*fESD, *fMCEvent);
            else fEfFd->Fill(*fESD);
  PostData(1,fEfFd);
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
