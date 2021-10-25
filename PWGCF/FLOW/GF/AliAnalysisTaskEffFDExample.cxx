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
  fEfFd(0)
{
};
AliAnalysisTaskEffFDExample::AliAnalysisTaskEffFDExample(const char *name, Bool_t IsMC):
  AliAnalysisTaskSE(name),
  fIsMC(IsMC),
  fMCEvent(0),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
  fEfFd(0)
{
  DefineOutput(1,AliEffFDContainer::Class());//AliEffFDContainer::Class());
};
AliAnalysisTaskEffFDExample::~AliAnalysisTaskEffFDExample() {
};
void AliAnalysisTaskEffFDExample::UserCreateOutputObjects(){
  OpenFile(1);
  //Setting up default pT and centrality bins
  const Int_t NCentBinsDefault=11;
  Double_t CentBinsDefault[NCentBinsDefault+1] = {0,5,10,20,30,40,50,60,70,80,90,100}; //Last bin to include V0M beyond anchor point
  const Int_t NPtBinsDefault = 25;
  Double_t PtBinsDefault[NPtBinsDefault+1] = {0.20, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
                     1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90,
                     2.00, 2.20, 2.40, 2.60, 2.80, 3.00};
  fEfFd = new AliEffFDContainer("EffAndFD","EffAndFD",fIsMC);
  fEfFd->SetCentralityBins(NCentBinsDefault,CentBinsDefault);
  fEfFd->SetPtBins(NPtBinsDefault,PtBinsDefault);
  fEfFd->SetEta(0.8);
  fEfFd->SetCentralityEstimator("V0M");

  /* //Example how to add custom cuts:
  AliESDtrackCuts *tc = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  fEfFd->AddCut(tc); */
  //Instead, filter bits 32, 64, 256, and 512 are implemented:
  fEfFd->AddCut(96);
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
