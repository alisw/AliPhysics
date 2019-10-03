#ifndef __CINT__
#include "AliAnalysisTaskCorrelation3p.h"
#endif
AliAnalysisTaskCorrelation3p* AddTaskThreePartGenerator (const char* name = "ThreePartGenerator",
						      const char* options = "",
						      const Double_t MinTriggerPt = 4.0,
						      const Double_t MaxTriggerPt = 8.0,
						      const Double_t MinAssociatedPt = 2.0,
						      const Double_t MaxAssociatedPt = 4.0,
						      const Double_t Acceptancecut = 0.8,
						      const Double_t MaxNumberOfTracks = 200,
						      UInt_t offlineTriggerMask = AliVEvent::kMB,
						      const Int_t MaxNEventsMix = 1000,
						      const Int_t MinNTracksMix = 100,
						      const Int_t NMBins = 1,
						      const Int_t NZBins = 1,
						      const Double_t Mbin0 = 0.,
						      const Double_t Mbin1 = 100.,
						      const Double_t Zbin0 = -10.,
						      const Double_t Zbin1 = 10.
 						    )
{
  //Add a task AliAnalysisTaskCorrelation3p to the analysis train in charged track analysis, for pp data 
  //Author: Paul Baetzing
  /* $Id$ */
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr){
    ::Error("AddTaskThreePartTracks", "No analysis manager to connect to");
    return NULL;
  }
    if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskThreePartTracks", "This task requires an input event handler");
    return NULL;
  }
  AliAnalysisTaskCorrelation3p* task = new AliAnalysisTaskCorrelation3p(Form("%sTask", name), options);

  task->SetTrigger(AliAnalysisTaskCorrelation3p::tracks);
  task->SetMinTriggerPt(MinTriggerPt);
  task->SetMaxTriggerPt(MaxTriggerPt);
  task->SetMinAssociatedPt(MinAssociatedPt);
  task->SetMaxAssociatedPt(MaxAssociatedPt);
  task->SetAcceptanceCut(Acceptancecut);
  task->SetMaxNumberOfTracks(MaxNumberOfTracks);
  task->SetPeriod(AliAnalysisTaskCorrelation3p::P10h);
  task->SetGenerate();
  //Mixing scheme:
  Double_t *Mbin = new Double_t[NMBins+1];
  Mbin[0] = Mbin0;
  Mbin[1] = Mbin1;
  TArrayD tMbin(NMBins+1, Mbin);
  Double_t *Zbin = new Double_t[NZBins+1];
  Zbin[0] = Zbin0;
  Zbin[1] = Zbin1;
  TArrayD tZbin(NZBins+1, Zbin);  
  task->SetMixingScheme(MaxNEventsMix,MinNTracksMix,tMbin,tZbin);

  task->SelectCollisionCandidates(offlineTriggerMask);

  
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  TString cname(Form("%sCoutput1", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
};
