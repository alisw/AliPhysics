#ifndef __CINT__
#include "AliAnalysisTaskCorrelation3p.h"
#endif
AliAnalysisTaskCorrelation3p* AddTaskThreePartTracks (const char* name = "ThreePartTracks",
						      const char* options = "",
						      const char* centrality = "V0M",
						      const Double_t MinTriggerPt = 4.0,
						      const Double_t MaxTriggerPt = 8.0,
						      const Double_t MinAssociatedPt = 2.0,
						      const Double_t MaxAssociatedPt = 4.0,
						      const Double_t Acceptancecut = 0.8,
						      const Double_t MaxNumberOfTracks = 200,
						      const char* period = "10b",
						      UInt_t offlineTriggerMask = AliVEvent::kMB,
						      const Int_t MaxNEventsMix = 1000,
						      const Int_t MinNTracksMix = 100,
						      const Int_t NMBins = 7,
						      const Int_t NZBins = 5,
						      const Double_t Mbin0 = 0.,
						      const Double_t Mbin1 = 5.,
						      const Double_t Mbin2 = 10.,
						      const Double_t Mbin3 = 20.,
						      const Double_t Mbin4 = 40.,
						      const Double_t Mbin5 = 60.,
						      const Double_t Mbin6 = 80.,
						      const Double_t Mbin7 = 90.,
						      const Double_t Zbin0 = -10.,
						      const Double_t Zbin1 = -5.,
						      const Double_t Zbin2 = -2.0.,
						      const Double_t Zbin3 = 2.0.,
						      const Double_t Zbin4 = 5.,
						      const Double_t Zbin5 = 10.
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

  task->SetCentralityEstimator(centrality);
  task->SetTrigger(AliAnalysisTaskCorrelation3p::tracks);
  task->SetMinTriggerPt(MinTriggerPt);
  task->SetMaxTriggerPt(MaxTriggerPt);
  task->SetMinAssociatedPt(MinAssociatedPt);
  task->SetMaxAssociatedPt(MaxAssociatedPt);
  task->SetAcceptanceCut(Acceptancecut);
  task->SetMaxNumberOfTracks(MaxNumberOfTracks);
  
  //Mixing scheme:
  Double_t *Mbin = new Double_t[NMBins+1];
  Mbin[0] = Mbin0;
  Mbin[1] = Mbin1;
  if(NMBins>1) Mbin[2] = Mbin2;
  if(NMBins>2) Mbin[3] = Mbin3;
  if(NMBins>3) Mbin[4] = Mbin4;
  if(NMBins>4) Mbin[5] = Mbin5;
  if(NMBins>5) Mbin[6] = Mbin6;
  if(NMBins>6) Mbin[7] = Mbin6;
  TArrayD tMbin(NMBins+1, Mbin);
  Double_t *Zbin = new Double_t[NZBins+1];
  Zbin[0] = Zbin0;
  Zbin[1] = Zbin1;
  if(NZBins>1) Zbin[2] = Zbin2;
  if(NZBins>2) Zbin[3] = Zbin3;
  if(NZBins>3) Zbin[4] = Zbin4;
  if(NZBins>4) Zbin[5] = Zbin5;
  TArrayD tZbin(NZBins+1, Zbin);  
  task->SetMixingScheme(MaxNEventsMix,MinNTracksMix,tMbin,tZbin);

  if( TString(period).Contains("10b") )
    task->SetPeriod(AliAnalysisTaskCorrelation3p::P10b);
  if( TString(period).Contains("10c") )
    task->SetPeriod(AliAnalysisTaskCorrelation3p::P10c);
  if( TString(period).Contains("10d") )
    task->SetPeriod(AliAnalysisTaskCorrelation3p::P10d);  
  if( TString(period).Contains("10e") )
    task->SetPeriod(AliAnalysisTaskCorrelation3p::P10e);
  if( TString(period).Contains("11a") )
    task->SetPeriod(AliAnalysisTaskCorrelation3p::P11a);
  task->SelectCollisionCandidates(offlineTriggerMask);

  
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  TString cname(Form("%sCoutput1", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
};
