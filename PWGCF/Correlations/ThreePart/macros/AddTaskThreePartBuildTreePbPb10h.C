#ifndef __CINT__
#include "AliAnalysisTaskBuildCorrTree.h"
#endif
AliAnalysisTaskBuildCorrTree* AddTaskThreePartBuildTreePbPb (const char* name = "ThreePartBuildTreePbPb",
						      const char* centrality = "V0M",
						      const Double_t MinPt = 1.0,
						      const Double_t MaxPt = 100.0,
						      const Double_t Acceptancecut = 0.9,
						      const char* period = "10h",
						      UInt_t offlineTriggerMask = AliVEvent::kMB,
						      const Int_t MaxNEventsMix = 10,
						      const Int_t MinNTracksMix = 2000,
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
  //Add a task AliAnalysisTaskBuildCorrTree to the analysis train in charged track analysis, for PbPb data 
  //Defaults to 11h data with MB trigger.
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
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  AliEPSelectionTask *eventplaneTask = new AliEPSelectionTask("EventplaneSelection");
  eventplaneTask->SelectCollisionCandidates(offlineTriggerMask);
  if (inputDataType == "AOD"){
    eventplaneTask->SetInput("AOD");
  }
  eventplaneTask->SetTrackType("TPC");
  eventplaneTask->SetUsePtWeight();
  eventplaneTask->SetUsePhiWeight();
  eventplaneTask->SetSaveTrackContribution();
  mgr->AddTask(eventplaneTask);
  mgr->ConnectInput(eventplaneTask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *coutputEP = mgr->CreateContainer("EPStat",TList::Class(), AliAnalysisManager::kOutputContainer,"EventStat_temp.root");
  mgr->ConnectOutput(eventplaneTask,1,coutputEP);
  
  
  
  const char* fname = Form("%s_%1.0f_%1.0f",name,MinPt,MaxPt);
  const char* tname = Form("%s_%1.0f_%1.0f_%1.0f_%1.0f",name,MinPt,MaxPt);
  AliAnalysisTaskBuildCorrTree* task = new AliAnalysisTaskBuildCorrTree(Form("%sTask", tname), "");

  task->SetCentralityEstimator(centrality);
  task->SetMinPt(MinPt);
  task->SetMaxPt(MaxPt);
  task->SetAcceptanceCut(Acceptancecut);
  task->SetCollisionType(AliAnalysisTaskBuildCorrTree::PbPb);
  

  
  //Mixing scheme:
  Double_t *Mbin = new Double_t[NMBins+1];
  Mbin[0] = Mbin0;
  Mbin[1] = Mbin1;
  if(NMBins>1) Mbin[2] = Mbin2;
  if(NMBins>2) Mbin[3] = Mbin3;
  if(NMBins>3) Mbin[4] = Mbin4;
  if(NMBins>4) Mbin[5] = Mbin5;
  if(NMBins>5) Mbin[6] = Mbin6;
  if(NMBins>6) Mbin[7] = Mbin7;

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
  task->SelectCollisionCandidates(offlineTriggerMask);

  
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  TString cname(Form("%sCoutput1", tname));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), fname));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("dstTree",TTree::Class(),AliAnalysisManager::kOutputContainer,"dstTree.root");
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);

    return task;
};
