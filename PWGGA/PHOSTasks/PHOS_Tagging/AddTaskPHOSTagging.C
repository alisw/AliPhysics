AliAnalysisTaskTaggedPhotons* AddTaskPHOSTagging (const char* name = "PHOSTagging",
					    const char* options = "",
					    UInt_t offlineTriggerMask = AliVEvent::kCentral,
					    Float_t timeCut = 100.e-9, //accept clusters with |t|<timeCut
                                            Bool_t ignorePHI7Events=kFALSE,
                                            Int_t centralityEstinator=1) //Centrality estimator, see cxx code for list (separate for pp and pPb
{
  //Add a task AliAnalysisTaskTaggedPhotons to the analysis train
  //Author: Dmitri Peresunko

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSTagging", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSTagging", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskTaggedPhotons* task = new AliAnalysisTaskTaggedPhotons(Form("%sTask%d", name,centralityEstinator));

  task->SelectCollisionCandidates(offlineTriggerMask);
 
  task->SetTimeCut(25.e-9) ;
  task->SetTrigger(ignorePHI7Events) ;
  task->SetCentralityEstimator(centralityEstinator) ; 
  
  Int_t binLimits[8]={5,10,15,20,30,50,100,300};
  TArrayI multBins(8,binLimits) ;
  task->SetMultiplicityBins(multBins) ;
  
 
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  TString cname(Form("%sCoutput1", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), THashList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
