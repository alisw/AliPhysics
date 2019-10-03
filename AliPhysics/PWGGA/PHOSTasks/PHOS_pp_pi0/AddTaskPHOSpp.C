AliAnalysisTaskPi0PP * AddTaskPHOSpp (const char* name = "PHOSpp",
					    Double_t tCut = 525.e-09)
{
  //Add a task AliAnalysisTaskPi0PP to the analysis train
  //Author: Pooja Pareek
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSGammaFlow", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSGammaFlow", "This task requires an input event handler");
    return NULL;
  }
  
  AliAnalysisTaskPi0PP * task = new AliAnalysisTaskPi0PP(Form("%s", name));
  
  
  TFile *fBadMap = TFile::Open("alien:///alice/cern.ch/user/k/kharlov/Analysis/AnalysisTask/BadMap_LHC10de_Majority.root");
  if(fBadMap->IsOpen()){
    printf("\n\n...Adding PHOS bad channel map \n") ;
    gROOT->cd();
    char key1[55] ;
    for(Int_t mod=1;mod<4; mod++){
      sprintf(key1,"PHOS_BadMap_mod%d",mod) ;
      TH2I * h = (TH2I*)fBadMap->Get(key1) ;
      if(h)
        task->SetPHOSBadMap(mod,h) ;
    }
    fBadMap->Close() ;
  }

  task->SelectCollisionCandidates(AliVEvent::kMB);

  
  
  // Set BC gap for LHC10e in seconds
  task->SetBCgap(tCut);

  // Set abs.recalibration for LHC11a
  task->SetRecalib(1, 0.9822696);
  task->SetRecalib(2, 0.9861288);
  task->SetRecalib(3, 1.0072);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  TString cname(Form("%s", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectInput(task , 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
