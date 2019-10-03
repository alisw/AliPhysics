AliAnalysisTaskPi0Flow* AddTaskPHOSHijingEff(const char* name = "PHOSHijingEff",
					    const char* options = "")
{
  //Add a task AliAnalysisTaskPi0Flow to the analysis train
  //Author: Henrik Qvigstad
  /* $Id: AddTaskPHOSPi0Flow.C 59900 2012-12-10 02:17:18Z hqvigsta $ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSHijingEff", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSHijingEff", "This task requires an input event handler");
    return NULL;
  }

  if (!mgr->GetMCtruthEventHandler()){
    ::Error("AddTaskPHOSHijingEff", "This task requires an MC input event handler");
    return NULL;
  }
    
  AliPHOSHijingEfficiency* task = new AliPHOSHijingEfficiency(Form("%sTask", name));

  TFile *fBadMap = TFile::Open("alien:///alice/cern.ch/user/p/prsnko/BadMaps/BadMap_LHC10h_period234.root");

  if(fBadMap->IsOpen()){
    printf("\n\n...Adding PHOS bad channel map \n") ;
    gROOT->cd();
    char key[55] ;
    for(Int_t mod=1;mod<4; mod++){
      sprintf(key,"PHOS_BadMap_mod%d",mod) ;
      TH2I * h = (TH2I*)fBadMap->Get(key) ;
      if(h)
	task->SetPHOSBadMap(mod,h) ;
    }
    fBadMap->Close() ;
  }
  
  
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  TString cname(Form("%sCoutput1", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
