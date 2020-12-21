AliAnalysisTaskBeauty *AddTaskBeauty(Bool_t applyeventw = kFALSE,TString file_momentum_smear="", TString versionsmearing="", Int_t processtype=0, UInt_t rndmseed=0, Double_t ptmax=1e20)
{
  AliAnalysisTaskBeauty* task = new  AliAnalysisTaskBeauty("");
  task->SetProcessType(processtype);
  task->SetSeed(rndmseed);
  task->SetPtCutHigh(ptmax);
  task->SetApplyEventw(applyeventw);

  // smearing
  gSystem->Exec(Form("alien_cp %s smearingfile.root",file_momentum_smear.Data()));
  TFile f("smearingfile.root");
  if (f.IsOpen() && ((TObjArray*)f.Get("ptSlices"))!=0x0) { // Old smearing file, only momentum.
    task->SetResolutionP((TObjArray*)f.Get("ptSlices"), kFALSE);
  }
  else { // New smearing file (or no file), from AliAnalysisTaskElectronEfficiency, postprocessed by LMeeAnaFW/Resolution/MakeResolutionArray.cxx.
    task->ReadResoFile(&f);
  }
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Printf("AliAnalysisTaskBeauty: No analysis manager to connect to.");
    return NULL;
  }
  
  if(!mgr->GetMCtruthEventHandler()){
    Printf("AliAnalysisTaskBeauty: This task requires an input MC event handler");
    return NULL;
  }
  
  mgr->AddTask(task);
  
  //Input and Output Slots:
  //AliAnalysisDataContainer *cinputSim = mgr->CreateContainer(inname,TChain::Class(), AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  //outputfile += ":KineSimulations";
  TString listname("lowee");
  if(processtype>0) listname += processtype;
  if(applyeventw) listname += "eventweight";
  listname += versionsmearing;
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(listname.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  
  return task;
}
