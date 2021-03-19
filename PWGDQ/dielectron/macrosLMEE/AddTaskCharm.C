AliAnalysisTaskCharm *AddTaskCharm(Bool_t applyeventw = kFALSE,Bool_t applyweight=kFALSE,Bool_t selectoneccbar=kFALSE,Bool_t selectcleanhistory=kFALSE,TString file_momentum_smear="", TString versionsmearing="", TString file_efficiency="", Int_t fTypeEff=0, Int_t processtype=0, UInt_t rndmseed=0, Double_t ptmax=8.0,TString file_cnm="",TString cnm="",Bool_t takeptofDCNM=kFALSE)
{
  AliAnalysisTaskCharm* task = new  AliAnalysisTaskCharm("");
  task->SetProcessType(processtype);
  task->SetSeed(rndmseed);
  task->SetPtCutHigh(ptmax);
  task->SetApplywm(applyweight);
  task->SetApplyEventw(applyeventw);
  task->Selectoneccbar(selectoneccbar);
  task->Selectcleanhistory(selectcleanhistory);
  task->TakeptOfDCNM(takeptofDCNM);

  // Smearing
  if(file_momentum_smear.Contains("alien")) {
    gSystem->Exec(Form("alien_cp %s .",file_momentum_smear.Data()));
    TObjArray* Strings = file_momentum_smear.Tokenize("/");
    TString namefile = Form("%s/%s",gSystem->pwd(),Strings->At(Strings->GetEntriesFast()-1)->GetName());
    printf("Resolution file is %s copied from %s\n",namefile.Data(),file_momentum_smear.Data());
    TFile f(namefile.Data());
    if (f.IsOpen() && ((TObjArray*)f.Get("ptSlices"))!=0x0) { // Old smearing file, only momentum.
      task->SetResolutionP((TObjArray*)f.Get("ptSlices"), kFALSE);
    }
    else { // New smearing file (or no file), from AliAnalysisTaskElectronEfficiency, postprocessed by LMeeAnaFW/Resolution/MakeResolutionArray.cxx.
      task->ReadResoFile(&f);
    }
  }
  // CNM
  if(file_cnm.Contains("alien")) {
    gSystem->Exec(Form("alien_cp %s cnmfile.root",file_cnm.Data()));
    TObjArray* Strings = file_cnm.Tokenize("/");
    TString namefile = Form("%s/%s",gSystem->pwd(),Strings->At(Strings->GetEntriesFast()-1)->GetName());
    printf("CNM file is %s copied from %s\n",namefile.Data(),file_cnm.Data());
    TFile fcnm(namefile.Data());
    if (fcnm.IsOpen()){
      if((TGraph*)fcnm.Get(cnm.Data())!=0x0) { // apply cnm scaling.
	task->ScaleByCNM(kTRUE,(TGraph*)fcnm.Get(cnm.Data()));
      }
    }
  }
  // Efficiency
  if(file_efficiency.Contains("alien")) {
    gSystem->Exec(Form("alien_cp %s .",file_efficiency.Data()));
    TObjArray* Strings = file_efficiency.Tokenize("/");
    TString namefile = Form("%s/%s",gSystem->pwd(),Strings->At(Strings->GetEntriesFast()-1)->GetName());
    printf("Efficiency file is %s copied from %s\n",namefile.Data(),file_efficiency.Data());
    TFile fefficiency(namefile.Data());
    if (fefficiency.IsOpen()){
      task->SetEffType(fTypeEff);
      task->ReadEffFile(&fefficiency);
    }
  }
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Printf("AliAnalysisTaskCharm: No analysis manager to connect to.");
    return NULL;
  }
  
  if(!mgr->GetMCtruthEventHandler()){
    Printf("AliAnalysisTaskCharm: This task requires an input MC event handler");
    return NULL;
  }
  
  mgr->AddTask(task);
  
  //Input and Output Slots:
  //AliAnalysisDataContainer *cinputSim = mgr->CreateContainer(inname,TChain::Class(), AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  //outputfile += ":KineSimulations";
  TString listname("lowee");
  if(applyeventw) listname += "ew";
  if(processtype>0) listname += processtype;
  if(applyweight) listname += "weight";
  if(selectoneccbar) listname += "oneccbar";
  if(selectcleanhistory) listname += "cleanhistory";
  if(takeptofDCNM) listname += "takeptDCNM";
  listname += versionsmearing;
  listname += cnm;

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(listname.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  
  return task;
}
