AliAnalysisTaskBeauty *AddTaskBeauty(Bool_t applyeventw = kFALSE,TString file_momentum_smear="", TString versionsmearing="", TString file_efficiency="", Int_t fTypeEff=0, Int_t processtype=0, UInt_t rndmseed=0, Double_t ptmax=1e20, Bool_t selectonebbbar=kFALSE, Float_t MinOpAng=0.0, TString file_raa="", TString hname_raa="")
{

  TString listname("lowee");
  if(processtype>0) listname += processtype;
  if(applyeventw) listname += "eventweight";
  listname += versionsmearing;
  
  AliAnalysisTaskBeauty* task = new  AliAnalysisTaskBeauty(Form("Task_%s",listname.Data()));
  task->SetProcessType(processtype);
  task->SetSeed(rndmseed);
  task->SetPtCutHigh(ptmax);
  task->SetMinOpAng(MinOpAng);
  task->SetApplyEventw(applyeventw);
  task->Selectonebbbar(selectonebbbar);

  // Smearing
  if(file_momentum_smear.Contains("alien")) {
    gSystem->Exec(Form("alien_cp %s file:./",file_momentum_smear.Data()));
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
  
  // Efficiency
  if(file_efficiency.Contains("alien")) {
    gSystem->Exec(Form("alien_cp %s file:./",file_efficiency.Data()));
    TObjArray* Strings = file_efficiency.Tokenize("/");
    TString namefile = Form("%s/%s",gSystem->pwd(),Strings->At(Strings->GetEntriesFast()-1)->GetName());
    printf("Efficiency file is %s copied from %s\n",namefile.Data(),file_efficiency.Data());
    TFile fefficiency(namefile.Data());
    if (fefficiency.IsOpen()){
      task->SetEffType(fTypeEff);
      task->ReadEffFile(&fefficiency);
    }
  }

  //RAA
  if(file_raa.Contains("alien")) {
    gSystem->Exec(Form("alien_cp %s file:./",file_raa.Data()));
    TObjArray* Strings = file_raa.Tokenize("/");
    TString namefile = Form("%s/%s",gSystem->pwd(),Strings->At(Strings->GetEntriesFast()-1)->GetName());
    printf("RAA file is %s copied from %s\n",namefile.Data(),file_raa.Data());
    TFile fraa(namefile.Data());
    if (fraa.IsOpen()){
      if(namefile.Contains("FIT")){
	// fit functions
	if((TF1*)fraa.Get(hname_raa.Data()) != 0x0) { // apply RAA weighting.
	  TF1 *h1RAA = (TF1*)fraa.Get(hname_raa.Data());
	  task->ScaleByRAA(kTRUE);
	  task->SetTF1RAA(h1RAA);
	}
      }
      else {
	// histo
	if((TH1F*)fraa.Get(hname_raa.Data()) != 0x0) { // apply RAA weighting.
	  TH1F *h1RAA = (TH1F*)fraa.Get(hname_raa.Data());
	  h1RAA->SetDirectory(0);
	  task->ScaleByRAA(kTRUE);
	  task->SetTH1FRAA(h1RAA);
	}
      }
    }
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
  
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(listname.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  
  return task;
}
