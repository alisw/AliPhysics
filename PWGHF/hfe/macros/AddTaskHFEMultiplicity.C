AliAnalysisTask* AddTaskHFEMultiplicity(TString suffixName = "",
					Bool_t PhysSelINT7 = kTRUE,
					Bool_t useTender   =  kTRUE,
					Bool_t ClsTypeEMC  =  kTRUE,
					Bool_t ClsTypeDCAL =  kTRUE,
					TString estimatorFilename="",
					Double_t refMult=61.2
					)
  
{ // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }
  
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":MyTask";
  
  if(PhysSelINT7){
    AliAnalysisTaskHFEMultiplicity* HFEtaskINT7 = new AliAnalysisTaskHFEMultiplicity("");
    mgr->AddTask(HFEtaskINT7); //HFEtask is my task
    HFEtaskINT7->SelectCollisionCandidates(AliVEvent::kINT7);
    HFEtaskINT7->SetTenderSwitch(useTender);
    
    HFEtaskINT7->SetClusterTypeEMC(ClsTypeEMC);
    HFEtaskINT7->SetClusterTypeDCAL(ClsTypeDCAL);

 
    TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
    if(!fileEstimator)  {
      AliFatal("File with multiplicity estimator not found\n");
      return;
    }
    HFEtaskINT7->SetReferenceMultiplicity(refMult);
    const Char_t* profilebasename="SPDmult";
    
    const Char_t* periodNames[2] = {"LHC16s", "LHC16r"};
    TProfile* multEstimatorAvg[2];
    for(Int_t ip=0; ip<2; ip++) {
      cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
      multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
      if (!multEstimatorAvg[ip]) {
	AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
	return;
      }
    }

    HFEtaskINT7->SetMultiProfileLHC16s(multEstimatorAvg[0]);
    HFEtaskINT7->SetMultiProfileLHC16r(multEstimatorAvg[1]);

    // Create containers for input/output
    TString finDirname         = "_INT7";
    TString outBasicname       = "EID";
   // TString profname       = "coutputProf";
        
    finDirname 	      +=   suffixName.Data();
    outBasicname      +=   finDirname.Data();
    //profname          +=   finDirname.Data();
        
        
        
    mgr->ConnectInput(HFEtaskINT7,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(HFEtaskINT7,1,mgr->CreateContainer(outBasicname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
   // mgr->ConnectOutput(HFEtaskINT7,2,mgr->CreateContainer(profname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  }
    
   if(!PhysSelINT7){
     if(ClsTypeEMC){
            
     AliAnalysisTaskHFEMultiplicity* HFEtaskEG1 = new AliAnalysisTaskHFEMultiplicity("");
     mgr->AddTask(HFEtaskEG1);
     HFEtaskEG1->SelectCollisionCandidates(AliVEvent::kEMCEGA);
     HFEtaskEG1->SetEMCalTriggerEG1(kTRUE);
     HFEtaskEG1->SetTenderSwitch(useTender);
     HFEtaskEG1->SetClusterTypeEMC(ClsTypeEMC);
     HFEtaskEG1->SetClusterTypeDCAL(ClsTypeDCAL);

    TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
    if(!fileEstimator)  {
      AliFatal("File with multiplicity estimator not found\n");
      return;
    }
    HFEtaskEG1->SetReferenceMultiplicity(refMult);
    const Char_t* profilebasename="SPDmultEG1";
    
    const Char_t* periodNames[2] = {"LHC16s", "LHC16r"};
    TProfile* multEstimatorAvgEG1[2];
    for(Int_t ip=0; ip<2; ip++) {
      cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
      multEstimatorAvgEG1[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
      if (!multEstimatorAvgEG1[ip]) {
	AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
	return;
      }
    }

    HFEtaskEG1->SetMultiProfileLHC16s(multEstimatorAvgEG1[0]);
    HFEtaskEG1->SetMultiProfileLHC16r(multEstimatorAvgEG1[1]);
            
            
     TString finDirname         = "_TrigGAEG1";
     TString outBasicname       = "EID";
            
     finDirname 	      +=   suffixName.Data();
     outBasicname      +=   finDirname.Data();
            
            
            
     mgr->ConnectInput(HFEtaskEG1,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(HFEtaskEG1,1,mgr->CreateContainer(outBasicname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
     // EMCal EGA EG2
            
     AliAnalysisTaskHFEMultiplicity* HFEtaskEG2 = new AliAnalysisTaskHFEMultiplicity("");
     mgr->AddTask(HFEtaskEG2);
     HFEtaskEG2->SelectCollisionCandidates(AliVEvent::kEMCEGA);
     HFEtaskEG2->SetEMCalTriggerEG2(kTRUE);
     HFEtaskEG2->SetTenderSwitch(useTender);
     HFEtaskEG2->SetClusterTypeEMC(ClsTypeEMC);
     HFEtaskEG2->SetClusterTypeDCAL(ClsTypeDCAL);

      TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
    if(!fileEstimator)  {
      AliFatal("File with multiplicity estimator not found\n");
      return;
    }
    HFEtaskEG2->SetReferenceMultiplicity(refMult);
    const Char_t* profilebasename="SPDmultEG2";
    
    const Char_t* periodNames[2] = {"LHC16s", "LHC16r"};
    TProfile* multEstimatorAvgEG2[2];
    for(Int_t ip=0; ip<2; ip++) {
      cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
      multEstimatorAvgEG2[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
      if (!multEstimatorAvgEG2[ip]) {
	AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
	return;
      }
    }

    HFEtaskEG2->SetMultiProfileLHC16s(multEstimatorAvgEG2[0]);
    HFEtaskEG2->SetMultiProfileLHC16r(multEstimatorAvgEG2[1]);
            
            
     TString finDirname         = "_TrigGAEG2";
     TString outBasicname       = "EID";
            
     finDirname 	      +=   suffixName.Data();
     outBasicname      +=   finDirname.Data();
            
            
            
     mgr->ConnectInput(HFEtaskEG2,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(HFEtaskEG2,1,mgr->CreateContainer(outBasicname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
            
     }
        
     if(ClsTypeDCAL){
     // DCal EGA DG1
     AliAnalysisTaskHFEMultiplicity* HFEtaskDG1 = new AliAnalysisTaskHFEMultiplicity("");
     mgr->AddTask(HFEtaskDG1);
     HFEtaskDG1->SelectCollisionCandidates(AliVEvent::kEMCEGA);
     HFEtaskDG1->SetEMCalTriggerDG1(kTRUE);
     HFEtaskDG1->SetTenderSwitch(useTender);
     HFEtaskDG1->SetClusterTypeEMC(ClsTypeEMC);
     HFEtaskDG1->SetClusterTypeDCAL(ClsTypeDCAL);
      
      TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
    if(!fileEstimator)  {
      AliFatal("File with multiplicity estimator not found\n");
      return;
    }
    HFEtaskDG1->SetReferenceMultiplicity(refMult);
    const Char_t* profilebasename="SPDmultDG1";
    
    const Char_t* periodNames[2] = {"LHC16s", "LHC16r"};
    TProfile* multEstimatorAvgDG1[2];
    for(Int_t ip=0; ip<2; ip++) {
      cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
      multEstimatorAvgDG1[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
      if (!multEstimatorAvgDG1[ip]) {
	AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
	return;
      }
    }

    HFEtaskDG1->SetMultiProfileLHC16s(multEstimatorAvgDG1[0]);
    HFEtaskDG1->SetMultiProfileLHC16r(multEstimatorAvgDG1[1]);
            
            
     TString finDirname         = "_TrigGADG1";
     TString outBasicname       = "EID";
            
     finDirname 	      +=   suffixName.Data();
     outBasicname      +=   finDirname.Data();
            
            
            
     mgr->ConnectInput(HFEtaskDG1,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(HFEtaskDG1,1,mgr->CreateContainer(outBasicname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
     // DCal EGA DG2

     AliAnalysisTaskHFEMultiplicity* HFEtaskDG2 = new AliAnalysisTaskHFEMultiplicity("");
     mgr->AddTask(HFEtaskDG2);
     HFEtaskDG2->SelectCollisionCandidates(AliVEvent::kEMCEGA);
     HFEtaskDG2->SetEMCalTriggerDG2(kTRUE);
     HFEtaskDG2->SetTenderSwitch(useTender);
     HFEtaskDG2->SetClusterTypeEMC(ClsTypeEMC);
     HFEtaskDG2->SetClusterTypeDCAL(ClsTypeDCAL);
     
     
      TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
    if(!fileEstimator)  {
      AliFatal("File with multiplicity estimator not found\n");
      return;
    }
    HFEtaskDG2->SetReferenceMultiplicity(refMult);
    const Char_t* profilebasename="SPDmultDG2";
    
    const Char_t* periodNames[2] = {"LHC16s", "LHC16r"};
    TProfile* multEstimatorAvgDG2[2];
    for(Int_t ip=0; ip<2; ip++) {
      cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
      multEstimatorAvgDG2[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
      if (!multEstimatorAvgDG2[ip]) {
	AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
	return;
      }
    }

    HFEtaskDG2->SetMultiProfileLHC16s(multEstimatorAvgDG2[0]);
    HFEtaskDG2->SetMultiProfileLHC16r(multEstimatorAvgDG2[1]);
            
            
     TString finDirname         = "_TrigGADG2";
     TString outBasicname       = "EID";
            
     finDirname 	      +=   suffixName.Data();
     outBasicname      +=   finDirname.Data();
            
            
            
     mgr->ConnectInput(HFEtaskDG2,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(HFEtaskDG2,1,mgr->CreateContainer(outBasicname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
            
     }
     }
    
  
    
  return NULL;
}
