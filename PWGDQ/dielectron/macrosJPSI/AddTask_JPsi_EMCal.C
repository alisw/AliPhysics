AliAnalysisTask* AddTask_JPsi_EMCal(

			Bool_t 	isMC 			= kFALSE, 
			Bool_t 	isAOD 			= kTRUE,
			char * period			= "16l",
			Int_t trigger_index=3,
            Int_t config=0,
            Bool_t     isTender             = kFALSE,
            TString estimator="profile_SPD_16l",
            TString estimatorV0="2dprofile_V0_16l",
            Bool_t local = kTRUE,
            Bool_t alienconf = kFALSE,
            Bool_t V0correction = kFALSE,
            Bool_t SPDcorrection = kFALSE,
            TString cfg = "Config_JPsi_EMCal"
			
)
{
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	
	if (!mgr) {
	::Error("AddTask_JPsi_EMCal", "No analysis manager to connect to.");
	return NULL;
	}
	
	if (!mgr->GetInputEventHandler()) {
	::Error("AddTask_JPsi_EMCal", "This task requires an input event handler");
	return NULL;
	}
	
    
    
    
    
	//_______________________
	//Config Task
    TString alienPath("alien:///alice/cern.ch/user/c/cjahnke/MacrosJPsi");
    TString alirootPath("$ALICE_PHYSICS/PWGDQ/dielectron/files");
    TString alirootPath_config("$ALICE_PHYSICS/PWGDQ/dielectron/macrosJPSI");
    
    TString configFile("");
    if(cfg.IsNull()) cfg="Config_JPsi_EMCal";
    
    // >>> local config
    if(local){
        configFile="Config_JPsi_EMCal.C";
    }
    
    // >>> aliroot config
    else if(!alienconf){
        configFile=alirootPath_config.Data();
        printf("Configuration file from aliphysics\n");
    }
    // >>> alien config
    else{
        if(!gSystem->Exec(Form("alien_cp %s/%s.C .",alienPath.Data(),cfg.Data()))) {
            configFile=gSystem->pwd();
            printf("Configuration file from alien cjahnke\n");
        }
        else {
            printf("ERROR: couldn't copy file %s/%s.C from grid \n", alienPath.Data(),cfg.Data() );
            return;
        }
    }
    // add config to path
    if(!local){
        configFile+="/";
        configFile+=cfg.Data();
        configFile+=".C";
    }
    
    // load dielectron configuration file (only once)
    if (!gROOT->GetListOfGlobalFunctions()->FindObject("Config_JPsi_EMCal")){
        gROOT->LoadMacro(configFile.Data());
    }

  
	//gROOT->LoadMacro("Config_JPsi_EMCal.C");
	AliAnalysisTask_JPsi_EMCal *task = Config_JPsi_EMCal(isMC,isAOD, period,trigger_index, config, isTender);
	
	//_______________________
	//Trigger
	if(!isMC)
	{
		if(trigger_index==0){
			task->SelectCollisionCandidates(AliVEvent::kINT7);
			
		}
		if(trigger_index==1){
			task->SelectCollisionCandidates(AliVEvent::kEMC7);
			
		}
		if(trigger_index==3){
			task->SelectCollisionCandidates(AliVEvent::kEMCEGA);
			
		}
		if(trigger_index==4){
			task->SelectCollisionCandidates(AliVEvent::kEMCEGA);
		}
		
		if(trigger_index==6){
			task->SelectCollisionCandidates(AliVEvent::kEMCEGA);
		}
		
		//For DCAL trigger -> Separated by trigger word, not physics selection!!!
		if(trigger_index==7){
			task->SelectCollisionCandidates(AliVEvent::kEMCEGA);
		}
		if(trigger_index==8){
			task->SelectCollisionCandidates(AliVEvent::kEMCEGA);
		}
	
		
    }
	
	
    
//===========================================================================================================================
//Including SPD corrections
//TString estimatorFilename="estimators.root";
   
    if(SPDcorrection){
        TString estimatorFilename("");
    
        if(local){
            estimatorFilename="profile_SPD_16l.root";
        }
        else if(!alienconf){
            estimatorFilename=alirootPath.Data();
        }
        else{
        
            if(!gSystem->Exec(Form("alien_cp %s/%s.root .",alienPath.Data(),estimator.Data()))) {
                estimatorFilename=gSystem->pwd();
                printf("Not sure if copy or not from alien!!!\n");
            }
            else {
                printf("ERROR: couldn't copy file %s/%s.root from grid \n", alienPath.Data(),estimator.Data() );
                return;
            }
        }
    
    
        if(!local){
            estimatorFilename+="/";
            estimatorFilename+=estimator.Data();
            estimatorFilename+=".root";
        }
    
        printf("SPD estimator name: %s\n",estimatorFilename.Data());
    
    
    
    
        TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
        if(!fileEstimator)  {
            AliFatal("File with multiplicity estimator not found ");
            return;
        }
        // task->SetReferenceMultiplicity(refMult);
        //const Char_t* profilebasename="SPDmult10";
        const Char_t* profilebasename="profile_SPD";
    
        // const Char_t* periodNames[2] = {"LHC16l", "LHC16k"};
        //const Char_t* periodNames[1] = {"LHC16l"};
        const Char_t* periodNames[1] = {"16l"};
        //TProfile* multEstimatorAvg[2];
        TProfile* multEstimatorAvg[1];
        for(Int_t ip=0; ip<1; ip++) {
            cout<< " Estimator used (test): "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
            multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
            if (!multEstimatorAvg[ip]) {
                AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
            return;
            }
        }
    
        task->SetMultiProfileLHC16l(multEstimatorAvg[0]);
        
    }//close if SPDcorrection
    
//===========================================================================================================================
//including V0 corrections
    
  // TString estimatorV0Filename="2dprofile_V0_16l.root";
    if(V0correction){
    
        TString estimatorV0Filename("");
    
        if(local){
            estimatorV0Filename="2dprofile_V0_16l.root";
        }
        else if(!alienconf){
            estimatorV0Filename=alirootPath.Data();
        }
        else{
        
            if(!gSystem->Exec(Form("alien_cp %s/%s.root .",alienPath.Data(),estimatorV0.Data()))) {
                estimatorV0Filename=gSystem->pwd();
                printf("Not sure if copy V0  or not from alien!!!\n");
            }
            else {
                printf("ERROR: couldn't copy file %s/%s.root from grid \n", alienPath.Data(),estimatorV0.Data() );
                return;
            }
        }
    
    
        if(!local){
            estimatorV0Filename+="/";
            estimatorV0Filename+=estimatorV0.Data();
            estimatorV0Filename+=".root";
        }
  
    
        printf("V0 estimator name: %s\n",estimatorV0Filename.Data());
    
   
        TFile* fileEstimatorV0=TFile::Open(estimatorV0Filename.Data());
        if(!fileEstimatorV0)  {
            AliFatal("File with multiplicity estimator for V0 not found\n");
            return;
        }
    
        TProfile2D* multEstimatorV0[1];
        for(Int_t iv=0; iv<1; iv++) {
            cout<< " Estimator V0 used (test): "<<Form("%s",estimatorV0Filename.Data())<<endl;
            multEstimatorV0[iv] = (TProfile2D*)(fileEstimatorV0->Get("2dprofile_V0_16l;1")->Clone("2dprofile_V0_16l_clone"));
            if (!multEstimatorV0[iv]) {
                AliFatal("Multiplicity estimator for V0 not found! Please check your estimator file");
                return;
            }
        }
        task->SetMultiProfileV0LHC16l(multEstimatorV0[0]);
    }//close if V0 correction
	
	mgr->AddTask(task);
    

	
	//Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("hist_%d", trigger_index),  TList::Class(),    AliAnalysisManager::kOutputContainer, "OutPutData.root");

	//Connect input/output
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput);
    
    mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("Multi_%d", trigger_index), TList::Class(), AliAnalysisManager::kOutputContainer, "OutPutData.root"));
	
	return task;
}
