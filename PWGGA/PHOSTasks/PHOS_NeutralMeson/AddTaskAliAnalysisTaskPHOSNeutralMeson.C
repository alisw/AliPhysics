void AddTaskAliAnalysisTaskPHOSNeutralMeson(
		TString TaskName, 
		TString ContName, 
		TString CollisionCandidates, 
		TString BadMapName = "defaultTenderBM", 
		Bool_t UsePHOSTender = true, 
		Bool_t ApplyBadMap_manually = false,
		Float_t DistToBadCell = 0.0) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAliFraNeutralMeson", "No analysis manager to connect to.");
    return NULL;
  }  
  
	// ****************************************************** 
	//  *****  SETTINGS ***** ***** ***** ***** ***** ***** 
	Bool_t usePHOSTender = UsePHOSTender;
	Bool_t applyBadMap_manually = ApplyBadMap_manually;
	Bool_t usePhysicsSelection = false;
	TString collisionCandidates = CollisionCandidates; //kPHI7, kINT7, kAny, kMB, kERT1
	TString badMapName = BadMapName;    				
															// BM13c_6_15_allMd_all_noMod2 = BM_d_13c_22_6_2015_allMod_allMod2Off.root (old name)
															// defaultTenderBM  (switches off SetForceUsingBadMap for Tender!)
															// BM_d_13c_22_6_2015_allMod_all.root    BM_d_13c_28_4_15
	Bool_t DoDistToBadCell = false;
	if (DistToBadCell > 0.00001) DoDistToBadCell = true;
																			
	// ***** ***** ***** ***** ***** ***** ***** ***** ***** 
	// ****************************************************** 
	
	TString pathToBadMap;
	
	if (badMapName !="defaultTenderBM") {
		gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/m/mhecker/BadMaps/%s.root .",badMapName.Data()));
		pathToBadMap = Form("%s/",gSystem->pwd());
		pathToBadMap += badMapName;
		pathToBadMap += ".root";
	}

	// **** PHOS TENDER  ***********
	if(usePHOSTender) {
		cout<<"######  PHOS TENDER ! ######"<<endl;

		
		// ~~~~~~~~~~~~~ Loading Tender for AODs ~~~~~~~~~~~~~~~~~~

		Int_t tenderPassNumber = 1;  // does not matter (just has to be >= 0). tender then finds the correct pass itself
		gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C");
		AliPHOSTenderTask *tenderPHOS = AddAODPHOSTender("PHOSTenderTask","PHOStender","",tenderPassNumber);
		AliPHOSTenderSupply *PHOSSupply = tenderPHOS->GetPHOSTenderSupply();
		
		if(badMapName !="defaultTenderBM") {
			printf("Setting Tender-BadMap to private BM \n");
			PHOSSupply->ForceUsingBadMap(pathToBadMap.Data());
		} else printf("Running Tender with its default BadMap \n");
	}
  
  // ********** PHYSICS SELECTION ******************
	if(usePhysicsSelection) {
		gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
		AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
		printf(" ********* AliPhysicsSelectionTask configured ********* ");
	}

  
	// Initialize Task
	AliAnalysisTaskPHOSNeutralMeson *task = new AliAnalysisTaskPHOSNeutralMeson(TaskName.Data());
	
	// Set variables for Outputfile Name
	if(usePHOSTender) { 
		task->SetSelectedSelectedTenderForOutputName("PHSTndr");
	} else {
		task->SetSelectedSelectedTenderForOutputName("noTndr");
	}
	
	if(usePhysicsSelection) {
		task->SetSelectedPhysSelectionForOutputName("AliPhySl");
	} else {
		task->SetSelectedPhysSelectionForOutputName("noPhySl");
	}
	
	task->SetUsedBadmapForOutputName(badMapName);
 
 
	// ******* SET BAD MAP ******************
	task->SetApplyBadMapManually(applyBadMap_manually);
	
	if(applyBadMap_manually) {
		if (badMapName == "defaultTenderBM") AliError("Cannot apply default tender bad map in task, now applying empty bad map. Specify own bad map to fix it.");
		else {
			TFile *fBadMap = TFile::Open(pathToBadMap.Data());
			
			if(fBadMap->IsOpen()){
				printf("\n\n...Adding PHOS bad channel map (MANUALY) \n") ;
				gROOT->cd();
				Char_t key[55] ;
				for(Int_t mod=1;mod<4; mod++){
					sprintf(key,"PHOS_BadMap_mod%d",mod) ;
					TH2I * h = (TH2I*)fBadMap->Get(key) ;
					if(h)
					task->SetPHOSBadMap(mod,h) ;
				}
				fBadMap->Close() ;
			}
		}
	}
	// ***** END OF SET BAD MAP ***************
  
	//~~~~~~~ Geometry ~~~~~~~~~~~~~~
	//task->SetEtaAccMinMax(-0.14, 0.14); 	    // Set Min and Max of eta for PHOS / EMCAL acceptance (MC)
	//task->SetPhiAccMinMax(-1.72568, -0.72635); // Set Min and Max of phi for PHOS / EMCAL acceptance (MC)
  
	//~~~~~~~~ MC related ~~~~~~~~~~~~~
	task->SetFillFeedDownHistos(false); 
  
	// ~~~~~~ Event Cuts ~~~~~~~~~~~~~
	task->SetUseOnlyCINT1events(false);
	task->SetDoZvertexCut(true);
	task->SetUseIsVertexSelected2013pA(true); 
	task->SetZvertexCut(10.0);										
	
	
	// ~~~~~~ Cluster Cuts ~~~~~~~~~~~
	task->SetClusterMinCells(3);
	task->SetClusterMinE(0.3); 	//GeV
	task->SetClusterMinM02(0.2);  
	task->SetDoTimingCut(true);	//for MC: false
	task->SetTimingCutMinMax(-0.10e-6, 0.10e-6);

	task->SetDoDistToBadCellCut(DoDistToBadCell);  
	task->SetDistToBadCell(DistToBadCell);  // [cm]  

	// ~~~~~~~~ Filling options ~~~~~~
	task->SetFillCellIdVsE(false);   	//Neccessary for QA (BadMap) (increases Outputfile-size)
	task->SetFillTimingHistos(true);  	//Neccessary to determine timing cut values.
	task->SetFillClusterHitmaps(true);	// To check the BadMaps and cluster distributions
	task->SetFillhMassPtModules(false); 	// For module-wise pi0 analysis
	task->SetFillhMassPtTiming(false);  	// For calculation of timing efficiency
    
	task->SetFillNTupelClusterEnergy(false);  //Neccesary to calculate module dependent energy recal factors
  
	// ~~~~ recalibration options ~~~~
	//task->SetDoClusterEnergyRecalibration(false);
	//task->SetRecalibrationOption("lhc12d_12h");
	task->SetRecalibrateModuleWise(false);
	task->SetModuleWiseRecalibrationFactors(1.0,1.0,1.0); 
	//task->SetModuleWiseRecalibrationFactors(0.9942,0.9822,1.0072); //pp
  
  
	// SelectCollisionCandidates according to collisionCandidates set in Settings (line 14)
	if(collisionCandidates == "kPHI7") {
		task->SelectCollisionCandidates(AliVEvent::kPHI7);
		task->SetSelectedCollCandidatesForOutputName("kPHI7");
	}
	
	if(collisionCandidates == "kINT7") {
		task->SelectCollisionCandidates(AliVEvent::kINT7);
		task->SetSelectedCollCandidatesForOutputName("kINT7");
	}
	
	if(collisionCandidates == "kAny") {
		task->SelectCollisionCandidates(AliVEvent::kAny);
		task->SetSelectedCollCandidatesForOutputName("kAny");
	}
	
	if(collisionCandidates == "kMB") {
		task->SelectCollisionCandidates(AliVEvent::kMB);
		task->SetSelectedCollCandidatesForOutputName("kMB");
	}
	
	if(collisionCandidates == "kERT1") {
		task->SelectCollisionCandidates(AliVEvent::kERT1);
		task->SetSelectedCollCandidatesForOutputName("kERT1");
	}


	// ********************* IMPORTANT: ****************************
	// The OutputfileName cannot be to long!
	// The Gridjobs crash in a ERROR_RE if it is too long
	// It has been successfully tested with a filename length of 135 symbols
	// *************************************************************
	// make OutPutFileName after making all Settings
	TString OutPutFileName = task->MakeOutputFileName();
	
	mgr->AddTask(task);
	AliAnalysisDataContainer *coutput0 = 
	mgr->CreateContainer(ContName.Data(),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         OutPutFileName);
	mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task,1,coutput0);
    
	RequestMemory(task,500*1024); // request 500mb memory for task
} 
