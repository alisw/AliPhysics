AliAnalysisTaskPHOSNeutralMeson* AddTaskAliAnalysisTaskPHOSNeutralMeson(
		TString taskName,
		TString collisionCandidates, 
		Bool_t usePhysicsSelection,
		Bool_t usePHOSTender,
		Int_t clusterMinCells, 
		Float_t clusterMinE,
		Float_t clusterMinM02,
		Float_t distToBadCellOnCellLevel,
		Float_t distToBadCell,
		Bool_t doTimingCut,
		Float_t timingCutMin,
		Float_t timingCutMax,
		Float_t zVertexCut,
		TString recalOption,
		Double_t recalFactorMod1,
		Double_t recalFactorMod2,
		Double_t recalFactorMod3,
		TString mPtHistoMode,
		Bool_t applyBadMap_manually,
		Bool_t useIsVertexSelected2013pA,
		TString badMapName = "defaultTenderBM",
		Bool_t fillMCHistos = kFALSE,
		Bool_t analyseAddedSignals = kFALSE,
		TString mcParticleToAnalyse = "pi0",
		Double_t expoPara1 = 0.0,
		Double_t expoPara2 = 0.0,
		Double_t expoPara3 = 0.0,
		Double_t expoPara4 = 0.0,
		TString additionalFileNameString = "",
		Bool_t fillHitMapsAddedS = kFALSE,
		Bool_t  fillDecayInfoAddedS = kFALSE, 
		Char_t *suffix = "") 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAliFraNeutralMeson", "No analysis manager to connect to.");
    return NULL;
  }  
  
	// ****************************************************** 
	//  *****  SETTINGS ***** ***** ***** ***** ***** ***** 
	// for possible bad maps see alien:///alice/cern.ch/user/m/mhecker/BadMaps/:
	// BM13c_6_15_allMd_all_noMod2 = BM_d_13c_22_6_2015_allMod_allMod2Off.root (old name)
	// defaultTenderBM  (switches off SetForceUsingBadMap for Tender!)
	// BM_d_13c_22_6_2015_allMod_all.root    BM_d_13c_28_4_15
	
	Bool_t DoDistToBadCell = false;
	if (distToBadCell > 0.00001) DoDistToBadCell = true;
	Bool_t DoDistToBadCellOnCellLevel = false;
	if (distToBadCellOnCellLevel > 0) DoDistToBadCellOnCellLevel = true;
																			
	// ***** ***** ***** ***** ***** ***** ***** ***** ***** 
	// ****************************************************** 
	
	TString ContName = "contPHOSNeutralMeson";
	TString combinedName;
	combinedName.Form("%s%s",ContName.Data(), suffix);
	
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
	AliAnalysisTaskPHOSNeutralMeson *task = new AliAnalysisTaskPHOSNeutralMeson(taskName.Data());
	
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
	//task->SetEtaAccMinMax(-0.14, 0.14); 	    
	//task->SetPhiAccMinMax(-1.72568, -0.72635); 
	task->SetEtaAccMinMax(-0.13, 0.13);				// Set Min and Max of eta for PHOS acceptance (MC)
	Double_t phiMin = 260 * (TMath::Pi()/180);
	Double_t phiMax = 320 * (TMath::Pi()/180);
	task->SetPhiAccMinMax(phiMin,phiMax); //60 degrees (in rad)
	
	
	//~~~~~~~~ MC related ~~~~~~~~~~~~~
	task->SetfFillMCHistos(fillMCHistos);
	task->SetAnalyseAddedSignals(analyseAddedSignals);
	//task->SetFillFeedDownHistos(true); 	
	task->SetAnalyseMCPi0OrEta(mcParticleToAnalyse);  //pi0 eta both 
	task->SetExponentialParameters(expoPara1, expoPara2, expoPara3, expoPara4); //for weighing added signals
  	task->SetAdditionalFileNameString(additionalFileNameString); //String appears at the end of the filname. 
  	task->SetFillClusterHitmapsAddedSignals(fillHitMapsAddedS);
  	task->SetFillDecayPhotonInfoAddedSig(fillDecayInfoAddedS);
	
	task->SetFillNDaughtersVsPtAndR(kFALSE);				//implement bool to list of arguments of AddTask if needed
	task->SetFillMPtForSingleOrMultContrClus(kFALSE);	//implement bool to list of arguments of AddTask if needed
	 
 
	// ~~~~~~ Event Cuts ~~~~~~~~~~~~~
	task->SetUseOnlyCINT1events(kFALSE);
	task->SetDoZvertexCut(kTRUE);
	task->SetUseIsVertexSelected2013pA(useIsVertexSelected2013pA); 
	task->SetZvertexCut(zVertexCut);	// 10.0									
	
	
	// ~~~~~~ Cluster Cuts ~~~~~~~~~~~
	task->SetClusterMinCells(clusterMinCells); //3
	task->SetClusterMinE(clusterMinE); 	//0.3 GeV
	task->SetClusterMinM02(clusterMinM02);  //0.2
	task->SetDoTimingCut(doTimingCut);	//for MC: false
	task->SetTimingCutMinMax(timingCutMin, timingCutMax); //-0.10e-6,0.10e-6

	task->SetDoDistToBadCellCut(DoDistToBadCell);  
	task->SetDistToBadCell(distToBadCell);  // [cm]  
	
	task->SetDoDistToBadCellCutOnCellLevel(DoDistToBadCellOnCellLevel);
	task->SetDistToBadCellOnCellLevel(distToBadCellOnCellLevel);

	// ~~~~~~~~ Filling options ~~~~~~
	task->SetFillCellIdVsE(false);   	//Neccessary for QA (BadMap) (increases Outputfile-size)
	task->SetFillTimingHistos(true);  	//Neccessary to determine timing cut values.
	task->SetFillClusterHitmaps(true);	// To check the BadMaps and cluster distributions

	if (mPtHistoMode != "nrmlMptHst") {
		if (mPtHistoMode == "hMPTvsMod") task->SetFillhMassPtModules(true); 	// For module-wise pi0 analysis
		else if (mPtHistoMode == "hMPTvsTim") task->SetFillhMassPtTiming(true);  	// For calculation of timing efficiency
		else if (mPtHistoMode == "hMPTnewAsy") task->SetFillNewAsymmClasses(true);
		else mPtHistoMode = "nrmlMptHst"; // if none of the implemented (above) modes is set, the string is set to default
	}
    
	task->SetFillNTupelClusterEnergy(false);  //Neccesary to calculate module dependent energy recal factors
  
	// ~~~~ recalibration options ~~~~
	//task->SetDoClusterEnergyRecalibration(false);
	task->SetRecalibrationOption(recalOption); // look in .cxx for possible options
	task->SetRecalibrateModuleWise(false);
	task->SetModuleWiseRecalibrationFactors(recalFactorMod1,recalFactorMod2,recalFactorMod3); 
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

	AliAnalysisDataContainer *coutput0 = mgr->CreateContainer(Form("%s:%s", ContName.Data(),OutPutFileName.Data()), TList::Class(),  AliAnalysisManager::kOutputContainer, Form("%s:dirNeutralMeson", mgr->GetCommonFileName()));
	mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task,1,coutput0);
	
	return task;
} 
