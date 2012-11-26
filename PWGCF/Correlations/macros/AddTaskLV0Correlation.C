AliLeadingV0Correlation* AddTaskLV0Correlation(TString fListName                = "LV0Correlation",
					       TString fCollisiontype           = "PbPb",
					       Bool_t  fAnalysisMC              = 0,
					       Double_t aPVzCut                  = 10,
					       UInt_t   aFilterBit               = 128,
					       Double_t aCutV0Radius             = 0.5,
					       Double_t aCutDCANegToPV           = 0.06,
					       Double_t aCutDCAPosToPV           = 0.06,
					       Double_t aCutDCAV0Daughters       = 1.0,
					       Double_t aCutV0CosPA              = 0.995,
					       Double_t aSpecialArmenterosCutK0s = 5
					       ) 
{
	// Get the current analysis manager.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {Error("AddTaskLV0Correlation.C", "No Analysis Manager ");return 0;}
	
	//PVz Binning for pool	
	Double_t pvzbinlimits[] = {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};
	Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
	
	//Cent(PbPb) or Mult(PP) Binning for pool	
	Double_t cent_mult_binlimits[] = {0,20,40,60,80,100}; 
	Int_t cent_mult_bin_numb = sizeof(cent_mult_binlimits)/sizeof(Double_t) - 1;

	//Correlation task
	AliLeadingV0Correlation *myTask = new AliLeadingV0Correlation(fListName.Data());
	
	myTask->SetCollidingSystem(fCollisiontype);
	myTask->SetMCAnalysis(fAnalysisMC);
	myTask->SetMaxNEventsInPool(2000);
	myTask->SetMinNTracksInPool(1000);
	myTask->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);
	myTask->SetPoolCentBinLimits(cent_mult_bin_numb,cent_mult_binlimits);
	myTask->SetPrimeryVertexCut(aPVzCut);
	myTask->SetEatCut(0.9);
	myTask->SetFilterBit(aFilterBit);
	myTask->SetTrigPtBinLimits(6.0,12.0);
	myTask->SetAssocPtBinLimits(0.15,6.0);
	myTask->SetCutV0Radius(aCutV0Radius);
	myTask->SetCutDCANToP(aCutDCANegToPV);
	myTask->SetCutDCAPToP(aCutDCAPosToPV);
	myTask->SetCutDCADaughters(aCutDCAV0Daughters);
	myTask->SetCutCPA(aCutV0CosPA);
	myTask->SetCutArmenterosK0s(aSpecialArmenterosCutK0s);
	
	mgr->AddTask(myTask);
	
	// Create containers for input/output
	TString outputFileName = AliAnalysisManager::GetCommonFileName();
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer(fListName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
	
	//____________________________________________//
	mgr->ConnectInput(myTask,0,cinput);	
	mgr->ConnectOutput(myTask,1,coutput);
	
return myTask;
}


