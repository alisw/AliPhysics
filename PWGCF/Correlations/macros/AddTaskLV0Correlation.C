AliLeadingV0Correlation* AddTaskLV0Correlation(TString fListName = "LV0Correlation",
					       TString  fCollisiontype           = "PbPb",
						   TString	fTriggerMask             = "kMB",
						   TString  fwithPID                 = "withPID",
					       Bool_t   fAnalysisMC              = 0,
					       Double_t aPVzCut                  = 10,
					       UInt_t   aFilterBit               = 128,
						   Double_t aRapidityCut             = 0.75,
						   Double_t aEtaCut                  = 0.9,
					       Double_t aCutV0Radius             = 0.5,
					       Double_t aCutDCANegToPV           = 0.06,
					       Double_t aCutDCAPosToPV           = 0.06,
					       Double_t aCutDCAV0Daughters       = 1.0,
					       Double_t aCutV0CosPA              = 0.995,
					       Double_t aSpecialArmenterosCutK0s = 5,
						   Double_t aCTauK0                  = 20,
						   Double_t aCTauLambda              = 30,
						   Int_t    fmaxEventsinPool         = 2000,
						   Int_t    fmaxTracksinPool         = 1000,
						   Double_t  aTrigLow                = 6,
						   Double_t  aTrigHigh               = 12) 
{
	// Get the current analysis manager.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {Error("AddTaskLV0Correlation.C", "No Analysis Manager ");return 0;}
	
	//PVz Binning for pool	
	Double_t pvzbinlimits[] = {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};
	Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
	
	//Cent(PbPb) or Mult(PP) Binning for pool	
	Double_t cent_mult_binlimits[] = {0,10,20,30,40,50,60,70,80,90,100}; 
	Int_t cent_mult_bin_numb = sizeof(cent_mult_binlimits)/sizeof(Double_t) - 1;

	//Correlation task
	AliLeadingV0Correlation *myTask = new AliLeadingV0Correlation(fListName.Data());
	
	myTask->SetCollidingSystem(fCollisiontype);
	myTask->SetTrigger(fTriggerMask);
	myTask->SetMCAnalysis(fAnalysisMC);
	myTask->SetMaxNEventsInPool(fmaxEventsinPool);
	myTask->SetMinNTracksInPool(fmaxTracksinPool);
	myTask->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);
	myTask->SetPoolCentBinLimits(cent_mult_bin_numb,cent_mult_binlimits);
	myTask->SetPrimeryVertexCut(aPVzCut);
	myTask->SetEatCut(aEtaCut);
	myTask->SetFilterBit(aFilterBit);
	myTask->SetTrigPtBinLimits(aTrigLow,aTrigHigh);
	myTask->SetAssocPtBinLimits(0.15,aTrigLow);
	myTask->SetCutV0Radius(aCutV0Radius);
	myTask->SetCutDCANToP(aCutDCANegToPV);
	myTask->SetCutDCAPToP(aCutDCAPosToPV);
	myTask->SetCutDCADaughters(aCutDCAV0Daughters);
	myTask->SetCutCPA(aCutV0CosPA);
	myTask->SetCutArmenterosK0s(aSpecialArmenterosCutK0s);
	myTask->SetCutCTauK0(aCTauK0);
	myTask->SetCutCTauLambda(aCTauLambda);
	myTask->SetUSEPID(fwithPID);
	myTask->SetCutRap(aRapidityCut);
	
	mgr->AddTask(myTask);
	
	// Create containers for input/output
	TString outputFileName = AliAnalysisManager::GetCommonFileName();
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer(fListName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
	
	//____________________________________________//
	mgr->ConnectInput(myTask,0,cinput);	
	mgr->ConnectOutput(myTask,1,coutput);
	
return myTask;
}


