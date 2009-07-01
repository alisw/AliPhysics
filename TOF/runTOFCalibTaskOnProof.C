void runTOFCalibTaskOnProof() {
	TStopwatch timer;
	timer.Start();
	
	printf("*** Open PROOF ***");

	gEnv->SetValue("XSec.GSI.DelegProxy","2");
	TProof::Open("alicecaf");
		
	gProof->UploadPackage("STEERBase.par");
	gProof->EnablePackage("STEERBase");
	gProof->UploadPackage("ESD.par");
	gProof->EnablePackage("ESD");
	gProof->UploadPackage("AOD.par");
	gProof->EnablePackage("AOD");
	gProof->UploadPackage("ANALYSIS.par");
	gProof->EnablePackage("ANALYSIS");
	gProof->UploadPackage("ANALYSISalice.par");
	gProof->EnablePackage("ANALYSISalice");
	
	gProof->Load("AliTOFArray.cxx++g");
	gProof->Load("AliTOFCalibTask.cxx++g");
	gROOT->LoadMacro("AddTaskTOFCalib.C");
	cout << "Loaded AddTaskTOFCalib macro "<< endl;
	
	gProof->ShowEnabledPackages();
	
	//ANALYSIS PART	
	//____________________________________________//
	// Make the analysis manager

	AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
	AliESDInputHandler* esdH = new AliESDInputHandler;
	
	esdH->SetInactiveBranches("FMD CaloCluster");
	mgr->SetInputEventHandler(esdH);  
	
	Bool_t domc = kFALSE;
	if (domc) {
		AliMCEventHandler *mcH = new AliMCEventHandler;
		mgr->SetMCtruthEventHandler(mcH);
	}

	//____________________________________________//
	// 1st TOFCalib task
	
	AliTOFCalibTask *taskTOFCalib = AddTaskTOFCalib();
	
	if (!mgr->InitAnalysis()) return;
	mgr->PrintStatus();
	mgr->StartAnalysis("proof","/COMMON/COMMON/LHC08c11_10TeV_0.5T",1000);
	
	timer.Stop();
	timer.Print();
}

