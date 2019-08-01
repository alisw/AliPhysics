/*
Author: Antonio Ortiz (aortizve@cern.ch, antonio.ortiz@nucleares.unam.mx)
Last update: 29/07/2019 
 */

AliAnalysisTaskUeRSpherocity* AddTaskUeRSpherocity(
		Bool_t AnalysisMC = kFALSE, 
		const Char_t* taskname = "UeRSpherocity" 
		)
{

	// Get the pointer to the existing analysis manager via the static
	//access method
	//=========================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskESA", "No analysis manager to connect to.");
		return NULL;
	}  

	// Check the analysis type using the event handlers connected to the
	// analysis manager The availability of MC handler can also be
	// checked here.
	// =========================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTaskESA", "This task requires an input event handler");
		return NULL;
	}  

	// Create spherocity percentile binning 
	//========================================================================
	const Int_t nSobins = 10;
	Float_t percSpheroBins[nSobins+1]={ 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0 };
	TH1D * hSpheroDef = new TH1D("hHelper","",nSobins,percSpheroBins);

	// Create the task and configure it 
	//========================================================================
	AliAnalysisTaskUeRSpherocity* taskESA = new AliAnalysisTaskUeRSpherocity("SpherocityTask");
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	taskESA->SetAnalysisType(type);
	//taskESA->SetAnalysisMC(AnalysisMC);
	taskESA->SetDebugLevel(0);
	//taskESA->SetHisto(hSpheroDef);
	mgr->AddTask(taskESA);

	mgr->ConnectInput(taskESA,0,mgr->GetCommonInputContainer());
	// same for the output
	mgr->ConnectOutput(taskESA,1,mgr->CreateContainer(Form("cList%s",taskname), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));

	return taskESA;

}
