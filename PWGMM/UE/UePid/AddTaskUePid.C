/*
Author: Antonio Ortiz (aortizve@cern.ch, antonio.ortiz@nucleares.unam.mx)
Last update: 11/09/2018 
 */

AliAnalysisTaskUePid* AddTaskUePid(
		Bool_t AnalysisMC = kFALSE, 
		const Char_t* taskname = "UePid",
		TString periodo = "k" 
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

	// Create multiplicity binning 
	//========================================================================
	const Int_t nMultbins = 35;
	Double_t Multbins[nMultbins+1]={
		1.00, 3.00, 5.00, 7.00, 9.00, 11.0, 13.0, 15.0, 17.0, 19.0,
		21.0, 23.0, 25.0, 27.0, 29.0, 31.0, 33.0, 35.0, 37.0, 39.0,
		41.0, 44.0, 47.0, 50.0, 53.0, 56.0, 59.0, 62.0, 65.0, 69.0,
		73.0, 77.0, 81.0, 86.0, 92.0, 100
	};
	TH1D * hMultDef = new TH1D("hHelperMult","",nMultbins,Multbins);


	// Create the task and configure it 
	//========================================================================
	AliAnalysisTaskUePid* taskUePid = new AliAnalysisTaskUePid("UePidTask");
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	taskUePid->SetAnalysisType(type);
	taskUePid->SetDebugLevel(0);
	taskUePid->SetPeriod(periodo);
	mgr->AddTask(taskUePid);

	mgr->ConnectInput(taskUePid,0,mgr->GetCommonInputContainer());
	// same for the output
	mgr->ConnectOutput(taskUePid,1,mgr->CreateContainer(Form("cList%s",taskname), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));

	return taskUePid;

}
