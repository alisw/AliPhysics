/*   This macro produces:  pT spectra in different multiplicity and Delta phi bins
     Aditya Nath Mishra Wigner RCP, Budapest, Hungary
     Please report bugs to: amishra@cern.ch / aditya.nath.mishra@wigner.hu
     last update: 16/06/2020

*/

AliAnalysisTaskUeSpectraDphi* AddTaskUeSpectraDphi(Bool_t AnalysisMC = kFALSE,
						   const Char_t* taskname = "UeSpectraDphi")
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

	// Create the task and configure it
	//========================================================================

	UInt_t trigSel = AliVEvent::kINT7;

	AliAnalysisTaskUeSpectraDphi* taskUeSpectraDphi = new AliAnalysisTaskUeSpectraDphi("UeSpectraDphi");
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	//	taskUeSpectraDphi->SetAnalysisType(type);
	taskUeSpectraDphi->SetTrigger(trigSel);
	taskUeSpectraDphi->SetAnalysisMC(AnalysisMC);
	taskUeSpectraDphi->SetDebugLevel(0);
	taskUeSpectraDphi->SetEtaCut(0.8);
	taskUeSpectraDphi->SetVtxCut(10.0);
	taskUeSpectraDphi->SetTPCNclCut(70);
	taskUeSpectraDphi->SetPileUpRej(kTRUE);
	mgr->AddTask(taskUeSpectraDphi);

	mgr->ConnectInput(taskUeSpectraDphi,0,mgr->GetCommonInputContainer());
	// same for the output
	mgr->ConnectOutput(taskUeSpectraDphi,1,mgr->CreateContainer(Form("cList%s",taskname), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));

	return taskUeSpectraDphi;

}
