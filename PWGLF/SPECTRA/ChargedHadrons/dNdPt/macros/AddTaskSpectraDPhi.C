/*
  Author: Aditya Nath Mishra Wigner RCP (amishra@cern.ch)
Last update: 05/05/2019 
 */

AliAnalysisTaskSpectraDPhi* AddTaskSpectraDPhi(
		Bool_t AnalysisMC = kFALSE, 
		const Char_t* taskname = "SpectraDPhi"
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
		Error("AddTaskSpectraDPhi", "This task requires an input event handler");
		return NULL;
	}  

	// Create the task and configure it 
	//========================================================================
	AliAnalysisTaskSpectraDPhi* taskSpectraDPhi = new AliAnalysisTaskSpectraDPhi("SpectraDPhi");
	if(!taskSpectraDPhi) return 0x0;
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	
	Bool_t is13TeV = kTRUE;
	UInt_t trigSel = AliVEvent::kINT7;
	
	taskSpectraDPhi->SetTrigger(trigSel);
	taskSpectraDPhi->SetAnalysisMC(AnalysisMC);
	taskSpectraDPhi->SetAnalysisType(type);
	taskSpectraDPhi->SetDebugLevel(0);
	taskSpectraDPhi->SetEtaCut(0.8);
	taskSpectraDPhi->SetVtxCut(10.0);
	taskSpectraDPhi->SetPtLeadMin(0.0);
	taskSpectraDPhi->SetPileUpRej(kTRUE);	
	//	taskSpectraDPhi->SetAveMultiInTrans(4.939);
	//task->SetAveRecMultiInTrans(5.003); // reco MC EPOS LHC
	//task->SetAveGenMultiInTrans(7.72); // true MC EPOS LHC
	mgr->AddTask(taskSpectraDPhi);


	mgr->ConnectInput(taskSpectraDPhi,0,mgr->GetCommonInputContainer());
	// same for the output
	mgr->ConnectOutput(taskSpectraDPhi,1,mgr->CreateContainer(Form("cList%s",taskname), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));
	
	return taskSpectraDPhi;
	
}
