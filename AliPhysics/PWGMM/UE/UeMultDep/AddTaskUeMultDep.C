/*
  Author: Aditya Nath Mishra ICN UNAM (amishra@cern.ch, aditya.mishra@correo.nucleares.unam.mx)
Last update: 07/03/2019 
 */

AliAnalysisTaskUeMultDep* AddTaskUeMultDep(
		Bool_t AnalysisMC = kFALSE, 
		const Char_t* taskname = "UeMultDep"
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
	const Int_t nMultbins = 200;
	Double_t Multbins[nMultbins+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
				       21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,
				       39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,
				       57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,
				       75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,
				       93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,
				       109,110,111,112,113,114,115,116,117,118,119,120,121,122,
				       123,124,125,126,127,128,129,130,131,132,133,134,135,136,
				       137,138,139,140,141,142,143,144,145,146,147,148,149,150,
				       151,152,153,154,155,156,157,158,159,160,161,162,163,164,
				       165,166,167,168,169,170,171,172,173,174,175,176,177,178,
				       179,180,181,182,183,184,185,186,187,188,189,190,191,192,
				       193,194,195,196,197,198,199,200};
		
	TH1D * hMultDef = new TH1D("hHelperMult","",nMultbins,Multbins);


	// Create the task and configure it 
	//========================================================================
	AliAnalysisTaskUeMultDep* taskUeMultDep = new AliAnalysisTaskUeMultDep("UeMultDep");
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	taskUeMultDep->SetAnalysisType(type);
	taskUeMultDep->SetDebugLevel(0);
	mgr->AddTask(taskUeMultDep);

	mgr->ConnectInput(taskUeMultDep,0,mgr->GetCommonInputContainer());
	// same for the output
	mgr->ConnectOutput(taskUeMultDep,1,mgr->CreateContainer(Form("cList%s",taskname), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));
	
	return taskUeMultDep;
	
}
