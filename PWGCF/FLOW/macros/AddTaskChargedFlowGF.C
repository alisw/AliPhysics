AliAnalysisTaskChargedFlowGF* AddTaskChargedFlowGF(
	Int_t			fFilterbit 		= 768,
	Double_t	fEtaCut 			= 1.0,
	Double_t	fVtxCut				= 10.0,
	Double_t	fMinPt				= 0.2,
	Double_t	fMaxPt				= 3.0,
	Int_t			fTPCclusters	= 70,
	Bool_t		fUseDCAzCut		= false,
	Double_t	fDCAz					= 1.0,
	Bool_t		fUseDCAxyCut	= false,
	Double_t	fDCAxy				= 0.2,
	Int_t			IsSample			= 10,
	Int_t			trigger				= 1,
	Bool_t		fNUE 					= false, 
	Bool_t		fNUA					= true,
		//..MC
	Bool_t		fIsMC 				= false,
		//....
	TString		fPeriod 			= "",
	TString		uniqueID 			= ""
)
{
  // Creates a pid task and adds it to the analysis manager
  
  // Get the pointer to the existing analysis manager via the static
  //access method
  //=========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskChargedFlowGF.C", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskChargedFlowGF", "This task requires an input event handler");
    return NULL;
  }  
  
  // Create the task and configure it 
  //========================================================================
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
  Double_t ptBins[] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0};
  //Double_t ptBins[] = {0.2, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
    
  AliAnalysisTaskChargedFlowGF* taskFlowEp = new AliAnalysisTaskChargedFlowGF("taskFlowEp");
  taskFlowEp->SetDebugLevel(3);
	taskFlowEp->SetFilterbit(fFilterbit);
  taskFlowEp->SetEtaCut(fEtaCut);
  taskFlowEp->SetVtxCut(fVtxCut);
  taskFlowEp->SetMinPt(fMinPt);
  taskFlowEp->SetMaxPt(fMaxPt);
	taskFlowEp->SetTPCclusters(fTPCclusters);
	taskFlowEp->SetUseDCAzCut(fUseDCAzCut);
	taskFlowEp->SetDCAzCut(fDCAz);
	taskFlowEp->SetUseDCAxyCut(fUseDCAxyCut);
	taskFlowEp->SetDCAxyCut(fDCAxy);
  taskFlowEp->SetIsSample(IsSample);
	taskFlowEp->SetTrigger(trigger);
	taskFlowEp->SetNUEFlag(fNUE);
	taskFlowEp->SetNUA(fNUA);
		//..MC
	taskFlowEp->SetMCFlag(fIsMC);
		//....
	taskFlowEp->SetPeriod(fPeriod);
  mgr->AddTask(taskFlowEp);
 
  
  // Create ONLY the output containers for the data produced by the
  // task.  Get and connect other common input/output containers via
  // the manager as below
   //=======================================================================
  //TString fileName = AliAnalysisManager::GetCommonFileName();
  //fileName+=suffixName;
  AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *cout_hist = mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("v4New.root:%s", uniqueID.Data()));
  if(fIsMC == true) AliAnalysisDataContainer *cout_histMC = mgr->CreateContainer(Form("outputMC_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("v4New.root:%s", uniqueID.Data()));
  mgr->ConnectInput (taskFlowEp, 0, cinput);
  mgr->ConnectOutput(taskFlowEp, 1, cout_hist);
  if(fIsMC == true) mgr->ConnectOutput(taskFlowEp, 2, cout_histMC);

  // Return task pointer at the end
  return taskFlowEp;
}
// 
// EOF
// 
