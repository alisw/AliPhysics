#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliAnalysisTaskAcceptanceCorr.h"
#endif

AliAnalysisTaskAcceptanceCorr* AddTaskAcceptanceCorr(
    Int_t           trigger                 = 0,
    Int_t           fSystFlag               = 0,
    TString         fPeriod                 = "LHC15o",
    TString         fNtrksName              = "Mult",
		TString		      uniqueID        	= ""
		)
{
        // The common parameters
	Double_t	fEtaCut 			= 0.8;

	// Creates a pid task and adds it to the analysis manager
	// Get the pointer to the existing analysis manager via the static
	//access method
	//=========================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskAcceptanceCorr.C", "No analysis manager to connect to."); return NULL;
	}

	// Check the analysis type using the event handlers connected to the
	// analysis manager The availability of MC handler can also be
	// checked here.
	// =========================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTaskNonLinearFlow.C", "This task requires an input event handler");
		return NULL;
	}

	// Create the task and configure it
	//========================================================================
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

	AliAnalysisTaskAcceptanceCorr* taskFlowEp = new AliAnalysisTaskAcceptanceCorr("taskFlowEp");
	taskFlowEp->SetDebugLevel(3);
	taskFlowEp->SetEtaCut(fEtaCut);
	taskFlowEp->SetMinPt(0.2);
	taskFlowEp->SetMaxPt(3.0);
		
	taskFlowEp->SetTrigger(trigger);
	taskFlowEp->SetNtrksName(fNtrksName);
  taskFlowEp->SetSystFlag(fSystFlag);

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
	mgr->ConnectInput (taskFlowEp, 0, cinput);
	AliAnalysisDataContainer *physics_hist = mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("AnalysisResults.root:%s", uniqueID.Data()));
	mgr->ConnectOutput(taskFlowEp, 1, physics_hist);

  int outSlotCounter=2;

	// Return task pointer at the end
	return taskFlowEp;
}
//
// EOF
//
