#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisPtN.h"

AliAnalysisPtN* AddTaskPtN(TString uniqueID = "")
{
    // Creates a pid task and adds it to the analysis manager
	// Get the pointer to the existing analysis manager via the static
	//access method
	//=========================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskPtN.C", "No analysis manager to connect to.");
		return NULL;
	}

	// Check the analysis type using the event handlers connected to the
	// analysis manager The availability of MC handler can also be
	// checked here.
	// =========================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTaskPtN.C", "This task requires an input event handler");
		return NULL;
	}

	// Create the task and configure it
	//========================================================================
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

	AliAnalysisPtN* task = new AliAnalysisPtN("PtN");
    if(!task) return 0x0;
    mgr->AddTask(task);

    // Create ONLY the output containers for the data produced by the
	// task.  Get and connect other common input/output containers via
	// the manager as below
	//=======================================================================
	//TString fileName = AliAnalysisManager::GetCommonFileName();
	//fileName+=suffixName;
	AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();

	TGrid::Connect("alien:");
    
    AliAnalysisDataContainer *NUE = mgr->CreateContainer(Form("NUE"), TList::Class(), AliAnalysisManager::kInputContainer);
    TFile* fileNUE = TFile::Open("alien:///alice/cern.ch/user/z/zyifan/weight/NUE/Efficiency_LHC20e3a_wSyst.root");
    TList* lstNUE = dynamic_cast<TList*>(fileNUE->Get("EffAndFD"));
    NUE->SetData(lstNUE);

	// AliAnalysisDataContainer *NUE = mgr->CreateContainer(Form("NUE"), TList::Class(), AliAnalysisManager::kInputContainer);
    // TFile* fileNUE = TFile::Open("/Users/macnbi/alice/Efficiency_LHC20e3a_wSyst.root");
    // TList* lstNUE = dynamic_cast<TList*>(fileNUE->Get("EffAndFD"));
    // NUE->SetData(lstNUE); 

	AliAnalysisDataContainer* cout = mgr->CreateContainer(Form("QA_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("AnalysisResults.root:%s", uniqueID.Data()));
	mgr->ConnectInput (task, 0, cinput);
	mgr->ConnectInput(task,1,NUE);
	mgr->ConnectOutput(task, 1, cout);

    return task;
}


