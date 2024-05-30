#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisPtN.h"

AliAnalysisPtN* AddTaskPtN(
	TString     fNUE = "LHC20j6a",
	TString     uniqueID        = " "
	)
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
 
    mgr->AddTask(task);

    // Create ONLY the output containers for the data produced by the
	// task.  Get and connect other common input/output containers via
	// the manager as below
	//=======================================================================
	//TString fileName = AliAnalysisManager::GetCommonFileName();
	//fileName+=suffixName;
	AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
	mgr->ConnectInput (task, 0, cinput);


	TGrid::Connect("alien:");

	TObjArray *AllContainers = mgr->GetContainers();
    
	if(!AllContainers->FindObject("NUE")) {

	TFile* fileNUE=nullptr;
    AliAnalysisDataContainer *NUE = mgr->CreateContainer(Form("NUE"), TList::Class(), AliAnalysisManager::kInputContainer);
    if (fNUE.EqualTo("LHC20e3a")){
	    fileNUE = TFile::Open("alien:///alice/cern.ch/user/z/zyifan/weight/NUE/Efficiency_LHC20e3a_wSyst.root");
    } else if (fNUE.EqualTo("LHC20j6a")){
		fileNUE = TFile::Open("alien:///alice/cern.ch/user/z/zyifan/weight/NUE/Efficiency_LHC20j6a_wSyst.root");
	}
	if(!fileNUE) return task;
	TList* lstNUE = dynamic_cast<TList*>(fileNUE->Get("EffAndFD"));
    NUE->SetData(lstNUE);
	mgr->ConnectInput(task,1,NUE);
	}else{
		mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("NUE"));
	}

	AliAnalysisDataContainer* cout = mgr->CreateContainer(Form("QA_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("AnalysisResults.root:%s", uniqueID.Data()));
	mgr->ConnectOutput(task,1, cout);

    return task;
}


