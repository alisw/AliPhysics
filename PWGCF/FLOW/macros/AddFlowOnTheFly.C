#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliAnalysisTaskFlowOnTheFly.h"
#endif

AliAnalysisTaskFlowOnTheFly* AddFlowOnTheFly(
		Double_t	fMinPt			= 0.2,
		Double_t	fMaxPt			= 3.0,
        Int_t           fSystFlag               = 0,
		Bool_t		fUseImpactXaxis	= false,
		TString		uniqueID        	= "Default"
    )
{
	TString name = "MyFlowOntheFlyTask";
     // The common parameters
	Double_t	fEtaCut 			= 0.8;
	Double_t	fVtxCut				= 10.0;


    // get the manager via the static access member. since it's static, you don't need
    // to create an instance of the class here to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }

    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyResults";      
	//fileName += uniqueID; // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskFlowOnTheFly* task = new AliAnalysisTaskFlowOnTheFly(name.Data());   
    if(!task) return 0x0;
	task->SetDebugLevel(3);
	task->SetEtaCut(fEtaCut);
	task->SetVtxCut(fVtxCut); // For systematics
	task->SetVtxCutDefault(fVtxCut);
	task->SetMinPt(fMinPt);
	task->SetMaxPt(fMaxPt);
    task->SetSystFlag(fSystFlag);
	task->SetfUseImpactXaxis(fUseImpactXaxis);
    
    // add your task to the manager
    mgr->AddTask(task);

    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid

    return task;
}
