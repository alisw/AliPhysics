
#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSpectraRT.h"
#include "AliAnalysisFilter.h"
#include "TInterpreter.h"
#include "TChain.h"
#include <TString.h>
#include <TList.h>
#endif

#include <string>


////class AliAnalysisDataContainer;

AliAnalysisTaskSpectraRT* AddTaskSpectraRT(
		bool AnalysisMC = kFALSE,
		bool MCClosure  = kFALSE,
		const char* Period  = "l",
		const double MeanMultTSdata = 7.266,
		const double MeanMultTSMCgen = 7.266,
		const double MeanMultTSMCrec = 7.266,
		bool IsTPCOnlyTrkCuts = kTRUE,
		const char* Container = "TPCOnly"		
		)   
{

	// get the manager via the static access member. since it's static, you don't need
	// an instance of the class to call the function
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

	AliAnalysisFilter* trackFilterGolden = new AliAnalysisFilter("trackFilter");
	AliESDtrackCuts* esdTrackCutsGolden = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
	trackFilterGolden->AddCuts(esdTrackCutsGolden);


	// by default, a file is open for writing. here, we get the filename
	TString fileName = AliAnalysisManager::GetCommonFileName();

	// now we create an instance of your task
	AliAnalysisTaskSpectraRT* task = new AliAnalysisTaskSpectraRT("taskHighPtDeDxpp");   
	if(!task) return 0x0;

	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	task->SetAnalysisType(type);
	task->SetAnalysisMC(AnalysisMC);
	task->SetMCClosure(MCClosure);
	task->SetPeriod(Period);
	task->SetMeanMultTSdata(MeanMultTSdata);
	task->SetMeanMultTSMCGen(MeanMultTSMCgen);
	task->SetMeanMultTSMCRec(MeanMultTSMCrec);
 	task->SetTrackCutsType(IsTPCOnlyTrkCuts);

	task->SetNcl(70);
	task->SetDebugLevel(0);
	task->SetEtaCut(0.8);
	task->SetTrackFilterGolden(trackFilterGolden);

	// add your task to the manager
	mgr->AddTask(task);
	// your task needs input: here we connect the manager to your task
	mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
	// same for the output
	mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("MyOutputContainer_PID_5GeV_%s",Container), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
	// in the end, this macro returns a pointer to your task. this will be convenient later on
	// when you will run your analysis in an analysis train on grid
	return task;
}

