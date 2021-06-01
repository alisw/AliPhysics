#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSpherocity.h"
#include "AliAnalysisFilter.h"
#include "TInterpreter.h"
#include "TChain.h"
#include <TString.h>
#include <TList.h>
#endif

#include <string>

AliAnalysisTaskSpherocity* AddTaskSpherocity(
		bool AnalysisMC = kFALSE,
		bool PostCalib = kFALSE,
		const char* period = "16l", 
		const char* multClass = "0_1",
		const bool IsV0M = kFALSE, 
		int TrkCutMode = 0
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

	// by default, a file is open for writing. here, we get the filename
	TString fileName = AliAnalysisManager::GetCommonFileName();
	fileName += ":Output";      // create a subfolder in the file

	// now we create an instance of your task
	AliAnalysisTaskSpherocity* task = new AliAnalysisTaskSpherocity("taskHighPtDeDxpp");   
	if(!task) return 0x0;

	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	task->SetAnalysisType(type);
	task->SetAnalysisMC(AnalysisMC);
	task->SetPeriod(period);
	task->SetEstimator(IsV0M);
	task->SetMinMult(0);
	task->SetNcl(70);
	task->SetDebugLevel(0);
	task->SetEtaCut(0.8);
	task->SetAnalysisTask(PostCalib);
	task->SetTrackCutsSystVars(TrkCutMode);

	task->SetJettyCutOff(0.680);
	task->SetJettyCutOff_0(0.487);
	task->SetJettyCutOff_1(0.577);
	task->SetJettyCutOff_2(0.624);
	task->SetIsotrCutOff(0.859);
	task->SetIsotrCutOff_0(0.942);
	task->SetIsotrCutOff_1(0.913);
	task->SetIsotrCutOff_2(0.892);

	std::string buf("MyOutputContainer");

	if(IsV0M)
	buf.append("_V0M");
	else
	buf.append("_Trks");
	printf("Name: %s\n",buf.c_str());

	// add your task to the manager
	mgr->AddTask(task);
	// your task needs input: here we connect the manager to your task
	mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
	// same for the output
	mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("%s_%s_TrkCutMode_%d",buf.c_str(),multClass,TrkCutMode), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
	// in the end, this macro returns a pointer to your task. this will be convenient later on
	// when you will run your analysis in an analysis train on grid
	return task;
}

