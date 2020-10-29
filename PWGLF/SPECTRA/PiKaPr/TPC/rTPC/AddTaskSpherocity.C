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
		float JettyValue = 0.5,
		float JettyValue_0 = 0.5,
		float JettyValue_1 = 0.5,
		float JettyValue_2 = 0.5,
		float IsotrValue = 0.7,
		float IsotrValue_0 = 0.7,
		float IsotrValue_1 = 0.7,
		float IsotrValue_2 = 0.7
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
	task->SetJettyCutOff(JettyValue);
	task->SetJettyCutOff_0(JettyValue_0);
	task->SetJettyCutOff_1(JettyValue_1);
	task->SetJettyCutOff_2(JettyValue_2);
	task->SetIsotrCutOff(IsotrValue);
	task->SetIsotrCutOff_0(IsotrValue_0);
	task->SetIsotrCutOff_1(IsotrValue_1);
	task->SetIsotrCutOff_2(IsotrValue_2);
	task->SetNcl(70);
	task->SetDebugLevel(0);
	task->SetEtaCut(0.8);
	task->SetAnalysisTask(PostCalib);

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
	mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("%s_%s",buf.c_str(),multClass), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
	// in the end, this macro returns a pointer to your task. this will be convenient later on
	// when you will run your analysis in an analysis train on grid
	return task;
}

