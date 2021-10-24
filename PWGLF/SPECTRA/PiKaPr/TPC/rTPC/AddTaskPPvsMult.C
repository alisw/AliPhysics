#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskPPvsMult.h"
#include "AliAnalysisFilter.h"
#include "TInterpreter.h"
#include "TChain.h"
#include <TString.h>
#include <TList.h>
#endif

AliAnalysisTaskPPvsMult* AddTaskPPvsMult(
		Bool_t AnalysisMC = kFALSE,
		////		Int_t system =1, // 0 for pp and 1 for Pb-Pb
		Bool_t PostCalib = kFALSE,
		////		Bool_t LowpT = kFALSE,
		////		Bool_t MakePid = kFALSE,
		const char* Period  = "16g",
		int Ncl = 70,
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

	AliAnalysisFilter* trackFilterTPC = new AliAnalysisFilter("trackFilterTPC");
	AliESDtrackCuts* esdTrackCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	trackFilterTPC->AddCuts(esdTrackCutsTPC);

	// by default, a file is open for writing. here, we get the filename
	TString fileName = AliAnalysisManager::GetCommonFileName();
	//fileName += Form(":%.2f-%.2f",minCent,maxCent);      // create a subfolder in the file
	fileName += ":Output";      // create a subfolder in the file


	// now we create an instance of your task
	AliAnalysisTaskPPvsMult* task = new AliAnalysisTaskPPvsMult("taskHighPtDeDxpp");   
	if(!task) return 0x0;

	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	//	task->SelectCollisionCandidates(AliVEvent::kAnyINT);
	task->SetAnalysisType(type);
	task->SetAnalysisMC(AnalysisMC);
	/////	task->SetAddLowPt(LowpT);
	task->SetPeriod(Period);

	//	if(system==1){
	//		task->SetAnalysisPbPb(kTRUE);
	//		task->SetMinCent(0.0);
	//		task->SetMaxCent(90.0);
	//		task->SetCentralityEstimator(centralityEstimator);
	//	}
	//	else
	//		task->SetAnalysisPbPb(kFALSE);

	task->SetNcl(Ncl);
	task->SetDebugLevel(0);
	task->SetEtaCut(0.8);
	//	task->SetVtxCut(10.0);
	//	task->SetTrigger(AliVEvent::kINT7);
	//	task->SetPileUpRej(ispileuprej);
	//Set Filtesr
	task->SetTrackFilterTPC(trackFilterTPC);
	//	task->SetStoreMcIn(AnalysisMC);     // def: kFALSE
	task->SetAnalysisTask(PostCalib);
	///	task->SetAnalysisPID(MakePid);
	task->SetTrackCutsSystVars(TrkCutMode);

	std::string buf("MyOutputContainer");

        if(Ncl==90) buf.append("_Ncl_90");

	// add your task to the manager
	mgr->AddTask(task);
	// your task needs input: here we connect the manager to your task
	mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
	// same for the output
	mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("%s_TrkCutMode_%d",buf.c_str(),TrkCutMode), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
	// in the end, this macro returns a pointer to your task. this will be convenient later on
	// when you will run your analysis in an analysis train on grid
	return task;
}
