#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#endif

AliAnalysisTaskEPCalib* AddTaskEPCalib(
      bool                  tpceston=false,
      bool                  tpcnua=false,
      bool                  tpcqmean=false,
      bool                  tpcshift=false,
      bool                  tpccalib=false,
      bool                  v0eston=false,
      bool                  v0gaineq=false,
      bool                  v0qmean=false,
      bool                  v0calib=false,
      bool                  v0QA=false,
      bool                  fillwNUA=false,
      int                     debug=0, // debug level controls amount of output statements
      double              Harmonic=2,
      TString             trigger="kMB",
      float                  ptmin=0.2, // minimum pt for Q-vector components
      float                  ptmax=2.0, // maximum pt for Q-vector components
      int                     cbinlo=0, // lower centrality bin for histogram array
      int                     cbinhg=8, // higher centrality bin for histogram array
      TString             period="LHC10h", // period
      float                  centcut=7.5, // centrality restriction for V0M and TRK
      TString             uniqueID=""
      )	
{	
	// Creates a pid task and adds it to the analysis manager
	// Get the pointer to the existing analysis manager via the static
	//access method
	//=========================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskEPCalib.C", "No analysis manager to connect to.");
		return NULL;
	}

	// Check the analysis type using the event handlers connected to the
	// analysis manager The availability of MC handler can also be
	// checked here.
	// =========================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTaskEPCalib.C", "This task requires an input event handler");
		return NULL;
	}
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  	// --- instantiate analysis task
  	AliAnalysisTaskEPCalib *task = new AliAnalysisTaskEPCalib("TaskEPCalib");
  	task->SetTPCEstOn(tpceston);
  	task->SetTPCNUAWeight(tpcnua);
  	task->SetFillTPCQMean(tpcqmean);
  	task->SetFillTPCShift(tpcshift);
  	task->SetTPCCalib(tpccalib);
  	task->SetVZEROEstOn(v0eston);
  	task->SetVZEROGainEq(v0gaineq);
  	task->SetFillVZEROQMean(v0qmean);
  	task->SetVZEROCalib(v0calib);
  	task->SetfQAV0(v0QA);
      task->SetfFillWNUA(fillwNUA);
	task->SetDebug(debug);
	task-> SetHarmonic(Harmonic);
	task->SetTrigger(trigger);
	task->SetPtMin(ptmin);
	task->SetPtMax(ptmax);
	task->SetCentBinLow(cbinlo);
	task->SetCentBinHigh(cbinhg);
	task->SetPeriod(period);
	task->SetCentCut(centcut);
	// task->SelectCollisionCandidates(AliVEvent::kINT7);
	mgr->AddTask(task);

	// Create ONLY the output containers for the data produced by the
	// task.  Get and connect other common input/output containers via
	// the manager as below
	//======================================================================
    	AliAnalysisDataContainer* cinput  = mgr->GetCommonInputContainer();
      const char* outputFileName = mgr->GetCommonFileName();
      AliAnalysisDataContainer* coutput = mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), 
                                                           AliAnalysisManager::kOutputContainer,                                                          
                                                           Form("%s:%s", outputFileName, uniqueID.Data()));
   	mgr->ConnectInput (task, 0, cinput);
  	mgr->ConnectOutput(task, 1, coutput);

	return task;
}	

