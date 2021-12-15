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
      bool                  v0eston=true,
      bool                  v0gaineq=false,
      bool                  v0qmean=false,
      bool                  v0calib=true,
      bool                  v0QA=true,
      int                     debug=0, // debug level controls amount of output statements
      double              Harmonic=3,
      int                     trigger=0,
      int                     filterBit=768, // AOD filter bit selection
      int                     nclscut=70, // ncls cut for all tracks 
      float                  chi2hg=4.0,
      float                  chi2lo=0.1,
      float                  dcacutz=3.2, // dcaz cut for all tracks
      float                  dcacutxy=2.4, // dcaxy cut for all tracks
      float                  ptmin=0.2, // minimum pt for Q-vector components
      float                  ptmax=2.0, // maximum pt for Q-vector components
      int                     cbinlo=0, // lower centrality bin for histogram array
      int                     cbinhg=8, // higher centrality bin for histogram array
      TString             period="LHC15o", // period
      TString             multComp="pileupByEDSTPC128", // multiplicity comparison
      float                  centcut=7.5 // centrality restriction for V0M and TRK
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
	task->SetDebug(debug);
	task-> SetHarmonic(Harmonic);
	task->SetTrigger(trigger);
	task->SetFilterBit(filterBit);
	task->SetNclsCut(nclscut);
  	task->SetChi2High(chi2hg);
  	task->SetChi2Low(chi2lo);
	task->SetDCAcutZ(dcacutz);
	task->SetDCAcutXY(dcacutxy);
	task->SetPtMin(ptmin);
	task->SetPtMax(ptmax);
	task->SetCentBinLow(cbinlo);
	task->SetCentBinHigh(cbinhg);
	task->SetPeriod(period);
	task->SetMultComp(multComp);
	task->SetCentCut(centcut);
	// task->SelectCollisionCandidates(AliVEvent::kINT7);
	mgr->AddTask(task);

	// Create ONLY the output containers for the data produced by the
	// task.  Get and connect other common input/output containers via
	// the manager as below
	//======================================================================
    	AliAnalysisDataContainer* cinput  = mgr->GetCommonInputContainer();
  	AliAnalysisDataContainer* coutput = mgr->CreateContainer("output", TList::Class(), 
                                                           AliAnalysisManager::kOutputContainer, 
                                                           mgr->GetCommonFileName());
   	mgr->ConnectInput (task, 0, cinput);
  	mgr->ConnectOutput(task, 1, coutput);
	return task;
}	

