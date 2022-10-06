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
      bool                  v0qmean18=false,
      bool                  v0calib=false,
      bool                  v0calib18=false,
      bool                  v0QA=false,
      bool                  fillwNUA=false,
      TString             period="LHC10h",
      int                     debug=0, // debug level controls amount of output statements
      double              Harmonic=2,
      TString             trigger="kMB",
      float                  ptmin=0.2, // minimum pt for Q-vector components
      float                  ptmax=2.0, // maximum pt for Q-vector components
      int                     cbinlo=0, // lower centrality bin for histogram array
      int                     cbinhg=8, // higher centrality bin for histogram array
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
  	AliAnalysisTaskEPCalib *task = new AliAnalysisTaskEPCalib("TaskEPCalib", period);
  	task->SetTPCEstOn(tpceston);
  	task->SetTPCNUAWeight(tpcnua);
  	task->SetFillTPCQMean(tpcqmean);
  	task->SetFillTPCShift(tpcshift);
  	task->SetTPCCalib(tpccalib);
  	task->SetVZEROEstOn(v0eston);
  	task->SetVZEROGainEq(v0gaineq);
  	task->SetFillVZEROQMean(v0qmean);
      task->SetFillVZEROQMean18(v0qmean18);
  	task->SetVZEROCalib(v0calib);
      task->SetVZEROCalib18(v0calib18);
  	task->SetfQAV0(v0QA);
      task->SetfFillWNUA(fillwNUA);
	task->SetDebug(debug);
	task-> SetHarmonic(Harmonic);
	task->SetTrigger(trigger);
	task->SetPtMin(ptmin);
	task->SetPtMax(ptmax);
	task->SetCentBinLow(cbinlo);
	task->SetCentBinHigh(cbinhg);
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

      if(period.EqualTo("LHC18q") || period.EqualTo("LHC18r") ) {
        
         Int_t inSlotCounter=1;
         TGrid::Connect("alien:");
         TObjArray *AllContainers = mgr->GetContainers();

          TFile *inV0Calib;
          if(!AllContainers->FindObject("V0Calib")) {
            if (period.EqualTo("LHC18q")) inV0Calib = TFile::Open("alien:///alice/cern.ch/user/m/mhaque/calib2021/CalibV0GainCorrectionLHC18q_Sept2021NoAvgQ.root");
            if (period.EqualTo("LHC18r")) inV0Calib = TFile::Open("alien:///alice/cern.ch/user/m/mhaque/calib2021/CalibV0GainCorrectionLHC18r_Sept2021NoAvgQ.root");

            AliAnalysisDataContainer *cin_V0Calib = mgr->CreateContainer(Form("inV0Calib_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kInputContainer);             
            TList* V0Calib_list = NULL;
            V0Calib_list = dynamic_cast<TList*>(inV0Calib->Get("fWgtsV0ZDC"));
            if (!V0Calib_list) printf("Read TList wrong!\n");
            cin_V0Calib->SetData(V0Calib_list);      
            mgr->ConnectInput(task,inSlotCounter,cin_V0Calib);
            inSlotCounter++;
          }
          else {
            mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("V0Calib"));
            inSlotCounter++;
            printf("V0Calib already loaded\n");
          }
      };

	return task;
}
	