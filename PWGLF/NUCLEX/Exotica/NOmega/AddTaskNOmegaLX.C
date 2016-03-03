#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskNOmegaLX.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliMCEventHandler.h"
#include "AliPID.h"
#endif

AliAnalysisTaskNOmegaLX *AddTaskNOmegaLX(TString finname="",
						 Bool_t readMC=kFALSE,
						 Bool_t additionalChecks=kFALSE,
						 TString suffix = ""){

	//------------------------------------------------------------------------------------------
	// version 2.10 (2016/02/11)
	//------------------------------------------------------------------------------------------

  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( finname.EqualTo("") ) {
    stdcuts=kTRUE;
  } else {
      filecuts=TFile::Open(finname.Data());
      if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
  AliFatal("Input file not found : check your cut object");
      }
  }

  // Creates, configures and attaches to the train the task for pi, K , p spectra
  // with ITS standalone tracks
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  ::Info("AddTaskNOmegaLX","Adding a new task with this settings readMC = %i",readMC);
 
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskNOmegaLX", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskNOmegaLX", "This task requires an input event handler");
    return NULL;
  }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    ::Error("AddTaskNOmegaLX", "This task requires to run on ESD");
    return NULL;
  }

  // Add MC handler (for kinematics)
  if(readMC){
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }


  // Create and configure the task

  TString taskname = "NOmegaLX";
	TString combinedName;
	combinedName.Form("NOmegaLX%s", suffix.Data());
  AliAnalysisTaskNOmegaLX *task = new AliAnalysisTaskNOmegaLX(combinedName);
  task->SetMC(readMC);
  mgr->AddTask(task);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":Dibaryon";
  
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("CutPar_%s",combinedName.Data()), TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("Vertex_%s",combinedName.Data()), TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
//  coutput2->SetSpecialOutput();
  mgr->ConnectOutput(task, 2, coutput2);


  return task;
}
