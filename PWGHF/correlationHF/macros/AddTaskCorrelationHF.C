//-*- Mode: C++ -*-
// $Id$

#ifndef __CINT__
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisCuts.h"
//#include "AliFlowTrackSimple.h"      // added as hint for hidden library dependency to libPWGflowBase
//#include "AliFlowCandidateTrack.h"   // added as hint for hidden library dependency to libPWGflowTasks
//#include "AliCFContainer.h"          // added as hint for hidden library dependency to libCORRFW
//#include "AliAODRecoDecayHF2Prong.h" // added as hint for hidden library dependency to libPWGHFvertexingHF
#include "AliAnalysisTaskDxHFEParticleSelection.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliLog.h"
#include "TObject.h"
#include "TClass.h"
#include "TDirectory.h"
using namespace std;
#endif

/// @file   AddTaskCorrelationHF.C
/// @author Matthias.Richter@ift.uib.no
/// @date   2012-05-09
/// @brief  Add the D0-HFE correlation task to the manager
///
void AddTaskCorrelationHF(TString analysisName="PWGHFcorrelationHF", TString ofile="PWGHFcorrelationHF.root")
{
  AliAnalysisManager *pManager = AliAnalysisManager::GetAnalysisManager();
  if (!pManager) {
    ::Error("AddTaskCorrelationHF", "No analysis manager to connect to.");
    return;
  }

  if (gDirectory) {
    const char* confObjectName="run_single_task_configuration";
    TObject* confObject=gDirectory->FindObject(confObjectName);
    if (confObject) {
      TString configuration=confObject->GetTitle();
      TObjArray* tokens=configuration.Tokenize(" ");
      if (tokens) {
	TIter next(tokens);
	TObject* token;
	while ((token=next())) {
	  TString argument=token->GetName();
	  if (argument.BeginsWith("file=")) {
	    argument.ReplaceAll("file=", "");
	    ofile=argument;
	  } else if (argument.BeginsWith("name=")) {
	    argument.ReplaceAll("name=", "");
	    analysisName=argument;
	  }
	}
	delete tokens;
      }
    }
  }

  ::Info("AddTaskCorrelationHF", Form("initializing analysis '%s', output file '%s'", analysisName.Data(), ofile.Data()));

  AliAnalysisTaskSE *pTask=new AliAnalysisTaskDxHFEParticleSelection;
  if (!pTask) {
    ::Error("AddTaskCorrelationHF", "failed to create task.");
    return;
  }
  pManager->AddTask(pTask);

  AliAnalysisDataContainer *pContainer=pManager->CreateContainer(analysisName, TObject::Class(), AliAnalysisManager::kOutputContainer, ofile);    
  pManager->ConnectInput(pTask,0,pManager->GetCommonInputContainer());
  pManager->ConnectOutput(pTask,1,pContainer);

  return;
}
