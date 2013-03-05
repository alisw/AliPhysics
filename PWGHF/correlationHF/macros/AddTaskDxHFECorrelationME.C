//-*- Mode: C++ -*-
// $Id$

#ifndef __CINT__
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisCuts.h"
//#include "AliFlowTrackSimple.h"      // added as hint for hidden library dependency to libPWGflowBase
//#include "AliFlowCandidateTrack.h"   // added as hint for hidden library dependency to libPWGflowTasks
//#include "AliCFContainer.h"          // added as hint for hidden library dependency to libCORRFW
//#include "AliAODRecoDecayHF2Prong.h" // added as hint for hidden library dependency to libPWGHFvertexingHF
#include "AliAnalysisTaskDxHFECorrelation.h"
#include "AliDxHFECorrelation.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliHFEcuts.h"
#include "AliLog.h"
#include "TObject.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFAssociatedTrackCuts.h"
using namespace std;
#endif


/// @file   AddTaskDxHFECorrelation.C
/// @author Matthias.Richter@ift.uib.no, Hege.Erdal@ift.uib.no
/// @date   2013-02-12
/// @brief  Add two instances of D0-HFE correlation task to the manager, SE and ME
///
int AddTaskDxHFECorrelationME(TString configuration="")
{
  TString path("AddTaskDxHFECorrelation.C");
  if (gSystem->AccessPathName(path)!=0) {
    // first try local macro, than AliRoot default path
    path="$ALICE_ROOT/PWGHF/vertexingHF/macros/AddTaskDxHFECorrelation.C";
  }
  gROOT->LoadMacro(path);

  TString taskOptions="";
  Bool_t nameset=kFALSE;
  Bool_t triggerset=kFALSE;
  // look for configuration arguments if nothing specified
  // in the function call
  if (configuration.IsNull() && gDirectory) {
    const char* confObjectName="run_single_task_configuration";
    TObject* confObject=gDirectory->FindObject(confObjectName);
    if (confObject) {
      configuration=confObject->GetTitle();
    }
  }

  // argument scan
  // TODO: currently the name and trigger type are fixed, this requirement
  // can be loosened later
  TString delimiter(" ");
  TStringToken configurationToken(configuration, delimiter);
  while (configurationToken.NextToken()) {
    TString argument=configurationToken;
    if (argument.BeginsWith("name=")) {
      //To remove various instances of name, and force directory to be DxHFE
      if(!nameset){
	taskOptions+=" name=DxHFE";
	nameset=kTRUE;
      }
    }
    //At the moment: Make sure trigger is D0
    else if (argument.BeginsWith("trigger=")) {
      if(!triggerset){
	taskOptions+=" trigger=D0";
	triggerset=kTRUE;
      }
    }
    else {
      // simply pass argument
      taskOptions+=" "+argument;
    }
  }
  if (!nameset) taskOptions+=" name=DxHFE";
  if (!taskOptions) taskOptions+=" trigger=D0";


  if(!AddTaskDxHFECorrelation(taskOptions)) {
    printf("Problem setting up the single event correlation task, returning\n");
    return 0;
  }

  taskOptions+=" event-mixing";
  if(!AddTaskDxHFECorrelation(taskOptions)) {
    printf("Problem setting up the mixed event correlation task, returning\n");
    return 0;
  }

  return 1;
}                

