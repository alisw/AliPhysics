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
  //First check to see if user wants to see help
  if (configuration.BeginsWith("help") || 
      configuration.BeginsWith("--help") || 
      configuration.BeginsWith("-h") || 
      configuration.BeginsWith("options") ) {
    cout <<"\n\n============================================" << endl;
    cout << "Keywords for AddTaskDxHFECorrelationME.C:\n"
	 << "file=                         - Filename to store output in\n"
	 << "name=                         - Name of analysis, will correspond to directory inside the file \n"
	 << "cutname=                      - Filename where information on event pool for event-mixing is stored (if use external file)\n"
	 << "runD0MassReference            - If you also want to include D2H task for D0selection (for comparison purposes)\n"
	 << "mc                            - Run on MC\n"
	 << "usekine                       - To run on kinematical level \n"
	 << "PbPb                          - To run on PbPb \n"
	 << "trigger=D/D0/electron         - Which particle to trigger on \n"
	 << "\nD0 settings: \n"
	 << "fillD0scheme=both/D0/D0bar    - Which fillsheme to use for D0\n"
	 << "\nelectron settings: \n"
	 << "useinvmasscut                 - If you want to use invariant mass cut (default is 100MeV/c)\n" 
	 << "invmasscut=                   - If you want to specify a different invariant mass cut \n"
	 << "extraname=                    - extraname for directory and list if you run several tasks at once\n"
	 << "tpcclusters=                  - How many TPC clusters to use on single track cuts for electrons (default=120)\n"
	 << "itsclusters=                  - How many itsclusters to be used in single track cuts for electrons (default=4) \n"
	 << "itsreq=                       - (kFirst,kAny,kNone) Which ITSpixel requirement you want to impose\n"
	 << "elmcreco=                     - (aftertrackcuts/aftertofpid/afterfullpid) Where you want to stop in track selection to look for electrons for mc \n\n";
    return;
  } 

  TString taskOptions="";
  Bool_t nameset=kFALSE;
  Bool_t triggerset=kFALSE;
  Int_t system=0;
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
	taskOptions+=" "+argument;
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
    else if(argument.BeginsWith("PbPb") ||
	    argument.BeginsWith("system=1")){
      system=1;
    }
    else {
      // simply pass argument
      taskOptions+=" "+argument;
    }
  }
  TString path;
  if(system!=0)
    path="AddTaskDxHFECorrelationPbPb.C";
  else
    path="AddTaskDxHFECorrelation.C";

  if (gSystem->AccessPathName(path)!=0) {
    // first try local macro, than AliRoot default path
    if(system!=0)
      path="$ALICE_ROOT/PWGHF/correlationHF/macros/AddTaskDxHFECorrelationPbPb.C";
    else
      path="$ALICE_ROOT/PWGHF/correlationHF/macros/AddTaskDxHFECorrelation.C";
  }

  gROOT->LoadMacro(path);

  if (!nameset) taskOptions+=" name=DxHFE";
  if (!taskOptions) taskOptions+=" trigger=D0";

  if(system!=0){
    if(!AddTaskDxHFECorrelationPbPb(taskOptions)) {
      printf("Problem setting up the single event correlation task, returning\n");
      return 0;
    }

    taskOptions+=" event-mixing";
    if(!AddTaskDxHFECorrelationPbPb(taskOptions)) {
      printf("Problem setting up the mixed event correlation task, returning\n");
      return 0;
    }
  }
  else{
    if(!AddTaskDxHFECorrelation(taskOptions)) {
      printf("Problem setting up the single event correlation task, returning\n");
      return 0;
    }

    taskOptions+=" event-mixing";
    if(!AddTaskDxHFECorrelation(taskOptions)) {
      printf("Problem setting up the mixed event correlation task, returning\n");
      return 0;
    }
  }
  return 1;
}                

