#ifndef __CINT__
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisCuts.h"
//#include "AliFlowTrackSimple.h"      // added as hint for hidden library dependency to libPWGflowBase
//#include "AliFlowCandidateTrack.h"   // added as hint for hidden library dependency to libPWGflowTasks
//#include "AliCFContainer.h"          // added as hint for hidden library dependency to libCORRFW
//#include "AliAODRecoDecayHF2Prong.h" // added as hint for hidden library dependency to libPWGHFvertexingHF
#include "AliHFAssociatedTrackCuts.h"
#include "AliAnalysisTaskDxHFEParticleSelection.h"
#include "AliAnalysisTaskDxHFECorrelation.h"
#include "AliDxHFEParticleSelection.h"
#include "AliDxHFEParticleSelectionEl.h"
#include "AliDxHFEParticleSelectionD0.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliHFEcuts.h"
#include "AliLog.h"
#include "TObject.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "AliRDHFCutsD0toKpi.h"
using namespace std;
#endif

/// @file   AddTaskDxHFE.C
/// @author Matthias.Richter@ift.uib.no, Hege.Erdal@ift.uib.no
/// @date   2013-02-12
/// @brief  Add the ParticleSelection task to the manager
///
int AddTaskDxHFE()
{

  TString ofilename;
  Int_t system=0;
  TString taskOptions;
  Int_t Particle=AliAnalysisTaskDxHFEParticleSelection::kD0;


  /////////////////////////////////////////////////////////////////////////////    ofilename=
  
  gROOT->LoadMacro("AddTaskDxHFECorrelation.C");
  gROOT->LoadMacro("AddTaskDxHFECorrelationME.C");
  gROOT->LoadMacro("AddTaskDxHFEParticleSelection.C");

  TString filename;
  TString configuration;
  if(gDirectory){
    const char* confObjectName="run_single_task_configuration";
    TObject* confObject=gDirectory->FindObject(confObjectName);
    if (confObject) {
      configuration=confObject->GetTitle();
      cout << endl <<"=======================================" << endl;
      cout << "configuration: " << configuration.Data() << endl;
      TObjArray* tokens=configuration.Tokenize(" ");
      if (tokens) {
	TIter next(tokens);
	TObject* token;
	while ((token=next())) {
	  TString argument=token->GetName();
	  if (argument.BeginsWith("file=")) {
	    filename=argument;
	  }
	  
	}	    
      }
      delete tokens;
    }
  }

  // //================================================================
  // //  D0e (ME) with cuts on 120clusters on electron, trigger=D0, fill both, plus add D2H invmasstask

  // taskOptions=filename+" name=DxHFE trigger=D0 fillD0scheme=both runD0MassReference";

  // if(gDirectory) gDirectory->Clear();
  // if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  // if(!AddTaskDxHFECorrelationME(taskOptions.Data()))
  //   return 0;



  //==========================================//
  //                                          //
  //                                          //
  //           Without Inv. Mass              //
  //                                          //
  //                                          //
  //                                          //
  //==========================================//



  //================================================================
  //  D0e (ME) with cuts 
 
    taskOptions=filename+" name=DxHFE extraname=ITS3kAnyCorr itsclusters=3 itsreq=kAny particle=electron p-Pb";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFECorrelation(taskOptions.Data()))
    return 0;

  taskOptions=filename+" name=DxHFE extraname=ITS4kFIRSTCorr itsclusters=4 itsreq=kFIRST particle=electron p-Pb";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFECorrelation(taskOptions.Data()))
    return 0;
  
  taskOptions=filename+" name=DxHFE extraname=ITS3kAny itsclusters=3 itsreq=kAny particle=electron storelastcutstep p-Pb";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
    return 0;

  taskOptions=filename+" name=DxHFE extraname=ITS4kAny itsclusters=4 itsreq=kAny particle=electron storelastcutstep p-Pb";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
    return 0;

  taskOptions=filename+" name=DxHFE extraname=ITS3kFIRST itsclusters=3 itsreq=kFIRST particle=electron storelastcutstep p-Pb";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
    return 0;

  taskOptions=filename+" name=DxHFE extraname=ITS4kFIRST itsclusters=4 itsreq=kFIRST particle=electron storelastcutstep p-Pb";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
    return 0;


  //==========================================//
  //                                          //
  //                                          //
  //             With Inv. Mass               //
  //                                          //
  //                                          //
  //                                          //
  //==========================================//
  /*
  taskOptions=filename+" name=DxHFE trigger=D0 fillD0scheme=both extraname=DefaultIM particle=electron useinvmass";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))
    return 0;


  taskOptions=filename+" name=DxHFE trigger=D0 fillD0scheme=both extraname=ITS4kAnyIM itsclusters=4 itsreq=kAny particle=electron useinvmass";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
    return 0;

  taskOptions=filename+" name=DxHFE trigger=D0 fillD0scheme=both extraname=ITS4kFIRSTIM itsclusters=4 itsreq=kFIRST particle=electron useinvmass";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
    return 0;

  taskOptions=filename+" name=DxHFE trigger=D0 fillD0scheme=both extraname=ITS4kAnyCAitsIM itsclusters=4 itsreq=kAny particle=electron elreco=afterhfeits useinvmass";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
  return 0;*/
    return 1;
}

