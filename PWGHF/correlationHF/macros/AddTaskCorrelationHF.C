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

/// @file   AddTaskCorrelationHF.C
/// @author Matthias.Richter@ift.uib.no
/// @date   2012-05-09
/// @brief  Add the D0-HFE correlation task to the manager
///
void AddTaskCorrelationHF(Bool_t bUseMC=kFALSE, TString analysisName="PWGHFcorrelationDxHF")
{
  AliAnalysisManager *pManager = AliAnalysisManager::GetAnalysisManager();
  if (!pManager) {
    ::Error("AddTaskCorrelationHF", "No analysis manager to connect to.");
    return;
  }

  TString ofilename;

  // look for configuration arguments
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
	    ofilename=argument;
	  } else if (argument.BeginsWith("name=")) {
	    argument.ReplaceAll("name=", "");
	    analysisName=argument;
	  }
	  if (argument.BeginsWith("mc")) {
	    bUseMC=kTRUE;
	  }
	}
	delete tokens;
      }
    }
  }

  // check for existence of PID task and add if not available
  const char* pidTaskName="PIDResponseTask";
  const char* pidTaskMacro="$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C";
  AliAnalysisTask* pidTask=pManager->GetTask(pidTaskName);
  if (!pidTask) {
    gROOT->LoadMacro(pidTaskMacro);
    TString pidFunction;
    pidFunction.Form("AddTaskPIDResponse(%d, %d)", bUseMC, kTRUE);
    gROOT->ProcessLine(pidFunction);
    if (pManager->GetTask(pidTaskName)==NULL) {
      ::Error("AddTaskCorrelationHF", Form("failed to add PID task '%s' from macro '%s'",
					   pidTaskName, pidTaskMacro));
      return;
    }
  } else {
    // TODO: would like to check if the PID task was set up
    // with consistent parameters, however there are no getters at the moment
    ::Info("AddTaskCorrelationHF", Form("PID task '%s' already existing", pidTaskName));
  }

  TString cutname="cutsD0Corr";

  if (ofilename.IsNull()) ofilename=AliAnalysisManager::GetCommonFileName();
  ofilename+=":"+analysisName;

  ::Info("AddTaskCorrelationHF", Form("\ninitializing analysis '%s'%s, output file '%s'\n", analysisName.Data(), bUseMC?" (using MC)":"", ofilename.Data()));
    
  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  RDHFD0toKpi->SetStandardCutsPP2010();

  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsTPCTOF","HFE Standard Cuts");
  hfecuts->CreateStandardCuts();
  ///______________________________________________________________________
  /// Standard Cuts defined by the HFE Group
  /*	
  **********STANDARD CUTS DESCRIPTION****************
  SetRequireProdVertex();
  fProdVtx[0] = 0;
  fProdVtx[1] = 3;
  fProdVtx[2] = 0;
  fProdVtx[3] = 3;
  fMinClustersTPC = 80;
  fMinClustersITS = 4;
  fMinTrackletsTRD = 0;
  SetRequireITSPixel();
  fCutITSPixel = AliHFEextraCuts::kFirst;
  fMaxChi2clusterITS = -1.;
  fMaxChi2clusterTPC = 4.;
  fMinClusterRatioTPC = 0.6;
  fPtRange[0] = 0.1;
  fPtRange[1] = 20.;
  SetRequireKineMCCuts();
  *************************************************** */
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);	
  hfecuts->SetMinNClustersTPC(120);	//Default = 80
  hfecuts->SetMinNClustersTPCPID(80);	//Default = 80
  hfecuts->SetMinRatioTPCclusters(0.6); 	//Default = 0.6
	
  ///ITS
  //hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny); 	//Cut on SPD
  //hfecuts->SetCutITSdrift(AliHFEextraCuts::kAny); 	//Cut on SDD
  //hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetMinNClustersITS(4); 					//Default = 4
	
  ///TOF
  //hfecuts->SetTOFPIDStep(kTRUE);
		
  ///Additional Cuts
  hfecuts->SetPtRange(0.30, 10.5);
  hfecuts->SetMaxImpactParam(1.,2.);
  hfecuts->SetVertexRange(10.);

  const char* taskName=AliAnalysisTaskDxHFECorrelation::Class()->GetName();
  if (pManager->GetTask(taskName)) {
    ::Warning("AddTaskCorrelationHF", Form("task '%s' already existing, skipping ...",
					   taskName));
    return;
  }
  
  AliAnalysisTaskDxHFECorrelation *pTask=new AliAnalysisTaskDxHFECorrelation();
  if (!pTask) {
    ::Error("AddTaskCorrelationHF", "failed to create task.");
    return;
  }
  pTask->SetFillOnlyD0D0bar(1); //0=both, 1=D0 only, 2=D0bar only
  pTask->SetCutsD0(RDHFD0toKpi);
  pTask->SetCutsHFE(hfecuts);
  pTask->SetUseMC(bUseMC);
  pManager->AddTask(pTask);

  TString listName="DxHFElist";

  //Should also add container with HFEcuts. PID?
  AliAnalysisDataContainer *pContainer=pManager->CreateContainer(listName, TList::Class(), AliAnalysisManager::kOutputContainer, ofilename.Data());    
  AliAnalysisDataContainer *pContainer2=pManager->CreateContainer(cutname,AliRDHFCutsD0toKpi::Class(),AliAnalysisManager::kOutputContainer, ofilename.Data()); //cuts

  pManager->ConnectInput(pTask,0,pManager->GetCommonInputContainer());
  pManager->ConnectOutput(pTask,1,pContainer);
  pManager->ConnectOutput(pTask,2,pContainer2);

  return;
}
