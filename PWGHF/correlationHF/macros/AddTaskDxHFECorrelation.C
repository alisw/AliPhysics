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
#include "AliHFAssociatedTrackCuts.h"
using namespace std;
#endif

const char* poolInfoName="PoolInfo";
AliAnalysisCuts* createDefaultPoolConfig();

/// @file   AddTaskDxHFECorrelation.C
/// @author Matthias.Richter@ift.uib.no
/// @date   2012-05-09
/// @brief  Add the D0-HFE correlation task to the manager
///
void AddTaskDxHFECorrelation(Bool_t bUseMC=kFALSE, TString analysisName="PWGHFcorrelationDxHF")
{
  AliAnalysisManager *pManager = AliAnalysisManager::GetAnalysisManager();
  if (!pManager) {
    ::Error("AddTaskDxHFECorrelation", "No analysis manager to connect to.");
    return;
  }

  TString ofilename;
  Int_t system=0;
  Bool_t bEventMixing=kFALSE;
  TString poolConfigFile="";
  TString taskOptions;

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
	  if (argument.BeginsWith("cutname=")) {
	    argument.ReplaceAll("cutname=", "");
	    poolConfigFile=argument;
	  }
	  if (argument.BeginsWith("mc")) {
	    bUseMC=kTRUE;
	    taskOptions+=" mc";
	  }
	  if (argument.BeginsWith("event-mixing") ||
	      argument.BeginsWith("mixing")/*deprecated, to be removed later*/) {
	    bEventMixing=kTRUE;
	    taskOptions+=" event-mixing";
	  }
	  if (argument.BeginsWith("PbPb")) {
	    system=1;
	    taskOptions+=" PbPb";
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
      ::Error("AddTaskDxHFECorrelation", Form("failed to add PID task '%s' from macro '%s'",
					   pidTaskName, pidTaskMacro));
      return;
    }
  } else {
    // TODO: would like to check if the PID task was set up
    // with consistent parameters, however there are no getters at the moment
    ::Info("AddTaskDxHFECorrelation", Form("PID task '%s' already existing", pidTaskName));
  }

  if (ofilename.IsNull()) ofilename=AliAnalysisManager::GetCommonFileName();
  ofilename+=":"+analysisName;

  ///______________________________________________________________________
  /// Cuts For D0

  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  RDHFD0toKpi->SetStandardCutsPP2010();

  ///______________________________________________________________________
  /// Cuts for HFE
  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsTPCTOF","HFE Standard Cuts");
  hfecuts->CreateStandardCuts();

  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound,AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetMinNClustersTPC(120);	//Default = 80
  hfecuts->SetMinNClustersTPCPID(80);	//Default = 80
  hfecuts->SetMinRatioTPCclusters(0.6); 	//Default = 0.6
	
  ///ITS
  //hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny); 	//Cut on SPD
  //hfecuts->SetCutITSdrift(AliHFEextraCuts::kAny); 	//Cut on SDD
  //hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetMinNClustersITS(4); 					//Default = 4
	
  ///TOF
  hfecuts->SetTOFPIDStep(kTRUE);
		
  ///Additional Cuts
  hfecuts->SetPtRange(0.30, 10.5);
  hfecuts->SetMaxImpactParam(1.,2.);
  hfecuts->SetVertexRange(10.);

  ///______________________________________________________________________
  /// Info for Pool
  // TODO: Don't think we need the MC part of AliHFCorrelator, needs to be checked
  AliAnalysisCuts* poolConfiguration=NULL;
  if (poolConfigFile.IsNull()) {
    // load the default configuration from below if no file is specified
    poolConfiguration=createDefaultPoolConfig();
  } else {
    // load configuration from file, and abort if something goes wrong
    TFile* filePoolConfiguration=TFile::Open(poolConfigFile.Data());
    if(!filePoolConfiguration){
      ::Error("AddTaskDxHFECorrelation", Form("Pool configuration object file %s not found, exiting", poolConfigFile.Data()));
      return;
    }
    TObject* pObj=filePoolConfiguration->Get(poolInfoName);
    if (!pObj) {
      ::Error("AddTaskDxHFECorrelation", Form("No Pool configuration object with name '%s' found in file %s, exiting", poolInfoName, poolConfigFile.Data()));
      return;
    }
    poolConfiguration = dynamic_cast<AliHFAssociatedTrackCuts*>(pObj);
    if (!poolConfiguration) {
      ::Error("AddTaskDxHFECorrelation", Form("Pool configuration object '%s' has inconsistent class type %s, exiting", poolInfoName, pObj->ClassName()));
      return;
    }
  }

  if(!poolConfiguration){
    ::Error("AddTaskDxHFECorrelation", Form("Pool configuration not found"));
    return;
  } 
  poolConfiguration->Print();

  const char* taskName=AliAnalysisTaskDxHFECorrelation::Class()->GetName();
  if (pManager->GetTask(taskName)) {
    ::Warning("AddTaskDxHFECorrelation", Form("task '%s' already existing, skipping ...",
					   taskName));
    return;
  }
  
  AliAnalysisTaskDxHFECorrelation *pTask=new AliAnalysisTaskDxHFECorrelation();
  if (!pTask) {
    ::Error("AddTaskDxHFECorrelation", "failed to create task.");
    return;
  }
  pTask->SetFillOnlyD0D0bar(0); //0=both, 1=D0 only, 2=D0bar only
  pTask->SetCutsD0(RDHFD0toKpi);
  pTask->SetCutsHFE(hfecuts);
  pTask->SetCuts(poolConfiguration);
  pTask->SetUseMC(bUseMC);
  pTask->SetUseEventMixing(bEventMixing);
  pTask->SetSystem(system);
  // TODO: the switches above can be consolidated and collected in the options
  pTask->SetOption(taskOptions);

  pManager->AddTask(pTask);

  TString listName="DxHFElist";
  TString cutnameD0="cutsD0Corr";
  TString cutnameEl="cutsElCorr";
  TString cutnamePool="PoolInfo";

  // The AnalysisManager handles the output file name in the following way:
  // The output file names are set by the function SetOutputFiles
  // If the file name given to the container begins with one of the initialized
  // file names, the data is stored in the corresponding file in a folder with
  // the full name specified to the container
  // E.g. output file has been set to "myanalysis", the container is created with
  // file name "myanalysis_A", data ends up in file "myanalysis" in folder
  // "myanalysis_A"
  // IMPORTANT: choosing a file name with a different stem at this point will
  // probably lead to an empty file.
  if(bEventMixing){ 
    ofilename+="ME";
    listName+="ME";
    cutnameD0+="ME";
    cutnameEl+="ME";
    cutnamePool+="ME";
  }

  if(bEventMixing) ::Info("AddTaskDxHFECorrelation", Form("\ninitializing analysis '%s'%s, output file '%s', Event Mixing Analysis\n", analysisName.Data(), bUseMC?" (using MC)":"", ofilename.Data()));
  if(!bEventMixing)  ::Info("AddTaskDxHFECorrelation", Form("\ninitializing analysis '%s'%s, output file '%s', Single Event Analysis\n", analysisName.Data(), bUseMC?" (using MC)":"", ofilename.Data()));


  AliAnalysisDataContainer *pContainer=pManager->CreateContainer(listName, TList::Class(), AliAnalysisManager::kOutputContainer, ofilename.Data());    
  AliAnalysisDataContainer *pContainer2=pManager->CreateContainer(cutnameD0,AliRDHFCutsD0toKpi::Class(),AliAnalysisManager::kOutputContainer, ofilename.Data()); //cuts D0
  AliAnalysisDataContainer *pContainer3=pManager->CreateContainer(cutnameEl,AliHFEcuts::Class(),AliAnalysisManager::kOutputContainer, ofilename.Data()); //cuts El
  AliAnalysisDataContainer *pContainer4=pManager->CreateContainer(cutnamePool,AliHFAssociatedTrackCuts::Class(),AliAnalysisManager::kOutputContainer, ofilename.Data()); //cuts El

  pManager->ConnectInput(pTask,0,pManager->GetCommonInputContainer());
  pManager->ConnectOutput(pTask,1,pContainer);
  pManager->ConnectOutput(pTask,2,pContainer2);
  pManager->ConnectOutput(pTask,3,pContainer3);
  pManager->ConnectOutput(pTask,4,pContainer4);

  return;
}

// Note: AliHFAssociatedTrackCuts keeps an instance of the external
// pointer, the arrays thus need to be global
// TODO: try a proper implementation of AliHFAssociatedTrackCuts later
const Int_t    nofMBins=5;
const Double_t MBins[nofMBins+1]={0,20,40,60,80,500};
const Double_t * MultiplicityBins = MBins;
const Int_t    nofZBins=5;
const Double_t ZBins[nofZBins+1]={-10,-5,-2.5,2.5,5,10};
const Double_t *ZVrtxBins = ZBins;

AliAnalysisCuts* createDefaultPoolConfig()
{
  AliHFAssociatedTrackCuts* HFCorrelationCuts=new AliHFAssociatedTrackCuts();
  HFCorrelationCuts->SetName("PoolInfo");
  HFCorrelationCuts->SetTitle("Info on Pool for EventMixing");

  // NEED to check this
  HFCorrelationCuts->SetMaxNEventsInPool(200);
  HFCorrelationCuts->SetMinNTracksInPool(100);
  HFCorrelationCuts->SetMinEventsToMix(8);
  HFCorrelationCuts->SetNofPoolBins(nofZBins,nofMBins); // Note: the arrays have dimension x+1
  HFCorrelationCuts->SetPoolBins(ZVrtxBins,MultiplicityBins);

  TString description = "Info on Pool for EventMixing";   
  HFCorrelationCuts->AddDescription(description);

  return HFCorrelationCuts;
}
