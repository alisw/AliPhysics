//-*- Mode: C++ -*-
// $Id: AddTaskDxHFECorrelation.C 60786 2013-02-08 18:16:19Z arossi $

#ifndef __CINT__
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisCuts.h"
//#include "AliFlowTrackSimple.h"      // added as hint for hidden library dependency to libPWGflowBase
//#include "AliFlowCandidateTrack.h"   // added as hint for hidden library dependency to libPWGflowTasks
//#include "AliCFContainer.h"          // added as hint for hidden library dependency to libCORRFW
//#include "AliAODRecoDecayHF2Prong.h" // added as hint for hidden library dependency to libPWGHFvertexingHF
#include "AliHFAssociatedTrackCuts.h"
#include "AliAnalysisTaskDxHFEParticleSelection.h"
#include "AliDxHFEParticleSelection.h"
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

const char* poolInfoName="PoolInfo";
AliAnalysisCuts* createDefaultPoolConfig();

/// @file   AddTaskDxHFEParticleSelection.C
/// @author Matthias.Richter@ift.uib.no, Hege.Erdal@ift.uib.no
/// @date   2013-02-12
/// @brief  Add the ParticleSelection task to the manager
///
int AddTaskDxHFEParticleSelection(AliHFEpid *fPID=NULL, TString extraName="", TString analysisName="PWGHFCJParticleSelection")
{
  AliAnalysisManager *pManager = AliAnalysisManager::GetAnalysisManager();
  if (!pManager) {
    ::Error("AddTaskDxHFEParticleSelection", "No analysis manager to connect to.");
    return;
  }

  Bool_t bUseMC=kFALSE;
  TString ofilename;
  Int_t system=0;
  TString taskOptions;
  Int_t Particle=AliAnalysisTaskDxHFEParticleSelection::kD0;

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
	    analysisName=argument+"PartSel";
	  }
	  if (argument.BeginsWith("mc")) {
	    bUseMC=kTRUE;
	    taskOptions+=" mc";
	  }
	  if (argument.BeginsWith("PbPb")) {
	    system=1;
	    taskOptions+=" system=PbPb";
	  }
	  if (argument.BeginsWith("fillD0scheme=")){
	    argument.ReplaceAll("fillD0scheme=","");
	    taskOptions+=" fillD0scheme="+argument;
	  }
	  if (argument.BeginsWith("particle=")) {
	    taskOptions+=" "+argument;
	    argument.ReplaceAll("particle=","");
	    if (argument.CompareTo("D0")==0){ 
	      Particle=AliAnalysisTaskDxHFEParticleSelection::kD0; 
	    }
	    else if (argument.CompareTo("electron")==0){ 
	      Particle=AliAnalysisTaskDxHFEParticleSelection::kElectron; 
	    }
	    
	  }
	}
	    
      }
      delete tokens;
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
      ::Error("AddTaskDxHFEParticleSelection", Form("failed to add PID task '%s' from macro '%s'",
					      pidTaskName, pidTaskMacro));
      return 0;
    }
  } else {
    // TODO: would like to check if the PID task was set up
    // with consistent parameters, however there are no getters at the moment
    ::Info("AddTaskDxHFEParticleSelection", Form("PID task '%s' already existing", pidTaskName));
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

  // ________________________________________________________________________
  // PID for HFE
  // PID for Only TOF
  AliHFEpid *fPIDOnlyTOF = new AliHFEpid("hfePidTOF");
  if(!fPIDOnlyTOF->GetNumberOfPIDdetectors()) { 
    fPIDOnlyTOF->AddDetector("TOF",0);
  }
  fPIDOnlyTOF->ConfigureTOF(3); // number of sigma TOF
  fPIDOnlyTOF->InitializePID();
  
  // PID object for TPC and TOF combined
  // Check if PID is set from outside (passed as argument)
  if(!fPID){
    ::Info("AddTaskDxHFEParticleSelection",Form("Setting up new combined PID object"));
    fPID = new AliHFEpid("hfePid");
    if(!fPID->GetNumberOfPIDdetectors()) { 
      fPID->AddDetector("TOF",0);
      fPID->AddDetector("TPC",1);
    }
    //Add settings for asymmetric cut on nSigma TPC
    const int paramSize=4;
    Double_t params[paramSize];
    memset(params, 0, sizeof(Double_t)*paramSize);
    params[0]=-1.;
    fPID->ConfigureTPCdefaultCut(NULL, params, 3.);
    fPID->InitializePID();
  }

  //=========================================================
  //Create TList of cut (and pid) objects for D0 or electron
  TList *Cutlist = new TList;
  if(Particle==AliAnalysisTaskDxHFEParticleSelection::kD0){
    Cutlist->SetName("cut objects D0");
    Cutlist->Add(RDHFD0toKpi);
  }
  else if(Particle==AliAnalysisTaskDxHFEParticleSelection::kElectron){
    Cutlist->SetName("cut objects HFE");
    Cutlist->Add(hfecuts);
    Cutlist->Add(fPID);
    Cutlist->Add(fPIDOnlyTOF);
  }

  //=======================Setting up the task=========================================  
  AliAnalysisTaskDxHFEParticleSelection *pTask=new AliAnalysisTaskDxHFEParticleSelection(taskOptions);
  if (!pTask) {
    ::Error("AddTaskDxHFEParticleSelection", "failed to create task.");
    return 0;
  }
  pTask->SetCutList(Cutlist);
  pManager->AddTask(pTask);

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

  TString listName="";

  TString cutname="";
  if(Particle==AliAnalysisTaskDxHFEParticleSelection::kD0){
    listName="D0list"+extraName;
    cutname="cutsD0Selection"+extraName;
  }
  else if(Particle==AliAnalysisTaskDxHFEParticleSelection::kElectron){
    listName="ElList"+extraName;
    cutname="cutsElectronSelection"+extraName;
  }

  ::Info("AddTaskDxHFEParticleSelection", Form("\ninitializing analysis '%s'%s, output file '%s'", analysisName.Data(), bUseMC?" (using MC)":"", ofilename.Data()));


  AliAnalysisDataContainer *pContainer=pManager->CreateContainer(listName, TList::Class(), AliAnalysisManager::kOutputContainer, ofilename.Data());    
  AliAnalysisDataContainer *pContainer2=pManager->CreateContainer(cutname,TList::Class(),AliAnalysisManager::kOutputContainer, ofilename.Data()); //cuts D0/El

  pManager->ConnectInput(pTask,0,pManager->GetCommonInputContainer());
  pManager->ConnectOutput(pTask,1,pContainer);
  pManager->ConnectOutput(pTask,2,pContainer2);


  return 1;
}

