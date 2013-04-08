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
#incldue "AliReducedParticle.h"
#include "AliHFCorrelator.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliHFEcuts.h"
#include "AliLog.h"
#include "TObject.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFEextraCuts.h"
using namespace std;
#endif

const char* poolInfoName="PoolInfo";
AliAnalysisCuts* createDefaultPoolConfig();
AliAnalysisCuts* createPbPbPoolConfig();

/// @file   AddTaskDxHFECorrelation.C
/// @author Matthias.Richter@ift.uib.no
/// @date   2012-05-09
/// @brief  Add the D0-HFE correlation task to the manager
///
int AddTaskDxHFECorrelation(TString configuration="", TString analysisName="PWGHFcorrelationDxHF")
{
  //First check to see if user wants to see help
  if (configuration.BeginsWith("help") || 
      configuration.BeginsWith("--help") || 
      configuration.BeginsWith("-h") || 
      configuration.BeginsWith("options") ) {
    cout <<"\n\n============================================" << endl;
    cout << "Keywords for AddTaskDxHFECorrelation.C:\n"
	 << "file=                         - Filename to store output in\n"
	 << "name=                         - Name of analysis, will correspond to directory inside the file \n"
	 << "cutname=                      - Filename where information on event pool for event-mixing is stored (if use external file)\n"
	 << "runD0MassReference            - If you also want to include D2H task for D0selection (for comparison purposes)\n"
	 << "mc                            - Run on MC\n"
	 << "usekine                       - To run on kinematical level \n"
	 << "event-mixing/mixing           - Whether to also run event-mixing (NB! Use AddTaskDxHFECorrelationME.C for eventmixing)\n"
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
  TString ofilename;
  Int_t system=0;
  Bool_t bUseMC=kFALSE;
  Bool_t bEventMixing=kFALSE;
  Bool_t bRunD0MassReference=kFALSE;
  TString poolConfigFile="";
  TString taskOptions;
  Int_t NrTPCclusters=120; // quick fix for problem sending hfe track cut object to addtask
  Int_t NrITSclusters=4;
  Int_t ITSreq=AliHFEextraCuts::kFirst;
  Int_t triggerParticle=AliDxHFECorrelation::kD;
  Bool_t bUseMCReco=kFALSE;
  Bool_t bUseKine=kFALSE;
  TString extraname="";

  cout << endl << "===============================================" << endl;
  cout << "Setting up Correlation task: " << configuration << endl;

  // look for configuration arguments if nothing specified
  // in the function call
  if (configuration.IsNull() && gDirectory) {
    const char* confObjectName="run_single_task_configuration";
    TObject* confObject=gDirectory->FindObject(confObjectName);
    if (confObject) {
      configuration=confObject->GetTitle();
    }
  }
  {// deprecated, but keep for formatting
    {// deprecated, but keep for formatting
      TObjArray* tokens=configuration.Tokenize(" ");
      if (tokens) {
	TIter next(tokens);
	TObject* token;
	while ((token=next())) {
	  TString argument=token->GetName();
	  if (argument.BeginsWith("file=")) {
	    argument.ReplaceAll("file=", "");
	    ofilename=argument;
	  }
	  else if (argument.BeginsWith("name=")) {
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
	  if(argument.BeginsWith("tpcclusters=")){
	    argument.ReplaceAll("tpcclusters=", "");
	    NrTPCclusters=argument.Atoi();
	    ::Info("AddTaskDxHFEParticleSelection",Form("Setting nr TPC clusters to %d",NrTPCclusters));
	  }
	  if (argument.BeginsWith("usekine") ||argument.BeginsWith("kine")) {
	    bUseKine=kTRUE;
	    taskOptions+=" usekine";
	  }
	  if (argument.BeginsWith("event-mixing") ||
	      argument.BeginsWith("mixing")/*deprecated, to be removed later*/) {
	    bEventMixing=kTRUE;
	    taskOptions+=" event-mixing";
	  }
	  if (argument.BeginsWith("PbPb")) {
	    system=1;
	    taskOptions+=" system=Pb-Pb";
	  }
	  if (argument.BeginsWith("fillD0scheme=")){
	    taskOptions+=" "+argument;
	  }
	  if(argument.BeginsWith("elmcreco")){
	    bUseMCReco=kTRUE;
	    taskOptions+=" "+argument;
	  }
	  if (argument.BeginsWith("trigger=")) {
	    taskOptions+=" "+argument;
	    argument.ReplaceAll("trigger=","");
	    if (argument.CompareTo("D0")==0) triggerParticle=AliDxHFECorrelation::kD;
	    else if (argument.CompareTo("D")==0) triggerParticle=AliDxHFECorrelation::kD;
	    else if (argument.CompareTo("electron")==0) triggerParticle=AliDxHFECorrelation::kElectron;
	  }
	  if (argument.CompareTo("runD0MassReference")==0){
	    bRunD0MassReference=kTRUE;
	  }
	  if(argument.BeginsWith("useinvmasscut"))
	    taskOptions+=" "+argument;
	  if(argument.BeginsWith("invmasscut="))
	    taskOptions+=" "+argument;

	  if(argument.BeginsWith("extraname=")){
	    argument.ReplaceAll("extraname=", "");
	    extraname=argument;
	  }
	  if(argument.BeginsWith("itsclusters=")){
	    argument.ReplaceAll("itsclusters=", "");
	    NrITSclusters=argument.Atoi();
	   }
	  if(argument.BeginsWith("itsreq=")){
	    argument.ReplaceAll("itsreq=", "");
	    if(argument.CompareTo("kFirst")==0) ITSreq=AliHFEextraCuts::kFirst;
	    else if(argument.CompareTo("kAny")==0) ITSreq=AliHFEextraCuts::kAny;
	    else if(argument.CompareTo("kNone")==0) ITSreq=AliHFEextraCuts::kNone;
	  }
	  
	}	
      }
      delete tokens;
    }
  }

  AliAnalysisManager *pManager = AliAnalysisManager::GetAnalysisManager();
  if (!pManager) {
    ::Error("AddTaskDxHFECorrelation", "No analysis manager to connect to.");
    return;
  }

  if(bUseMCReco && bUseKine) {
    ::Fatal("AddTaskDxHFECorrelation","CAN'T SET BOTH usekine AND elmcreco AT THE SAME TIME");
    return;
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
      return 0;
    }
  } else {
    // TODO: would like to check if the PID task was set up
    // with consistent parameters, however there are no getters at the moment
    ::Info("AddTaskDxHFECorrelation", Form("PID task '%s' already existing", pidTaskName));
  }

  // optionally add D0Mass task for reference analysis
  if (bRunD0MassReference && !bEventMixing) {

    TString path("AddTaskD0Mass.C");
    if (gSystem->AccessPathName(path)!=0) {
      // first try local macro, than AliRoot default path
      path="$ALICE_ROOT/PWGHF/vertexingHF/macros/AddTaskD0Mass.C";
    }
    if (gSystem->AccessPathName(path)==0) {
      cout << "Setting up D0Mass reference task " << path << endl;
    } else {
      cout << "Can not find D0Mass reference task " << path << endl;
    }
    gROOT->LoadMacro(path);

    const char* taskName=AliAnalysisTaskSED0Mass::Class()->GetName();
    if (pManager->GetTask(taskName)) {
      ::Warning("AddTaskDxHFECorrelation", Form("task '%s' already existing, skipping ...",
						taskName));
    }
    else{
      //flag, readMC,filldistr,cutonDistr, system, flagD0D0bar,minC,maxC,finDirname,finname, finObjname,flagAOD049,FillMassPt, FillImptPar
      AliAnalysisTaskSED0Mass *d0massTask = AddTaskD0Mass(0,bUseMC,kTRUE,kTRUE, 0, 0, 0, 0, "", "","D0toKpiCuts", kFALSE, false, false);
    }
  }

  if(triggerParticle==AliDxHFECorrelation::kElectron)
    analysisName="HFExD";
  if (ofilename.IsNull()) ofilename=AliAnalysisManager::GetCommonFileName();
  ofilename+=":"+analysisName;

  ///______________________________________________________________________
  /// Cuts For D0

  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  // TODO: we might want to move this to separate functions if more data
  // sets are going to be handled
  if (system==0) {
  RDHFD0toKpi->SetStandardCutsPP2010();
  } else {
  // TODO: think about p-Pb
  RDHFD0toKpi->SetStandardCutsPbPb2011();

  // For centrality 0-10%, add centrality flattening
  //NB! NEED FOR THE MOMENT THE FILE!
  TFile *fFlat=TFile::Open("CentrDistrBins005.root","READ");
  TCanvas *c=fFlat->Get("cintegral");
  TH1F *hfl=(TH1F*)c->FindObject("hint");
  RDHFD0toKpi->SetHistoForCentralityFlattening(hfl,0.,10.,0.,0);
  //  RDHFD0toKpi->SetUseCentrality(AliRDHFCuts::kCentV0M);

  RDHFD0toKpi->SetMinCentrality(0.);// 40.*1.01
  RDHFD0toKpi->SetMaxCentrality(10.);// 80.*1.01
  }

  ///______________________________________________________________________
  /// Cuts for HFE
  TString hfeCutsName;
  if (system==0) hfeCutsName="HFE Standard Cuts";
  else hfeCutsName="HFE Cuts PbPb";
  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsTPCTOF", hfeCutsName);
  hfecuts->CreateStandardCuts();

  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound,AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetMinNClustersTPC(NrTPCclusters);	//Default = 80
  hfecuts->SetMinNClustersTPCPID(80);	//Default = 80
  hfecuts->SetMinRatioTPCclusters(0.6); 	//Default = 0.6
	
  ///ITS
  hfecuts->SetCutITSpixel(ITSreq);        	//Cut on SPD
  //hfecuts->SetCutITSdrift(AliHFEextraCuts::kAny); 	//Cut on SDD
  //hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetMinNClustersITS(NrITSclusters); //Default = 4
	
  ///TOF
  hfecuts->SetTOFPIDStep(kTRUE);
		
  ///Additional Cuts
  hfecuts->SetPtRange(0.30, 10.5);
  hfecuts->SetMaxImpactParam(1.,2.);
  hfecuts->SetVertexRange(10.);

  //
  // PID for HFE
  // PID for Only TOF
  AliHFEpid *fPIDOnlyTOF = new AliHFEpid("hfePidTOF");
  if(!fPIDOnlyTOF->GetNumberOfPIDdetectors()) { 
    fPIDOnlyTOF->AddDetector("TOF",0);
  }
  fPIDOnlyTOF->ConfigureTOF(3); // number of sigma TOF
  fPIDOnlyTOF->InitializePID();
  
  // PID object for TPC and TOF combined
  AliHFEpid *fPID = new AliHFEpid("hfePid");
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

  //Create TList of HFE pid and track cuts
  TList *listHFE = new TList;
  listHFE->SetName("cut objects HFE");
  listHFE->Add(hfecuts);
  listHFE->Add(fPID);
  listHFE->Add(fPIDOnlyTOF);


  ///______________________________________________________________________
  /// Info for Pool
  // TODO: Don't think we need the MC part of AliHFCorrelator, needs to be checked
  AliAnalysisCuts* poolConfiguration=NULL;
  if (poolConfigFile.IsNull()) {
    // load the default configuration from below if no file is specified
    if (system==0) poolConfiguration=createDefaultPoolConfig();
    else poolConfiguration=createPbPbPoolConfig();
  } else {
    // load configuration from file, and abort if something goes wrong
    TFile* filePoolConfiguration=TFile::Open(poolConfigFile.Data());
    if(!filePoolConfiguration){
      ::Error("AddTaskDxHFECorrelation", Form("Pool configuration object file %s not found, exiting", poolConfigFile.Data()));
      return 0;
    }
    TObject* pObj=filePoolConfiguration->Get(poolInfoName);
    if (!pObj) {
      ::Error("AddTaskDxHFECorrelation", Form("No Pool configuration object with name '%s' found in file %s, exiting", poolInfoName, poolConfigFile.Data()));
      return 0;
    }
    poolConfiguration = dynamic_cast<AliHFAssociatedTrackCuts*>(pObj);
    if (!poolConfiguration) {
      ::Error("AddTaskDxHFECorrelation", Form("Pool configuration object '%s' has inconsistent class type %s, exiting", poolInfoName, pObj->ClassName()));
      return 0;
    }
  }

  if(!poolConfiguration){
    ::Fatal("AddTaskDxHFECorrelation", Form("Pool configuration not found"));
    return 0;
  } 
  poolConfiguration->Print();

  //Taken out, causes problem when adding more than one task
  /*const char* taskName=AliAnalysisTaskDxHFECorrelation::Class()->GetName();
    if (pManager->GetTask(taskName)) {
    ::Warning("AddTaskDxHFECorrelation", Form("task '%s' already existing, skipping ...",
    taskName));
    return 0;
    }*/
  
  AliAnalysisTaskDxHFECorrelation *pTask=new AliAnalysisTaskDxHFECorrelation(taskOptions);
  if (!pTask) {
    ::Error("AddTaskDxHFECorrelation", "failed to create task.");
    return 0;
  }
  //TODO: Could also consider putting RDHFD0toKpi in a list (ParticleSelectionD0 allows it)
  pTask->SetCutsD0(RDHFD0toKpi);
  pTask->SetCutsHFE(listHFE);
  pTask->SetCuts(poolConfiguration);

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

  TString listName;
  TString cutnameD0="cutsD0Corr";
  TString cutnameEl="cutsElCorr";
  TString cutnamePool="PoolInfo";
  if(triggerParticle==AliDxHFECorrelation::kElectron){
    cutnameD0+="Eltrigg";
    cutnameEl+="Eltrigg";
    cutnamePool+="Eltrigg";
    listName="HFExDlist";
  }
  else listName="DxHFElist";

  listName+=extraname;
  cutnameD0+=extraname;
  cutnameEl+=extraname;
  cutnamePool+=extraname;


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
  AliAnalysisDataContainer *pContainer3=pManager->CreateContainer(cutnameEl,TList::Class(),AliAnalysisManager::kOutputContainer, ofilename.Data()); //cuts El
  AliAnalysisDataContainer *pContainer4=pManager->CreateContainer(cutnamePool,AliHFAssociatedTrackCuts::Class(),AliAnalysisManager::kOutputContainer, ofilename.Data()); // contains event pool info

  pManager->ConnectInput(pTask,0,pManager->GetCommonInputContainer());
  pManager->ConnectOutput(pTask,1,pContainer);
  pManager->ConnectOutput(pTask,2,pContainer2);
  pManager->ConnectOutput(pTask,3,pContainer3);
  pManager->ConnectOutput(pTask,4,pContainer4);

  return 1;
}

// old signature kept for backward compatibility
int AddTaskDxHFECorrelation(Bool_t bUseMC, TString analysisName)
{
  TString arguments(bUseMC?"mc":"");
  if (!analysisName.IsNull()) {
    arguments+=Form(" name=%s", analysisName.Data())
  }
  AddTaskDxHFECorrelation(arguments)
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

AliAnalysisCuts* createPbPbPoolConfig()
{
  AliHFAssociatedTrackCuts* HFCorrelationCuts=new AliHFAssociatedTrackCuts();
  HFCorrelationCuts->SetName("PoolInfo");
  HFCorrelationCuts->SetTitle("Info on Pool for EventMixing");

  // NEED to check this
  HFCorrelationCuts->SetMaxNEventsInPool(250);
  HFCorrelationCuts->SetMinNTracksInPool(80);
  HFCorrelationCuts->SetMinEventsToMix(5);
  HFCorrelationCuts->SetNofPoolBins(nofZBins,nofMBins); // Note: the arrays have dimension x+1
  HFCorrelationCuts->SetPoolBins(ZVrtxBins,MultiplicityBins);

  TString description = "Info on Pool for EventMixing";   
  HFCorrelationCuts->AddDescription(description);

  return HFCorrelationCuts;
}
