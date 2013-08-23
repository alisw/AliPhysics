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
#include "AliVEvent.h"
using namespace std;
#endif

const char* poolInfoName="PoolInfo";
AliAnalysisCuts* createDefaultPoolConfig();
AliAnalysisCuts* createPbPbPoolConfig();

/// @file   AddTaskDxHFEParticleSelection.C
/// @author Matthias.Richter@ift.uib.no, Hege.Erdal@ift.uib.no
/// @date   2013-02-12
/// @brief  Add the ParticleSelection task to the manager
///
int AddTaskDxHFEParticleSelection(TString configuration="",TString analysisName="PWGHFCJParticleSelection")
{

  //First check to see if user wants to see help
  if (configuration.BeginsWith("help") || 
      configuration.BeginsWith("--help") || 
      configuration.BeginsWith("-h") || 
      configuration.BeginsWith("options") ) {
    cout <<"\n\n============================================" << endl;
    cout << "Keywords for AddTaskDxHFEParticleSelection.C:\n"
	 << "file=                         - Filename to store output in\n"
	 << "name=                         - Name of analysis, will correspond to directory inside the file \n"
	 << "cutname=                      - Filename where information on event pool for event-mixing is stored (if use external file)\n"
	 << "runD0MassReference            - If you also want to include D2H task for D0selection (for comparison purposes)\n"
	 << "mc                            - Run on MC\n"
	 << "PbPb                          - Run on PbPbn"
	 << "usekine                       - To run on kinematical level \n"
	 << "particle=D0/electron          - Which particle to run analysis on \n"
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

  AliAnalysisManager *pManager = AliAnalysisManager::GetAnalysisManager();
  if (!pManager) {
    ::Error("AddTaskDxHFEParticleSelection", "No analysis manager to connect to.");
    return;
  }
  TString cutFilename="";
  Bool_t bUseMC=kFALSE;
  TString ofilename;
  Int_t system=0; 
 TString poolConfigFile="";
  TString taskOptions;
  Bool_t bUseKine=kFALSE;
  Bool_t bUseMCReco=kFALSE;
  Bool_t bUseEMCAL=kFALSE;
  Int_t NrTPCclusters=120; // quick fix for problems sending track cut objects in some instances to task
  Int_t NrITSclusters=4; // quick fix for problem sending hfe track cut object to addtask
  Int_t ITSreq=AliHFEextraCuts::kFirst;
  ULong64_t triggerMask=AliVEvent::kAnyINT;
  Int_t Particle=AliAnalysisTaskDxHFEParticleSelection::kD0;
  TString extraname="";

  // look for configuration arguments
  cout << endl << "===============================================" << endl;
  cout << "Setting up Particle Selection task: " << configuration << endl;

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
	  } else if (argument.BeginsWith("name=")) {
	    argument.ReplaceAll("name=", "");
	    analysisName=" "+argument+"PartSel";
	  }
	  if (argument.BeginsWith("cutFilename=")) { //--------------------//
	    argument.ReplaceAll("cutFilename=", ""); //   Move this to     //
	    cutFilename=argument;                    //     cutname?       //
	  }                                          //--------------------//
	  if (argument.BeginsWith("mc")) {
	    bUseMC=kTRUE;
	    taskOptions+=" mc";
	  }
	  if (argument.BeginsWith("EMCALPID")) {
	    bUseEMCAL=kTRUE;
	    taskOptions+=" EMCALPID";
	  }
	  if (argument.BeginsWith("PbPb") || argument.BeginsWith("Pb-Pb")) {
	    system=1;
	    taskOptions+=" system=Pb-Pb";
	    cout << "Use PbPb" << endl;
	  }
	  if (argument.BeginsWith("system=p-Pb") ||
	      argument.BeginsWith("pPb") ||
	      argument.BeginsWith("p-Pb") ||
	      argument.BeginsWith("system=2")) {
	    system=2;
	    taskOptions+=" system=p-Pb";
	  }

	  if(argument.BeginsWith("tpcclusters=")){
	    argument.ReplaceAll("tpcclusters=", "");
	    NrTPCclusters=argument.Atoi();
	    ::Info("AddTaskDxHFEParticleSelection",Form("Setting nr TPC clusters to %d",NrTPCclusters));
	  }
	  if (argument.BeginsWith("fillD0scheme=")){
	    argument.ReplaceAll("fillD0scheme=","");
	    taskOptions+=" fillD0scheme="+argument;
	  }
	  if(argument.BeginsWith("elmcreco")){
	    bUseMCReco=kTRUE;
	    taskOptions+=" "+argument;
	  }
	  if (argument.BeginsWith("usekine") ||argument.BeginsWith("kine")) {
	    bUseKine=kTRUE;
	    taskOptions+=" usekine";
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
	  if(argument.BeginsWith("useinvmasscut"))
	    taskOptions+=" "+argument;
	  if(argument.BeginsWith("twoselectedinvmasscut"))
	    taskOptions+=" "+argument;
	  if(argument.BeginsWith("invmasscut="))
	    taskOptions+=" "+argument;
	  if(argument.BeginsWith("impactparamcut"))
	    taskOptions+=" "+argument;
	  if(argument.BeginsWith("etacut"))
	    taskOptions+=" "+argument;
	  if(argument.BeginsWith("ElSelection="))
	    taskOptions+=" "+argument;
	  if(argument.BeginsWith("storelastcutstep"))
	    taskOptions+=" "+argument;
	  if(argument.BeginsWith("notusefilterbit")){
	    taskOptions+=" "+argument;
	  }
	  if(argument.BeginsWith("filterbit=")){
	    taskOptions+=" "+argument;
	  }
	  if(argument.BeginsWith("maxPtCombinedPID="))
	    taskOptions+=" "+argument;
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
	  if(argument.BeginsWith("triggermask=")){
	    argument.ReplaceAll("triggermask=", "");
	    if(argument.CompareTo("kINT7")==0) triggerMask=AliVEvent::kINT7;
	    else if(argument.CompareTo("kEMC1")==0) triggerMask=AliVEvent::kEMC1;
	    else if(argument.CompareTo("kEMC7")==0) triggerMask=AliVEvent::kEMC7;
	    else if(argument.CompareTo("kEMC8")==0) triggerMask=AliVEvent::kEMC8;
	    else if(argument.CompareTo("allEMCAL")==0) triggerMask=(AliVEvent::kEMC1|AliVEvent::kEMC7|AliVEvent::kEMC8);
	  }

	  if(argument.BeginsWith("extraname=")){
	    argument.ReplaceAll("extraname=", "");
	    extraname=argument;
	  }
	}
	    
      }
      delete tokens;
    }
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

  if(cutFilename=="")
{
  ///______________________________________________________________________
  /// Cuts For D0
  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  if (system==0) {
    RDHFD0toKpi->SetStandardCutsPP2010();
    RDHFD0toKpi->SetTriggerMask(triggerMask); //pPb
  }
  else if (system==1) {
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
  else if (system==2) {
    RDHFD0toKpi->SetStandardCutsPP2010();
    RDHFD0toKpi->SetTriggerMask(AliVEvent::kINT7); //pPb
    RDHFD0toKpi->SetTriggerClass(""); //pPb  
  } 
  else {
    //Warning, no system set.
  }

  ///______________________________________________________________________
  /// Cuts for HFE
  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsTPCTOF","HFE Standard Cuts");
  hfecuts->CreateStandardCuts();

  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound,AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetMinNClustersTPC(NrTPCclusters);	//Default = 80
  hfecuts->SetMinNClustersTPCPID(80);	//Default = 80
  hfecuts->SetMinRatioTPCclusters(0.6); 	//Default = 0.6
	
  ///ITS
  hfecuts->SetCutITSpixel(ITSreq); 	//Cut on SPD
  //hfecuts->SetCutITSdrift(AliHFEextraCuts::kAny); 	//Cut on SDD
  //hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetMinNClustersITS(NrITSclusters);		//Default = 4
    
  ///TOF
  hfecuts->SetTOFPIDStep(kTRUE);
  hfecuts->SetEtaRange(-0.8,0.8);
		
  ///Additional Cuts
  hfecuts->SetPtRange(0.3, 10);
  hfecuts->SetMaxImpactParam(1.,2.);
  hfecuts->SetVertexRange(10.);

  // ________________________________________________________________________
  // PID for HFE
  // PID for Only TOF
  
  //Moved out of TPC+TOF PID, to make it easier available
  const int paramSize=4;
  Double_t params[paramSize];
  memset(params, 0, sizeof(Double_t)*paramSize);
  params[0]=-1.;

  if(!bUseEMCAL)
    {
      AliHFEpid *fPIDOnlyTOF = new AliHFEpid("hfePidTOF");
      if(!fPIDOnlyTOF->GetNumberOfPIDdetectors()) { 
	fPIDOnlyTOF->AddDetector("TOF",0);
      }
      fPIDOnlyTOF->ConfigureTOF(3); // number of sigma TOF
      fPIDOnlyTOF->InitializePID();
      
      // PID object for TPC and TOF combined
      // Check if PID is set from outside (passed as argument)
      ::Info("AddTaskDxHFEParticleSelection",Form("Setting up new combined PID object"));
      AliHFEpid* fPID = new AliHFEpid("hfePid");
      if(!fPID->GetNumberOfPIDdetectors()) { 
	fPID->AddDetector("TOF",0);
	fPID->AddDetector("TPC",1);
      }
      
      //Add settings for asymmetric cut on nSigma TPC
      //      const int paramSize=4;
      //      Double_t params[paramSize];
      //      memset(params, 0, sizeof(Double_t)*paramSize);
      //      params[0]=-1.;
      fPID->ConfigureTPCdefaultCut(NULL, params, 3.);
      fPID->InitializePID();
    }
  // PID for Only TPC
  AliHFEpid *fPIDOnlyTPC = new AliHFEpid("hfePidTPC");
  if(!fPIDOnlyTPC->GetNumberOfPIDdetectors()) { 
    fPIDOnlyTPC->AddDetector("TPC",0);
  }
  fPIDOnlyTPC->ConfigureTPCdefaultCut(NULL, params, 3.);
  fPIDOnlyTPC->InitializePID();
  
  if(bUseEMCAL)
    {
      printf("Using EMCAL\n");
      // PID for EMCAL
      AliHFEpid *fPIDEMCAL = new AliHFEpid("hfePidEMCAL");
      if(!fPIDEMCAL->GetNumberOfPIDdetectors()) {
	fPIDEMCAL->AddDetector("EMCAL", 1);
      }
      fPIDEMCAL->SortDetectors();
      fPIDEMCAL->InitializePID();
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
    if(!bUseEMCAL){
	Cutlist->Add(fPID);
	Cutlist->Add(fPIDOnlyTOF);
      }
    Cutlist->Add(fPIDOnlyTPC);
    if(bUseEMCAL) Cutlist->Add(fPIDEMCAL);
  }
 }
  else //if there is a cutfile
    {
      TFile *filecuts;
      TString finname="Cutlist.root";
      filecuts=TFile::Open(finname.Data());
      TString fRDHFcutsObj="D0toKpiCutsStandard";
      AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
      RDHFD0toKpi = (AliRDHFCutsD0toKpi*)filecuts->Get(fRDHFcutsObj.Data());
      printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n TEsting cutfile\n\n\n\n\n\n\n\n\n\n");
      RDHFD0toKpi->PrintAll();
      printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n TEsting cutfile done\n\n\n\n\n\n\n\n\n\n");

      
      ///______________________________________________________________________
      /// Cuts for HFE
      AliHFEcuts *hfecuts = new AliHFEcuts();//("hfeCutsTPCTOF","HFE Standard Cuts");
      TString fHFEcutsObj="hfeCutsTPCTOF";
      hfecuts=(AliHFEcuts*)filecuts->Get(fHFEcutsObj.Data());

      // ________________________________________________________________________
      // PID for HFE
      // PID for Only TOF
      if(!bUseEMCAL)
	{
	  AliHFEpid *fPIDOnlyTOF = new AliHFEpid("hfePidTOF");
	  TString fHFEpidTOFobj="hfePidTOF";
	  fPIDOnlyTOF=(AliHFEpid*)filecuts->Get(fHFEpidTOFobj.Data());
	  
	  // PID object for TPC and TOF combined
	  // Check if PID is set from outside (passed as argument)
	  
	  AliHFEpid* fPID = new AliHFEpid("hfePid");
	  TString fHFEpidobj="hfePid";
	  fPID=(AliHFEpid*)filecuts->Get(fHFEpidobj.Data());
	}
      AliHFEpid *fPIDOnlyTPC = new AliHFEpid("hfePidTPC");
      TString fHFEpidTPCobj="hfePidTPC";
      fPIDOnlyTPC=(AliHFEpid*)filecuts->Get(fHFEpidTPCobj.Data());
      
      //  if(bUseEMCAL)
      //{
      //[FIX] implement
      //}
	  
  //=========================================================
      //Create TList of cut (and pid) objects for D0 or electron
      TList *Cutlist = new TList;
      /*      if(Particle==AliAnalysisTaskDxHFEParticleSelection::kD0){
	Cutlist->SetName("cut objects D0");
	Cutlist->Add(RDHFD0toKpi);
      }
      else if(Particle==AliAnalysisTaskDxHFEParticleSelection::kElectron){
      */Cutlist->SetName("cut objects HFE");
	Cutlist->Add(hfecuts);

	if(!bUseEMCAL)	
	  {
	    Cutlist->Add(fPID);
	    Cutlist->Add(fPIDOnlyTOF);
	  }
	if(bUseEMCAL){
	  Cutlist->Add(fPIDEMCAL);
	}
	Cutlist->Add(fPIDOnlyTPC);
	// }
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
    listName="D0list"+extraname;
    cutname="cutsD0Selection"+extraname;
  }
  else if(Particle==AliAnalysisTaskDxHFEParticleSelection::kElectron){
    listName="ElList"+extraname;
    cutname="cutsElectronSelection"+extraname;
  }

  ::Info("AddTaskDxHFEParticleSelection", Form("\ninitializing analysis '%s'%s, output file '%s'", analysisName.Data(), bUseMC?" (using MC)":"", ofilename.Data()));


  AliAnalysisDataContainer *pContainer=pManager->CreateContainer(listName, TList::Class(), AliAnalysisManager::kOutputContainer, ofilename.Data());    
  AliAnalysisDataContainer *pContainer2=pManager->CreateContainer(cutname,TList::Class(),AliAnalysisManager::kOutputContainer, ofilename.Data()); //cuts D0/El

  pManager->ConnectInput(pTask,0,pManager->GetCommonInputContainer());
  pManager->ConnectOutput(pTask,1,pContainer);
  pManager->ConnectOutput(pTask,2,pContainer2);


  return 1;
}

