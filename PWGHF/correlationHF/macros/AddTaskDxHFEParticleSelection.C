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
	 << "particle=D0/electron          - Which particle to run analysis on \n"
	 << "name=                         - Name of analysis, will correspond to directory inside the file \n"
	 << "cutname=                      - Filename where information on event pool for event-mixing is stored (if use external file)\n"
	 << "extraname=                    - extraname for directory and list if you run several tasks at once\n"
	 << "runD0MassReference            - If you also want to include D2H task for D0selection (for comparison purposes)\n"
	 << "mc                            - Run on MC\n"
	 << "usekine                       - To run on kinematical level \n"
	 << "triggermask=                  - Which triggers to use\n"
	 << "PbPb/Pb-Pb/system=1           - To use Pb-Pb collision system\n"
	 << "pPb/p-Pb/system=2/system=p-Pb - To use p-Pb collision system\n"
	 << "event-mixing/mixing           - Whether to also run event-mixing (NB! Use AddTaskDxHFECorrelationME.C for eventmixing)\n"
	 << "useTrackEff                   - If you want to use tracking efficiency (need to attach efficiency maps\n"
	 << "TrackEffName=                 - The file where the efficiency map is stored\n"
	 << "useD0Eff                      - If you want to use tracking efficiency (need to attach efficiency maps\n"
	 << "D0EffName=                    - The file where the efficiency map for D0 is stored\n"
	 << "\nD0 settings: \n"
	 << "fillD0scheme=both/D0/D0bar    - Which fillsheme to use for D0\n";
    cout << "\nELECTRON SETTINGS: \n"
	 << "useinvmasscut                 - If you want to use invariant mass cut (default is 100MeV/c)\n" 
	 << "twoselectedinvmasscut         - If you want to use invariant mass selection with stricter cut on partner\n"
	 << "invmasscut=                   - If you want to specify a different invariant mass cut \n"
         << "ElSelection=                  - To only store a certain source of electron candidates (hadron/nonHFE/HFE/Onlyc/Onlyb)\n"  
         << "storelastcutstep              - Store the last cutstep (so store all track with final cut info)\n"  
         << "filterbit=                    - which filterbit to use (default is 0)\n"  
         << "notusefilterbit               - not use filterbits at all\n"  
	 << "maxPtCombinedPID=             - Max Pt to use TOF PID, above only use TPC PID\n"
         << "EMCALPID                      - Use also Emcal for PID, together with TPC\n"  
	 << "tpcclusters=                  - How many TPC clusters to use on single track cuts for electrons (default=120)\n"
	 << "itsclusters=                  - How many itsclusters to be used in single track cuts for electrons (default=4) \n"
	 << "itsreq=                       - (kFirst,kAny,kNone) Which ITSpixel requirement you want to impose\n"
         << "impactparamcut=               - To set impact parameter cut in radial direction (used to select out all tracks)\n"  
         << "etacut                        - To set eta cut (used to select out all tracks) \n"  
	 << "elmcreco=                     - (aftertrackcuts/aftertofpid/afterfullpid) Where you want to stop in track selection to look for electrons for mc \n"
	 << "elreco=                       - Where you want to stop in track selection to look for electrons for (see code for different selections  \n\n";
    return;
  } 

  AliAnalysisManager *pManager = AliAnalysisManager::GetAnalysisManager();
  if (!pManager) {
    ::Error("AddTaskDxHFEParticleSelection", "No analysis manager to connect to.");
    return;
  }
  TString cutFilename="";
  Bool_t bUseMC=kFALSE;
  Bool_t bTuneOnData=kFALSE;
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
  TString cutFilenameD0="";
  TString cutFilenameEl="";

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
	    continue;
	  }
	  if (argument.BeginsWith("name=")) {
	    argument.ReplaceAll("name=", "");
	    analysisName=" "+argument+"PartSel";
	    continue;
	  }
	  if (argument.BeginsWith("cutFilenameD0=")) { 
	    argument.ReplaceAll("cutFilenameD0=", ""); 
	    cutFilenameD0=argument;                    
 	    continue;
 	  } 
	  //electron cutfile not yet enabled
	  /*
	    if (argument.BeginsWith("cutFilenameEl=")) { 
	    argument.ReplaceAll("cutFilenameEl=", ""); 
	    cutFilenameEl=argument;                    
	    continue;
	    } 
	  */  

	  if (argument.BeginsWith("mc")) {
	    bUseMC=kTRUE;
	    bTuneOnData=kTRUE; //Consider separating these, but for now this is default when using mc 
	    taskOptions+=" mc";
	    continue;
	  }
	  if (argument.BeginsWith("EMCALPID")) {
	    bUseEMCAL=kTRUE;
	    taskOptions+=" EMCALPID";
	    continue;
	  }
	  if (argument.BeginsWith("PbPb") || argument.BeginsWith("Pb-Pb")) {
	    system=1;
	    taskOptions+=" system=Pb-Pb";
	    cout << "Use PbPb" << endl;
	    continue;
	  }
	  if (argument.BeginsWith("system=p-Pb") ||
	      argument.BeginsWith("pPb") ||
	      argument.BeginsWith("p-Pb") ||
	      argument.BeginsWith("system=2")) {
	    system=2;
	    taskOptions+=" system=p-Pb";
	    continue;
	  }

	  if(argument.BeginsWith("tpcclusters=")){
	    argument.ReplaceAll("tpcclusters=", "");
	    NrTPCclusters=argument.Atoi();
	    ::Info("AddTaskDxHFEParticleSelection",Form("Setting nr TPC clusters to %d",NrTPCclusters));
	    continue;
	  }
	  if(argument.BeginsWith("elmcreco")){
	    bUseMCReco=kTRUE;
	    taskOptions+=" "+argument;
	    continue;
	  }
	  if (argument.BeginsWith("usekine") ||argument.BeginsWith("kine")) {
	    bUseKine=kTRUE;
	    taskOptions+=" usekine";
	    continue;
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
	    continue;
	  }
	  if(argument.BeginsWith("itsclusters=")){
	    argument.ReplaceAll("itsclusters=", "");
	    NrITSclusters=argument.Atoi();
	    continue;
	  }
	  if(argument.BeginsWith("itsreq=")){
	    argument.ReplaceAll("itsreq=", "");
	    if(argument.CompareTo("kFirst")==0) ITSreq=AliHFEextraCuts::kFirst;
	    else if(argument.CompareTo("kAny")==0) ITSreq=AliHFEextraCuts::kAny;
	    else if(argument.CompareTo("kNone")==0) ITSreq=AliHFEextraCuts::kNone;
	    continue;
	  }
	  if(argument.BeginsWith("triggermask=")){
	    argument.ReplaceAll("triggermask=", "");
	    if(argument.CompareTo("kINT7")==0) triggerMask=AliVEvent::kINT7;
	    else if(argument.CompareTo("kEMC1")==0) triggerMask=AliVEvent::kEMC1;
	    else if(argument.CompareTo("kEMC7")==0) triggerMask=AliVEvent::kEMC7;
	    else if(argument.CompareTo("kEMC8")==0) triggerMask=AliVEvent::kEMC8;
	    else if(argument.CompareTo("allEMCAL")==0) triggerMask=(AliVEvent::kEMC1|AliVEvent::kEMC7|AliVEvent::kEMC8);
	    else if(argument.CompareTo("kEMCEGA")==0) (triggerMask=AliVEvent::kEMCEGA);
	    else if(argument.CompareTo("kAnyInt")==0) (triggerMask=AliVEvent::kAnyINT);
	    else if(argument.CompareTo("kMB")==0) (triggerMask=AliVEvent::kMB);
	    else ::Fatal("AddTaskDxHFECorrelation",Form("trigger argument not known %s",argument));
	    continue;
	  }

	  if(argument.BeginsWith("extraname=")){
	    argument.ReplaceAll("extraname=", "");
	    extraname=argument;
	    continue;
	  }
	  cout << "Adding argument " << argument << endl;
	  taskOptions+=" "+argument;
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
    // isMC, autoMCeds, tuneOnData, recoPass, cachePID, detResponse, useTPCEtaCorrection, useTPCMultiplicityCorrection, recoDataPass //Warning, multcorrection seems to not work for LHC13b2_efix
    pidFunction.Form("AddTaskPIDResponse(%d, %d, %d, %d, %d, %s, %d, %d, -1)", bUseMC, kTRUE, bTuneOnData, 2, kFALSE,"\"\"" ,kTRUE, !bUseMC); //PIL IMPORTANT! "!bUseMC" is only like this at the moment since lhc13b2_efix doesnt support tpcmultcorrection...
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

  if(cutFilenameD0=="")
    {
      ///______________________________________________________________________
      /// Cuts For D0
      AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
      if (system==0) {
	RDHFD0toKpi->SetStandardCutsPP2010();
	RDHFD0toKpi->SetTriggerMask(triggerMask);
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
	//New cuts (from cutfile used to make D0 eff map)
	//____________________________________________________
	//Set Centrality
	Float_t minc=0,maxc=100;
	
	// Cuts for D0
	AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
	RDHFD0toKpi->SetName("D0toKpiCuts");
	RDHFD0toKpi->SetTitle("Cuts for D0 analysis");
	
	// PILE UP REJECTION
	RDHFD0toKpi->SetOptPileup(1);  	   //per DATI (spegni per MC)		
	RDHFD0toKpi->ConfigurePileupCuts(5,0.8);  //per DATI (spegni per MC)
	
	//Event cuts
	RDHFD0toKpi->SetMinVtxContr(1);
	RDHFD0toKpi->SetMaxVtxZ(10.);
	
	//Trigger selection
	RDHFD0toKpi->SetTriggerClass("");
	RDHFD0toKpi->SetTriggerMask(triggerMask);
	
	//Quality tracks for daughters
	AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
	esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
	esdTrackCuts->SetRequireTPCRefit(kTRUE);
	esdTrackCuts->SetRequireITSRefit(kTRUE);
	//esdTrackCuts->SetMinNClustersITS(4); // default is 5
	//esdTrackCuts->SetMinNClustersTPC(120);
	esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); // default is kBoth, otherwise kAny
	esdTrackCuts->SetMinDCAToVertexXY(0.);
	esdTrackCuts->SetEtaRange(-0.8,0.8);
	esdTrackCuts->SetPtRange(0.3,1.e10);
	
	RDHFD0toKpi->AddTrackCuts(esdTrackCuts);
	
	//D0 selection topological cuts
	const Int_t nptbins =14;
	const Double_t ptmax = 9999.;
	const Int_t nvars=11;
	Float_t ptbins[nptbins+1];
	ptbins[0]=0.;
	ptbins[1]=0.5;	
	ptbins[2]=1.;
	ptbins[3]=2.;
	ptbins[4]=3.;
	ptbins[5]=4.;
	ptbins[6]=5.;
	ptbins[7]=6.;
	ptbins[8]=7.;
	ptbins[9]=8.;
	ptbins[10]=12.;
	ptbins[11]=16.;
	ptbins[12]=20.;
	ptbins[13]=24.;
	ptbins[14]=ptmax;
	
	RDHFD0toKpi->SetGlobalIndex(nvars,nptbins);
	RDHFD0toKpi->SetPtBins(nptbins+1,ptbins);
	
	Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.400,350.*1E-4,0.8,0.5,0.5,1000.*1E-4,1000.*1E-4,-0.000325,0.80,0.,3.2},/* pt<0.5*/
							{0.400,350.*1E-4,0.8,0.5,0.5,1000.*1E-4,1000.*1E-4,-0.000325,0.80,0.,3.2},/* 0.5<pt<1*/
							{0.400,300.*1E-4,0.8,0.4,0.4,1000.*1E-4,1000.*1E-4,-35000.*1E-8,0.90,0.,0.},/* 1<pt<2 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-30000.*1E-8,0.90,0.,0.},/* 2<pt<3 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-30000.*1E-8,0.90,0.,0.},/* 3<pt<4 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-15000.*1E-8,0.90,0.,0.},/* 4<pt<5 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-10000.*1E-8,0.90,0.,0.},/* 5<pt<6 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.85,0.,0.},/* 6<pt<7 */
							{0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.85,0.,0.},/* 7<pt<8 */
							{0.400,300.*1E-4,0.9,0.7,0.7,1000.*1E-4,1000.*1E-4,-5000.*1E-8,0.85,0.,0.},/* 8<pt<12 */
							{0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,10000.*1E-8,0.85,0.,0.},/* 12<pt<16 */
							{0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,10000.*1E-8,0.85,0.,0.},/* 16<pt<20 */
							{0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,10000.*1E-8,0.85,0.,0.},/* 20<pt<24 */
							{0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,10000.*1E-8,0.85,0.,0.}};/* pt>24 */
	
	//CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
	Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
	for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
	
	for (Int_t ibin=0;ibin<nptbins;ibin++){
	  for (Int_t ivar = 0; ivar<nvars; ivar++){
	    cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];      
	  }
	}
	
	RDHFD0toKpi->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);
	RDHFD0toKpi->SetUseSpecialCuts(kTRUE);
	RDHFD0toKpi->SetRemoveDaughtersFromPrim(kTRUE);
	
	for(Int_t iv=0;iv<nvars;iv++) delete [] cutsMatrixTransposeStand[iv];
	delete [] cutsMatrixTransposeStand;
	cutsMatrixTransposeStand=NULL;
	
	//D0 pid settings
	Bool_t pidflag=kTRUE;
	RDHFD0toKpi->SetUsePID(pidflag);
	if(pidflag) cout<<"PID is used"<<endl;
	else cout<<"PID is not used"<<endl;
	
	AliAODPidHF* pidObj=new AliAODPidHF();
	Int_t mode=1;
	const Int_t nlims=2;
	Double_t plims[nlims]={0.6,0.8}; //TPC limits in momentum [GeV/c]
	Bool_t compat=kTRUE; //effective only for this mode
	Bool_t asym=kTRUE;
	Double_t sigmas[5]={2.,1.,0.,3.,0.}; //to be checked and to be modified with new implementation of setters by Rossella
	pidObj->SetAsym(asym);// if you want to use the asymmetric bands in TPC
	pidObj->SetMatch(mode);
	pidObj->SetPLimit(plims,nlims);
	pidObj->SetSigma(sigmas);
	pidObj->SetCompat(compat);
	pidObj->SetPCompatTOF(2.);
	pidObj->SetSigmaForTPCCompat(3.);
	pidObj->SetSigmaForTOFCompat(3.);
	pidObj->SetTPC(kTRUE);
	pidObj->SetTOF(kTRUE);
	pidObj->SetOldPid(kFALSE);
	RDHFD0toKpi->SetPidHF(pidObj);
	RDHFD0toKpi->SetUsePID(kTRUE);
	RDHFD0toKpi->SetUseDefaultPID(kFALSE); //to use the AliAODPidHF
	RDHFD0toKpi->SetLowPt(kFALSE);
	RDHFD0toKpi->SetMaximumPforPID(999.);
	
	//activate pileup rejection (for pp)
	//  RDHFD0toKpi->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);
	
	TString cent="";
	//[FIXME] needed for pPb?
	//centrality selection (Pb-Pb)
	RDHFD0toKpi->SetMinCentrality(minc);
	RDHFD0toKpi->SetMaxCentrality(maxc);
	cent=Form("%.0f%.0f",minc,maxc);
	RDHFD0toKpi->SetUseCentrality(AliRDHFCuts::kCentV0A); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
      } 
      else {
	//Warning, no system set.
      }
    
    }else{ //If there is a cutfile for D0
    cout<<"Getting D0 cut object from: "<<cutFilenameD0<<endl;
    TFile *filecuts;
    //TString finname="Cutlist.root";
    filecuts=TFile::Open(cutFilenameD0.Data());
    TString fRDHFcutsObj="D0toKpiCuts";
    AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
    RDHFD0toKpi = (AliRDHFCutsD0toKpi*)filecuts->Get(fRDHFcutsObj.Data());
    RDHFD0toKpi->PrintAll();
  } 
  
  if(cutFilenameEl=="")
    { 
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
      params[0]=-3.; //[FIXME] Should be -1 for 2010pp, loosened up for pPb

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
      }
    }
  else //if there is a cutfile for electrons
    {
      TFile *filecuts;
      TString finname="Cutlist.root";
      filecuts=TFile::Open(finname.Data());
      TString fRDHFcutsObj="D0toKpiCutsStandard";
      AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
      RDHFD0toKpi = (AliRDHFCutsD0toKpi*)filecuts->Get(fRDHFcutsObj.Data());
      RDHFD0toKpi->PrintAll();
      
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

      //=========================================================
      //Create TList of cut (and pid) objects for D0 or electron
      TList *Cutlist = new TList;
      Cutlist->SetName("cut objects HFE");
      Cutlist->Add(hfecuts);

      if(!bUseEMCAL)	
	{
	  Cutlist->Add(fPID);
	  Cutlist->Add(fPIDOnlyTOF);
	}
      Cutlist->Add(fPIDOnlyTPC);
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

