#if ! defined __CINT__ || defined __MAKECINT__
 #include <TStopwatch.h>
#include <TString.h>
#include <TROOT.h>
#include <TChain.h>
#include <TSystem.h>
#include <AliAnalysisManager.h>
#include <AliRecoInputHandler.h>
#include <AliMCEventHandler.h>
#include "AnalysisInput/AliAnalysisTaskSEFT2Simulation.h"
#include <AliAODHandler.h>
#include <AliAnalysisTaskPIDResponse.h>
#include <AliCDBManager.h>
#include <AliAnalysisTaskESDfilter.h>
#endif

//______________________________________________________________________________
void run_mcgenFT2(Int_t nEvents =0, Int_t run_number=0, Bool_t ft2OutputTree = 0, Int_t seed =0)
{
	TStopwatch timer;
	timer.Start();
	//
	TStopwatch timerp1;
	timerp1.Start();
	//
	printf("########################################\n");
	printf("###     Generating MC Simulation     ###\n");
	// FT2 setup
	//==============================================================================
	Int_t debugLevel		= 0;
	TString taskname		= "AliAnalysisTaskSEFT2Simulation";
	//==============================================================================
	//===================== Below you find setters for the FT2 =====================
	//============ Please modify only when you know what you are doing =============
	//==================== Currently tuned for ITSu TDR Design: ====================
	//========================= (0,1,1,1,1,1,0,1,1,1.,1., ==========================
	//========= ft2_parameterization_LHC13d19.root,h-p_CrossSections.root) =========
	//====================== Anchored on: LHC13d91 (->LHC10h) ======================
	Int_t  streamLevel							= 0;
	Bool_t tuneOnDataOrMC						= kTRUE;
	Bool_t simMat										= kTRUE;
	Bool_t usePID										= kTRUE;
	Bool_t useKalman								= kTRUE;
	Bool_t allowDecay								= kTRUE;
	Bool_t allowAbsorption					= kFALSE;
	Bool_t allowConversions					= kTRUE;
	Bool_t standaloneTune						= kTRUE;
	Double_t maxStepTGeo						= 1.;
	Double_t dNdY										= 1.;
	TString tpcParameterizationFile = "ft2_parameterization_LHC13d19.root";
	TString crossSectionFile				= "h-p_CrossSections.root";
	TString outputFileName					= "jstiller_FT2_TestAnalysis.root";
	//==============================================================================
	// ### Settings that make sense when using the Alien plugin
	//==============================================================================
	Int_t       iCollision         = 1;       // 0=pp, 1=Pb-Pb
	Int_t       run_flag           = 1100;		// year (2011 = 1100)
	//==============================================================================
	Bool_t      doCDBconnect        = 1;				//												| 1
	Bool_t      usePhysicsSelection = kFALSE;	// use physics selection		| T
	Bool_t      useTender           = kFALSE;	// use tender wagon					| F
	Bool_t      useKFILTER          = kTRUE;	// use Kinematics filter		| T
	Bool_t      useTR               = kFALSE;	// use track references			| T
	// ### Analysis modules to be included. Some may not be yet fully implemented.
	//==============================================================================
	Int_t       iAODhandler         = 1;      // Analysis produces an AOD or dAOD's													| 1
	Int_t       iESDfilter          = 1;      // ESD to AOD filter (barrel + muon tracks)										| 1
	Int_t       iMUONcopyAOD        = 0;      // Task that copies only muon events in a separate AOD (PWG3)	| 0
	Int_t       iJETAN              = 0;      // Jet analysis (PWG4)																				| 1 (local mode: 0)
	Int_t       iJETANdelta         = 0;      // Jet delta AODs																							| 1 (local mode: 0)
	Int_t       iPWGHFvertexing     = 1;      // Vertexing HF task (PWG3)																		| 1
	Int_t       iPWGDQJPSIfilter    = 0;      // JPSI filtering (PWG3)																			| 0
	Int_t       iPWGHFd2h           = 0;      // D0->2 hadrons (PWG3)																				| 1
	Int_t       iPIDResponse        = 1;      // PID response																								| 1
	Int_t       iPWGLFForward       = 0;      // Forward mult task (PWGLF)																	| 0
	
	TString configPWGHFd2h = (iCollision==0)?"$ALICE_PHYSICS/../src/PWGHF/vertexingHF/ConfigVertexingHF.C"
	:"$ALICE_PHYSICS/../src/PWGHF/vertexingHF/upgrade/ConfigVertexingHF_ITSUpgrade_wPID.C";
	
	LoadAllLibraries();

       	//if (!TGrid::Connect("alien://")) return;
	
	
//	gROOT->ProcessLine(".L fastGen.C");
//	fastGen(nEvents,seed);
	gROOT->ProcessLine(Form(".x sim.C(%i,%i)",nEvents,run_number));
	
	gROOT->CloseFiles();
	
	printf("########################################\n");
	printf("###     Timing for MC Generation     ###\n");
	timerp1.Print();
	printf("########################################\n");
//	return;
	TStopwatch timerp2;
	timerp2.Start();

	// Make the analysis manager and connect event handlers
	AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "Production train");
	
	// Input
	AliRecoInputHandler *inputHandler = new AliRecoInputHandler("handler","handler for D2H");
	mgr->SetInputEventHandler(inputHandler);
	
	AliMCEventHandler* mcHandler = new AliMCEventHandler();
	mgr->SetMCtruthEventHandler(mcHandler);
	mcHandler->SetReadTR(kFALSE);
	mgr->SetDebugLevel(debugLevel);
	
	
	printf("########################################\n");
	printf("###               FT2                ###\n");
	
	gROOT->ProcessLine(".L FT2.cxx++");
	gROOT->LoadMacro("AddTaskFT2Simulation.C");
	gROOT->LoadMacro("AliAnalysisTaskSEFT2Simulation.cxx+");
	
	AliAnalysisTaskSEFT2Simulation *TaskTrackpFTSCovariance = (AliAnalysisTaskSEFT2Simulation*)AddTaskFT2Simulation(ft2OutputTree,tuneOnDataOrMC,simMat,tpcParameterizationFile.Data(),crossSectionFile.Data(),usePID,maxStepTGeo,dNdY,useKalman,allowDecay,allowAbsorption,allowConversions,outputFileName.Data(),run_number,standaloneTune,debugLevel,streamLevel,taskname.Data());
	
	printf("###              Done!               ###\n");
	printf("########################################\n");
	
	TChain *chain = CreateChain();
	if (!chain) return;
	
	mgr->SetSkipTerminate(kTRUE);
	if (mgr->InitAnalysis()) {
		mgr->PrintStatus();
		mgr->StartAnalysis("local", chain);
	}
	printf("########################################\n");
	printf("###    Timing for ESD Production     ###\n");
	timerp2.Print();
	printf("########################################\n");

	gROOT->CloseFiles();
//	return;

	TStopwatch timerp3;
	timerp3.Start();
	//---
	// Set temporary merging directory to current one
	gSystem->Setenv("TMPDIR", gSystem->pwd());
	// Set temporary compilation directory to current one
	gSystem->SetBuildDir(gSystem->pwd(), kTRUE);
	printf("########################################\n");
	printf("###     RUNNING FILTERING TRAIN      ###\n");
	printf("### Configuring analysis train for:  ###\n");
	if (usePhysicsSelection)printf("### Physics selection                ###\n");
	if (useTender)	printf("### TENDER                           ###\n");
	if (iESDfilter)	printf("### ESD filter                       ###\n");
	if (iMUONcopyAOD)	printf("### MUON copy AOD                    ###\n");
	if (iJETAN)	printf("### Jet analysis                     ###\n");
	if (iJETANdelta)	printf("### Jet delta AODs                   ###\n");
	if (iPWGHFvertexing)	printf("### PWGHF vertexing                  ###\n");
	if (iPWGDQJPSIfilter)	printf("### PWGDQ j/psi filter               ###\n");
	if (iPWGHFd2h)	printf("### PWGHF D0->2 hadrons QA           ###\n");
	printf("########################################\n");
	
	AliAnalysisManager *mgr2  = new AliAnalysisManager("Analysis Train", "Production train");
	mgr2->SetDebugLevel(debugLevel);
	AliESDInputHandler *esdHandler = new AliESDInputHandler();
	mgr2->SetInputEventHandler(esdHandler);
	AliMCEventHandler* mcHandler2 = new AliMCEventHandler();
	mgr2->SetMCtruthEventHandler(mcHandler2);
	mcHandler2->SetPreReadMode(1);
	mcHandler2->SetReadTR(useTR);
	
	AliAODHandler* aodHandler   = new AliAODHandler();
	aodHandler->SetOutputFileName("AliAOD.root");
	mgr2->SetOutputEventHandler(aodHandler);
	
	
	if (iPIDResponse) {
		gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
		AliAnalysisTaskPIDResponse *PIDResponse = AddTaskPIDResponse(kTRUE);//,kTRUE,kTRUE);
	}
	// CDB connection
	//
	
	if (doCDBconnect && !useTender) {
		gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
		AliTaskCDBconnect *taskCDB = AddTaskCDBconnect("alien://folder=/alice/data/2010/OCDB",run_number);
		if (!taskCDB) return;
		AliCDBManager *cdb = AliCDBManager::Instance();
		cdb->SetDefaultStorage("alien://folder=/alice/data/2010/OCDB");
		cdb->SetSpecificStorage("ITS/Align/Data","alien://folder=/alice/simulation/LS1_upgrade/Ideal");
		//	taskCDB->Setrun_number(run_number);
	}
	
	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C");
	AliAnalysisTaskESDfilter *taskesdfilter =
	AddTaskESDFilter(useKFILTER,			//									| 1
									 iMUONcopyAOD,			// write Muon AOD					| F
									 kFALSE,						// write dimuon AOD				| F
									 kFALSE,						// usePhysicsSelection		| F
									 kFALSE,						// centrality OBSOLETE		| F
									 kTRUE,							// enable TPC only tracks	| T
									 kFALSE,						// disable cascades				| F
									 kFALSE,						// disable kinks					| F
									 run_flag);					// run flag (YY00)				| 1100
	
	if (iPWGHFvertexing) {
		gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskVertexingHF.C");
		if (!iPWGHFd2h) TFile::Cp(gSystem->ExpandPathName(configPWGHFd2h.Data()), "file:ConfigVertexingHF.C");
		AliAnalysisTaskSEVertexingHF *taskvertexingHF = AddTaskVertexingHF();
		if (!taskvertexingHF) ::Warning("AnalysisTrainNew", "AliAnalysisTaskSEVertexingHF cannot run for this train conditions - EXCLUDED");
		else mgr2->RegisterExtraFile("AliAOD.VertexingHF.root");
		taskvertexingHF->SelectCollisionCandidates(0);
		AliLog::SetClassDebugLevel("AliAnalysisVertexingHF",debugLevel);
	}
	
	TChain *chainEsd = CreateChainESD();
	if (!chainEsd) return;
	
	mgr2->SetSkipTerminate(kTRUE);
	if (mgr2->InitAnalysis()) {
		mgr2->PrintStatus();
		mgr2->StartAnalysis("local", chainEsd);
	}
	printf("########################################\n");
	printf("###     Timing for Vertexing HF      ###\n");
	timerp3.Print();
	printf("########################################\n");

	
	gROOT->CloseFiles();
	printf("########################################\n");
	printf("###      Timing for Full Chain       ###\n");
	timer.Print();
	printf("########################################\n");

}
//______________________________________________________________________________
TChain *CreateChain(){
	// Create the input chain
	TChain *chain = new TChain("TE");
	if (gSystem->AccessPathName("./galice.root"))
		::Error("AnalysisTrainNew.C::CreateChain", "File: galice.root not in dir");
	else
		chain->Add("galice.root");
	if (chain->GetNtrees()) return chain;
	return NULL;
}
//______________________________________________________________________________
TChain *CreateChainESD()
{
	// Create the input chain
	TChain *chainEsd = new TChain("esdTree");
	if (gSystem->AccessPathName("./jstiller_AliESDs.root"))
		::Error("AnalysisTrainNew.C::CreateChain", "File: jstiller_AliESDs.root not in ./data dir");
	else
		chainEsd->Add("jstiller_AliESDs.root");
	if (chainEsd->GetNtrees()) return chainEsd;
	return NULL;
}
//______________________________________________________________________________
void LoadAllLibraries(){
	
	printf("########################################\n");
	printf("###   Loading all paths requested!   ###\n");
	//
	gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks -I$ALICE_ROOT/PWGPP -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS/Upgrade -I$ALICE_ROOT/ITS/UPGRADE -I$ALICE_ROOT/STEER/CDB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/PWGPP -g");
	//
	printf("###               Done!              ###\n");
	printf("### Loading all libraries requested! ###\n");
	//
	TString loadMacroPath="$ALICE_PHYSICS/PWGHF/vertexingHF/macros/";
	//
	TString loadLibraries="LoadLibraries.C"; loadLibraries.Prepend(loadMacroPath.Data());
	gROOT->LoadMacro(loadLibraries.Data());
	gSystem->Load("libGui.so");
	gSystem->Load("libCore.so");
	gSystem->Load("libTree.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libMinuit.so");
	gSystem->Load("libProof.so");
	gSystem->Load("libmicrocern.so");
	gSystem->Load("liblhapdf.so");
	gSystem->Load("libpythia6_4_21");   // Pythia 6.4
	gSystem->Load("libpythia6.so");
	gSystem->Load("libEG.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libVMC.so");
	gSystem->Load("libEGPythia6.so");
	gSystem->Load("libSTEERBase.so");
	gSystem->Load("libESD.so");
	gSystem->Load("libRAWDatabase.so");
	gSystem->Load("libRAWDatarec.so");
	gSystem->Load("libAOD.so");
	gSystem->Load("libANALYSIS.so");
	gSystem->Load("libANALYSISalice.so");
	gSystem->Load("libCDB.so");
	gSystem->Load("libSTEER.so");
	gSystem->Load("libRAWDatasim.so");
	gSystem->Load("libFASTSIM.so");
	gSystem->Load("libEVGEN.so");
	gSystem->Load("libAliPythia6.so");
	gSystem->Load("libSTAT.so");
	gSystem->Load("libhijing.so");
	gSystem->Load("libTHijing.so");
	gSystem->Load("libSTRUCT.so");
	gSystem->Load("libPHOSUtils.so");
	gSystem->Load("libPHOSbase.so");
	gSystem->Load("libPHOSsim.so");
	gSystem->Load("libPHOSrec.so");
	gSystem->Load("libMUONcore.so");
	gSystem->Load("libMUONmapping.so");
	gSystem->Load("libMUONgeometry.so");
	gSystem->Load("libMUONcalib.so");
	gSystem->Load("libMUONraw.so");
	gSystem->Load("libMUONtrigger.so");
	gSystem->Load("libMUONbase.so");
	gSystem->Load("libMUONsim.so");
	gSystem->Load("libMUONrec.so");
	gSystem->Load("libMUONevaluation.so");
	gSystem->Load("libFMDbase.so");
	gSystem->Load("libFMDsim.so");
	gSystem->Load("libFMDrec.so");
	gSystem->Load("libPMDbase.so");
	gSystem->Load("libPMDsim.so");
	gSystem->Load("libPMDrec.so");
	gSystem->Load("libHMPIDbase.so");
	gSystem->Load("libHMPIDsim.so");
	gSystem->Load("libHMPIDrec.so");
	gSystem->Load("libT0base.so");
	gSystem->Load("libT0sim.so");
	gSystem->Load("libT0rec.so");
	gSystem->Load("libZDCbase.so");
	gSystem->Load("libZDCsim.so");
	gSystem->Load("libZDCrec.so");
	gSystem->Load("libACORDEbase.so");
	gSystem->Load("libACORDErec.so");
	gSystem->Load("libACORDEsim.so");
	gSystem->Load("libVZERObase.so");
	gSystem->Load("libVZEROrec.so");
	gSystem->Load("libVZEROsim.so");
	gSystem->Load("libEMCALraw.so");
	gSystem->Load("libEMCALUtils.so");
	gSystem->Load("libEMCALbase.so");
	gSystem->Load("libEMCALsim.so");
	gSystem->Load("libEMCALrec.so");
	gSystem->Load("libTPCbase.so");
	gSystem->Load("libTPCrec.so");
	gSystem->Load("libTPCsim.so");
	gSystem->Load("libITSbase.so");
	gSystem->Load("libITSsim.so");
	gSystem->Load("libITSrec.so");
	gSystem->Load("libTRDbase.so");
	gSystem->Load("libTRDsim.so");
	gSystem->Load("libTRDrec.so");
	gSystem->Load("libTOFbase.so");
	gSystem->Load("libTOFsim.so");
	gSystem->Load("libTOFrec.so");
	gSystem->Load("libHLTbase.so");
	gSystem->Load("libHLTinterface.so");
	gSystem->Load("libHLTsim.so");
	gSystem->Load("libHLTrec.so");
	gSystem->Load("libTENDER.so");
	gSystem->Load("libTENDERSupplies.so");
	gSystem->Load("libPWGPP.so");
	gSystem->Load("libITSUpgradeBase.so");
	gSystem->Load("libITSUpgradeSim.so");
	gSystem->Load("libITSUpgradeRec.so");
	gSystem->Load("libESDfilter.so");
	gSystem->Load("libHepMC.so");
	gSystem->Load("libTauola.so");
	gSystem->Load("libPhotos.so");
	gSystem->Load("libEvtGen.so");
	gSystem->Load("libEvtGenExternal.so");
	gSystem->Load("libTEvtGen.so");
	//
	printf("###               Done!              ###\n");
	printf("########################################\n");
	
}
//______________________________________________________________________________
Int_t main(){
	
	run_mcgenFT2();
	
}
