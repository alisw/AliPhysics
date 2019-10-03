#include <iostream>
#include <cstdio>
#include <ctime>
#include <time.h>
#include "TGrid.h"

AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers, 
                                     const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, 
				    Int_t workerTTL, Bool_t isMC, Int_t nFiles);
                                    
//______________________________________________________________________________
void runAODfiltering(
         Bool_t         useGrid             = kTRUE,
         UInt_t         numLocalFiles       = 1,
         TString        runPeriod           = "LHC16g1",
	 Int_t          iCollision          = 1,           // pp:0 , PbPb:1, pPb:2
         Bool_t         isMC                = kTRUE,
	 Bool_t         applyPS             = kTRUE,   
	 Int_t          PassForTuneOnData   = 1,           // recoPass (for MC to tune on which data pass)                                   ; default  2
	 Int_t          CustomSplinesPass   = -1,          // for data & MC: recoDataPass (to use custom splines object from this data pass) ; default -1
	 Bool_t         debug               = 0
     )
{

  // settings in case of grid calculation
  const char*    gridMode            = "test";                        // set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  const char*    pattern             = "*/AliESDs.root";      // file pattern (here one can specify subdirs like passX etc.) 
  const char*    gridDir             = "/alice/sim/2016/LHC16g1/";    // dir on alien, where the files live 
  const char*    runNumbers          = "246272";                      // considered run numbers 
  const char*    uniqueName          = "Refiltering";                 // name of the task 
  Int_t          nFiles              = 1;                             // number of files for local test
  Int_t          maxFilesPerWorker   = 20;                            // larger ->less jobs (can submit max 1500 jobs)
  Int_t          workerTTL           = 70000;                         // time to live [s] for each job


  gSystem->SetFPEMask();
  gSystem->Setenv("ETRAIN_ROOT", ".");
  gSystem->Setenv("ETRAIN_PERIOD", runPeriod);

  // Load necessary libraries
  LoadLibs();

  // Create analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager(uniqueName);

  // Check type of input and create handler for it
  TString localFiles("-1");
  localFiles = Form("/home/meike/analysis/data/fileLists/files_esd_%s.txt", runPeriod.Data());
  if(!useGrid) cout << "Using " << localFiles.Data() << " as input file list.\n";

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
  AliESDInputHandler* esdH = AddESDHandler();
  mgr->SetInputEventHandler(esdH);

  AliAODHandler* aodHandler   = new AliAODHandler();
  aodHandler->SetOutputFileName("AliAOD.root");
  mgr->SetOutputEventHandler(aodHandler);

  // Create MC handler, if MC is demanded
  if (isMC)
  {
    AliMCEventHandler* mcH = new AliMCEventHandler();
    mcH->SetPreReadMode(AliMCEventHandler::kLmPreRead);
    mcH->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(mcH); 
  }

  //========= Add CDBManager ====
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect("raw://",0);
    if (!taskCDB) return;

  // ============= add physics selection task========================
  if(applyPS){
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    mgr->RegisterExtraFile("event_stat.root");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);
    mgr->AddStatisticsTask(AliVEvent::kAny);
  }

  //============ Centrality/Multiplicity Selection==============
    if(iCollision){
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
      AliMultSelectionTask *taskCentrality = AddTaskMultSelection(); 
      if(isMC){ 
      taskCentrality->SetAlternateOADBforEstimators("LHC15o-DefaultMC-HIJING");
      }
    }

  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    AliAnalysisTask* taskPID = 0x0;
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    taskPID = AddTaskPIDResponse(isMC, kTRUE, isMC, PassForTuneOnData, kFALSE, "", kTRUE, kTRUE, CustomSplinesPass);      
  }    
   
                                                                                         
//========= Add Production Tasks ====
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskFilteredTree.C");
  AddTaskFilteredTree("FilterEvents_Trees.root");

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_ConversionAODProduction.C");
  AliAnalysisTask *taskconv = AddTask_ConversionAODProduction(iCollision, isMC, runPeriod.Data());
  mgr->RegisterExtraFile("AliAODGammaConversion.root");

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C");
  printf("Registering delta AOD file\n");
  mgr->RegisterExtraFile("AliAOD.Muons.root");
  mgr->RegisterExtraFile("AliAOD.Dimuons.root");
  AliAnalysisTaskESDfilter *taskesdfilter =
               AddTaskESDFilter(isMC,            // useKFILTER
                                kFALSE,         // write Muon AOD
                                kFALSE,         // write dimuon AOD
                                kFALSE,         // usePhysicsSelection
                                kFALSE,         // centrality OBSOLETE
                                kTRUE,          // enable TPS only tracks
                                kFALSE,         // disable cascades
                                kFALSE,         // disable kinks
                                1500,           // run flag (YY00) // The first 2 digits are the year, the second
                                                                    //2 digits are used to distinguish sub-periods (if needed)
                                3,
                                kFALSE,                  // useV0Filter
				kFALSE,                 // add SPD information to muon AOD only for pp
				kFALSE,                // isMuonCaloPass
				kFALSE);                // addPCMv0s, default kTRUE

  AliEMCALGeometry::GetInstance("","");

  //================================================================================================

  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();
  
  if (useGrid) 
    {  // GRID CALCULATION

      AliAnalysisGrid *plugin = CreateAlienHandler(uniqueName, gridDir, gridMode, runNumbers, pattern, maxFilesPerWorker, workerTTL, isMC, nFiles);
      mgr->SetGridHandler(plugin);
      
      // start analysis
      if(debug){
	mgr->SetUseProgressBar(kFALSE);
	mgr->SetDebugLevel(AliLog::kDebug);
      }
      else{
	mgr->SetUseProgressBar(1,100);
	mgr->SetDebugLevel(0);  
      }
      cout << "Starting GRID Analysis..." << endl;
      mgr->StartAnalysis("grid");
    }
  else
    {  // LOCAL CALCULATION
      TChain* chain = 0;
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
      chain = CreateESDChain(localFiles.Data(), numLocalFiles);
      
      // start analysis
      if(debug){
	mgr->SetUseProgressBar(kFALSE);
	mgr->SetDebugLevel(AliLog::kDebug);
      }
      else{
	mgr->SetUseProgressBar(1,10);
	mgr->SetDebugLevel(0);  
      }
      cout << "Starting LOCAL Analysis..." << endl;
      //TStopwatch *watch = new TStopwatch();
      mgr->StartAnalysis("local", chain);
      //watch->Stop();
      //cout << "*************************" << endl;
      //cout << "spent real time: " << (Int_t)(watch->RealTime()) << "s"<< endl;
      //cout << "spent cpu time: " << (Int_t)(watch->CpuTime()) << "s" << endl;
      //cout << "*************************" << endl;
    }
}

//______________________________________________________________________________
void LoadLibs()
{

  // Load common libraries (better too many than too few)
	gSystem->Load("libCore.so");  
	gSystem->Load("libTree.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libVMC.so");
	gSystem->Load("libGui.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libMinuit");
	gSystem->Load("libMinuit2");
	gSystem->Load("libSTEERBase.so");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libOADB");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");  
	gSystem->Load("libPWGGAGammaConv.so");
	gSystem->Load("libEve.so");   
	gSystem->Load("libCDB.so");
	gSystem->Load("libProof.so");
	gSystem->Load("libRAWDatabase.so");
	gSystem->Load("libSTEER.so");
	gSystem->Load("libEVGEN.so");
	gSystem->Load("libESDfilter.so");	
	gSystem->Load("libTender.so");
	gSystem->Load("libTOFbase.so");
	gSystem->Load("libTOFsim.so");
	gSystem->Load("libTOFrec.so");
	gSystem->Load("libTRDbase.so");
	gSystem->Load("libVZERObase.so");
	gSystem->Load("libVZEROrec.so");
	//gSystem->Load("libTenderSupplies.so");
	gSystem->Load("libXMLParser.so");
	gSystem->Load("libSTAT.so");
	gSystem->Load("libRAWDatarec.so");
	gSystem->Load("libANALYSIScalib.so");
	gSystem->Load("libCORRFW.so");
	gSystem->Load("libPWGUDbase.so");
	gSystem->Load("libTPCbase.so");
	gSystem->Load("libTPCrec.so");
	gSystem->Load("libTPCcalib.so");
	gSystem->Load("libTRDrec.so");
	gSystem->Load("libITSbase.so");
	gSystem->Load("libITSrec.so");
	gSystem->Load("libHMPIDbase.so");
	gSystem->Load("libPWGmuon.so");
	gSystem->Load("libPWGPP.so");
	gSystem->Load("libPWGHFbase.so");
	gSystem->Load("libPWGDQdielectron.so");
	gSystem->Load("libPWGHFhfe.so");
	gSystem->Load("libPHOSUtils.so");
	gSystem->Load("libPHOSbase.so");
	gSystem->Load("libPHOSpi0Calib.so");
	gSystem->Load("libPHOSrec.so");
	gSystem->Load("libPWGGAPHOSTasks.so");
	gSystem->Load("libPHOSshuttle.so");
	gSystem->Load("libPHOSsim.so");
	gSystem->Load("libPWGCaloTrackCorrBase.so");
	gSystem->Load("libEMCALUtils.so");
	gSystem->Load("libEMCALraw.so");
	gSystem->Load("libEMCALbase.so");
	gSystem->Load("libEMCALrec.so");
	gSystem->Load("libPWGTools.so");
	//gSystem->Load("libPWGEMCAL.so");
	gSystem->Load("libPWGGAUtils.so");
	gSystem->Load("libPWGGAEMCALTasks.so");
	gSystem->Load("libPWGGAGammaConv.so");
	gSystem->Load("libPWGGAPHOSTasks.so");
	gSystem->Load("libPWGCFCorrelationsBase.so");
	gSystem->Load("libPWGCFCorrelationsDPhi.so");
  
}

AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers, 
                                     const char* pattern, Int_t maxFilesPerWorker, 
				    Int_t workerTTL, Bool_t isMC, Int_t nFiles)
{


	AliAnalysisAlien *plugin = new AliAnalysisAlien();

	/*************** versions etc. *****************/
	//if(!strcmp(gridMode, "full")) plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include"); // when own aliphysics is used via .par file
	//else plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
	plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

	//plugin->SetAdditionalLibs("libTree.so libVMC.so libGui.so libXMLParser.so libMinuit.so libMinuit2.so libProof.so libGeom.so libPhysics.so libSTEERBase.so libESD.so libAOD.so libOADB.so libANALYSIS.so libANALYSISalice.so libCDB.so libRAWDatabase.so libSTEER.so libEVGEN.so libESDfilter.so libCORRFW.so libTOFbase.so libTOFsim.so libTOFrec.so libRAWDatarec.so libTPCbase.so libTPCrec.so libTPCcalib.so libITSbase.so libITSrec.so libTRDbase.so libTRDrec.so libSTAT.so libHMPIDbase.so libVZERObase.so libVZEROrec.so libTENDERSupplies libPWGmuon.so libPWGPP.so libPWGHFbase.so libPWGDQdielectron.so libPWGHFhfe.so libPHOSUtils.so libPHOSbase.so libPHOSpi0Calib.so libPHOSrec.so libPWGGAPHOSTasks.so libPHOSshuttle.so libPHOSsim.so libPWGCaloTrackCorrBase.so libEMCALUtils.so libEMCALraw.so libEMCALbase.so libEMCALrec.so libPWGTools.so libPWGGAUtils.so libPWGGAEMCALTasks.so libPWGGAPHOSTasks.so libPWGCFCorrelationsBase.so libPWGCFCorrelationsDPhi.so libPWGGAGammaConv.so");
	plugin->SetAdditionalLibs("libTree.so libVMC.so libGui.so libXMLParser.so libMinuit.so libMinuit2.so libProof.so libGeom.so libPhysics.so libSTEERBase.so libESD.so libAOD.so libOADB.so libANALYSIS.so libANALYSISalice.so libCDB.so libRAWDatabase.so libSTEER.so libEVGEN.so libESDfilter.so libCORRFW.so libTOFbase.so libTOFsim.so libTOFrec.so libRAWDatarec.so libTPCbase.so libTPCrec.so libTPCcalib.so libITSbase.so libITSrec.so libTRDbase.so libTRDrec.so libSTAT.so libHMPIDbase.so libVZERObase.so libVZEROrec.so libPWGmuon.so libPWGPP.so libPWGHFbase.so libPWGDQdielectron.so libPWGHFhfe.so libPHOSUtils.so libPHOSbase.so libPHOSpi0Calib.so libPHOSrec.so libPWGGAPHOSTasks.so libPHOSshuttle.so libPHOSsim.so libPWGCaloTrackCorrBase.so libEMCALUtils.so libEMCALraw.so libEMCALbase.so libEMCALrec.so libPWGTools.so libPWGGAUtils.so libPWGGAEMCALTasks.so libPWGGAPHOSTasks.so libPWGCFCorrelationsBase.so libPWGCFCorrelationsDPhi.so libPWGGAGammaConv.so");

	plugin->SetAliPhysicsVersion("vAN-20170531-1");

	if(!strcmp(gridMode, "full")){
	  cout << "Will use .par files!" << endl;
	  plugin->SetupPar("ESDfilter.par");  
	  plugin->EnablePackage("ESDfilter.par");   
	  plugin->SetupPar("PWGGAGammaConv.par");  // Compile the par file // The current directory should have precedence in LD_LIBRARY_PATH
	  plugin->EnablePackage("PWGGAGammaConv.par");  // par file must be in current directory
	}      
       

	/************ names of files and folders ***********************/
        TDatime currentTime;
	TString tmpName(uniqueName);
	// Only add current date and time when not in terminate mode! In this case the exact name has to be supplied by the user
	if(strcmp(gridMode, "terminate"))
	{
		tmpName += "_";
		tmpName += currentTime.GetDate();
		tmpName += "_";
		tmpName += currentTime.GetTime();
	}
	TString macroName = Form("%s.C", tmpName.Data());
	TString execName = Form("%s.sh", tmpName.Data());
	TString jdlName = Form("%s.jdl", tmpName.Data());
	plugin->SetAnalysisMacro(macroName.Data());
	plugin->SetExecutable(execName.Data());
	plugin->SetJDLName(jdlName.Data());
	plugin->SetGridWorkingDir(Form("work/%s",tmpName.Data()));
	plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/tmpName/<GridOutputDir>/<runNo>/<i>/		

	/******************** settings ********************************/
	plugin->SetRunMode(gridMode);  
	plugin->SetGridDataDir(gridDir); 
	plugin->SetDataPattern(pattern); 
	if (!isMC) plugin->SetRunPrefix("000");
	plugin->AddRunList(runNumbers);
	plugin->SetNtestFiles(nFiles);  
	plugin->SetSplitMaxInputFileNumber(maxFilesPerWorker);
	plugin->SetTTL(workerTTL);
	plugin->SetMaxMergeStages(1);
	plugin->SetMergeViaJDL(1);
	plugin->SetOutputToRunNo(kTRUE);
	plugin->SetKeepLogs(kTRUE);
	return plugin;
}
