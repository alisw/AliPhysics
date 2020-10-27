#ifdef __ECLIPSE_IDE
//  few includes and external declarations just for the IDE
#include "Riostream.h"
#include "TSystem.h"
#include "TChain.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliAnalysisTaskPIDqa.h"
#include "AliPhysicsSelectionTask.h"
#include "AliCentralitySelectionTask.h"
#include "AliTaskCDBconnect.h"
#include "AliAnalysisTaskPIDCombined.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMultSelectionTask.h"
#include "runAnalysis.H"


#include "AliQnCorrectionsCutsSet.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsHistos.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowTrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include <AliVParticle.h>
#include <AliVTrack.h>
#include "AliVVertex.h"
#include <AliLog.h>
#include <AliPID.h>
#include <AliPIDResponse.h>
#include <AliITSPIDResponse.h>
#include <AliTPCPIDResponse.h>
#include <AliTRDPIDResponse.h>
#include <AliTOFPIDResponse.h>
#include "AliTriggerAnalysis.h"
#include "AliOADBContainer.h"
#include "AliEventplane.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskPhiSA.h"

//----- must include-------                                                                                                                    
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliPhysicsSelection.h"
#include "AliFlowEventSimple.h"





#endif // ifdef __ECLIPSE_IDE declaration and includes for the ECLIPSE IDE

class AliAnalysisGrid;
TString addincpath = "-I$ALICE_PHYSICS/include -I$ALICE_ROOT/include";
Int_t RunNo        = 297590;
//______________________________________________________________________________
void runPhiFlow(		
		//local or grid
		//const char* runtype = "local",
		const char* runtype = "grid", 
		//Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). 
		//Full & Test work for proof
		const char *gridmode ="full", 
		//const char *gridmode = "terminate", 
		//const char *gridmode = "test", 
		//1 = MCEvent handler is on (MC truth), 0 = MCEvent handler is off (MC reconstructed/real data)
		const bool bMCtruth = 0, 
		//1 = looking at MC truth or reconstructed, 0 = looking at real data
		const bool bMCphyssel = 0, 
		const char *analysislevel="ESD",
		//const Bool_t useShift = kTRUE
		const Bool_t useShift = kFALSE,
		const Float_t pairrapidity=0.5,
		const char *systematiccut = "RClsCrRowsOvFCls0.9"

		)
{
  TStopwatch timer;
  timer.Start();
  // check run type
  if(runtype != "local" && runtype != "grid"){
    Printf("\n\tIncorrect run option, check first argument of run macro");
    Printf("\tint runtype = local or grid\n");
    return;
  }
  Printf("%s analysis chosen",runtype);



  gSystem->Load("libCore.so");        
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");   
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libEventMixing.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGLFresonances.so");
  gSystem->Load("libPWGPPevcharQn.so");
  gSystem->Load("libPWGPPevcharQnInterface.so");


  
  // sets name of grid generated macros
  const char *taskname = Form("NewTpcTofPidTaskPhiWt%d",RunNo);

  // add aliroot indlude path
  gSystem->AddIncludePath(addincpath.Data());
  
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");


  AliAnalysisGrid      *alienHandler   =   NULL;
  
  AliAnalysisManager *mgr = new AliAnalysisManager(taskname);

  // analysis manager
  if(runtype == "grid"){
    AliAnalysisGrid *plugin = CreateAlienHandler(taskname, gridmode, RunNo); 
    mgr->SetGridHandler(plugin);
  }


  if(analysislevel == "ESD")
    {
      if(runtype == "local"){
	gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/Correlations/macros/dphicorrelations/CreateESDChain.C");
	TChain* chain = CreateESDChain("filelist.txt", 1);
      }
      AliVEventHandler* esdH = new AliESDInputHandler();
      mgr->SetInputEventHandler(esdH);
      
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(0);
      mgr->AddTask(physSelTask);
      if(!physSelTask) { Printf("no physSelTask"); return; }
      /*       
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"); // Needed for LHC2015o
      AliMultSelectionTask *taskM = AddTaskMultSelection(kTRUE);            // kFALSE == User mode, kTRUE == Calibration mode
      taskM->SetSelectedTriggerClass(AliVEvent::kINT7); 
      mgr->AddTask(taskM);
      */
      //taskCentrality->SetPass(2);    
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
      AliEPSelectionTask *taskEP = AddTaskEventplane();
    }
  
  if(analysislevel == "AOD")
    {
      if(runtype == "local"){
	TChain* chain = new TChain("aodTree");
	chain->Add("/wrk15a/ajay/Data/AliAOD.root");
      } 
     
      //Add AOD handler
      AliVEventHandler *aodH = new AliAODInputHandler();
      mgr->SetInputEventHandler(aodH);
    }
  // mc event handler
  if(bMCtruth) {
    AliMCEventHandler* mchandler = new AliMCEventHandler();
    // Not reading track references
    mchandler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mchandler);
  }   
  
  // === PID RESPONSE==================================================
  gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPIDResponse.C");
  AddTaskPIDResponse(kFALSE);


  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/extra/RESOSA/AddTaskPhiSA.C");
  AddTaskPhiSA(RunNo, analysislevel, kFALSE, bMCtruth, runtype, pairrapidity, systematiccut, 1);
  
  // enable debug printouts
  mgr->SetDebugLevel(2);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  // start analysis
  Printf("Starting Analysis....");
  
  if(runtype == "local")
    {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/Correlations/macros/dphicorrelations/CreateESDChain.C");
      TChain* chain = CreateESDChain("filelist.txt", 1);
    mgr->StartAnalysis(runtype,chain);
    }
  if(runtype == "grid"){
    Printf("Starting Analysis....%s mode\n",runtype);
    mgr->StartAnalysis(runtype);
  }
  timer.Stop();
  timer.Print();
}

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode, Int_t RunNo)
{
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    
    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode(gridmode);
 
    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    plugin->SetROOTVersion("v5-34-30-alice10-72");
    plugin->SetAliROOTVersion("v5-09-55-2");
    plugin->SetAliPhysicsVersion("vAN-20200910-1");

    plugin->AddIncludePath(addincpath.Data());
    // plugin->SetGridDataDir("/alice/data/2010/LHC10h");
    plugin->SetGridDataDir("/alice/data/2018/LHC18r");
    // Data pattern for reconstructed data
    plugin->SetDataPattern("pass3/*ESDs.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
    
    //plugin->SetGridDataDir("/alice/data/2011/LHC11h_2");
    //plugin->SetDataPattern("*ESDs/pass2/*ESDs.root");

    plugin->SetRunPrefix("000");   // real data
    plugin->AddRunNumber(RunNo);
    //plugin->AddRunNumber(246984);
    
        
    plugin->SetNrunsPerMaster(1);
    plugin->SetOutputToRunNo();
    //comment out the next line when using the "terminate" option, unless
    //you want separate merged files for each run
    plugin->SetMergeViaJDL();

    //plugin->SetGridWorkingDir("2018_SA/Phi/systematics/data/LHC10h/rFnClp9");
    plugin->SetGridWorkingDir("2015_SA/Phi/systematics/data/LHC15o/rFnClp9/donelocally/callcheck20");
    //Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("out"); // In this case will be $HOME/taskname/out
    //Declare the analysis source files names separated by blancs. To be compiled runtime
    //using ACLiC on the worker nodes.
    plugin->SetAnalysisSource("$ALICE_PHYSICS/PWGLF/RESONANCES/extra/RESOSA/AliAnalysisTaskPhiSA.cxx");
    //Declare all libraries (other than the default ones for the framework. These will be
    //plugin->SetAdditionalLibs("AliAnalysisTaskPhiSA.h AliAnalysisTaskPhiSA.cxx");
    // Add aditional AliRoot libraries
    plugin->SetAdditionalLibs("libPWGPPevchar.so libPWGPPevcharQn.so libPWGPPevcharQnInterface.so libGui.so libProof.so libMinuit.so libXMLParser.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libCORRFW.so libTOFbase.so libPWGmuon.so libPWGflowBase.so libPWGflowTasks.so AliAnalysisTaskPhiSA.h AliAnalysisTaskPhiSA.cxx");
    //plugin->SetAdditionalLibs("libPWGPPevchar.so");
    //loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.


    
    //Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    plugin->SetAnalysisMacro(Form("%s.C",taskname));
    //Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
    plugin->SetSplitMaxInputFileNumber(100);
    plugin->SetExecutableCommand("aliroot -q -b");
    plugin->SetExecutableArgs(">& std.log");
    //Optionally modify the executable name (default analysis.sh)
    plugin->SetExecutable(Form("%s.sh",taskname));
    //set number of test files to use in "test" mode
    plugin->SetNtestFiles(2);
    //Optionally resubmit threshold.
    plugin->SetMasterResubmitThreshold(90);
    //Optionally set time to live (default 30000 sec)
    plugin->SetTTL(30000);
    //Optionally set input format (default xml-single)
    plugin->SetInputFormat("xml-single");
    //Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName(Form("%s.jdl",taskname));
    //Optionally modify job price (default 1)
    plugin->SetPrice(1);      
    //Optionally modify split mode (default 'se')    
    plugin->SetSplitMode("se");
    return plugin;
}

