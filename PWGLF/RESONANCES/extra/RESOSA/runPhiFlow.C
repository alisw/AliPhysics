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
//Int_t RunNo        = 297590;
Int_t RunNo        = 296550;
//______________________________________________________________________________
void runPhiFlow(		
		//local or grid
		//const char* runtype = "local",
		const char* runtype = "grid", 
		//Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). 
		//Full & Test work for proof
		//const char *gridmode ="full", 
		//const char *gridmode = "terminate", 
		const char *gridmode = "test", 
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
	//TChain* chain = CreateESDChain("$ALICE_PHYSICS/PWGLF/RESONANCES/RESOSA/filelist.txt", 1);
	TChain* chain = CreateESDChain("alien:///alice/cern.ch/user/p/prottay/checktxt/filelist.txt", 1);
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
  //gROOT->LoadMacro("AddTaskPhiSA.C");
  AddTaskPhiSA(RunNo, analysislevel, kFALSE, bMCtruth, runtype, pairrapidity, systematiccut, 0);
  
  // enable debug printouts
  mgr->SetDebugLevel(2);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  // start analysis
  Printf("Starting Analysis....");
  
  if(runtype == "local")
    {
      //gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/Correlations/macros/dphicorrelations/CreateESDChain.C");
      //TChain* chain = CreateESDChain("$ALICE_PHYSICS/PWGLF/RESONANCES/RESOSA/filelist.txt", 1);
      //TChain* chain = CreateESDChain("filelist.txt", 1);
      mgr->StartAnalysis(runtype,chain);
    }
  if(runtype == "grid"){
    Printf("Starting Analysis....%s mode\n",runtype);
    mgr->StartAnalysis(runtype);
  }
  timer.Stop();
  timer.Print();
}
////////////////////////////////////
//AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
//task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral); 
////////////////////////////////////

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode, Int_t RunNo)
{
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    
    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode(gridmode);
    /*
    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    plugin->SetROOTVersion("v5-34-30-alice8-5");
    plugin->SetAliROOTVersion("v5-09-21-1");
    plugin->SetAliPhysicsVersion("vAN-20180118-1");
    */
    /*
    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    plugin->SetROOTVersion("v5-34-30-alice10-21");
    plugin->SetAliROOTVersion("v5-09-39-1");
    plugin->SetAliPhysicsVersion("vAN-20181022-1");
    */

    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    plugin->SetROOTVersion("v5-34-30-alice10-72");
    plugin->SetAliROOTVersion("v5-09-55-2");
    plugin->SetAliPhysicsVersion("vAN-20200910-1");

    plugin->AddIncludePath(addincpath.Data());
    // Method 1: Create automatically XML collections using alien 'find' command.
    // Define production directory LFN
    // plugin->SetGridDataDir("/alice/data/2010/LHC10h");
    plugin->SetGridDataDir("/alice/data/2018/LHC18q");
    // Data pattern for reconstructed data
    plugin->SetDataPattern("pass3/*ESDs.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
    
    //plugin->SetGridDataDir("/alice/data/2011/LHC11h_2");
    //plugin->SetDataPattern("*ESDs/pass2/*ESDs.root");

    plugin->SetRunPrefix("000");   // real data
    plugin->AddRunNumber(RunNo);
    //plugin->AddRunNumber(296550);
    //plugin->AddRunNumber(246984);
    /*
    plugin->AddRunNumber(137638); 
    plugin->AddRunNumber(137608); 
    plugin->AddRunNumber(137595); 
    plugin->AddRunNumber(137549); 
    plugin->AddRunNumber(137546); 
    plugin->AddRunNumber(137544); 
    plugin->AddRunNumber(137541);
    plugin->AddRunNumber(137539); 
    plugin->AddRunNumber(137531); 
    plugin->AddRunNumber(137530); 
    plugin->AddRunNumber(137443); 
    plugin->AddRunNumber(137441); 
    plugin->AddRunNumber(137440); 
    plugin->AddRunNumber(137439); 
    plugin->AddRunNumber(137434); 
    plugin->AddRunNumber(137432); 
    plugin->AddRunNumber(137431); 
    plugin->AddRunNumber(137366); 
    plugin->AddRunNumber(137243);
    plugin->AddRunNumber(137236); 
    plugin->AddRunNumber(137235); 
    plugin->AddRunNumber(137232); 
    plugin->AddRunNumber(137231);
    plugin->AddRunNumber(137162); 
    plugin->AddRunNumber(137161);
    */
    /*plugin->AddRunNumber(139510);
    plugin->AddRunNumber(139507);
    plugin->AddRunNumber(139505);
    plugin->AddRunNumber(139503);
    plugin->AddRunNumber(139465);
    plugin->AddRunNumber(139438);
    plugin->AddRunNumber(139437);
    plugin->AddRunNumber(139360);
    plugin->AddRunNumber(139329);
    plugin->AddRunNumber(139328);
    plugin->AddRunNumber(139314);
    plugin->AddRunNumber(139310);
    plugin->AddRunNumber(139309);
    plugin->AddRunNumber(139173);
    plugin->AddRunNumber(139107);
    plugin->AddRunNumber(139105);
    plugin->AddRunNumber(139038);
    plugin->AddRunNumber(139037);
    plugin->AddRunNumber(139036);
    plugin->AddRunNumber(139029);
    plugin->AddRunNumber(139028);
    plugin->AddRunNumber(138872);
    plugin->AddRunNumber(138871);
    plugin->AddRunNumber(138870);
    */
    //...then add run numbers to be considered
    //plugin->AddRunNumber(169238);
    //plugin->AddRunNumber(169167);
    //plugin->AddRunNumber(169145);
    
        
    plugin->SetNrunsPerMaster(1);
    plugin->SetOutputToRunNo();
    //comment out the next line when using the "terminate" option, unless
    //you want separate merged files for each run
    plugin->SetMergeViaJDL();

    //plugin->SetGridWorkingDir("2018_SA/Phi/systematics/data/LHC10h/rFnClp9");
    plugin->SetGridWorkingDir("2015_SA/Phi/systematics/data/LHC15o/rFnClp9/donelocally/checktxt2");
    //Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("out"); // In this case will be $HOME/taskname/out
    //Declare the analysis source files names separated by blancs. To be compiled runtime
    //using ACLiC on the worker nodes.
    plugin->SetAnalysisSource("$ALICE_PHYSICS/PWGLF/RESONANCES/extra/RESOSA/AliAnalysisTaskPhiSA.cxx");
    //plugin->SetAnalysisSource("AliAnalysisTaskPhiSA.cxx");
    //Declare all libraries (other than the default ones for the framework. These will be
    //plugin->SetAdditionalLibs("AliAnalysisTaskPhiSA.h AliAnalysisTaskPhiSA.cxx");
    // Add aditional AliRoot libraries
    plugin->SetAdditionalLibs("libPWGPPevchar.so libPWGPPevcharQn.so libPWGPPevcharQnInterface.so libGui.so libProof.so libMinuit.so libXMLParser.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libCORRFW.so libTOFbase.so libPWGmuon.so libPWGflowBase.so libPWGflowTasks.so AliAnalysisTaskPhiSA.h AliAnalysisTaskPhiSA.cxx");
    //plugin->SetAdditionalLibs("libPWGPPevchar.so");
    //loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.


    //Declare the output file names separated by blancs.
    //(can be like: file.root or file.root@ALICE::Niham::File)
    //To only save certain files, use SetDefaultOutputs(kFALSE), and then
    //SetOutputFiles("AnalysisResults.local.ESD.root") to choose which files to save
    //plugin->SetDefaultOutputs(kFALSE);
    //plugin->SetDefaultOutputs(kFALSE);

    //plugin->SetOutputFiles("AnalysisResultsWithPhiPtWeights.grid.ESD.root");//with phi wt
    //plugin->SetOutputFiles("AnalysisResultsWithPhiPtWtShifting.grid.ESD.root");//WO phi wt
    //plugin->SetOutputArchive("");//subhash

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


//////////////////////| LHC11h |///(good run list)//////////////
//170204, 170193, 170163, 170159, 170155, 170091, 170089, 170088, 170085, 170084,
//170083, 170081, 170040, 170027, 169965, 169923, 169859, 169858, 169855, 169846,
//169838, 169837, 169835, 169553, 169417, 169415, 169411, 169238, 169167, 169160,
//169156, 169148, 169145, 169144, 169138, 169094, 169091, 169035, 168992, 168988,
//168826, 168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362,
//168361, 168342, 168341, 168325, 168322, 168311, 168310, 168115, 168108, 168107,
//168105, 168076, 168069, 167988, 167987, 167985, 
////////////////////////////////////////////////////////////////


//////////////////////| LHC10h |//////////////////////////////////////
//139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310,
//139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870,
//138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582,
//138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201,
//138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693,
//137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541,
//137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137366, 137243,
//137236, 137235, 137232, 137231,137162, 137161
/////////////////////////////////////////////////////////////////////
