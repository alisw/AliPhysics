#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "AliAnalysisGrid.h"
#include "TSystem.h"
#include "TROOT.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisGrid.h"
#include "AliVEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisAlien.h"
#include "AliMCEventHandler.h"
#include "AliPhysicsSelectionTask.h"
#include "TRegexp.h"
#include "AliTriggerAnalysis.h"
#include "TChain.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliAODHandler.h"
#include "AliAnalysisTaskNanoAODFilter.h"
#include "AliESEHelpers.h"

#endif
void LoadLibs();

class AliAnalysisGrid;
class AliESETrkCut;
class AliESEEvtCut;

AliESETrkCut * TrkCuts() {

  AliESETrkCut * trk = new AliESETrkCut;

  AliSpectraAODTrackCuts  * trcuts = new AliSpectraAODTrackCuts("TrackCuts");  
  trcuts->SetDCA(100000);
  trcuts->SetTrackBits(1);
  trcuts->SetPt(50);
  trcuts->SetPtTOFMatching(0.6);   
  trcuts->SetEta(-0.8,0.8);
  trcuts->SetMinTPCcls(70);
  trcuts->PrintCuts();

  trk->SetTrackCuts(trcuts);
  trk->Init();

  return trk;

}

AliESEEvtCut * EvtCuts(Int_t mc) {

  AliESEEvtCut * evt = new AliESEEvtCut;

  AliSpectraAODEventCuts * evcuts = new AliSpectraAODEventCuts("EventCuts");
  evcuts->SetQVectorCut(0,100);
  evcuts->SetCentralityCutMax(100);  
  evcuts->SetCentralityCutMin(0);
  if(mc>0)evcuts->SetIsMC(kTRUE);
  TFile * fCalib = new TFile("./calibV0New.root");
  evcuts->SetCalibFile(fCalib);
  evcuts->SetIsLHC10h(kTRUE);
  evcuts->PrintCuts();
  //  evcuts->SetEventSelectionBit(AliVEvent::kAny);

  evt->SetEventCuts(evcuts);
  evt->Init();

  return evt;

}

//______________________________________________________________________________
void runGridESE(
		 const char *gridmode = "test", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). Full & T		 
		 const char * taskname = "DevNanoAOD",
		 const int iMCtruth = 2, 
		 const char * addTaskString = ".x AddTaskNanoAODFilter.C(%d,1)" // 
		 )
{
  LoadLibs();
  // analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager("NanoAOD Filter", "NanoAOD filter for nanoAOD production");

  AliAnalysisGrid *plugin = CreateAlienHandler(taskname, gridmode); 
  plugin->SetFileForTestMode("files.txt"); // file should contain path name to a local directory containg *ESDs.root or AOD etc
  mgr->SetGridHandler(plugin);

    
  AliAODInputHandler* iH = new AliAODInputHandler();
  mgr->SetInputEventHandler(iH);

  // Define aod output handler
  AliAODHandler* aodOutputHandler = new AliAODHandler();
  aodOutputHandler->SetOutputFileName("AliAOD.NanoAOD.root");
  mgr->SetOutputEventHandler(aodOutputHandler);
  
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse(iMCtruth);
  taskPID->SetUseTPCEtaCorrection(kTRUE); 

  // create task
  cout << "Macro: "<< addTaskString << " " << Form(addTaskString, iMCtruth) << endl;

  AliAnalysisTaskNanoAODFilter * task = (AliAnalysisTaskNanoAODFilter*) gROOT->ProcessLine(Form(addTaskString, iMCtruth));

  // Set Track event and vertex cuts here!
  task->SetVarListTrack("pt,theta,phi,cstNSigmaTPCPi,cstNSigmaTPCKa,cstNSigmaTPCPr,cstNSigmaTOFPi,cstNSigmaTOFKa,cstNSigmaTOFPr,cstBayesTPCPi,cstBayesTPCKa,cstBayesTPCPr,cstBayesTOFPi,cstBayesTOFKa,cstBayesTOFPr");
  task->SetVarListHeader("cstCentr,cstQVec");
  AliESETrkCut * trkCuts = TrkCuts();
  AliESEEvtCut * evtCuts = EvtCuts(iMCtruth);
  evtCuts->SetTrackCuts(trkCuts->GetTrackCuts());
  AliAnalysisESESetter * setter  = new AliAnalysisESESetter;
  setter->SetEventCuts(evtCuts->GetEventCuts());

  task->SetTrkCuts(trkCuts);
  task->SetEvtCuts(evtCuts);
  task->AddSetter(setter);

  //task->SelectCollisionCandidates(AliVEvent::kMB);// FIXME
  // enable debug printouts
  mgr->SetDebugLevel(10);
  //    mgr->SetNSysInfo(100);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  // start analysis
  // Always read the same file:
  TChain * chain = new TChain("aodTree");
  chain->Add("./AliAOD.root");

  Printf("Starting Analysis....");
  mgr->StartAnalysis("grid", chain,123456789);

}

//______________________________________________________________________________

AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode)
{
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode(gridmode);

    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    plugin->SetROOTVersion("v5-34-08-5");
    plugin->SetAliROOTVersion("v5-05-73-AN");

    // Declare input data to be processed.
    //    plugin->SetCheckCopy(kFALSE);

    // Method 1: Create automatically XML collections using alien 'find' command.
    // Define production directory LFN
    plugin->SetGridDataDir("/alice/sim/LHC11a10a_bis");
    // On real reconstructed data:
    // plugin->SetGridDataDir("/alice/data/2009/LHC09d");
    // Set data search pattern
    //plugin->SetDataPattern("*ESDs.root"); // THIS CHOOSES ALL PASSES
    // Data pattern for reconstructed data
    plugin->SetDataPattern("*AOD120/*/AliAOD.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
    //    plugin->SetDataPattern("ESDs/pass2/AOD038/*AliAOD.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
    //    plugin->SetRunPrefix("000");   // real data
    // ...then add run numbers to be considered
    Int_t runlist[1]={137366};  
    for (Int_t ind=0; ind<1; ind++) {
//     plugin->AddRunNumber(138275);
      plugin->AddRunNumber(runlist[ind]);
    }
    //plugin->SetRunRange(114917,115322);
    plugin->SetNrunsPerMaster(10); // 1
    plugin->SetOutputToRunNo();
    // comment out the next line when using the "terminate" option, unless
    // you want separate merged files for each run
    plugin->SetMergeViaJDL();

    // Method 2: Declare existing data files (raw collections, xml collections, root file)
    // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
    // XML collections added via this method can be combined with the first method if
    // the content is compatible (using or not tags)
    //   plugin->AddDataFile("tag.xml");
    //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");

    // Define alien work directory where all files will be copied. Relative to alien $HOME.
    plugin->SetGridWorkingDir(taskname);

    // Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("out"); // In this case will be $HOME/taskname/out

    // Declare the analysis source files names separated by blancs. To be compiled runtime
    // using ACLiC on the worker nodes.
    //    plugin->SetAnalysisSource("AliAnalysisTaskEx01.cxx");

    // Declare all libraries (other than the default ones for the framework. These will be
    // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    plugin->SetAdditionalLibs("libPWGLFspectra.so libPWGDevNanoAOD.so");

    // Declare the output file names separated by blancs.
    // (can be like: file.root or file.root@ALICE::Niham::File)
    // To only save certain files, use SetDefaultOutputs(kFALSE), and then
    // SetOutputFiles("list.root other.filename") to choose which files to save
    plugin->SetDefaultOutputs();
    //plugin->SetOutputFiles("list.root");

    // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    plugin->SetAnalysisMacro(Form("%s.C",taskname));

    // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
    plugin->SetSplitMaxInputFileNumber(100);

    // Optionally modify the executable name (default analysis.sh)
    plugin->SetExecutable(Form("%s.sh",taskname));

    // set number of test files to use in "test" mode
    plugin->SetNtestFiles(10);

    // Optionally resubmit threshold.
    plugin->SetMasterResubmitThreshold(90);

    // Optionally set time to live (default 30000 sec)
    plugin->SetTTL(30000);

    // Optionally set input format (default xml-single)
    plugin->SetInputFormat("xml-single");

    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName(Form("%s.jdl",taskname));

    // Optionally modify job price (default 1)
    plugin->SetPrice(1);      

    // Optionally modify split mode (default 'se')    
    plugin->SetSplitMode("se");
    
    return plugin;
}



void LoadLibs() {
  gSystem->Load("libCore");  
  gSystem->Load("libGeom");
  gSystem->Load("libPhysics");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libProof");
  gSystem->Load("libMatrix");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  //  return;
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");

  //  gSystem->Load("libNanoAOD");
  gSystem->Load("libPWGLFspectra");
  gSystem->Load("libPWGDevNanoAOD");

}
