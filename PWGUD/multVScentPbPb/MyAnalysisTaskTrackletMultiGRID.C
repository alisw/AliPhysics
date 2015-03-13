Bool_t needRecPoints = kFALSE;
Bool_t InputHandlerSetup(TString format = "esd", Bool_t useKine = kTRUE);
AliAnalysisGrid PreparePlugin(TString mode, //can be "full", "test", "offline", "submit" or "terminate"
			      Int_t   runNumber,
			      TString dataDir,
			      TString dataPattern,
			      TString outFile,
			      TString runPrefix,
			      TString outDir="testOutput",
			      TString rootVer="v5-34-08-6",
			      TString alirootVer="v5-05-Rev-30");
//====================================================================
int nCentBins = 5;
const float centBins[] = {0,20,40,60,80,100};
UInt_t trigSel = AliVEvent::kINT7;
//
TString rootVer="v5-34-08-6";
TString alirootVer="v5-05-Rev-30";

//====================================================================
void MyAnalysisTaskTrackletMultiGRID
(
 TString mode="full",//test",
 Int_t run=195483,
 Int_t nEvents=-1,
 TString dataset="/alice/sim/2014/LHC14i2",
 TString dataPattern="*ESDs.root", //"*ESDs/pass4/*ESDs.root"
 TString outDir   = "mlt",
 Float_t etaMin     =-2.3,          // min eta range to fill in histos
 Float_t etaMax     = 2.3,          // max eta range to fill in histos
 Float_t zMin       = -18,         // process events with Z vertex min
 Float_t zMax       =  18,         //                     max positions
 const char* useCentVar = "V0A",          // centrality variable to use
 TString outFName = "trmult.root",
 //
 Float_t cutSigNStd  = 1.5,        // cut on weighed distance used to extract signal
 Float_t cutSigDPhiS = -1,        // cut on dPhi-phiBent used to extract signal (if negative -> dphi*sqrt(cutSigNStd)
 Bool_t  useMC  = kFALSE,          // fill MC info
 //
 Bool_t doRec  = kFALSE,          // fill data histos from new reco
 Bool_t doInj  = kFALSE,          // create Inj. bg
 Bool_t doRot  = kFALSE,          // create Rot. bg
 // 
 // specific parameters for reconstruction
 float  phiRot      = 3.14159e+00, // angle for bg. generation with rotation
 float  injScale    = 1.,//0.7,    // inject injScale*Ncl(Lr1/Lr2) hits
 Bool_t scaleDTheta = kTRUE,       // scale dTheta by 1/sin^2(theta) in trackleting
 float  nStdDev     = 25.,         // number of st.dev. for tracklet cut to keep
 float  dphi        = 0.06,        // dphi window (sigma of tracklet cut)
 float  dtht        = 0.025,       // dtheta .... (if negative, abs will be used with additional cut on |dthetaX|, apart from w.distance
 float  phishift    = 0.0045,      // bending shift
 Bool_t remOvl      = kTRUE,       
 float  ovlPhiCut   = 0.005, 
 float  ovlZetaCut  = 0.05,
 Bool_t checkReconstructables = kFALSE, //kTRUE, // fill histos for reconstructable (needs useMC and doRec) 
 //
 )
{
  //  
  if (cutSigDPhiS<0) cutSigDPhiS = TMath::Sqrt(cutSigNStd)*dphi;
  //
  needRecPoints = doRec || doInj || doRot;
    //
  printf("Start Analysis for %s,  Event Cuts: %.1f<eta<%.1f, %.2f<Zv<%.2f\n",
	 dataset.Data(),etaMin,etaMax,zMin,zMax);
  printf("Centrality variable: %s\n",useCentVar);
  printf("Tracklet cuts: dPhi:%.3f dTheta:%.3f phiShift:%.4f | Keep %.1f NstDev\n"
	 "Scale dTheta: %s | Signal Selection: NstDev:%.1f, dPhiS: %.3f\n", 
	 dphi,dtht,phishift,nStdDev,scaleDTheta ? "ON":"OFF",
	 cutSigNStd,cutSigDPhiS);
  //
  if (dataset.Contains("sim")) useMC = kTRUE; // always override if simulation
  printf("UseMC: %s \n",useMC ? "ON":"OFF");
  printf("Operations:  \n"
	 "Reco:%s (RemOvl:%s phi:%.3f zeta:%.3f)\n"
	 "Inj:%s (scale: %.2f)\n"
	 "Rot:%s (phi: %.4f)\n",
	 doRec ? "ON":"OFF",remOvl? "ON":"OFF",ovlPhiCut,ovlZetaCut,
	 doInj ? "ON":"OFF",injScale,
	 doRot ? "ON":"OFF",phiRot);
  //
  if (nEvents<0) nEvents = int(1e9);
  TString prefRun = "";
  if (!useMC) prefRun = "000";
  //
  outDir += run;
  if (useMC) outDir+="MC";
  //
  gROOT->ProcessLine(".include /home/shahoian/alicesw/aliroot/v5-05-Release/src/ITS");
  //
  gROOT->LoadMacro("AliITSMultRecBg.cxx++g");
  gROOT->LoadMacro("AliTrackletTaskMulti.cxx++g");
  // ALICE stuff
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("mgr");
  //
  AliAnalysisGrid* plugin = PreparePlugin(mode,run,dataset,dataPattern,outFName,
					  prefRun,outDir,rootVer,alirootVer);
  mgr->SetGridHandler(plugin);
  //
  InputHandlerSetup("esd",useMC);
  //
  printf("Loading Centrality task\n");
  gROOT->LoadMacro("$ALICE_ROOT/../src/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
  if (useMC) taskCentrality->SetMCInput();
  //
  //================================================================================
  printf("Requesting physics selection in %s mode\n",useMC ? "MC":"Data");
  gROOT->LoadMacro("$ALICE_ROOT/../src/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelectionTask = AddTaskPhysicsSelection(useMC);
  //==================================================================================
  //
  AliTrackletTaskMulti *mltTask = new AliTrackletTaskMulti("AliTrackletTaskMulti");
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist", TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outFName.Data());
  mgr->AddTask(mltTask);
  mgr->ConnectInput(mltTask, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(mltTask,1,coutput1);
  //
  mltTask->SetUseCentralityVar(useCentVar);
  mltTask->SetCentPercentiles(centBins, nCentBins);
  mltTask->SetTriggerSelection(trigSel);
  //
  mltTask->SetDoNormalReco(doRec);
  mltTask->SetDoInjection(doInj);
  mltTask->SetDoRotation(doRot);
  //
  mltTask->SetUseMC(useMC);
  mltTask->SetCheckReconstructables(checkReconstructables);
  //
  mltTask->SetEtaMin(etaMin);
  mltTask->SetEtaMax(etaMax);
  mltTask->SetZVertexMin(zMin);
  mltTask->SetZVertexMax(zMax);
  //
  mltTask->SetDPhiSCut(cutSigDPhiS);
  mltTask->SetNStdCut(cutSigNStd);
  //
  mltTask->SetScaleDThetaBySin2T(scaleDTheta);
  mltTask->SetNStdDev(nStdDev);
  mltTask->SetPhiWindow(dphi);
  mltTask->SetThetaWindow(dtht);
  mltTask->SetPhiShift(phishift);
  mltTask->SetPhiOverlapCut(ovlPhiCut);
  mltTask->SetZetaOverlapCut(ovlZetaCut);
  mltTask->SetPhiRot(phiRot);
  mltTask->SetInjScale(injScale);
  mltTask->SetRemoveOverlaps(remOvl);
  //
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid",nEvents);
  //
}


//________________________________________________________________
Bool_t InputHandlerSetup(TString format, Bool_t useKine)
{
  format.ToLower();
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cin = mgr->GetCommonInputContainer();
  if (cin) return;
  if (!format.CompareTo("esd"))  {
    AliESDInputHandler *esdInputHandler = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!esdInputHandler) {
      Info("CustomAnalysisTaskInputSetup", "Creating esdInputHandler ...");
      if (needRecPoints)
	esdInputHandler = new AliESDInputHandlerRP();
      else 
	esdInputHandler = new AliESDInputHandler();
      //
      mgr->SetInputEventHandler(esdInputHandler);
    }
    else if (needRecPoints && !esdInputHandler->InheritsFrom("AliESDInputHandlerRP")) {
      printf("Error: AliESDInputHandlerRP is needed\n");
      return kFALSE;
    }
    //
    if (useKine) {
      AliMCEventHandler* mcInputHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!mcInputHandler) {
        Info("CustomAnalysisTaskInputSetup", "Creating mcInputHandler ...");
        AliMCEventHandler* mcInputHandler = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(mcInputHandler);
      }
      mcInputHandler->SetPreReadMode(AliMCEventHandler::kLmPreRead);
    }    
  }
  else if (!format.CompareTo("aod")) {
    printf("This task requires ESDs\n"); return kFALSE;
    AliAODInputHandler *aodInputHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!aodInputHandler) {
      Info("CustomAnalysisTaskInputSetup", "Creating aodInputHandler ...");
      aodInputHandler = new AliAODInputHandler();
      mgr->SetInputEventHandler(aodInputHandler);
    }
  }
  else  {
    Info("Wrong input format!!! Only ESD and AOD are supported. Skipping Task ...");
    return kFALSE;
  }

  return kTRUE;
}


AliAnalysisGrid* PreparePlugin(TString mode, //can be "full", "test", "offline", "submit" or "terminate"
			       Int_t   runNumber,
			       TString dataDir,
			       TString dataPattern,
			       TString outFname,
			       TString runPrefix,
			       TString outDir,
			       TString rootVer,
			       TString alirootVer
			       )
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  //  if (!AliAnalysisGrid::CreateToken()) return NULL;
  // Load common libraries
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");   
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include $ALICE_ROOT/ITS");
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I./");
  //
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetRunMode(mode.Data());
  //
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion(rootVer.Data());
  plugin->SetAliROOTVersion(alirootVer.Data());
  //
  plugin->SetGridDataDir(dataDir.Data()); 
  plugin->SetDataPattern(dataPattern.Data());
  plugin->SetRunPrefix(runPrefix);
  plugin->AddRunNumber(runNumber);

  //  plugin->SetCheckCopy(kFALSE); //tmp
  //
  plugin->SetGridWorkingDir("mltTask");
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir(outDir.Data()); // In this case will be $HOME/work/output
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  plugin->SetAnalysisSource("AliITSMultRecBg.cxx AliTrackletTaskMulti.cxx");
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalLibs("AliITSMultRecBg.cxx AliITSMultRecBg.h AliTrackletTaskMulti.cxx AliTrackletTaskMulti.h");
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  plugin->SetDefaultOutputs(kFALSE);
  printf("*****SET %s\n",outFname.Data());
  plugin->SetOutputFiles(outFname.Data());
  // Optionally define the files to be archived.
  plugin->SetOutputArchive("log_archive.zip:stdout,stderr,*.root");
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("MacroMultTask.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(50);
  // Optionally modify the executable name (default analysis.sh)
  plugin->SetExecutable("ExecutableMultTask.sh");

  plugin->SetExecutableCommand("aliroot -b -q");
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  //   plugin->SetMaxInitFailed(5);
  // Optionally resubmit threshold.
  //   plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(28000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("JdlMultTask.jdl");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");

  // merging via
  plugin->SetMergeViaJDL(kTRUE);
  //plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeStages(5); 
  //
  
  plugin->SetUser("shahoian");

  return plugin;
  //
}
