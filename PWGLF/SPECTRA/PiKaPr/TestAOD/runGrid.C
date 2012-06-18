class  AliAnalysisManager;
class  AliAnalysisAlien;

void runGrid(TString mode="test",Int_t mc=1,Int_t day=15,Int_t month=6, Int_t year=2012) 
{
  AliLog::SetGlobalDebugLevel(100);
  // Load common libraries
  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
  
  // Load common libraries
  gSystem->Load("libCore.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so"); 
  gSystem->Load("libGui.so");
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libCDB.so");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libProof.so");
  gSystem->Load("libRAWDatabase.so");
  gSystem->Load("libSTEER.so");
  //__________________________________________________________________________
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gSystem->Load("libPWGLFspectra.so");
  
  // gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD/HistogramNames.h");
  // gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD/AliSpectraAODTrackCuts.cxx+g");
  // gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD/AliSpectraAODEventCuts.cxx+g");
  // gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD/AliSpectraAODHistoManager.cxx+g");
  // gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD/AliSpectraAODPID.cxx+g");
  // gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD/AliAnalysisTaskSpectraAOD.cxx+g");
  // Use AliRoot includes to compile our task
  
  // Create and configure the alien handler plugin
  AliAnalysisGrid *alienHandler = CreateAlienHandler(mode,mc,day,month,year);  
  if (!alienHandler) return;
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);
  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  
  // Add PID task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse(mc);
  // Printf("OADB PATH:::::%s",taskPID->GetOADBPath());
  // taskPID->SetOADBPath("alien:///alice/cern.ch/user/a/akalweit/ForLeornado/OADB");
  // Printf("OADB PATH:::::%s",taskPID->GetOADBPath());

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
  AliVZEROEPSelectionTask *selTask = AddTaskVZEROEPSelection();
  
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD/AddTaskSpectraAOD.C");
  // AliAnalysisTaskSpectraAOD* taskAOD=AddTaskSpectraAOD(mc,
  // 						       0,//min Cent
  // 						       100,//CentCutMax
  // 						       0,//QvecCutMin
  // 						       100,//QvecCutMax
  // 						       -0.8,//EtaMin
  // 						       0.8,//EtaMax
  // 						       3.,//Nsigmapid
  // 						       5.,//pt
  // 						       5.,//p
  // 						       .5,//y
  // 						       .6,//ptTofMatch
  // 						       1024,//trkbit
  // 						       1,//trkbitQVector
  // 						       kFALSE,//UseCentPatchAOD049
  // 						       100000,//DCA
  // 						       70//minNclsTPC
  // 						       );
  AliAnalysisTaskSpectraAOD* taskAOD=AddTaskSpectraAOD();
  
  mgr->SetDebugLevel(2);

  //mgr->Init();
  if (!mgr->InitAnalysis())return;
  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid");
}


AliAnalysisGrid* CreateAlienHandler(TString mode="test",Int_t mc=1,Int_t day=0,Int_t month=0, Int_t year=2012)
{
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TOF -I$ALICE_ROOT/PWGLF");
  plugin->SetAdditionalLibs("libPWGLFspectra.so");
  
//plugin->SetAdditionalLibs("Histograms.h HistogramNames.h AliSpectraAODHistoManager.cxx AliSpectraAODHistoManager.h AliSpectraAODTrackCuts.cxx AliSpectraAODTrackCuts.h AliSpectraAODEventCuts.cxx AliSpectraAODEventCuts.h AliSpectraAODPID.cxx AliSpectraAODPID.h AliAnalysisTaskSpectraAOD.cxx AliAnalysisTaskSpectraAOD.h");
  //plugin->SetAnalysisSource("Histograms.h HistogramNames.h HistogramNames.cxx AliSpectraAODHistoManager.cxx AliSpectraAODTrackCuts.cxx AliSpectraAODEventCuts.cxx AliSpectraAODPID.cxx AliAnalysisTaskSpectraAOD.cxx");
  plugin->SetOverwriteMode();
  plugin->SetExecutableCommand("aliroot -q -b");  
  plugin->SetRunMode(mode.Data());
  plugin->SetNtestFiles(1);
  //Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-33-02b");
  plugin->SetAliROOTVersion("v5-03-28-AN");
  // Declare input data to be processed.
  
  if(mc)
    {
      plugin->SetGridDataDir("/alice/sim/LHC11a10a_bis");
      plugin->SetDataPattern("AOD090/*AliAOD.root");
      plugin->SetRunPrefix(""); 
      plugin->SetAnalysisMacro("TaskAODPbPbMC.C");
      plugin->SetExecutable("TaskAODPbPbMC.sh");
      plugin->SetJDLName("TaskAODPbPbMC.jdl");
      plugin->SetGridWorkingDir(Form("/AODPbPb%d%d%d/mc/",day,month,year));
    }  
  else
    {
      plugin->SetGridDataDir("/alice/data/2010/LHC10h");
      plugin->SetDataPattern("ESDs/pass2/AOD086/*AliAOD.root");
      plugin->SetRunPrefix("000"); 
      plugin->SetAnalysisMacro("TaskAODPbPbdata.C");
      plugin->SetExecutable("TaskAODPbPbdata.sh");
      plugin->SetJDLName("TaskAODPbPbdata.jdl");
      //plugin->SetSplitMaxInputFileNumber(500);
      plugin->SetGridWorkingDir(Form("/AODPbPb%d%d%d/data/",day,month,year));
    }   
  FILE* listruns=fopen("RunListGrid-AOD086.txt","r");
  
  Int_t irun;
  while(!feof(listruns))
    {
      fscanf(listruns,"%d\n",&irun);
      plugin->AddRunNumber(irun);
    }
  // Declare alien output directory. Relative to working directory.
  //plugin->SetGridOutputDir("/alice/cern.ch/user/l/lmilano/AODAnalysis/AOD086TrackBit10/mc1/output/000/Stage_1"); // In this case will be $HOME/work/output
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  //plugin->SetNrunsPerMaster(3); // 
  
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  //plugin->SetDefaultOutputs(0);
  //plugin->SetOutputFiles("AnalysisResults.root.root");
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeStages(2);
  
  plugin->SetTTL(100000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  //plugin->SetSplitMaxInputFileNumber();
  plugin->SetSplitMode("se");
  return plugin;
}
