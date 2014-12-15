class  AliAnalysisManager;
class  AliAnalysisAlien;

void runGrid(TString mode="test",Bool_t mc=0,Int_t day=19,Int_t month=7, Int_t year=2012) 
{
  //to be used with Aliroot > v5-03-32-AN
  AliLog::SetGlobalDebugLevel(100);
  // Load common libraries
  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
  
  // Load common libraries
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libGui");
  gSystem->Load("libXMLParser");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libCDB");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libProof");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libSTEER");
  //__________________________________________________________________________
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gSystem->Load("libPWGLFspectra");
  
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
  //LOOP OVER SELECTION
  //                            0    1    2    3    4    5    6    7    8    9
  Double_t Nsigmapid=3.;
  Double_t pt=10000.;
  Double_t p=10000.;
  Double_t y=.5;
  Double_t ptTofMatch=.6;
  UInt_t trkbit=1024;
  UInt_t trkbitQVector=1;
  Bool_t UseCentPatchAOD049=kFALSE;
  Double_t DCA=100000;
  UInt_t minNclsTPC=70;
  Int_t nrebin=30;
  TString opt="";
  Int_t centcut=30;
  
  for(Int_t icut=0;icut<1;icut++){
    //if(icut!=0)continue;
    //MB
    AliAnalysisTaskSpectraAOD *taskAOD =AddTaskSpectraAOD(mc,centcut,centcut+1,0,100,-0.8,0.8,Nsigmapid,pt,p,y,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,nrebin,opt);
    taskAOD->GetOutputSlot(1)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(1)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(2)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(2)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(3)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(3)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(4)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(4)->GetContainer()->GetName(),taskAOD->GetName()));
    //taskAOD->GetEventCuts()->SetMultiplicityCut(250,450);
    taskAOD->GetEventCuts()->PrintCuts();
    //smallq
    AliAnalysisTaskSpectraAOD *taskAOD =AddTaskSpectraAOD(mc,centcut,centcut+1,0,0.4,-0.8,0.8,Nsigmapid,pt,p,y,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,nrebin,opt);
    taskAOD->GetOutputSlot(1)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(1)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(2)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(2)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(3)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(3)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(4)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(4)->GetContainer()->GetName(),taskAOD->GetName()));
    //taskAOD->GetEventCuts()->SetMultiplicityCut(250,450);
    taskAOD->GetEventCuts()->PrintCuts();
    taskAOD->GetEventCuts()->SetCentFromV0(kTRUE);
    //largeq
    AliAnalysisTaskSpectraAOD *taskAOD =AddTaskSpectraAOD(mc,centcut,centcut+1,1.7,100,-0.8,0.8,Nsigmapid,pt,p,y,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,nrebin,opt);
    taskAOD->GetOutputSlot(1)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(1)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(2)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(2)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(3)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(3)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(4)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(4)->GetContainer()->GetName(),taskAOD->GetName()));
    //taskAOD->GetEventCuts()->SetMultiplicityCut(250,450);
    taskAOD->GetEventCuts()->PrintCuts();
    //increase centrality
    centcut++;
    }
  
  mgr->SetDebugLevel(2);
  
  //mgr->Init();
  if (!mgr->InitAnalysis())return;
  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid");
}


AliAnalysisGrid* CreateAlienHandler(TString mode="test",Bool_t mc=1,Int_t day=0,Int_t month=0, Int_t year=2012)
{
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TOF");
  plugin->SetAdditionalLibs("libPWGLFspectra.so");
  // gROOT->LoadMacro("AliSpectraAODTrackCuts.cxx+g"); 
  // gROOT->LoadMacro("AliSpectraAODEventCuts.cxx+g"); 
  // gROOT->LoadMacro("AliSpectraAODHistoManager.cxx+g"); 
  // gROOT->LoadMacro("AliSpectraAODPID.cxx+g"); 
  // gROOT->LoadMacro("AliAnalysisTaskSpectraAOD.cxx+g"); 
  
  //plugin->SetAnalysisSource("AliSpectraAODTrackCuts.cxx AliSpectraAODEventCuts.cxx AliSpectraAODHistoManager.cxx AliSpectraAODPID.cxx AliAnalysisTaskSpectraAOD.cxx");
  //plugin->SetAdditionalLibs("Histograms.h HistogramNames.h HistogramNames.cxx AliSpectraAODTrackCuts.cxx AliSpectraAODEventCuts.cxx AliSpectraAODHistoManager.cxx AliSpectraAODPID.cxx AliAnalysisTaskSpectraAOD.cxx AliSpectraAODTrackCuts.h AliSpectraAODEventCuts.h AliSpectraAODHistoManager.h AliSpectraAODPID.h AliAnalysisTaskSpectraAOD.h");

  plugin->SetOverwriteMode();
  plugin->SetExecutableCommand("aliroot -q -b");  
  plugin->SetRunMode(mode.Data());
  plugin->SetNtestFiles(1);
  //Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-33-02b");
  plugin->SetAliROOTVersion("v5-03-35-AN");
  // Declare input data to be processed.
  if(mc)
    {
      // plugin->SetGridDataDir("/alice/sim/LHC11a10a_bis");
      // plugin->SetDataPattern("AOD090/*AliAOD.root");
      // plugin->SetRunPrefix(""); 
      // plugin->SetAnalysisMacro("TaskAODPbPbMC.C");
      // plugin->SetExecutable("TaskAODPbPbMC.sh");
      // plugin->SetJDLName("TaskAODPbPbMC.jdl");
      // //plugin->SetSplitMaxInputFileNumber(500);
      // plugin->SetGridWorkingDir(Form("/AODPbPb%d%d%d/mc/",day,month,year));
      
      plugin->SetGridDataDir(" /alice/sim/2012/LHC12a11e/");
      plugin->SetDataPattern("AOD081/*AliAOD.root");
      plugin->SetRunPrefix(""); 
      plugin->SetAnalysisMacro("TaskAODPbPbMC.C");
      plugin->SetExecutable("TaskAODPbPbMC.sh");
      plugin->SetJDLName("TaskAODPbPbMC.jdl");
      //plugin->SetSplitMaxInputFileNumber(200);
      plugin->SetGridWorkingDir(Form("/AODPbPb%d%d%d/mcAMPTNoPhysSel_RangeXLarger/",day,month,year));
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
      plugin->SetGridWorkingDir(Form("/AODPbPb%d%d%d/dataFinalPlot3040/",day,month,year));
    }   
  FILE* listruns=fopen("RunListGrid-AOD086.txt","r");
  Int_t irun;
  while(!feof(listruns))
    {
      fscanf(listruns,"%d\n",&irun);
      plugin->AddRunNumber(irun);
    }
  // plugin->AddRunNumber(139437);
  // plugin->AddRunNumber(139038);
  // plugin->AddRunNumber(138653);
  // plugin->AddRunNumber(138534);
  // plugin->AddRunNumber(137686);
 
  // Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  //plugin->SetGridOutputDir("/alice/cern.ch/user/l/lmilano/AODPbPb472012/data/output/000/");
  plugin->SetNrunsPerMaster(24); // 
  
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  //plugin->SetDefaultOutputs(0);
  //plugin->SetOutputFiles("AnalysisResults.root.root");
  plugin->SetMergeViaJDL(1);
  plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeStages(2);
  //plugin->SetMaxMergeFiles(100);
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
