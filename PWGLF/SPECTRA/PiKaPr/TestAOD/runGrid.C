class  AliAnalysisManager;
class  AliAnalysisAlien;

void runGrid(TString mode="test",Bool_t mc=1,Int_t day=15,Int_t month=6, Int_t year=2012) 
{
  //to be used with Aliroot > v5-03-32-AN
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
  Double_t CentCutMin[4]= {     0,  30,  30,  30};
  Double_t CentCutMax[4]= {     5,  40,  40,  40};
  Double_t QvecCutMin[4]={      0,   0,   0, 1.5};
  Double_t QvecCutMax[4]={    100, 100, 0.4, 100};
  Double_t EtaMin[4]={       -0.8,-0.8,-0.8,-0.8};
  Double_t EtaMax[4]={        0.8, 0.8, 0.8, 0.8};
  Double_t Nsigmapid=3.;
  Double_t pt=10.;
  Double_t p=10.;
  Double_t y=.5;
  Double_t ptTofMatch=.6;
  UInt_t trkbit=1024;
  UInt_t trkbitQVector=1;
  Bool_t UseCentPatchAOD049=kFALSE;
  Double_t DCA=100000;
  UInt_t minNclsTPC=70;
  Int_t nrebin=0;
  TString opt="";
  
  for(Int_t icut=0;icut<4;icut++){
    //if(icut!=0)continue;
    AliAnalysisTaskSpectraAOD *taskAOD =AddTaskSpectraAOD(mc,CentCutMin[icut],CentCutMax[icut],QvecCutMin[icut],QvecCutMax[icut],EtaMin[icut],EtaMax[icut],Nsigmapid,pt,p,y,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,nrebin,opt);
    taskAOD->GetOutputSlot(1)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(1)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(2)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(2)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(3)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(3)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(4)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(4)->GetContainer()->GetName(),taskAOD->GetName()));
   
    // __R_ADDTASK__->GetOutputSlot(1)->GetContainer()->SetName(Form("%s_%s",__R_ADDTASK__->GetOutputSlot(1)->GetContainer()->GetName(),__R_ADDTASK__->GetName()));
    // __R_ADDTASK__->GetOutputSlot(2)->GetContainer()->SetName(Form("%s_%s",__R_ADDTASK__->GetOutputSlot(2)->GetContainer()->GetName(),__R_ADDTASK__->GetName()));
    // __R_ADDTASK__->GetOutputSlot(3)->GetContainer()->SetName(Form("%s_%s",__R_ADDTASK__->GetOutputSlot(3)->GetContainer()->GetName(),__R_ADDTASK__->GetName()));
    // __R_ADDTASK__->GetOutputSlot(4)->GetContainer()->SetName(Form("%s_%s",__R_ADDTASK__->GetOutputSlot(4)->GetContainer()->GetName(),__R_ADDTASK__->GetName()));
    
    // taskAOD->GetEventCuts()->SetMultiplicityCut(800,2000);
    // taskAOD->GetEventCuts()->SetVertexCut(-8,10);
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
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TOF -I$ALICE_ROOT/PWGLF");
  plugin->SetAdditionalLibs("libPWGLFspectra.so");
  plugin->SetOverwriteMode();
  plugin->SetExecutableCommand("aliroot -q -b");  
  plugin->SetRunMode(mode.Data());
  plugin->SetNtestFiles(1);
  //Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-33-02b");
  plugin->SetAliROOTVersion("v5-03-32-AN");
  // Declare input data to be processed.
  if(mc)
    {
      plugin->SetGridDataDir("/alice/sim/LHC11a10a_bis");
      plugin->SetDataPattern("AOD090/*AliAOD.root");
      plugin->SetRunPrefix(""); 
      plugin->SetAnalysisMacro("TaskAODPbPbMC.C");
      plugin->SetExecutable("TaskAODPbPbMC.sh");
      plugin->SetJDLName("TaskAODPbPbMC.jdl");
      //plugin->SetSplitMaxInputFileNumber(500);
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
