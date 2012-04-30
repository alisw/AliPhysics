class  AliAnalysisManager;
class  AliAnalysisAlien;

void runGrid(TString mode="test",Int_t mc=0,Int_t sub=0,Int_t hi=1,TString fname="AODAnalysis") 
{
  AliLog::SetGlobalDebugLevel(100);
  // Load common libraries
   gEnv->SetValue("XSec.GSI.DelegProxy", "2");
  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libOADB.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gROOT->LoadMacro("AliSpectraAODTrackCuts.cxx+g");
  gROOT->LoadMacro("AliSpectraAODEventCuts.cxx+g");
  gROOT->LoadMacro("AliSpectraAODHistoManager.cxx+g");
  gROOT->LoadMacro("AliAnalysisTaskSpectraAOD.cxx+g");
  // Use AliRoot includes to compile our task
  // Create and configure the alien handler plugin
  
  AliAnalysisGrid *alienHandler = CreateAlienHandler(mode,mc,sub,fname);  
  if (!alienHandler) return;
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);
  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  
  gROOT->LoadMacro("AliSpectraAODTrackCuts.cxx+g");
  gROOT->LoadMacro("AliSpectraAODEventCuts.cxx+g");
  gROOT->LoadMacro("AliSpectraAODHistoManager.cxx+g");
  gROOT->LoadMacro("AliAnalysisTaskSpectraAOD.cxx+g");

  // Physics selection
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  //  AliPhysicsSelectionTask * physicsSelectionTask = AddTaskPhysicsSelection(isMC,kTRUE,0);
  AliPhysicsSelectionTask * physicsSelectionTask = AddTaskPhysicsSelection(mc,kTRUE,0);
  //  AliPhysicsSelectionTask * physicsSelectionTask = AddTaskPhysicsSelection(isMC,kFALSE,0);
  if(!physicsSelectionTask ) { Printf("no physSelTask"); return; }
  
  // Add PID task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  Bool_t isMC = kFALSE;
  AliAnalysisTask * taskPID = AddTaskPIDResponse(mc);
  mgr->AddTask(taskPID);
  
  // gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  // AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
  // if(hi)
  //   {
  //     gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  //     AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
  //   }	
  // if(mc)
  //   {
  //     AliMCEventHandler *mch = new AliMCEventHandler();
  //     mch->Init("");
  //     mgr->SetMCtruthEventHandler(mch);
  //     //     //physSelTask->GetPhysicsSelection()->SetAnalyzeMC(); 
  //     //     //taskCentrality->SetMCInput(); 
  //   }
  // AliAODInputHandler *inputHandler = new AliAODInputHandler("handler","handler for D2H");
  // inputHandler->AddFriend("./AliAOD.root");
  // mgr->SetInputEventHandler(inputHandler);

  //setting the analysis
  Int_t iCut=0;
  Double_t CentCut[2]={0,5};
  Double_t qVecCut[2]={0,100};
  
  AliAnalysisTaskSpectraAOD *task = new AliAnalysisTaskSpectraAOD("TaskAODExercise");
  mgr->AddTask(task);
  //physics selection
  task->SelectCollisionCandidates();     
  // Set the cuts
  AliSpectraAODEventCuts * vcuts = new AliSpectraAODEventCuts("Event Cuts");
  AliSpectraAODTrackCuts  * tcuts = new AliSpectraAODTrackCuts("Track Cuts");
  //tcuts->SetTrackType(5); //AOD 046 & 047. Standard Cut with loose DCA
  tcuts->SetTrackType(6); //AOD 046 & 047. Standard Cut with tight DCA
  tcuts->SetEta(.8);
  tcuts->SetPt(5);
  tcuts->SetPtTOFMatching(0.6);   
  tcuts->SetQvecMin(qVecCut[0]);   
  tcuts->SetQvecMax(qVecCut[1]);    
  vcuts->SetCentralityCutMax(CentCut[0]);  
  vcuts->SetCentralityCutMin(CentCut[1]);
  task->SetEventCuts(vcuts);
  task->SetTrackCuts(tcuts);
  task->SetNSigmaForIdentification(5.); // FIXME
  task->SetYCut(.5);
  vcuts->PrintCuts();
  tcuts->PrintCuts();
  if (mc)
    {
      task->SetIsMC(kTRUE);
      AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
      AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer(Form("chistpt%d",iCut), AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer, 
								  Form("Pt.AOD.1._MC_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCut[0],CentCut[1],qVecCut[0],qVecCut[1]));
      AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer(Form("cvcutpt%d",iCut), AliSpectraAODEventCuts::Class(),    AliAnalysisManager::kOutputContainer, 
								  Form("Pt.AOD.1._MC_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCut[0],CentCut[1],qVecCut[0],qVecCut[1]));
      AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer(Form("ctcutpt%d",iCut), AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, 
								  Form("Pt.AOD.1._MC_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCut[0],CentCut[1],qVecCut[0],qVecCut[1]));
    }else{	 
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer(Form("chistpt%d",iCut), AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer, 
								Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCut[0],CentCut[1],qVecCut[0],qVecCut[1]));
    AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer(Form("cvcutpt%d",iCut), AliSpectraAODEventCuts::Class(),    AliAnalysisManager::kOutputContainer, 
								Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCut[0],CentCut[1],qVecCut[0],qVecCut[1]));
    AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer(Form("ctcutpt%d",iCut), AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, 
								Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCut[0],CentCut[1],qVecCut[0],qVecCut[1]));
  }
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputpt1);
  mgr->ConnectOutput(task, 2, coutputpt2);
  mgr->ConnectOutput(task, 3, coutputpt3);
  
  mgr->SetDebugLevel(2);
  
  
  


  

  //mgr->Init();
  if (!mgr->InitAnalysis())return;
  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid");
}


AliAnalysisGrid* CreateAlienHandler(TString mode="test",Int_t mc=0,Int_t sub=0,TString fname){
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();
  plugin->SetExecutableCommand("aliroot -q -b");  
  plugin->SetRunMode(mode.Data());
  plugin->SetNtestFiles(1);
  //Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-33-02b");
  plugin->SetAliROOTVersion("v5-03-18-AN");
  // Declare input data to be processed.
  if(mc)
    {
      plugin->SetGridDataDir("/alice/sim/LHC11a10a");
      plugin->SetDataPattern("AOD048/*AliAOD.root");
      plugin->SetFriendChainName("./AliAOD.root");
      plugin->SetRunPrefix(""); 
      plugin->SetAnalysisMacro(Form("TaskAODPbPbMC%d.C",sub));
      plugin->SetExecutable(Form("TaskAODPbPbMC%d.sh",sub));
      plugin->SetJDLName(Form("TaskAODPbPbMC%d.jdl",sub));
      //plugin->SetSplitMaxInputFileNumber(500);
      plugin->SetGridWorkingDir(Form("%s/mc%d/",fname.Data(),sub));
    }  
  else
    {
      plugin->SetGridDataDir("/alice/data/2010/LHC10h");
      plugin->SetDataPattern("ESDs/pass2/AOD049/*AliAOD.root");
      plugin->SetFriendChainName("./AliAOD.root");
      plugin->SetRunPrefix("000"); 
      plugin->SetAnalysisMacro(Form("TaskAODPbPbMC%d.C",sub));
      plugin->SetExecutable(Form("TaskAODPbPbMC%d.sh",sub));
      plugin->SetJDLName(Form("TaskAODPbPb%d.jdl",sub));
      //plugin->SetSplitMaxInputFileNumber(500);
      plugin->SetGridWorkingDir(Form("%s/data%d/",fname.Data(),sub));
    }  
  
  Int_t irun;
  FILE* listruns=fopen("RunListGrid.txt","r");
  while(!feof(listruns)){
    fscanf(listruns,"%d\n",&irun);
    plugin->AddRunNumber(irun);
  }
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  plugin->SetNrunsPerMaster(60); // 
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TOF");
  plugin->SetAnalysisSource("Histograms.h HistogramNames.h AliSpectraAODHistoManager.cxx+ AliSpectraAODEventCuts.cxx+ AliSpectraAODTrackCuts.cxx+ AliAnalysisTaskSpectraAOD.cxx+");
  plugin->SetAdditionalLibs("Histograms.h HistogramNames.h AliSpectraAODHistoManager.cxx AliSpectraAODHistoManager.h AliSpectraAODEventCuts.cxx AliSpectraAODEventCuts.h AliSpectraAODTrackCuts.cxx AliSpectraAODTrackCuts.h AliAnalysisTaskSpectraAOD.cxx AliAnalysisTaskSpectraAOD.h");
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
  delete runs;
}
