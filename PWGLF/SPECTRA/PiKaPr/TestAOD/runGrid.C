class  AliAnalysisManager;
class  AliAnalysisAlien;

void runGrid(TString mode="test",Int_t mc=1,Int_t sub=0,Int_t hi=1,TString fname="AODAnalysis_1June2012") 
{
  //0 is AOD048-049 in this case you can choos FilterBit5 (loose DCA) or 6 (tight DCA)!!!!!!!!!
  //1 is AOD086-090
  AliLog::SetGlobalDebugLevel(100);
  // Load common libraries
  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
  
  gSystem->Load("libTree.so");
  //gSystem->Load("libGeom.so");
  //gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libOADB.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");

  // gROOT->LoadMacro("$ALICE_ROOT/SPECTRA/PiKaPr/TestAOD/AliSpectraAODTrackCuts.cxx+g");
  // gROOT->LoadMacro("$ALICE_ROOT/SPECTRA/PiKaPr/TestAOD/AliSpectraAODEventCuts.cxx+g");
  // gROOT->LoadMacro("$ALICE_ROOT/SPECTRA/PiKaPr/TestAOD/AliSpectraAODHistoManager.cxx+g");
  // gROOT->LoadMacro("$ALICE_ROOT/SPECTRA/PiKaPr/TestAOD/AliSpectraAODPID.cxx+g");
  // gROOT->LoadMacro("$ALICE_ROOT/SPECTRA/PiKaPr/TestAOD/AliAnalysisTaskSpectraAOD.cxx+g");
  gROOT->LoadMacro("AliSpectraAODTrackCuts.cxx+g");
  gROOT->LoadMacro("AliSpectraAODEventCuts.cxx+g");
  gROOT->LoadMacro("AliSpectraAODHistoManager.cxx+g");
  gROOT->LoadMacro("AliSpectraAODPID.cxx+g");
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
  
  // Physics selection
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask * physicsSelectionTask = AddTaskPhysicsSelection(mc,kTRUE,0);
  if(!physicsSelectionTask ) { Printf("no physSelTask"); return; }
  
  // Add PID task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  Bool_t isMC = kFALSE;
  AliAnalysisTask * taskPID = AddTaskPIDResponse(mc);
  mgr->AddTask(taskPID);
  
  //setting the analysis
  Int_t iCut=0;
  //Double_t CentCut[2]={0,100};
  Double_t CentCut[2]={0,100};
  Double_t qVecCut[2]={0,100};
   
  //PID object
  using namespace AliSpectraNameSpace;
  AliSpectraAODPID *pid = new AliSpectraAODPID(kNSigmaTPCTOF); 
  pid->SetNSigmaCut(3.);
  
  AliAnalysisTaskSpectraAOD *task = new AliAnalysisTaskSpectraAOD("TaskAODExercise");
  mgr->AddTask(task);
  //physics selection
  task->SelectCollisionCandidates();     
  
  // Set the cuts
  AliSpectraAODEventCuts * vcuts = new AliSpectraAODEventCuts("Event Cuts");
  AliSpectraAODTrackCuts  * tcuts = new AliSpectraAODTrackCuts("Track Cuts");
  if(sub==0){
    tcuts->SetTrackType(5); //AOD 046 & 047. Standard Cuts with loose DCA
    //tcuts->SetTrackType(6); //AOD 046 & 047. Standard Cuts with tight DCA
  }
  if(sub==1)tcuts->SetTrackType(10); //AOD 086 & 090. Standard Raa cuts
  
  // set pid object
  task->SetPID(pid);
  tcuts->SetEta(.8);
  tcuts->SetPt(5);
  tcuts->SetY(.5);
  tcuts->SetPtTOFMatching(0.6);   
  tcuts->SetQvecMin(qVecCut[0]);   
  tcuts->SetQvecMax(qVecCut[1]);    
  vcuts->SetCentralityCutMin(CentCut[0]);
  vcuts->SetCentralityCutMax(CentCut[1]);  
  task->SetEventCuts(vcuts);
  task->SetTrackCuts(tcuts);
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
      AliAnalysisDataContainer *coutputpt4 = mgr->CreateContainer(Form("cpidpt%d",iCut), AliSpectraAODPID::Class(),     AliAnalysisManager::kOutputContainer, 
								  Form("Pt.AOD.1._MC_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCut[0],CentCut[1],qVecCut[0],qVecCut[1]));
    }else{	 
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer(Form("chistpt%d",iCut), AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer, 
								Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCut[0],CentCut[1],qVecCut[0],qVecCut[1]));
    AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer(Form("cvcutpt%d",iCut), AliSpectraAODEventCuts::Class(),    AliAnalysisManager::kOutputContainer, 
								Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCut[0],CentCut[1],qVecCut[0],qVecCut[1]));
    AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer(Form("ctcutpt%d",iCut), AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, 
								Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCut[0],CentCut[1],qVecCut[0],qVecCut[1]));
    AliAnalysisDataContainer *coutputpt4 = mgr->CreateContainer(Form("cpidpt%d",iCut), AliSpectraAODPID::Class(),     AliAnalysisManager::kOutputContainer, 
								Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCut[0],CentCut[1],qVecCut[0],qVecCut[1]));
  }
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputpt1);
  mgr->ConnectOutput(task, 2, coutputpt2);
  mgr->ConnectOutput(task, 3, coutputpt3);
  mgr->ConnectOutput(task, 4, coutputpt4);
  
  mgr->SetDebugLevel(2);
  
  //mgr->Init();
  if (!mgr->InitAnalysis())return;
  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid");
}


AliAnalysisGrid* CreateAlienHandler(TString mode="test",Int_t mc=0,Int_t sub=0,TString fname){
  
  
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TOF");
  plugin->SetAdditionalLibs("Histograms.h HistogramNames.h AliSpectraAODHistoManager.cxx AliSpectraAODHistoManager.h AliSpectraAODTrackCuts.cxx AliSpectraAODTrackCuts.h AliSpectraAODEventCuts.cxx AliSpectraAODEventCuts.h AliSpectraAODPID.cxx AliSpectraAODPID.h AliAnalysisTaskSpectraAOD.cxx AliAnalysisTaskSpectraAOD.h");
  plugin->SetAnalysisSource("Histograms.h HistogramNames.h AliSpectraAODHistoManager.cxx AliSpectraAODTrackCuts.cxx AliSpectraAODEventCuts.cxx AliSpectraAODPID.cxx AliAnalysisTaskSpectraAOD.cxx");
  plugin->SetOverwriteMode();
  plugin->SetExecutableCommand("aliroot -q -b");  
  plugin->SetRunMode(mode.Data());
  plugin->SetNtestFiles(1);
  //Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-33-02b");
  plugin->SetAliROOTVersion("v5-04-25-AN");
  // Declare input data to be processed.
  if(sub==0){
  if(mc)
    {
      plugin->SetGridDataDir("/alice/sim/LHC11a10a");
      plugin->SetDataPattern("AOD048/*AliAOD.root");
      plugin->SetRunPrefix(""); 
      plugin->SetAnalysisMacro(Form("TaskAOD046PbPbMC%d.C",sub));
      plugin->SetExecutable(Form("TaskAOD046PbPbMC%d.sh",sub));
      plugin->SetJDLName(Form("TaskAOD046PbPbMC%d.jdl",sub));
      //plugin->SetSplitMaxInputFileNumber(500);
      plugin->SetGridWorkingDir(Form("%s/AOD048-049Filter5/mc%d/",fname.Data(),sub));
    }  
  else
    {
      plugin->SetGridDataDir("/alice/data/2010/LHC10h");
      plugin->SetDataPattern("ESDs/pass2/AOD049/*AliAOD.root");
      plugin->SetRunPrefix("000"); 
      plugin->SetAnalysisMacro(Form("TaskAOD046PbPbdata%d.C",sub));
      plugin->SetExecutable(Form("TaskAOD046PbPbdata%d.sh",sub));
      plugin->SetJDLName(Form("TaskAOD046PbPbdata%d.jdl",sub));
      //plugin->SetSplitMaxInputFileNumber(500);
      plugin->SetGridWorkingDir(Form("%s/AOD048-049Filter5/data%d/",fname.Data(),sub));
    }  
  FILE* listruns=fopen("RunListGrid-AOD046.txt","r");
  }
  if(sub==1){
    if(mc)
      {
	plugin->SetGridDataDir("/alice/sim/LHC11a10a_bis");
	plugin->SetDataPattern("AOD090/*AliAOD.root");
	plugin->SetRunPrefix(""); 
	plugin->SetAnalysisMacro(Form("TaskAOD086PbPbMC%d.C",sub));
	plugin->SetExecutable(Form("TaskAOD086PbPbMC%d.sh",sub));
	plugin->SetJDLName(Form("TaskAOD086PbPbMC%d.jdl",sub));
	//plugin->SetSplitMaxInputFileNumber(500);
	plugin->SetGridWorkingDir(Form("%s/AOD086-090Filter10/mc%d/",fname.Data(),sub));
      }  
    else
      {
	plugin->SetGridDataDir("/alice/data/2010/LHC10h");
	plugin->SetDataPattern("ESDs/pass2/AOD086/*AliAOD.root");
	plugin->SetRunPrefix("000"); 
	plugin->SetAnalysisMacro(Form("TaskAOD086PbPbdata%d.C",sub));
	plugin->SetExecutable(Form("TaskAOD086PbPbdata%d.sh",sub));
	plugin->SetJDLName(Form("TaskAOD086PbPbdata%d.jdl",sub));
	//plugin->SetSplitMaxInputFileNumber(500);
	plugin->SetGridWorkingDir(Form("%s/AOD086-090Filter10/data%d/",fname.Data(),sub));
          }  
    FILE* listruns=fopen("RunListGrid-AOD086.txt","r");
  }
  Int_t irun;
  while(!feof(listruns)){
    fscanf(listruns,"%d\n",&irun);
    plugin->AddRunNumber(irun);
  }
  // Declare alien output directory. Relative to working directory.
  //plugin->SetGridOutputDir("/alice/cern.ch/user/l/lmilano/AODAnalysis/AOD086TrackBit10/mc1/output/000/Stage_1"); // In this case will be $HOME/work/output
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  plugin->SetNrunsPerMaster(60); // 
  

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
