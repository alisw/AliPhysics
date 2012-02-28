//char *gridMode="";
char *gridMode = "full"; // "" for local
char *workingDir = "takuv2c123456_rev3_2_mc";
TString commonOutputFileName = Form("Results%s",workingDir);

class AliAnalysisAlien;
class AliAnalysisGrid;
void runTakuMC() {


  gSystem->Setenv("alien_CLOSE_SE", "ALICE::GSI::SE");
  cout << "alien_CLOSE_SE: " << gSystem->Getenv("alien_CLOSE_SE") << endl;
  TGrid::Connect("alien://");


  TStopwatch timer; timer.Start();
  LoadLibraries();
  if(gridMode!="")
    AliAnalysisGrid *alienHandler = CreateAlienHandler();
  else {
    TChain *chain = new TChain("esdTree");
    for(int i=1; i!=10; ++i)
      chain->Add( Form("/home/gunji/softwares/dielectron/data/esd138740_mc/%d0/AliESDs.root",i) );
  }
  AliAnalysisManager *mgr = new AliAnalysisManager("DielectronAnalysisManager");
  AliESDInputHandler *esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);
  if(gridMode!="")
    mgr->SetGridHandler(alienHandler);


  AliMCEventHandler *MC = new AliMCEventHandler;
  mgr->SetMCtruthEventHandler(MC);

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection(kTRUE);
  physicsSelTask->GetPhysicsSelection()->SetAnalyzeMC();


  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality =AddTaskCentrality();
  ///taskCentrality->SetPass(2);
  taskCentrality->SetMCInput();

  gROOT->LoadMacro("AliDielectronDebugTreeTaku.cxx++");
  gROOT->LoadMacro("AliDielectronHistosTaku.cxx++");
  gROOT->LoadMacro("AliDielectronTaku.cxx++");
  gROOT->LoadMacro("AliAnalysisTaskMultiDielectronNewTaku.cxx++");

  LoadAddLibraries();
  
  TString configFile("./ConfigJpsi2eeDataTaku.C");
  //Bool_t hasMC=(mgr->GetMCtruthEventHandler()!=0x0);
  //if (hasMC){
  //configFile="$ALICE_ROOT/PWG3/dielectron/macros/ConfigJpsi2eeEff.C";
  //}
  gROOT->LoadMacro(configFile.Data());

  gROOT->LoadMacro("AddTaskDielectronTaku.C");
  

  
  AddTaskDielectronTaku( 0, 99, commonOutputFileName.Data(), "CENT1" );
  /*
  AddTaskDielectronTaku( 10, 20, commonOutputFileName.Data(), "CENT2" );
  AddTaskDielectronTaku( 20, 30, commonOutputFileName.Data(), "CENT3" );
  AddTaskDielectronTaku( 30, 40, commonOutputFileName.Data(), "CENT4" );
  AddTaskDielectronTaku( 40, 50, commonOutputFileName.Data(), "CENT5" );  
  AddTaskDielectronTaku( 50, 60, commonOutputFileName.Data(), "CENT6" );
  AddTaskDielectronTaku( 60, 70, commonOutputFileName.Data(), "CENT7" );
  AddTaskDielectronTaku( 70, 80, commonOutputFileName.Data(), "CENT8" );
  AddTaskDielectronTaku( 80, 99, commonOutputFileName.Data(), "CENT9" );
  */

  mgr->SetDebugLevel(2);
  mgr->SetDebugLevel(0);
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(gridMode!="")
    mgr->StartAnalysis("grid");
  else
    mgr->StartAnalysis("local",chain);
  timer.Stop();
  timer.Print();  

}

AliAnalysisGrid* CreateAlienHandler() {
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetRunMode(gridMode);
  plugin->SetNtestFiles(1);
  plugin->SetAPIVersion("V1.1x");

  plugin->SetROOTVersion("v5-28-00d");
  plugin->SetAliROOTVersion("v4-21-25-AN");
  //  plugin->SetRunPrefix("000");
  gROOT->LoadMacro("AddRunsPbPbMC.C");
  AddRunsPbPbMC(plugin);
  plugin->SetOutputToRunNo();  
  plugin->SetGridWorkingDir(workingDir);
  plugin->SetGridOutputDir("output");
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWG2/FLOW/AliFlowCommon -I$ALICE_ROOT/PWG2/FLOW/AliFlowTasks -I$ALICE_ROOT/PWG3/dielectron/ -g");
  plugin->SetAnalysisSource("AliDielectronHistosTaku.cxx AliDielectronDebugTreeTaku.cxx AliDielectronTaku.cxx AliAnalysisTaskMultiDielectronNewTaku.cxx");
  //  plugin->SetAdditionalLibs("libCORRFW.so libPWG3base.so libPWG3dielectron.so libPWG3hfe.so libTENDER.so libTENDERSupplies.so  AliAnalysisTaskMultiDielectronNew.h AliAnalysisTaskMultiDielectronNew.cxx");
  plugin->SetAdditionalLibs("libCORRFW.so libPWG3base.so libPWG3dielectron.so libPWG3hfe.so libPWG2flowCommon.so libPWG2flowTasks.so AliDielectronHistosTaku.h AliDielectronHistosTaku.cxx AliDielectronDebugTreeTaku.h AliDielectronDebugTreeTaku.cxx AliDielectronTaku.h AliDielectronTaku.cxx AliAnalysisTaskMultiDielectronNewTaku.h AliAnalysisTaskMultiDielectronNewTaku.cxx");
  plugin->SetDefaultOutputs(kTRUE);
  plugin->SetAnalysisMacro(Form("%s.C",workingDir));
  plugin->SetExecutable(Form("%s.sh",workingDir));
  plugin->SetSplitMaxInputFileNumber(100);
  plugin->SetOverwriteMode(kTRUE);
  plugin->SetMaxInitFailed(20);
  plugin->SetTTL(90000);
  plugin->SetOutputToRunNo(kTRUE);
  plugin->SetMasterResubmitThreshold(90);
  plugin->SetInputFormat("xml-single");
  plugin->SetJDLName(Form("%s.jdl",workingDir));
  plugin->SetPrice(1);      
  plugin->SetSplitMode("se");
  plugin->SetKeepLogs(kTRUE);
  plugin->SetExecutableCommand("aliroot -b -q ");
  return plugin;
}

void LoadLibraries()
{
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/build/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWG2/FLOW/AliFlowCommon -I$ALICE_ROOT/PWG2/FLOW/AliFlowTasks -I$ALICE_ROOT/PWG3/dielectron/ -g");
  gSystem->Load("libCore");// no
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");// no
  gSystem->Load("libPhysics");
  gSystem->Load("libXMLParser");
  gSystem->Load("libProof");
  gSystem->Load("libMinuit");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libSTEER");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTOFbase");
  gSystem->Load("libTOFrec");
  gSystem->Load("libT0base");
  gSystem->Load("libT0rec");
  gSystem->Load("libPWG2flowCommon");
  gSystem->Load("libPWG2flowTasks");

  gSystem->Load("libTENDER");
  gSystem->Load("libTENDERSupplies");



  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG3base.so");
  gSystem->Load("libPWG3dielectron.so");
  gSystem->Load("libPWG3hfe.so");

}

void LoadAddLibraries(){
  gSystem->Load("./AliDielectronHistosTaku_cxx.so");
  gSystem->Load("./AliDielectronDebugTreeTaku_cxx.so");
  gSystem->Load("./AliDielectronTaku_cxx.so");


}
