 /********************************************
   Macro to launch TOF QA task on ESDs data

Author: fbellini@cern.ch
Last update: 18 october 2012
*********************************************/
class AliAnalysisGrid;
TString analysisMode = "grid"; // "local" or "grid" - needs to be "grid" for plugin test mode
Bool_t useAlienPlugin=kTRUE; //use kFALSE for local analysis

//input dataset
Int_t runList[2]={186084,186162};
Int_t runNmin=0;
Int_t runNmax=2;
Int_t year = 2012;
TString prod = "LHC12d";
TString myRecPass="cpass1";
TString myQAfileSuffix="_Barrel";
Bool_t isMC = kFALSE;
TString gridUser="fbellini";
Int_t gridNtestFiles = 1;

TString prefix="_";

//do not change! to be set by SetupIO()
TString myGridDataDir="";
TString myDataPattern="";
TString myWorkDir;
TString myOutDir;
TString myJDLname;
TString myExecutableName;
TString myMacroName;

//----------------------------------------------------------------------
void SetupIO(TString filesPrefix = "")
{
  //Setup I/O paths and file names
  myGridDataDir = "/alice";

  myWorkDir = "QA_pp2012";
  myOutDir = Form("%s", prod.Data()) ;
  myJDLname = Form("job_Run%s", myWorkDir.Data());
  myExecutableName = Form("Run%s", myWorkDir.Data());
  myMacroName = Form("Run%s", myWorkDir.Data());

  if (isMC) myDataPattern="*AliESDs.root";
  else myDataPattern=Form("%s/*AliESDs%s.root",myRecPass.Data(),myQAfileSuffix.Data());
  myOutDir.Append(Form("_%s",myRecPass.Data()));
  myJDLname.Append(Form("_%s",myRecPass.Data()));
  myExecutableName.Append(Form("_%s",myRecPass.Data()));
  myMacroName.Append(Form("_%s",myRecPass.Data()));

  if (isMC) {
    myGridDataDir.Append("/sim");
    myWorkDir.Append("_MC");
    myJDLname.Append("_MC");
    myExecutableName.Append("_MC");
    myMacroName.Append("_MC");
  } else {
    myGridDataDir.Append(Form("/data/%i",year));
  }
  myGridDataDir.Append(Form("/%s",prod.Data()));
  myOutDir.Append(filesPrefix.Data());
  myExecutableName.Append(Form("%s.sh",filesPrefix.Data()));
  myMacroName.Append(Form("%s.C",filesPrefix.Data()));
  myJDLname.Append(Form("_run%i-%i%s.jdl",runNmin,runNmax,filesPrefix.Data()));
  
  Printf("=========================================  Setup I/O:");
  Printf("myGridDataDir = %s", myGridDataDir.Data());
  Printf("myDataPattern = %s", myDataPattern.Data());
  Printf("myWorkDir = %s", myWorkDir.Data());
  Printf("myOutDir = %s", myOutDir.Data());
  Printf("myJDLname = %s", myJDLname.Data());
  Printf("myExecutableName = %s", myExecutableName.Data());
  Printf("myMacroName = %s", myMacroName.Data());
  Printf("=======================================================\n");
}
//----------------------------------------------------------------------
void SetGridUser(TString username){
  if (username)
    gridUser = username.Data();
  return;
}
//----------------------------------------------------------------------
void SetGridNtestFiles(Int_t nfiles = 1){
  if (nfiles<1) gridNtestFiles=1;
  else gridNtestFiles = nfiles;
  return;
}
//----------------------------------------------------------------------
void LoadLibraries()
{
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/PWGPP -I$ALICE_ROOT/PWGPP/TRD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTENDER");
  gSystem->Load("libPWGPP.so");      
}
//----------------------------------------------------------------------
void RunAnalysisTOFqaGrid(TString pluginmode="test", Int_t ntestfiles = 10, TString filesPrefix = "", TString gridUser="fbellini") 
{
  LoadLibraries();  
  SetGridUser(gridUser.Data());
  SetGridNtestFiles(ntestfiles);
  SetupIO(filesPrefix.Data());
  if(analysisMode=="grid") {
    TGrid::Connect("alien://");
  } 

  Long64_t nentries=100000, firstentry=0;
  Bool_t readMC = kFALSE;
 
  if(useAlienPlugin) {  
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode);
    if(!alienHandler) return;
  }
  
  // Prepare input chain
  TChain *chainESD = 0;
  if(!useAlienPlugin) {
    chainESD=new TChain("esdTree");
    chainESD->Add("AliESDs.root"); //used to test locally - modify with local path files
  }
  
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(1);
  if(useAlienPlugin) mgr->SetGridHandler(alienHandler);
  AliInputEventHandler *esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  // Analysis tasks (wagons of the train)   
  //
  TString taskName;
  
  //Wagon for physics event selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
  AliPhysicsSelection* physSel = physSelTask->GetPhysicsSelection();
  if (isMC)physSelTask->GetPhysicsSelection()->SetAnalyzeMC();
  //physSel->AddBackgroundIdentification(new AliBackgroundSelection());
   
  gROOT->LoadMacro("$ALICE_ROOT/PWGPP/TOF/AddTaskTOFQA.C");
  AliAnalysisTaskTOFqa *TOFqa = (AliAnalysisTaskTOFqa*) AddTaskTOFQA(kTRUE);
  
  if(readMC) {
    AliMCEventHandler  *mcH = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcH); 
  }
  //
  // Run the analysis
  //    
  if(chainESD) printf("CHAIN HAS %d ENTRIES\n",(Int_t)chainESD->GetEntries());
  
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(analysisMode=="grid" && !useAlienPlugin) analysisMode="local";
  mgr->StartAnalysis(analysisMode.Data(),chainESD,nentries,firstentry);
  
  return;
}

//_____________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(TString pluginmode="full")
{

  AliAnalysisAlien *plugin = new AliAnalysisAlien();  
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(pluginmode.Data());
  plugin->SetUser(gridUser.Data()); //MODIFICA
  plugin->SetNtestFiles(gridNtestFiles);
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-02"); //MODIFICA
  plugin->SetAliROOTVersion("v5-03-70-AN"); //MODIFICA
  //Set user grid output dir
  plugin->SetGridWorkingDir(myWorkDir.Data()); 
  plugin->SetGridOutputDir(myOutDir.Data()); 
  plugin->SetCheckCopy(kTRUE);
  plugin->SetGridDataDir(myGridDataDir.Data());
  plugin->SetDataPattern(myDataPattern.Data());
  if (!isMC) plugin->SetRunPrefix("000");
  for (Int_t irun=runNmin;irun<runNmax;irun++){
    plugin->AddRunNumber((Int_t )runList[irun]);
  }
  plugin->SetOutputToRunNo(1);
  plugin->SetNrunsPerMaster(1);
  plugin->SetExecutableCommand("aliroot -b -q");
  plugin->SetOutputToRunNo(1);
  plugin->SetNrunsPerMaster(1);

  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/PWGPP -I$ALICE_ROOT/PWGPP/TRD");   
  plugin->SetAdditionalLibs("libANALYSIS.so libANALYSISalice.so libCORRFW.so libTENDER.so libPWGPP.so ");//libTRDbase.so libTRDrec.so
  plugin->SetDefaultOutputs(kTRUE);
  plugin->SetAnalysisMacro("AnalysisTOFqaSigned.C"); //MODIFICA se vuoi
  plugin->SetExecutable("analysisTOFqaSigned.sh"); //MODIFICA se vuoi
  plugin->SetSplitMaxInputFileNumber(50);
  plugin->SetMaxInitFailed(15);
  plugin->SetTTL(80000);
  plugin->SetInputFormat("xml-single");
  plugin->SetJDLName("TaskAnalysisTOFqaSigned.jdl"); //MODIFICA se vuoi
  plugin->SetSplitMode("se");
  return plugin;
}

//-----------------------------------------------------------------------------
TChain *CreateESDChain(TString esdpath=".",Int_t ifirst=-1,Int_t ilast=-1) 
{
  TChain *chainESD = new TChain("esdTree");
  if(ifirst<0) {
    chainESD->Add("AliESDs.root");
  } else {
    for(Int_t i=ifirst; i<=ilast; i++) {
      TString esdfile=esdpath; esdfile+=i; esdfile.Append("/AliESDs.root");
      chainESD->Add(esdfile.Data());
    }
  }
  
  return chainESD;
}
