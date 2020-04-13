TChain* CreateChainXML(const char *xmlfile);
TChain* CreateChainTXT(const char *txtfile);
Bool_t InputHandlerSetup(TString formatFr = "AliESDfriends.root");

Bool_t barrelFlag = kFALSE;


//====================================================================
void runGloAlgTask
(
 TString data="wn.xml",
 // TString data="algColl.txt",
 Int_t nEvents=-1,
 UInt_t trigSel = AliVEvent::kAny
 )
{
  //
  if (!gGrid) {
    TGrid::Connect("alien://");
    if (!gGrid) {printf("Cannot connect\n");exit(1);}
  }
  //
  gROOT->ProcessLine(".include $ALICE_ROOT/include $ALICE_PHYSICS/include ./");
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I./");
  
  if (gClassTable->GetID("AliAlgSteer")<0) {
    gSystem->Load("libALIGN.so");
  }
  //
  if (gClassTable->GetID("AliGloAlgTask")<0) {
    gROOT->ProcessLine(".L AliGloAlgTask.cxx++g");
  }
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("mgr");
  //
  TChain *chain = data.EndsWith(".xml") ? CreateChainXML(data) : CreateChainTXT(data);
  TString formatFr = "AliESDfriends.root";
  if (barrelFlag) formatFr = "AliESDfriends_Barrel.root";
  printf("Deduced friend name : %s\n",formatFr.Data());
  //
  InputHandlerSetup(formatFr.Data());
  //
  //================================================================================
  //  printf("Requesting physics selection in %s mode\n",useMC ? "MC":"Data");
  //  gROOT->LoadMacro("$ALICE_ROOT/../src/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  //  AliPhysicsSelectionTask* physicsSelectionTask = AddTaskPhysicsSelection(useMC);
  //==================================================================================
  //
  AliGloAlgTask *algTask = new AliGloAlgTask("alg");
  //-------------------------------------------
  algTask->SetTriggerSelection(trigSel);
  algTask->SetConfMacroName("alignConf.C");
  //  algTask->SetConfMacroName("pedeF/alignConf.C");
  algTask->SetIniParFileName("millepede.res");
  //  algTask->SetApplyMPSolAlignment(kTRUE);
  //-------------------------------------------
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("clist", TList::Class(),AliAnalysisManager::kOutputContainer,"mpStatOut.root");
  mgr->AddTask(algTask);
  //
  mgr->ConnectInput(algTask, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(algTask,1,coutput1);
  //
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  // Start analysis in grid.
  if (nEvents<0) nEvents = chain->GetEntries();
  mgr->StartAnalysis("localfile",chain,nEvents);
  //
}

//________________________________________________________________
Bool_t InputHandlerSetup(TString esdFName)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cin = mgr->GetCommonInputContainer();
  if (cin) return;
  AliESDInputHandler *esdInputHandler = 
    dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdInputHandler) {
    Info("CustomAnalysisTaskInputSetup", "Creating esdInputHandler ...");
    esdInputHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdInputHandler);
  }
  esdInputHandler->SetFriendFileName(esdFName.Data());
  //  esdInputHandler->SetActiveBranches("ESDfriend*");
  esdInputHandler->SetReadFriends(kTRUE);
  return kTRUE;
}


//________________________________________________________________________________
TChain* CreateChainXML(const char *xmlfile)
{
// Create a chain using url's from xml file
   TString filename;
   Int_t run = 0;
   TString treename = "esdTree";
   printf("***************************************\n");
   printf("    Getting chain of trees %s\n", treename.Data());
   printf("***************************************\n");
   TGridCollection *coll = gGrid->OpenCollection(xmlfile);
   if (!coll) {
      ::Error("CreateChain", "Cannot create an AliEn collection from %s", xmlfile);
      return NULL;
   }
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   TChain *chain = new TChain(treename);
   coll->Reset();
   while (coll->Next()) {
      filename = coll->GetTURL();
      if (filename.EndsWith("Barrel.root")) barrelFlag = kTRUE;
      if (mgr) {
         Int_t nrun = AliAnalysisManager::GetRunFromAlienPath(filename);
         if (nrun && nrun != run) {
            printf("### Run number detected from chain: %d\n", nrun);
            mgr->SetRunFromPath(nrun);
            run = nrun;
         }
      }
      chain->Add(filename);
   }
   if (!chain->GetNtrees()) {
      ::Error("CreateChain", "No tree found from collection %s", xmlfile);
      return NULL;
   }
   printf("Created chain with %d entries in %d trees from %s\n",chain->GetEntries(),chain->GetNtrees(),xmlfile);
   return chain;
}

//________________________________________________________________________________
TChain* CreateChainTXT(const char* inpData)
{
  const char* chName="esdTree";
  TChain* chain = new TChain(chName);
  //
  TString inpDtStr = inpData;
  if (inpDtStr.EndsWith(".root")) {
    chain->AddFile(inpData);
  }
  else {
    //
    ifstream inpf(inpData);
    if (!inpf.good()) {
      printf("Failed on input filename %s\n",inpData);
      return kFALSE;
    }
    //
    TString flName;
    flName.ReadLine(inpf);
    while ( !flName.IsNull() ) {
      flName = flName.Strip(TString::kBoth,' ');
      if (flName.BeginsWith("//") || flName.BeginsWith("#")) {flName.ReadLine(inpf); continue;}
      flName = flName.Strip(TString::kBoth,',');
      flName = flName.Strip(TString::kBoth,'"');
      if (flName.EndsWith("Barrel.root")) barrelFlag = kTRUE;
      printf("Adding %s\n",flName.Data());
      chain->AddFile(flName.Data());
      flName.ReadLine(inpf);
    }
  }
  //
  int n = chain->GetEntries();
  if (n<1) {
    printf("Obtained chain is empty\n");
    return kFALSE;
  }
  printf("Created chain with %d entries in %d trees from %s\n",chain->GetEntries(),chain->GetNtrees(),inpData);
  return chain;
}
