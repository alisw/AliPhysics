void runBatch(int Nevents=20000000, const char* outfilename="AnalysisResults.root",  bool batchmode=kTRUE, const char* collectionfile="collection.xml") {
  TStopwatch timer;
  timer.Start();

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");

  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  //____________________________________________________//
  //_____________Setting up required packages___________//
  //____________________________________________________//
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  //ANALYSIS PART
  // gROOT->LoadMacro("AliAnalysisCheck.cxx+g");

gROOT->LoadMacro("AliAnalysisTaskParticleEfficiency.cxx+g");

  //gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  


  TChain *chain = new TChain("aodTree");

  const char *collectionfile="wn.xml";
  ifstream *istr = new ifstream(collectionfile);

  char fname[2000];
  char pname[2000];
  while (!istr->eof()) {
    fname[0] = '\0';
    (*istr) >> fname;
    if (strlen(fname) > 10) {
      sprintf(pname, "alien://%s", fname);
      chain->Add(pname);
    }
    }


  //chain->Add("../../data/AOD/PbPb/AOD95/AliAOD.root");
  //chain->Add("/opt/alice/workdir/TestConfig/data/Pythia/AOD/1/AliAOD.root");
  //chain->Add("/opt/alice/workdir/TestConfig/data/AOD/PbPb/PbPb.LHC10h/AOD86/1/AliAOD.root");
   //chain->Add("/opt/alice/workdir/TestConfig/data/AOD/PbPb/PbPb.LHC10h/AOD86/2/AliAOD.root");
   // chain->Add("/opt/alice/workdir/TestConfig/data/AOD/PbPb/PbPb.LHC10h/AOD86/3/AliAOD.root");
   //chain->Add("/opt/alice/workdir/TestConfig/data/AOD/PbPb/PbPb.LHC10h/AOD86/4/AliAOD.root");
  // chain->Add("../../data/AOD/PbPb/New.PbPb.2012.01/AliAOD.root");
  //chain->Add("../../data/AOD/PbPb/New.PbPb.2012.01.2/AliAOD.root");
// else {
 

//   chain = CreateChainFromCollection(collectionfile, "esdTree");
    //  }

   
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

  AliAODInputHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
  //____________________________________________//
  // 1st Pt task


  //AddTaskPIDResponse
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  Bool_t isMC=kTRUE; Bool_t tuneOnData = kTRUE; // kTRUE in case of MC
  AddTaskPIDResponse(isMC,kTRUE,tuneOnData); 

  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  //AddTaskPIDqa();

  AliAnalysisTaskParticleEfficiency *myTask = new AliAnalysisTaskParticleEfficiency("MyTask");
  if(!myTask) exit(-1);
  mgr->AddTask(myTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("MyList", TList::Class(),AliAnalysisManager::kOutputContainer,outfilename);
 
  //____________________________________________//
  mgr->ConnectInput(myTask,0,cinput);	
  mgr->ConnectOutput(myTask,1,coutput2);	

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain, Nevents);

  timer.Stop();
  timer.Print();
}

//______________________________________________________________________________
TChain *CreateChainFromCollection(const char* xmlfile, const char *treeName="esdTree")
{
// Create a chain from an alien collection.
   TAlienCollection * myCollection  = TAlienCollection::Open(xmlfile);

   if (!myCollection) {
      ::Error("CreateChainSingle", "Cannot create an AliEn collection from %s", xmlfile) ;
      return NULL ;
  }

  TChain* chain = new TChain(treeName);
  myCollection->Reset() ;
  while ( myCollection->Next() ) chain->Add(myCollection->GetTURL("")) ;
  chain->ls();
  return chain;
}
