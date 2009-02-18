void RunLocalFMDAnalysisSE(const Char_t* filename= "AliESDs.root", const Char_t* cdbPath="local://$ALICE_ROOT/OCDB", const Char_t* outFile = "fmd_analysis.root") {

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libFMDanalysis");
  
  //TGrid::Connect("alien://",0,0,"t"); 
 
  //TChain* chain = CreateChainSingle(collectionName);  
  //TFile::Open(filename);
  TChain* chain = new TChain("esdTree");
  chain->Add(filename);
  
  ///chain->Add("/home/canute/ALICE/Simulations/TestOfAnalysis/AliESDs.root");
  
  //(TChain*)gFile->Get("esdTree");
  if (!chain) return;
  
  //
  // Make the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "A test setup for the analysis train");
  // ESD input handler
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdHandler);
  
  AliMCEventHandler *mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  
  AliAODHandler* aodHandler   = new AliAODHandler();
  mgr->SetOutputEventHandler(aodHandler);
  aodHandler->SetOutputFileName("AliAODs.root");
  
  
  AliAnalysisDataContainer *cin_esd = mgr->GetCommonInputContainer();
  // Output AOD container. Pointer to AOD put here only by filter task.
  // This container is managed by the AOD handler
  AliAnalysisDataContainer *cout_aod = mgr->GetCommonOutputContainer();
  
  AliFMDAnalysisTaskSE *fmdana = new AliFMDAnalysisTaskSE("FMDAnalysis");
  mgr->AddTask(fmdana);
  // Output histograms list for jet analysis                       
  AliAnalysisDataContainer *cout_fmd = mgr->CreateContainer("BackgroundCorrected", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer, outFile);
  // Dummy AOD output container for jet analysis (no client yet)
  AliAnalysisDataContainer *c_aodfmd = mgr->CreateContainer("cAODfmd", 
							    TTree::Class(),
							    AliAnalysisManager::kExchangeContainer);
  // Connect to data containers
  mgr->ConnectInput  (fmdana,     0, cin_esd  );
  mgr->ConnectOutput (fmdana,     0, c_aodfmd );
  mgr->ConnectOutput (fmdana,     1, cout_fmd );
       
  
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(cdbPath);
  cdb->SetSpecificStorage("FMD/*","local://$ALICE_ROOT");
  cdb->SetRun(0);
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  pars->Init();
  if (AliGeomManager::GetGeometry() == NULL)
    AliGeomManager::LoadGeometry();
  
  AliFMDGeometry* geo = AliFMDGeometry::Instance();
  geo->Init();
  geo->InitTransformations();
  
  TStopwatch timer;
  timer.Start();
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
  }   
  timer.Stop();
  timer.Print();
}
//______________________________________________________________________________
TChain* CreateChainSingle(const char* xmlfile, const char *treeName="esdTree")
{
   printf("*******************************\n");
   printf("*** Getting the ESD Chain   ***\n");
   printf("*******************************\n");
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
