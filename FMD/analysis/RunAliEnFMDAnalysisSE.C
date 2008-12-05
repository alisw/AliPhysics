void RunAliEnFMDAnalysisSE(const Char_t* collectionName="collection.xml", const Char_t* cdbPath="local://$ALICE_ROOT") {

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libFMDanalysis");
  
  TGrid::Connect("alien://",0,0,"t"); 
  TChain* chain = CreateChainSingle(collectionName);  
  
  if (!chain) return;
  
  //
  // Make the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "A test setup for the analysis train");
  // ESD input handler
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdHandler);

  AliAODHandler* aodHandler   = new AliAODHandler();
  mgr->SetOutputEventHandler(aodHandler);
  aodHandler->SetOutputFileName("AliAODs.root");
  
  
  AliAnalysisDataContainer *cin_esd = mgr->CreateContainer("cESD",TChain::Class(), 
							   AliAnalysisManager::kInputContainer);
  // Output AOD container. Pointer to AOD put here only by filter task.
  // This container is managed by the AOD handler
  AliAnalysisDataContainer *cout_aod = mgr->CreateContainer("cAOD", TTree::Class(),
							    AliAnalysisManager::kOutputContainer, "default");
  
  AliFMDAnalysisTaskSE *fmdana = new AliFMDAnalysisTaskSE("FMDAnalysis");
  mgr->AddTask(fmdana);
  // Output histograms list for jet analysis                       
  AliAnalysisDataContainer *cout_fmd = mgr->CreateContainer("BackgroundCorrected", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer, "fmd_analysis.root");
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
    mgr->StartAnalysis("local",chain, 1000);
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
