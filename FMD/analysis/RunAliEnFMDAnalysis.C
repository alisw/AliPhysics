void RunAliEnFMDAnalysis(const Char_t* collectionfile = "collection.xml",
			 const Char_t* cdbPath        = "local://$ALICE_ROOT",
			 const Char_t* outFile        = "fmd_analysis.root"){
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libFMDanalysis");
  
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
  
  cout<<"Creating task for analysis"<<endl;
  AliFMDAnalysisTaskESDReader *FMDana0 = new AliFMDAnalysisTaskESDReader("reader");
  FMDana0->SetDebugLevel(10);
  AliFMDAnalysisTaskSharing *FMDana1 = new AliFMDAnalysisTaskSharing("sharing");
  AliFMDAnalysisTaskDensity *FMDana2 = new AliFMDAnalysisTaskDensity("density");
  AliFMDAnalysisTaskBackgroundCorrection *FMDana3 = new AliFMDAnalysisTaskBackgroundCorrection("background");
  AliFMDAnalysisTaskDndeta *FMDana4 = new AliFMDAnalysisTaskDndeta("dNdeta");
  
  cout<<"Creating the manager"<<endl;
  AliAnalysisManager* mgr = new AliAnalysisManager("fmd_analysis","fmd_analysis");
  mgr->AddTask(FMDana0); 
  mgr->AddTask(FMDana1); 
  mgr->AddTask(FMDana2);
  mgr->AddTask(FMDana3);
  mgr->AddTask(FMDana4);
  
  
  AliAnalysisDataContainer* cin_esd = mgr->CreateContainer("esdTree",TTree::Class(),AliAnalysisManager::kInputContainer,"AliESDs.root");
  AliAnalysisDataContainer* cexchangevertex = mgr->CreateContainer("esdvertex",AliESDVertex::Class(),AliAnalysisManager::kExchangeContainer);
  AliAnalysisDataContainer* cexchange0 = mgr->CreateContainer("exchangeESDFMD0",AliESDEvent::Class(),AliAnalysisManager::kExchangeContainer);
  AliAnalysisDataContainer* cdiag1 = mgr->CreateContainer("diagSharing1",AliESDEvent::Class(),AliAnalysisManager::kExchangeContainer);
  AliAnalysisDataContainer* cdiag2 = mgr->CreateContainer("diagSharing2",TList::Class(),AliAnalysisManager::kOutputContainer,"edists.root");
  AliAnalysisDataContainer* cexchange1 = mgr->CreateContainer("exchangeESDFMD1",AliESDFMD::Class(),AliAnalysisManager::kExchangeContainer);
  AliAnalysisDataContainer* cexchange2 = mgr->CreateContainer("listOfhists",TList::Class(),AliAnalysisManager::kExchangeContainer);
  AliAnalysisDataContainer* cvertex = mgr->CreateContainer("vertex",TObjString::Class(),AliAnalysisManager::kExchangeContainer);
  AliAnalysisDataContainer* cexchange3 = mgr->CreateContainer("BackgroundCorrectedperevent",TList::Class(),AliAnalysisManager::kOutputContainer,"testOut.root");
  AliAnalysisDataContainer* coutput = mgr->CreateContainer("BackgroundCorrected",TList::Class(),AliAnalysisManager::kOutputContainer,outFile);
  
  mgr->ConnectInput(FMDana0, 0 , cin_esd);   
  mgr->ConnectOutput(FMDana0, 0 , cexchange0);
  
  mgr->ConnectInput(FMDana1, 0 , cexchange0);   

  mgr->ConnectOutput(FMDana1, 0 , cexchange1);  
  mgr->ConnectOutput(FMDana1, 1 , cexchangevertex);   
  mgr->ConnectOutput(FMDana1, 2 , cdiag1);
  mgr->ConnectOutput(FMDana1, 3 , cdiag2);
  
  
  mgr->ConnectInput(FMDana2, 0 , cexchange1);   
  mgr->ConnectInput(FMDana2, 1 , cexchangevertex);   
  mgr->ConnectOutput(FMDana2, 0 , cexchange2);
  
  mgr->ConnectInput(FMDana3, 0 , cexchange2);   
  mgr->ConnectOutput(FMDana3, 0 , cexchange3);
  mgr->ConnectOutput(FMDana3, 1 , cvertex);
  
  mgr->ConnectInput(FMDana4, 0 , cexchange3);   
  mgr->ConnectInput(FMDana4, 1 , cvertex);   
  mgr->ConnectOutput(FMDana4, 0 , coutput);
  
  TGrid::Connect("alien://",0,0,"t");
  
  
  TChain* chain = new TChain("esdTree","esdTree");
  
  TAlienCollection* coll =  TAlienCollection::Open(collectionfile);  
  coll->Reset();
  Int_t nFiles = 0;
  while(coll->Next() && nFiles<2) {
    cout<<coll->GetTURL("")<<endl;
    TString test(coll->GetTURL(""));
    chain->Add(coll->GetTURL(""));
    
    nFiles++;
  }
  
  mgr->InitAnalysis();
  mgr->PrintStatus();
  TStopwatch timer;
  timer.Start();
  cout<<"Executing analysis"<<endl;
  mgr->StartAnalysis("local",chain);
  timer.Stop();
  timer.Print();
}
