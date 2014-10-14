
void runTRDqaAnalysis(const char *chainName, int limit = 0) {
  //
  // runs the analysis train
  // parameters: 
  // chainName -- a name of a file with a list of ESDs
  // limit -- number of files to be processed
  //
  //
  
  /**/
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  /**/

  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTRDqaAnalysis.so");
  
  // Setup chain
  TChain *chain = new TChain("esdTree");
  //chain->SetBranchStatus("*",0);
  //chain->SetBranchStatus("*fTracks*",1);
  //chain->SetBranchStatus("ESDfriend*",1);
  //chain->SetBranchStatus("ESDfriend.*",1);

  int nfiles = 0;
  fstream  coll(chainName, ios_base::in);
  TString line;
  while (line.ReadLine(coll)) {
    cout << line.Data() << endl;
    chain->Add(line.Data()); 
    nfiles++;
    if (limit && nfiles > limit) break;
  } 
  
  // Create an analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("qaTasks", "title");
  AliAnalysisTask *tasks[5];
  AliAnalysisDataContainer *out[5];

  tasks[0] = new AliTRDqaElectronSpectra("trdElectronSpectra");
  tasks[1] = new AliTRDqaESDFriends("trdESDFriends");
  tasks[2] = new AliTRDqaEnergyDeposit("trdEnergyDeposit");
  tasks[3] = new AliTRDqaBasic("trdBasic");
  tasks[4] = new AliTRDqaJPsi("trdJPsi");
  
  AliAnalysisDataContainer *cinput = 
    mgr->CreateContainer("inputESD", TTree::Class(), AliAnalysisManager::kInputContainer);

  out[0] = mgr->CreateContainer("oES", TObjArray::Class(), AliAnalysisManager::kOutputContainer);
  out[1] = mgr->CreateContainer("oEF", TObjArray::Class(), AliAnalysisManager::kOutputContainer);
  out[2] = mgr->CreateContainer("oED", TObjArray::Class(), AliAnalysisManager::kOutputContainer);
  out[3] = mgr->CreateContainer("oBS", TObjArray::Class(), AliAnalysisManager::kOutputContainer);
  out[4] = mgr->CreateContainer("oEJ", TObjArray::Class(), AliAnalysisManager::kOutputContainer);  

  // register
  for(int i=0; i<4; i++) {
    mgr->AddTask(tasks[i]);
    mgr->ConnectInput(tasks[i],0,cinput);
    mgr->ConnectOutput(tasks[i],0,out[i]);
  }

  // Connect input data
  cout << "connect to data" << endl;
  cinput->SetData(chain);
  Long_t t0 = gSystem->Now();

  TStopwatch sw;
  sw.Start();
  
  cout << "Initializing" << endl;
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }
  
  sw.Stop();
  sw.Print();

  Long_t t1 = gSystem->Now();
  double time = 1e-3 * (t1-t0);
  
  //cout << "Size   = " << mgr->GetNBytes()*1e-6 << " MB" << endl;
  //cout << "Time   = " << time << " s" << endl;
  //cout << "Speed  = " << 1e-6*mgr->GetNBytes()/time << " MB/s" << endl;
  //cout << "Events = " << mgr->GetNEvents()/time << " Events/s" << endl;   

  //gSystem->Exit(0);
  

}
