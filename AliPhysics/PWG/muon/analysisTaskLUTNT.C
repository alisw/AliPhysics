//________________________________________________________________________
void analysisTaskLUTNT() {

  TChain* chain = new TChain("esdTree");

  TGridCollection* coll = gGrid->OpenCollection("wn.xml");

  TGridResult* result = coll->GetGridResult("",0,0);
  Int_t nFiles = 0;
  for(Int_t i = 0; i < result->GetEntries(); i++) {
    printf("TURL = %s \n",result->GetKey(i,"turl"));
    chain->Add(result->GetKey(i,"turl"));
    nFiles++;
    //if (nFiles == 10) break;
  }

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisManager");

  //____________________________________________//
  // ntuple task
  AliAnalysisTaskLUT *task = new AliAnalysisTaskLUT("TaskLUT");
  mgr->AddTask(task);

  // Create containers for input/output

  // input
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("cchain",TChain::Class(),AliAnalysisManager::kInputContainer);

  Char_t text[256];
  sprintf(text,"Ntuple.LUT.root");
  printf("Analysis output in %s \n",text);

  // output
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("cntuple", TNtuple::Class(),AliAnalysisManager::kOutputContainer,text);

  //____________________________________________//
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,0,coutput);

  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();

    TStopwatch timer;
    timer.Start();

    mgr->StartAnalysis("local",chain);

    timer.Stop();
    timer.Print();

  }
}                         
                      
