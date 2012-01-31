void run(const char* esdList = "esds.list", Bool_t bUseMCInfo = kFALSE, Int_t nFiles=100, Int_t nEvents=500, Int_t firstEvent =0) {
  TStopwatch timer;
  timer.Start();
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT -I$ALICE_ROOT/TRD");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("libPWG0dep.so");
  gSystem->Load("libPWG0selectors.so");
 
  Bool_t dodNdPtCutAnalysis = kTRUE;   // AlidNdPtTask (data/MCtruth)

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  mgr->SetDebugLevel(0);

  AliESDInputHandlerRP* esdH = new AliESDInputHandlerRP();
  esdH->SetReadFriends(0);
  mgr->SetInputEventHandler(esdH);  

  if(bUseMCInfo) {
    // Enable MC event handler
    AliMCEventHandler* handler = new AliMCEventHandler;
    //handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  // Chain 
  // TGrid::Connect("alien://");
  
  //TChain* chain = new TChain("esdTree");
  //chain->AddFile("inFile");

  // Create input chain
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain(esdList,nFiles);
  if(!chain) {
    printf("ERROR: chain cannot be created\n");
    return;
  }
  chain->Lookup();

  //
  // Wagons to run 
  //
  // dNdPtCutAnalysis (Jacek Otwinowski)
  //
  if (dodNdPtCutAnalysis) {
      // 
      // Optionally MC information can be used by setting the 1st argument to true
      //
      gROOT->LoadMacro("$ALICE_ROOT/PWG0/dNdPt/macros/AddTask_dNdPtCutAnalysisPbPb.C");
      AddTask_dNdPtCutAnalysisPbPb();
  }
  
  // Init
  if (!mgr->InitAnalysis()) 
      mgr->PrintStatus();
      mgr->PrintStatus();
  // Run on dataset
  //
  
  mgr->StartAnalysis("local",chain,nEvents, firstEvent);

  timer.Stop();
  timer.Print();
}

