//####################################################//
//###           To execute the macro:              ###//
//###      $ALIEN/globus/bin/grid-proxy-init       ###//
//###   $ALIEN/api/bin/alien-token-init pchrista   ###//
//###            root.exe -q DemoTags             ###//
//####################################################//

Int_t demoJETAN() {
  TStopwatch timer;
  timer.Start();

  const char* pararchiveJETAN = "JETAN";
  const char* pararchiveESD   = "ESD";
  //////////////////////////////////////////
  // Libraries required to load
  //////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // Setup PAR File
  if (pararchiveJETAN) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchiveJETAN);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchiveJETAN);

    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("**************************************\n");
      printf("*** Building PAR archive for JETAN ***\n");
      printf("**************************************\n");

      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
        Error("batchSelector","Cannot Build the PAR Archive! - Abort!");
        return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
    
    gSystem->ChangeDirectory("../");
  }


  if (pararchiveESD) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchiveESD);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchiveESD);

    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("**************************************\n");
      printf("*** Building PAR archive for JETAN ***\n");
      printf("**************************************\n");

      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
        Error("batchSelector","Cannot Build the PAR Archive! - Abort!");
        return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
    
    gSystem->ChangeDirectory("../");
  }

  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so"); 
  gSystem->Load("libJETAN.so");
  
  

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://"); 

  //////////////////////////////////////////////////////////////////
  // Create A tag analysis object and impose some selection criteria
  AliTagAnalysis *TagAna = new AliTagAnalysis(); 
  // create an EventTagCut object
  AliEventTagCuts *EvCuts = new AliEventTagCuts();
  AliRunTagCuts   *RuCuts = new AliRunTagCuts();
  //EvCuts->SetNChargedAbove1GeVRange(1, 1000);
  //EvCuts->SetMultiplicityRange(11,120);
  //EvCuts->SetNPionRange(2,10000);

  printf("*******************************\n");
  printf("*** Querying the tags       ***\n");
  printf("*******************************\n");
  
  //grid tags
  TAlienCollection* coll = TAlienCollection::Open("tag100.xml");
  TGridResult* TagResult = coll->GetGridResult("");
  TagAna->ChainGridTags(TagResult);

  //////////////////////////////////////////////////////////////////
  //Get the chain
  printf("*******************************\n");
  printf("*** Getting the Chain       ***\n");
  printf("*******************************\n");
  TChain* analysischain = 0x0;
  analysischain = TagAna->QueryTags(RuCuts, EvCuts);

  /////////////////////////////////////////////////////////////////
  // Run the Analysis Selector
  const char *selectorfile = "AliJetSelector.C";
  printf("*******************************\n");
  printf("*** Run Analysis Selector %s\n",selectorfile);
  printf("*******************************\n");
  analysischain->ls();
  
  analysischain->Process(selectorfile);

  timer.Stop();
  timer.Print();

  return 0;
}
