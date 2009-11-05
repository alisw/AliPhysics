void SetupPar(char* pararchivename);
TChain * CreateXMLChain(char* xmlfile);

void runElectronTask(const char *treelist = 0x0){
  if(!treelist){
    printf("Error: No ESD list specified\n");
    return;
  }
  if(gSystem->Getenv("ALICE_ROOT")){
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    gSystem->Load("libPWG3hfe");
  }
  else{
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
    SetupPar("CORRFW");
    SetupPar("PWG3hfe");
  }
 // AliLog::SetGlobalLogLevel(AliLog::kError);
  
  // Make the ESD chain
  TString treename = treelist;
  TChain *esdchain = 0x0;
  if(treename.EndsWith(".xml"))
    esdchain = CreateXMLChain(treelist);
  else{
    gROOT->LoadMacro("CreateESDChain.C");
    esdchain = CreateESDChain(treelist, -1);
  }
  //esdchain->SetBranchStatus("*", 0);
  esdchain->SetBranchStatus("Calo*", 0);
  esdchain->SetBranchStatus("*FMD*", 1);
  esdchain->SetBranchStatus("Tracks", 1);
  
  // Start the Analysis Manager and Create Handlers
  AliAnalysisManager *electronAnalysis = new AliAnalysisManager("Single Electron Analysis");
  electronAnalysis->SetInputEventHandler(new AliESDInputHandler);
  electronAnalysis->SetMCtruthEventHandler(new AliMCEventHandler);
  AliHFEcuts *hfecuts = new AliHFEcuts;
  hfecuts->CreateStandardCuts();
  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE("Heavy Flavour Analysis");
  task->SetHFECuts(hfecuts);
  task->SetPIDStrategy(4);
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
  task->SetQAOn(AliAnalysisTaskHFE::kMCqa);
  task->SetSecVtxOn();
  electronAnalysis->AddTask(task);
  task->ConnectInput(0, electronAnalysis->GetCommonInputContainer());
  task->ConnectOutput(0, electronAnalysis->CreateContainer("nEvents", TH1I::Class(), AliAnalysisManager::kOutputContainer, "PWG3hfe.root"));
  task->ConnectOutput(1, electronAnalysis->CreateContainer("Results", TList::Class(), AliAnalysisManager::kOutputContainer, "PWG3hfe.root"));
  task->ConnectOutput(2, electronAnalysis->CreateContainer("QA", TList::Class(), AliAnalysisManager::kOutputContainer, "PWG3hfe.root"));
  
  // Run the Analysis
  if(electronAnalysis->InitAnalysis()){  
    TStopwatch timer;
    timer.Start();
    electronAnalysis->StartAnalysis("local", esdchain);
    timer.Stop();
    timer.Print();
  }
}

//____________________________________________
TChain * CreateXMLChain(char* xmlfile)
{

  //TChain *chain = 0x0;
  const char *chainname="esdTree";
  TChain* chain = new TChain(chainname);
 
  TString input =xmlfile;
  cout<<" the input is::"<< xmlfile<<endl;

  char kXML  [1000];
  if(gSystem->Getenv("XML") )
    kXML = gSystem->Getenv("XML");
  else
    sprintf(kXML, xmlfile) ; 

  //    sprintf(kXML, "collection.xml") ; 
  
  cout<<"XML file "<<kXML<<endl;
  
  if (!TFile::Open(kXML)) {
    printf("No collection file with name -- %s -- was found\n",kXML);
    return ;
  }
  gSystem->Load("libNetx.so") ;
  gSystem->Load("libRAliEn.so");
  TGrid::Connect("alien://") ;



  //  TGridCollection * collection =  (TGridCollection*)gROOT->ProcessLine(Form("TAlienCollection::Open(\"%s\", 0)", kXML));
  TGridCollection * collection = (TGridCollection*) TAlienCollection::Open(kXML);
  if (! collection) {
    AliError(Form("%s not found", kXML)) ; 
    return kFALSE ; 
  }
  //collection->CheckIfOnline();

  TGridResult* result = collection->GetGridResult("",0 ,0);

  //  TList* analysisfilelist = result->GetFileInfoList();
  
  // Makes the ESD chain 
  printf("*** Getting the Chain       ***\n");

  for (Int_t index = 0; index < result->GetEntries(); index++) {
    TString alienURL = result->GetKey(index, "turl") ; 
    cout << "================== " << alienURL << endl ; 
    chain->Add(alienURL) ; 
    //alienURL.ReplaceAll("AliESDs.root",kXSFileName);
    // chainxs->Add(alienURL) ; 
  }


  //  chain->AddFileInfoList(analysisfilelist);
  if (chain) chain->ls();
  chain->SetBranchStatus("*Calo*",0);
  chain->SetBranchStatus("*FMD*",0);
  return chain;
}

//______________________________________________________________________________
void SetupPar(char* pararchivename)
{
  if (pararchivename) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    TString ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchivename);
    
    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");
      
      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
	Error("runProcess","Cannot Build the PAR Archive! - Abort!");
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
    
    gSystem->ChangeDirectory(ocwd.Data());
    printf("Current dir: %s\n", ocwd.Data());
  }
}
