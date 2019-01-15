void SetupPar(char* pararchivename);
void SwitchBranches(TChain *c);
TChain * CreateXMLChain(char* xmlfile);
TChain * CreateAODChain(const char * name = 0x0);

void runElectronTask(const char *treelist = 0x0, Bool_t hasMC = kTRUE){
  if(!treelist){
    printf("Error: No ESD list specified\n");
    return;
  }
  if(gSystem->Getenv("ALICE_ROOT")){
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
  }
  else{
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
    SetupPar("CORRFW");
    SetupPar("Util");
  }
  //SetupPar("HFE");
  gSystem->Load("libHFE");
//  gROOT->LoadMacro("AliAnalysisTaskHFE.cxx++");
 // AliLog::SetGlobalLogLevel(AliLog::kError);
  
  // Make the ESD chain
  TString treename = treelist;
  TChain *inputChain = 0x0;
  Bool_t isAOD = kFALSE;
  if(treename.EndsWith(".xml")){
    esdchain = CreateXMLChain(treelist);
    SwitchBranches(inputChain);
  }else{
    inputChain = CreateAODChain(treelist);
    if(inputChain) isAOD = kTRUE;
    else{
      gROOT->LoadMacro("CreateESDChain.C");
      inputChain = CreateESDChain(treelist, -1);
      SwitchBranches(inputChain);
    }
  }
 
  // Start the Analysis Manager and Create Handlers
  AliAnalysisManager *hfeAnalysis = new AliAnalysisManager("Single Electron Analysis");
  if(!isAOD) hfeAnalysis->SetInputEventHandler(new AliESDInputHandler);
  else hfeAnalysis->SetInputEventHandler(new AliAODInputHandler);
  if(hasMC && !isAOD) hfeAnalysis->SetMCtruthEventHandler(new AliMCEventHandler);
  AliHFEcuts *hfecuts = new AliHFEcuts;
  hfecuts->CreateStandardCuts();
  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE;
  if(hasMC) task->SetHasMCData();
  if(isAOD) task->SetAODAnalysis();
  task->SetHFECuts(hfecuts);
  task->SetPIDStrategy(4);
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
  task->SetQAOn(AliAnalysisTaskHFE::kMCqa);
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kSecVtx);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kIsElecBackGround);
  hfeAnalysis->AddTask(task);
  task->ConnectInput(0, hfeAnalysis->GetCommonInputContainer());
  task->ConnectOutput(1, hfeAnalysis->CreateContainer("nEvents", TH1I::Class(), AliAnalysisManager::kOutputContainer, "HFEtask.root"));
  task->ConnectOutput(2, hfeAnalysis->CreateContainer("Results", TList::Class(), AliAnalysisManager::kOutputContainer, "HFEtask.root"));
  task->ConnectOutput(3, hfeAnalysis->CreateContainer("QA", TList::Class(), AliAnalysisManager::kOutputContainer, "HFEtask.root"));
  //task->ConnectOutput(3, hfeAnalysis->CreateContainer("QA-new", TList::Class(), AliAnalysisManager::kOutputContainer, "HFEtask.root"));  
  // Run the Analysis
  if(hfeAnalysis->InitAnalysis()){  
    TStopwatch timer;
    timer.Start();
    hfeAnalysis->StartAnalysis("local", inputChain);
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
  gSystem->Load("libNetx") ;
  gSystem->Load("libRAliEn");
  TGrid::Connect("alien://") ;



  TGridCollection * collection = gGrid->OpenCollection(kXML);
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
TChain *CreateAODChain(const char * name){
  if(!name) return NULL;
  ifstream filelist(name);
  string filename;
  TChain *c = new TChain("aodTree");
  Bool_t isAODlist = kTRUE;
  while(getline(filelist,filename)){
    if(!strstr(filename.c_str(), "AliAOD")){
      isAODlist = kFALSE;
      break;
    }
    c->Add(filename.c_str());
  }
  if(!isAODlist){
    printf("No AOD anlysis");
    delete c;
    return NULL;
  }
  return c;
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

void SwitchBranches(TChain *c){
  //esdchain->SetBranchStatus("*", 0);
  c->SetBranchStatus("Calo*", 0);
  c->SetBranchStatus("*FMD*", 1);
  c->SetBranchStatus("Tracks", 1);
}
