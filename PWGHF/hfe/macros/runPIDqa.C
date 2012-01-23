void SetupPar(char* pararchivename);
void SwitchBranches(TChain *c);
TChain * CreateXMLChain(char* xmlfile);
TChain * CreateAODChain(const char * name = 0x0);

void runPIDqa(const char *treelist = 0x0, Bool_t hasMC = kTRUE, Int_t nFiles = 5, Int_t nSkip = 0){
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
  gSystem->Load("libPWG3hfeDevel.so");
  AliLog::SetClassDebugLevel("AliHFEV0pid", 0);
  AliLog::SetClassDebugLevel("AliHFEpidQA", 0);
  //  gROOT->LoadMacro("AliAnalysisTaskHFE.cxx++");
  AliLog::SetGlobalLogLevel(AliLog::kWarning);

  AliHFEtools::SetLogLevel(0);
  
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
      gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
      inputChain = CreateESDChain(treelist, nFiles, nSkip);
      SwitchBranches(inputChain);
    }
  }
 
  // Start the Analysis Manager and Create Handlers
  AliAnalysisManager *hfeAnalysis = new AliAnalysisManager("Single Electron Analysis");
  AliESDInputHandler *esdInputHandler;
  if(!isAOD) {
    hfeAnalysis->SetInputEventHandler(esdInputHandler = new AliESDInputHandler);
    // otherwise your proof analysis crashes if no ESD friends are available
    //esdInputHandler->SetReadFriends(kFALSE);   
  }
  else hfeAnalysis->SetInputEventHandler(new AliAODInputHandler);
  if(hasMC && !isAOD) hfeAnalysis->SetMCtruthEventHandler(new AliMCEventHandler);

  // test the NN reference file
  TFile *f = new TFile("NNref_v2_data.root", "READ");
  if(!f){
    AliError("NN reference file not found");
    return;
  }

  AliAnalysisTaskHFEpidQA *qaTask = new AliAnalysisTaskHFEpidQA("HFEpidQA");
  qaTask->SetV0pidQA();
  qaTask->SetNNref(f);
  hfeAnalysis->AddTask(qaTask);
  qaTask->ConnectInput(0, hfeAnalysis->GetCommonInputContainer());
  qaTask->ConnectOutput(1, hfeAnalysis->CreateContainer("PIDResults", TList::Class(), AliAnalysisManager::kOutputContainer, "PIDqa.root"));
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
    printf("No AOD anlysis\n");
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
