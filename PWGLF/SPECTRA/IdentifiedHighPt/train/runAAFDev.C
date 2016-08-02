/*
  How to debug:
  Check that AliRoot and ROOT versions are OK!
  
  to run:
  export X509_CERT_DIR=$ALIEN_ROOT/globus/share/certificates
  aliroot
  .L runAAFDev.C
  
  runAAF(2, "local aod PbPb MC", "test_aod_mc.txt", 4)
  
*/

UInt_t kTriggerInt = AliVEvent::kMB;
Float_t minCent = 0.0;
Float_t maxCent = 100.0;

void runAAF(Int_t nFilesMax, char* type = "local", char* textFileName = "esd.txt", Int_t task = 2);
TChain* CreateChainCAF(Int_t nFilesMax, TFileCollection* coll, char* treeName);
TChain* CreateChainLocal(Int_t nFilesMax, char* filename, char* treeName);
const Char_t* GetTreeName(Bool_t esdAna);

//________________________________________________________________
const Char_t* GetTreeName(Bool_t esdAna)
{
  if(esdAna)
    return "esdTree";
  else
    return "aodTree";    
}

//________________________________________________________________
void runAAF(Int_t nFilesMax, char* type, char* textFileName, Int_t task)
{
  Int_t analysisMC  = 0;
  Bool_t debug      = kTRUE;
  Bool_t readTR     = kFALSE;
  Bool_t esdAna     = kTRUE;  // run ESD analysis (kTRUE) or AOD analysis (kFALSE)
  Int_t runtype = -1;
  
  if(strstr(type,"aod") != NULL) {
    esdAna = kFALSE;; // AOD analysis
    cout << "Performing AOD analysis" << endl;
  } else 
    cout << "Performing ESD analysis" << endl;
  
  const Char_t* treeName = GetTreeName(esdAna);
  
  if(strstr(type,"PbPb") != NULL) runtype = 2; // PbPb
  if(strstr(type,"pp") != NULL) runtype = 3; // pp
  if(strstr(type,"900") != NULL){ // 900GeV pp
    if(runtype != -1){
      printf("conflicting run types\n");
      return;
    }
    runtype = 1;
  }
  if(runtype == -1) runtype = 0; // 7TeV pp
  
  char *test;
  test = strstr(type,"MC");
  if(test != NULL){
    analysisMC = kTRUE;
    cout << "Test: " << test << endl;
    if(sscanf(test,"MC%d",&analysisMC)){
      //  cout << "**************" << analysisMC << endl;
      if(!analysisMC) analysisMC = 1;
    }
  }
  cout << "MonteCarlo " << analysisMC << endl;
  
  const char* alirootver = "v5-02-17-AN";
  Int_t mode = -1; // mode 1: PROOF, 2: GRID, 0: LOCAL
  printf("===================================================================\n");
  if(strstr(type,"proof") != NULL){ 
    if(mode != -1){
      printf("===== CONFLICTING TYPES =====\n");
      return;
    }
    mode = 1; 
    printf("===============   RUNNING ANALYSIS IN CAF MODE   ================\n");
  }
  if(strstr(type,"grid") != NULL){
    if(mode != -1){
      printf("===== CONFLICTING TYPES =====\n");
      return;
    }
    printf("===============  RUNNING ANALYSIS IN GRID MODE   ================\n");
    mode = 2; 
  }
  if(strstr(type,"local") != NULL || mode<0){ 
    if(mode != -1){
      printf("===== CONFLICTING TYPES =====\n");
      return;
    }
    printf("===============  RUNNING ANALYSIS IN LOCAL MODE  ================\n");
    mode = 0;
  }
  
  printf("===================================================================\n");
  printf("===================================================================\n");
  if (analysisMC) printf(":: use MC    TRUE\n");
  else            printf(":: use MC    FALSE\n");
  if (readTR)     printf(":: read TR   TRUE\n");
  else            printf(":: read TR   FALSE\n");
  if (debug)      printf(":: debugging TRUE\n");
  else            printf(":: debugging FALSE\n");
  
  Char_t taskname[128];
  
  switch(task){
  case 1:
    sprintf(taskname,"HighPtDeDxCosmic");
    break;
  case 2:
    sprintf(taskname,"HighPtDeDx");
    break;
  case 3:
    sprintf(taskname,"HighPtDeDxV0");
    break;
  case 4:{
    sprintf(taskname,"EffDecomposition");
  }
    break;
  default:
    printf("Unknown task\n");
    return;
  }
  
  Char_t nameouputfiles[1280]={0};
  if(runtype ==2){
    for(Int_t i=0;i<6;++i){
      //Char_t tmpname[128]={0};
      sprintf(nameouputfiles,"%s %s_Tree_%1.0f_%1.0f.root",nameouputfiles,taskname,minCent[i],maxCent[i]);
      
    }
  }
  if(runtype ==3){
    sprintf(nameouputfiles,"%s %s_Tree.root",taskname);
  }
  


  cout<<"Files to be stored:"<<nameouputfiles<<endl;


  // Load common libraries   
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libCDB.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  // tender

  gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_ROOT")));


  cout << "mode " << mode << ", type " << type << endl;

  if (mode == 1){
    cout << "Connecting to CAF..." << endl;
    gEnv->SetValue("XSec.GSI.DelegProxy","2");

    //const char* AAF = "pchristi@skaf.saske.sk";    
    const char* AAF = "aortizve@alice-caf.cern.ch";    

    if(strstr(type,"reset") != NULL)
       TProof::Reset(AAF,kFALSE);
    //       TProof::Reset("pchristi@alice-caf.cern.ch",kTRUE);
    //    TProof::Open(AAF);
    TProof::Open(AAF, "workers=10");

    TList *list = new TList();
    list->Add(new TNamed("ALIROOT_EXTRA_LIBS", "CDB"));

    gProof->EnablePackage(Form("VO_ALICE@AliRoot::%s",alirootver), list);

    cout << "Connected to " << AAF << endl;
  }

 
  // Load the analysis macro (old version); should be removed when the library will be decided
  Char_t loadtask[128];

  if(task==4) {
    sprintf(loadtask,"AliAnaTaskV0%s.cxx++", taskname);
  }else{
    sprintf(loadtask,"AliAnalysisTask%s.cxx++", taskname);
  }
  
  if(mode == 1){ // PROOF
    
    // switch(task){
    // case 1: // ch fluct
    //   break;
    // case 2: // high pt dedx
    //gProof->Load("DebugClasses.C++g");
    //   break;
    // case 3: // high pt v0s
    //   gProof->Load("DebugClasses.C++g");
    //   break;
    // default:
    //   printf("Unknown task\n");
    //   return;
    // }

    gProof->Load(loadtask);
  }else{

    if(task==2) {
      
      cout << "Analysis task is loaded with AliRoot" << endl;
    } else {

      gROOT->LoadMacro(loadtask);
    }
  }

  // Make the analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager("PID histos", "testing analysis");


  // Dataset 
  switch(mode){
  case 1: // PROOF
    TFileCollection* proofColl = gProof->GetDataSet(textFileName);
    TFileCollection* stagedColl = proofColl->GetStagedSubset();
    TChain* chain = CreateChainCAF(nFilesMax, stagedColl, treeName);
//    TChain* chain = CreateChainCAF(nFilesMax, textFileName, treeName, analysisMC);
    break;
  case 2: // GRID
    //gSystem->Setenv("alien_CLOSE_SE", "ALICE::GSI::SE2");
    //    gSystem->Setenv("alien_CLOSE_SE", "ALICE::NIHAM::FILE");
    TGrid::Connect("alien://");
    //gSystem->Setenv("alien_CLOSE_SE", "ALICE::GSI::SE2");
    //    gSystem->Setenv("alien_CLOSE_SE", "ALICE::NIHAM::FILE");
    
    Char_t gridmode[64];
    sprintf(gridmode, "full");
    if(strstr(type,"term") != NULL) 
      sprintf(gridmode, "terminate");
    if(strstr(type,"off") != NULL) 
      sprintf(gridmode, "offline");
    if(strstr(type,"sub") != NULL) 
      sprintf(gridmode, "submit");
    if(strstr(type,"test") != NULL) 
      sprintf(gridmode, "test");


    gROOT->LoadMacro("CreateAlienHandler.C");
    //    AliAnalysisGrid *alienHandler = CreateAlienHandler(nFilesMax, analysisMC, runtype, taskname, gridmode);  
    AliAnalysisGrid *alienHandler = CreateAlienHandler(nFilesMax, analysisMC, esdAna, taskname, nameouputfiles, gridmode, textFileName, alirootver, task);  
    if (!alienHandler) return; 

    // DOES NOT WORK BECAUSE THERE ARE NO GETTERS?
    // // Here we can add extra files to the plugin
    // switch(task){
    // case 1: // ch fluct
    //   break;
    // case 2: // high pt dedx
    //   alienHandler->SetAnalysisSource(Form("DebugClasses.C %s", alienHandler->GetAnalysisSource()));
    //   alienHandler->SetAdditionalLibs(Form("DebugClasses.C %s", alienHandler->GetAdditionalLibs()));
    //   break;
    // case 3: // high pt v0s
    //   alienHandler->SetAnalysisSource(Form("DebugClasses.C %s", alienHandler->GetAnalysisSource()));
    //   alienHandler->SetAdditionalLibs(Form("DebugClasses.C %s", alienHandler->GetAdditionalLibs()));
    //   break;
    // default:
    //   printf("Unknown task\n");
    //   return;
    // }

    // Connect plugin to the analysis manager
    mgr->SetGridHandler(alienHandler); 
    break;
  case 0: // LOCAL
    // Process data - chain
    AliXRDPROOFtoolkit tool;
    TChain* chain = tool.MakeChain(textFileName,treeName, 0, 100);
    chain->Lookup();
    break;
  default:
    printf("Unknown mode");
    return;
  }
  
  // ESD input handler
  if(esdAna) {

    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdHandler);
  } else {

    AliAODInputHandler* aodHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodHandler);
  }

  // Monte Carlo handler
  if (analysisMC) {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    if(esdAna)
      mgr->SetMCtruthEventHandler(mcHandler);
    mcHandler->SetReadTR(readTR); 
  }

  // Debug if needed
  if (debug) 
    mgr->SetDebugLevel(3);


  // ######### Centrality task ###############  

  cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
  cout<<"esdAna="<<esdAna<<"  runtype="<<runtype<<endl;
    
  // ######### PHYSICS SELECTION ###############
  if (esdAna) { // PC: 17/8-15: this does not work for AODs..... 

    cout << "Adding task physics selection" << endl;
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
    if (analysisMC) {
      AliPhysicsSelection* physSel = physSelTask->GetPhysicsSelection();
      physSel->SetAnalyzeMC();
    }
  }

  if(esdAna && runtype==2) { // only for ESD and PbPb

    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality(); 
    //taskCentrality->SetPass(2);
    if (analysisMC)
      taskCentrality->SetMCInput(); 
  }

  // ######### PID task ###############

  if(task == 2) {
    cout << "Loading macro in AliRoot!" << endl;
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/SPECTRA/IdentifiedHighPt/train/macros/AddTask.C");      
    cout << "Adding task!" << endl;
    AliAnalysisTask* taskPid = AddTask(analysisMC, taskname, runtype, kTriggerInt, minCent, maxCent);
  } else {
    cout << "Adding AliAnaTaskV0EffDecomposition" << endl;
    gROOT->LoadMacro("AddTask.C");  
    AliAnalysisTask* taskPid = AddTask(taskname);
  }
  
  // Run the analysis  
  if (mgr->InitAnalysis()){
    mgr->PrintStatus();
    switch(mode){
    case 1: // PROOF
      mgr->StartAnalysis("proof",chain);
      break;
    case 2: // GRID
      mgr->StartAnalysis("grid");
      break;
    case 0:
      mgr->StartAnalysis("local",chain);
      break;
    default:
      printf("Unknown mode\n");
      return;
    }
  } 
}  

//__________________________________________________________
TChain* CreateChainCAF(Int_t nFilesMax, TFileCollection* coll, char* treeName)
{
  TIter iter(coll->GetList());
  
  TChain* target = new TChain(treeName);
  
  Int_t nFiles = 0;
  
  TFileInfo* fileInfo = 0;
  while ((fileInfo = dynamic_cast<TFileInfo*> (iter())) && (nFiles<nFilesMax || nFilesMax == 0)){
    if (fileInfo->GetFirstUrl()) {
      target->Add(fileInfo->GetFirstUrl()->GetUrl());
      nFiles++;
    }
  }
  
  Printf("Added %d files to chain", target->GetListOfFiles()->GetEntries());
  
  return target;
}

//__________________________________________________________
TChain* CreateChainLocal(Int_t nFilesMax, char* filename, char* treeName)
{
  // If file name ends in .root the chain will be created directly from that
  // fileName, otherwise it will assume that it is a text file and create it
  // from that.
  //
  // nFilesMax is only used for the text files
  // nFilesMax=0 means add all files

  TChain* chain = new TChain(treeName); 
  
  // Open the input stream
  ifstream in;
  in.open(filename);

  Int_t nFiles = 0;

  // Read the input list of files and add them to the chain
  TString file;
  while(in.good() && (nFiles<nFilesMax || nFilesMax<=0)) {
    in >> file;
    if (!file.Contains("root")) 
      continue; // protection

    nFiles++;
    chain->Add(file.Data());
  }

  in.close();

  return chain;
}