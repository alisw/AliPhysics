/*
  How to debug:
  Check that AliRoot and ROOT versions are OK!
  
  to run:
  aliroot
  .L runAAF.C
  runAAF(1, "dataESDs.txt", 2)

*/

Int_t mode        = 0; // mode 1: PROOF, 2: GRID, 0: LOCAL
char *gridmode    = "test";//full, terminate, test, reset
Bool_t debug      = kTRUE;
Bool_t readTR     = kFALSE;// always false
Int_t runtype     = 3;//3 pp, 2 PbPb, 4 pPb
Bool_t esdAna     = kTRUE; //FALSE if AOD
Bool_t analysisMC = kFALSE;
Bool_t ispileuprej= kFALSE;
char *centralityEstimator = "V0A";
UInt_t kTriggerInt = AliVEvent::kMB;//pPb: kINT7, other kMB

void runAAF(Int_t nFilesMax, char* textFileName = "esd.txt", Int_t task = 2);
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
void runAAF(Int_t nFilesMax, char* textFileName, Int_t task)
{  
  const Char_t* treeName = GetTreeName(esdAna); 
  cout << "MonteCarlo " << analysisMC << endl;  
  const char* alirootver = "v5-04-51-AN";
  
  
  printf("===================================================================\n");
  printf("===================================================================\n");
  if (analysisMC) printf(":: use MC    TRUE\n");
  else            printf(":: use MC    FALSE\n");
  if (readTR)     printf(":: read TR   TRUE\n");
  else            printf(":: read TR   FALSE\n");
  if (debug)      printf(":: debugging TRUE\n");
  else            printf(":: debugging FALSE\n");
  
  Char_t taskname[128];
  sprintf(taskname,"TransverseEventShapeTask");
  Char_t nameouputfiles[1280]={0};
  
  if(runtype ==3 || runtype ==4){
    sprintf(nameouputfiles,"%s %s.root",taskname);
  }

  cout<<"Files to be stored:  "<<nameouputfiles<<endl;
  // Load common libraries   
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libSTEER.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libCDB.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWGLFspectra.so"); 
  
  gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_ROOT")));
  gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_PHYSICS")));

  if (mode == 1){
    cout << "Connecting to CAF..." << endl;
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    const char* AAF = "cuautle@alice-caf.cern.ch";    
    
    if(strstr(gridmode,"reset") != NULL)
      TProof::Reset(AAF,kFALSE);
    TProof::Open(AAF, "workers=10");
    
    TList *list = new TList();
    list->Add(new TNamed("ALIROOT_EXTRA_LIBS", "CDB"));    
    gProof->EnablePackage(Form("VO_ALICE@AliRoot::%s",alirootver), list);    
    cout << "Connected to " << AAF << endl;
  }
    
  // Load the analysis macro (old version); should be removed when the library will be decided
  Char_t loadtask[128];
  sprintf(loadtask,"AliAna%s.cxx+", taskname);
  
  if(mode == 1){ // PROOF
    gProof->Load(loadtask);
  }else{
    gROOT->LoadMacro(loadtask);
  }
  
  // Make the analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager("PID histos", "testing analysis");
  
  
  // Dataset 
  switch(mode){
  case 1: // PROOF
    TFileCollection* proofColl = gProof->GetDataSet(textFileName);
    TFileCollection* stagedColl = proofColl->GetStagedSubset();
    TChain* chain = CreateChainCAF(nFilesMax, stagedColl, treeName);
    break;
  case 2: // GRID
    TGrid::Connect("alien://");    
    gROOT->LoadMacro("CreateAlienHandler.C");
    //    AliAnalysisGrid *alienHandler = CreateAlienHandler(nFilesMax, analysisMC, runtype, taskname, gridmode);  
    if(task==4)
      AliAnalysisGrid *alienHandler = CreateAlienHandler(nFilesMax, analysisMC, esdAna, taskname1, nameouputfiles, gridmode, textFileName, alirootver, task);  
    else
      AliAnalysisGrid *alienHandler = CreateAlienHandler(nFilesMax, analysisMC, esdAna, taskname, nameouputfiles, gridmode, textFileName, alirootver, task);  
    
    if (!alienHandler) return; 
    
    
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
    esdHandler->SetNeedField();
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
  if(esdAna){
    gROOT->LoadMacro("$(ALICE_PHYSICS)/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
    if (analysisMC) {
      AliPhysicsSelection* physSel = physSelTask->GetPhysicsSelection();
      physSel->SetAnalyzeMC();
    }
  }
  
    if(esdAna && runtype==2 || esdAna && runtype==4  || esdAna && runtype == 3)  { // only for ESD and pPb pp    
      //gROOT->LoadMacro("$(ALICE_PHYSICS)/OADB/macros/AddTaskCentrality.C");
      //AliCentralitySelectionTask *taskCentrality = AddTaskCentrality(); 

      //    taskCentrality->SetPass(2);
  if (analysisMC)
    taskCentrality->SetMCInput(); 
}

  
  // ######### ESA task ###############
  gROOT->LoadMacro(Form("Add%s.C", taskname));  
  AliAnaTransverseEventShapeTask* taskESA = AddTask(analysisMC);
 
//cout<<"!!!!!!!!!!!!!!!!!!!!   GetTest()="<<taskESA->GetTest()<<endl;
 //return; 
  // Run the analysis  
  if (mgr->InitAnalysis()){
    mgr->PrintStatus();
    switch(mode){
    case 1: // PROOF
      mgr->StartAnalysis("proof",chain);
      break;
    case 2: // GRID
      mgr->StartAnalysis("grid");
//cout<<"!!!!!!!!!!!!!!!!!!!!   GetTest()="<<taskESA->GetTest()<<endl;
      break;
    case 0:
      mgr->StartAnalysis("local",chain);
//cout<<"!!!!!!!!!!!!!!!!!!!!   GetTest()="<<taskESA->GetTest()<<endl;
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



