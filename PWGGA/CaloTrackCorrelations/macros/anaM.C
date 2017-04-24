/// \file anaM.C
/// \ingroup CaloTrackCorrMacros
/// \brief Example of execution macro for mixing frame analysis
///
/// Example macro to do analysis with the
/// analysis classes in CaloTrackCorrelations within the mixing frame
/// in local or grid modes.
///
/// Not being used for long time, keep it in case it is useful in future
///
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

//-------------------------------------------------
enum anaModes {mLocal, mGRID};
//mLocal    = 0: Analyze locally files in your computer
//mGRID     = 3: Analyze files on GRID

//---------------------------------------------------------------------------
//Settings to read locally several files, only for "mLocal" mode
//The different values are default, they can be set with environmental 
//variables: INDIR, PATTERN, NFILES, respectively
char * kInDir = "/Users/ymao/group/ana/7TeV/corr";
char * kPattern = ""; // Data are in files kInDir/kPattern+i
Int_t kFile = 1; // Number of files
//---------------------------------------------------------------------------
//Collection file for grid analysis
char * kXML = "collection.xml";

//---------------------------------------------------------------------------

const Bool_t kMC = kFALSE; //With real data kMC = kFALSE
TString kInputData = "ESD";//ESD, AOD, MC
TString kTreeName ;
const TString calorimeter = "EMCAL" ;
const Bool_t kUsePAR = kFALSE; //set to kFALSE for libraries
//const Bool_t kUsePAR = kTRUE; //set to kFALSE for libraries
const Bool_t kDoESDFilter = kFALSE;  //filter the tracks from the esd

Int_t mode = mGRID;

//________________________
/// Main execution method.
//________________________
void anaM()
{
  // Main
  //--------------------------------------------------------------------
  // Load analysis libraries
  // Look at the method below, 
  // change whatever you need for your analysis case
  // ------------------------------------------------------------------
  LoadLibraries() ;
  
  //-------------------------------------------------------------------------------------------------
  // Create chain from ESD and from cross sections files, look below for options.
  //------------------------------------------------------------------------------------------------- 
  if(kInputData == "ESD") kTreeName = "esdTree" ;
  else if(kInputData == "AOD") kTreeName = "aodTree" ;
  else if (kInputData == "MC") kTreeName = "TE" ;
  else {
    cout<<"Wrong  data type "<<kInputData<<endl;
    break;
  }
  
  TChain * chain   = new TChain(kTreeName) ;
  
  CreateChain(mode, chain);//, chainxs);
  cout<<"Chain created"<<endl;
  
  if( chain )
  {
    AliLog::SetGlobalLogLevel(AliLog::kError);//Minimum prints on screen
    
    //--------------------------------------
    // Make the analysis manager
    //-------------------------------------
    AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");
    // MC handler
    if( (kMC && (kInputData == "ESD")) || kInputData == "MC"){
      AliMCEventHandler* mcHandler = new AliMCEventHandler();
      mcHandler->SetReadTR(kFALSE);//Do not search TrackRef file
      mgr->SetMCtruthEventHandler(mcHandler);
      if( kInputData == "MC") mgr->SetInputEventHandler(NULL);
    }
    
//    // AOD output handler
//    AliAODHandler* aodoutHandler   = new AliAODHandler();
//    aodoutHandler->SetOutputFileName("AliAOD.root");
//    mgr->SetOutputEventHandler(aodoutHandler);
    
    //input
    Int_t maxiterations = 1;
    AliEventPoolLoop* pool = new AliEventPoolLoop(maxiterations);
    pool->SetChain(chain);
    Int_t eventsInPool = 10;
    AliMultiEventInputHandler *inpHandler = NULL ; 
    if(kInputData == "ESD"){
      // ESD handler
      printf("ESD MultiInput \n");
      inpHandler = new AliMultiEventInputHandler(eventsInPool, 0);
    }
    if(kInputData == "AOD"){
      // AOD handler
      inpHandler = new AliMultiEventInputHandler(eventsInPool, 1);	   	  
    }
    mgr->SetInputEventHandler(inpHandler);
    cout<<"Input handler "<<mgr->GetInputEventHandler()<<endl;
    mgr->SetEventPool(pool);
    inpHandler->SetEventPool(pool);
    
    //mgr->SetDebugLevel(-1); // For debugging, do not uncomment if you want no messages.
    
    // select triigger events for physics run
    
//    if(!kMC){
//      gROOT->LoadMacro("AddTaskPhysicsSelection.C");
//      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
//      mgr->AddTask(physSelTask);    
//    }
    
    //-------------------------------------------------------------------------
    //Define task, put here any other task that you want to use.
    //-------------------------------------------------------------------------
    
    //correlation analysis
    gROOT->LoadMacro("AddTaskCaloTrackCorrM.C");
    
    AliAnalysisTaskCaloTrackCorrelationM *taskEMCAL = AddTaskCaloTrackCorrM(kInputData,"EMCAL",kFALSE);

    mgr->AddTask(taskEMCAL);
    
    AliAnalysisTaskCaloTrackCorrelationM *taskPHOS  = AddTaskCaloTrackCorrM(kInputData,"PHOS", kFALSE);

    mgr->AddTask(taskPHOS);
    
    //gROOT->LoadMacro("AddTaskChargeCorr.C");
    AliAnalysisTaskCaloTrackCorrelationM *taskCharge  = AddTaskCaloTrackCorrM(kInputData, "CTS",kFALSE);
//    if(!kMC)
//      taskCharge->SelectCollisionCandidates();
    mgr->AddTask(taskCharge);
    
     //-----------------------
    // Run the analysis
    //-----------------------    
    //mgr->ResetAnalysis();
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("mix",chain);
    
    cout <<" Analysis ended sucessfully "<< endl ;
  }
  else cout << "Chain was not produced ! "<<endl;
  
}

//_____________________________
/// Load analysis libraries. Out of date.
//_____________________________
void  LoadLibraries()
{
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  if(kUsePAR){
    //--------------------------------------------------------
    //If you want to use root and par files from aliroot
    //--------------------------------------------------------  
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
    SetupPar("PHOSUtils");
    SetupPar("EMCALUtils");
    
    SetupPar("PWGCaloTrackCorrBase");
    SetupPar("PWGGACaloTrackCorrelations");
  }
  else{
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPHOSUtils");
    gSystem->Load("libEMCALUtils");
    gSystem->Load("libPWGCaloTrackCorrBase");
    gSystem->Load("libPWGGACaloTrackCorrelations");
  }
}

//_________________________________
/// Par files compilation method.
//_________________________________
void SetupPar(char* pararchivename)
{
  //Load par files, create analysis libraries
  //For testing, if par file already decompressed and modified
  //classes then do not decompress.
  
  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ; 
  TString parpar(Form("%s.par", pararchivename)) ; 
  if ( gSystem->AccessPathName(pararchivename) ) {  
    TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }
  
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(pararchivename);
  
  // check for BUILD.sh and execute
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
    printf("*******************************\n");
    printf("*** Building PAR archive    ***\n");
    cout<<pararchivename<<endl;
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
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    gROOT->Macro("PROOF-INF/SETUP.C");
  }
  
  gSystem->ChangeDirectory(ocwd.Data());
  printf("Current dir: %s\n", ocwd.Data());
}


//_____________________________________________________________________
/// Fills chain with data files paths.
//_____________________________________________________________________
void CreateChain(const anaModes mode, TChain * chain){//, TChain * chainxs){
  //Fills chain with data
  
  TString datafileName="";
  if(kInputData == "ESD") datafileName = "AliESDs.root" ;
  else if(kInputData == "AOD") datafileName = "AliAOD.root" ;
  else if(kInputData == "MC")  datafileName = "galice.root" ;
  
  TString ocwd = gSystem->WorkingDirectory();
  
  //---------------------------------------
  // Local files analysis
  //---------------------------------------
  if(mode == mLocal)
  {
    //If you want to add several ESD files sitting in a common directory INDIR
    //Specify as environmental variables the directory (INDIR), the number of files 
    //to analyze (NFILES) and the pattern name of the directories with files (PATTERN)
    
    cout<<"INDIR : "<<kInDir<<endl;
    cout<<"NFILES : "<<kFile<<endl;
    cout<<"PATTERN: " <<kPattern<<endl;
    
    
    //Loop on ESD files, add them to chain
    TString FileName ;      
    for (Int_t iFile = 0 ; iFile < kFile ; iFile++) {
      FileName = Form("%s/%s%d/%s", kInDir,kPattern,iFile,datafileName.Data()) ; 
      //cout << "FileName: " << FileName <<endl ;
      TFile * dataFile = 0 ; 
      //Check if file exists and add it, if not skip it
      if ( dataFile = TFile::Open(FileName.Data())) {
        if ( dataFile->Get(kTreeName) ) { 
          Int_t nEventsPerFile = ((TTree*) dataFile->Get(kTreeName)) ->GetEntries();
          printf(" ++++ Adding %s, with %d events \n", FileName.Data(), nEventsPerFile) ;
          chain->AddFile(FileName);
        }
      }
    }    
    printf("number of entries # %lld \n", chain->GetEntries()) ; 	
  }// local files analysis
  
  //------------------------------
  //GRID xml files
  //-----------------------------
  else if(mode == mGRID)
  {
    //Load necessary libraries and connect to the GRID
    gSystem->Load("libNetx") ;
    gSystem->Load("libRAliEn");
    TGrid::Connect("alien://") ;
    
    //Feed Grid with collection file
    //TGridCollection * collection =  (TGridCollection*)gROOT->ProcessLine(Form("TAlienCollection::Open(\"%s\", 0)", kXML));
    TGridCollection * collection = (TGridCollection*) TAlienCollection::Open(kXML);
    if (! collection) {
      AliError(Form("%s not found", kXML)) ; 
      return kFALSE ; 
    }
    TGridResult* result = collection->GetGridResult("",0 ,0);
    
    // Makes the ESD chain 
    printf("*** Getting the Chain       ***\n");
    Int_t nEventsPerFile = 0;
    for (Int_t index = 0; index < result->GetEntries(); index++) {
      TString alienURL = result->GetKey(index, "turl") ; 
      cout << "================== " << alienURL << endl ; 
      chain->Add(alienURL) ; 
      
    }
  }// xml analysis
  
  gSystem->ChangeDirectory(ocwd.Data());
}
