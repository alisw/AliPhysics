/// \file anaQA.C
/// \ingroup CaloTrackCorrMacrosQA
/// \brief Example of execution macro for calorimeter QA analysis
///
/// Example macro to do for calorimeter QA analysis
/// in local, grid or plugin modes.
///
/// Pay attention to the options and definitions
/// set in the lines below
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

//---------------------------------------------------------------------------
/// Different analysis modes
enum anaModes
{
    mLocal  = 0, /// Analyze locally files in your computer.
    mPROOF  = 1, /// Analyze files on GRID with Plugin
    mPlugin = 2, /// Analyze files on GRID with Plugin
    mGRID   = 3  /// Analyze files on GRID, jobs launched from aliensh
};

//---------------------------------------------------------------------------
// Settings to read locally several files, only for "mLocal" mode
// The different values are default, they can be set with environmental
// variables: INDIR, PATTERN, NFILES, respectivelly

char * kInDir   = "/user/data/files/";  /// Global,  path to data files
char * kPattern = ""; /// Data are in files kInDir/kPattern+i
Int_t  kFile    = 2;  /// Number of files to analyze in local mode.

//---------------------------------------------------------------------------
// Collection file for grid analysis

char * kXML = "collection.xml"; /// Global name for the xml collection file with data on grid

const Bool_t kMC = kFALSE; /// With real data kMC = kFALSE
const TString kInputData = "ESD"; /// ESD, AOD, MC
TString kTreeName = "esdTree";

//___________________________
/// Main execution method. It:
/// * 1) loads the needed libraries in method LoadLibraries
/// * 2) The analysis frame is initialized via de analysis manager
/// * 3) Different general analysis are initialized: Physics selection, centrality etc.
/// * 4) Specialized analysis are initialized: AliAnalysistaskCounter,  AliAnalysisTaskCaloTrackCorrelations with QA analysis in it.
/// * 5) The output/input containers are passed to the analysis manager
/// * 6) The analysis is executed
///
/// \param mode: analysis mode defined in enum anaModes
//___________________________
void anaQA(Int_t mode=mLocal)
{
  // Main

  //--------------------------------------------------------------------
  // Load analysis libraries
  // Look at the method below, 
  // change whatever you need for your analysis case
  // ------------------------------------------------------------------
  LoadLibraries(mode) ;
  //    TGeoManager::Import("geometry.root") ; //need file "geometry.root" in local dir!!!!

  //-------------------------------------------------------------------------------------------------
  // Create chain from ESD and from cross sections files, look below for options.
  //-------------------------------------------------------------------------------------------------
  if      (kInputData == "ESD") kTreeName = "esdTree" ;
  else if (kInputData == "AOD") kTreeName = "aodTree" ;
  else if (kInputData == "MC" ) kTreeName = "TE" ;
  else
  {
    cout<<"Wrong  data type "<<kInputData<<endl;
    break;
  }

  TChain *chain       = new TChain(kTreeName) ;
  CreateChain(mode, chain);  

  if(chain)
  {
    AliLog::SetGlobalLogLevel(AliLog::kError);//Minimum prints on screen
    
    //--------------------------------------
    // Make the analysis manager
    //-------------------------------------
    AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");
    // MC handler
    if((kMC || kInputData == "MC") && kInputData!="AOD"){
      AliMCEventHandler* mcHandler = new AliMCEventHandler();
      mcHandler->SetReadTR(kFALSE);//Do not search TrackRef file
      mgr->SetMCtruthEventHandler(mcHandler);
      if( kInputData == "MC") mgr->SetInputEventHandler(NULL);
    }
    
    //input
    if(kInputData == "ESD")
    {
      // ESD handler
      AliESDInputHandler *esdHandler = new AliESDInputHandler();
      mgr->SetInputEventHandler(esdHandler);
      cout<<"ESD handler "<<mgr->GetInputEventHandler()<<endl;
    }
    if(kInputData == "AOD")
    {
      // AOD handler
      AliAODInputHandler *aodHandler = new AliAODInputHandler();
      mgr->SetInputEventHandler(aodHandler);
      cout<<"AOD handler "<<mgr->GetInputEventHandler()<<endl;
    }
    
     //mgr->SetDebugLevel(-1); // For debugging, do not uncomment if you want no messages.

    // AOD output handler, needed for physics selection
    cout<<"Init output handler"<<endl;
    AliAODHandler* aodoutHandler   = new AliAODHandler();
    aodoutHandler->SetOutputFileName("aod.root");
    ////aodoutHandler->SetCreateNonStandardAOD();
    mgr->SetOutputEventHandler(aodoutHandler);
    
    //-------------------------------------------------------------------------
    // Define task, put here any other task that you want to use.
    //-------------------------------------------------------------------------
    
    TString outputFile = AliAnalysisManager::GetCommonFileName(); 
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    
    if(kInputData=="ESD" && !kMC)
    {
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
    }
    
    //Counting events tasks
    if(!kMC)
    {
      AliAnalysisTaskCounter * counterMB = new AliAnalysisTaskCounter("CounterMB");
      counterMB->SelectCollisionCandidates(AliVEvent::kMB);
      
      AliAnalysisDataContainer *coutputMB =
      mgr->CreateContainer("counterMB", TList::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());
      mgr->AddTask(counterMB);
      mgr->ConnectInput  (counterMB, 0, cinput1  );
      mgr->ConnectOutput (counterMB, 1, coutputMB);
      
      AliAnalysisTaskCounter * counterEMC = new AliAnalysisTaskCounter("CounterEMC");
      counterEMC->SelectCollisionCandidates(AliVEvent::kEMC7);
      
      AliAnalysisDataContainer *coutputEMC =
      mgr->CreateContainer("counterEMC", TList::Class(), AliAnalysisManager::kOutputContainer,  outputFile.Data());
      mgr->AddTask(counterEMC);
      mgr->ConnectInput  (counterEMC, 0, cinput1   );
      mgr->ConnectOutput (counterEMC, 1, coutputEMC);
      
      AliAnalysisTaskCounter * counterINT = new AliAnalysisTaskCounter("CounterINT");
      counterINT->SelectCollisionCandidates(AliVEvent::kINT7);
      
      AliAnalysisDataContainer *coutputINT =
      mgr->CreateContainer("counterINT7", TList::Class(), AliAnalysisManager::kOutputContainer,  outputFile.Data());
      mgr->AddTask(counterINT);
      mgr->ConnectInput  (counterINT, 0, cinput1   );
      mgr->ConnectOutput (counterINT, 1, coutputINT);
    }
    else
    {
      AliAnalysisDataContainer *coutput =
      mgr->CreateContainer("counter", TList::Class(), AliAnalysisManager::kOutputContainer,  outputFile.Data());
      mgr->AddTask(counter);
      mgr->ConnectInput  (counter, 0, cinput1);
      mgr->ConnectOutput (counter, 1, coutput);
    }
    
    // QA task
    
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/QA/AddTaskCalorimeterQA.C");
    if(!kMC)
    {
      AliAnalysisTaskCaloTrackCorrelation *taskQAEMC = AddTaskCalorimeterQA("EMC",kMC,"",2012);
      taskQAEMC->SelectCollisionCandidates(AliVEvent::kEMC7);
      AliAnalysisTaskCaloTrackCorrelation *taskQAINT = AddTaskCalorimeterQA("default",kMC,"",2012);
      taskQAINT->SelectCollisionCandidates(AliVEvent::kINT7);
    }
    else
    {
      AliAnalysisTaskCaloTrackCorrelation *taskQA = AddTaskCalorimeterQA("default",kMC,"",2012);
    }
    
    //-----------------------
    // Run the analysis
    //-----------------------    
    TString smode = "";
    if (mode==mLocal || mode == mLocalCAF) 
      smode = "local";
    else if (mode==mPROOF) 
      smode = "proof";
    else if (mode==mGRID) 
      smode = "local";
    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis(smode.Data(),chain);

cout <<" Analysis ended sucessfully "<< endl ;

  }
  else cout << "Chain was not produced ! "<<endl;
  
}

//_____________________________
/// Load analysis libraries.
//_____________________________
void  LoadLibraries(const anaModes mode)
{
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libMatrix");
  gSystem->Load("libPhysics");
  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (mode==mLocal || mode == mLocalCAF || mode == mGRID) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
 
    gSystem->Load("libSTEERBase");
    gSystem->Load("libProof");
    gSystem->Load("libOADB");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libESDfilter");

    gSystem->Load("libPHOSUtils");
    gSystem->Load("libEMCALUtils");
    
    gSystem->Load("libTender");
    gSystem->Load("libTenderSupplies");
    
    gSystem->Load("libCORRFW");
    gSystem->Load("libPWGTools");
    
    gSystem->Load("libPWGCaloTrackCorrBase");
    gSystem->Load("libPWGGACaloTrackCorrelations");

 
    //--------------------------------------------------------
    //If you want to use root and par files from aliroot
    //--------------------------------------------------------  
//     SetupPar("STEERBase");
//     SetupPar("ESD");
//     SetupPar("AOD");
//     SetupPar("ANALYSIS");
//     SetupPar("ANALYSISalice");
//     //If your analysis needs PHOS geometry uncomment following lines
//     SetupPar("PHOSUtils");
//     SetupPar("EMCALUtils");
//
//     SetupPar("PWGCaloTrackCorrBase");
//     SetupPar("PWGGACaloTrackCorrelations");
  }

  //---------------------------------------------------------
  // <<<<<<<<<< PROOF mode >>>>>>>>>>>>
  //---------------------------------------------------------
  else if (mode==mPROOF)
  {
    //
    // Connect to proof
    // Put appropriate username here
    // TProof::Reset("proof://mgheata@lxb6046.cern.ch"); 
    TProof::Open("proof://mgheata@lxb6046.cern.ch");
    
    //    gProof->ClearPackages();
    //    gProof->ClearPackage("ESD");
    //    gProof->ClearPackage("AOD");
    //    gProof->ClearPackage("ANALYSIS");   
    //    gProof->ClearPackage("PWG4PartCorrBase");
    //    gProof->ClearPackage("PWG4PartCorrDep");
    
    // Enable the STEERBase Package
    gProof->UploadPackage("STEERBase.par");
    gProof->EnablePackage("STEERBase");
    // Enable the ESD Package
    gProof->UploadPackage("ESD.par");
    gProof->EnablePackage("ESD");
    // Enable the AOD Package
    gProof->UploadPackage("AOD.par");
    gProof->EnablePackage("AOD");
    // Enable the Analysis Package
    gProof->UploadPackage("ANALYSIS.par");
    gProof->EnablePackage("ANALYSIS");
	// Enable the PHOS geometry Package
    //gProof->UploadPackage("PHOSUtils.par");
    //gProof->EnablePackage("PHOSUtils");
    // Enable PartCorr analysis
    gProof->UploadPackage("PWG4PartCorrBase.par");
    gProof->EnablePackage("PWG4PartCorrBase");
	gProof->UploadPackage("PWG4PartCorrDep.par");
    gProof->EnablePackage("PWG4PartCorrDep");    
    gProof->ShowEnabledPackages();
  }
}

//_________________________________
/// Load par files, create analysis libraries
/// For testing, if par file already decompressed and modified
/// classes then do not decompress.
//_________________________________
void SetupPar(char* pararchivename)
{
  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ;
    
  TString parpar(Form("%s.par", pararchivename)) ;
    
  if ( gSystem->AccessPathName(parpar.Data()) )
  {
    gSystem->ChangeDirectory(gSystem->Getenv("ALICE_PHYSICS")) ;
    TString processline(Form(".! make %s", parpar.Data())) ; 
    gROOT->ProcessLine(processline.Data()) ;
    gSystem->ChangeDirectory(cdir) ; 
    processline = Form(".! mv $ALICE_PHYSICS/%s .", parpar.Data()) ;
    gROOT->ProcessLine(processline.Data()) ;
  }
    
  if ( gSystem->AccessPathName(pararchivename) )
  {
    TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }
  
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(pararchivename);
  
  // check for BUILD.sh and execute
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh"))
  {
    printf("*******************************\n");
    printf("*** Building PAR archive    ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    
    if (gSystem->Exec("PROOF-INF/BUILD.sh"))
    {
      Error("runProcess","Cannot Build the PAR Archive! - Abort!");
      return -1;
    }
  }
  // check for SETUP.C and execute
  if (!gSystem->AccessPathName("PROOF-INF/SETUP.C"))
  {
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
void CreateChain(const anaModes mode, TChain * chain)
{
  TString ocwd = gSystem->WorkingDirectory();
  
  //-----------------------------------------------------------
  // Analysis of CAF data locally and with PROOF
  //-----------------------------------------------------------
  if(mode ==mPROOF || mode ==mLocalCAF)
  {
    // Chain from CAF
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/CreateESDChain.C");
    // The second parameter is the number of input files in the chain
    chain = CreateESDChain("ESD12001.txt", 5);  
  }
  
  //---------------------------------------
  // Local files analysis
  //---------------------------------------
  else if(mode == mLocal)
  {
    //If you want to add several ESD files sitting in a common directory INDIR
    //Specify as environmental variables the directory (INDIR), the number of files 
    //to analyze (NFILES) and the pattern name of the directories with files (PATTERN)

    if(gSystem->Getenv("INDIR"))  
      kInDir = gSystem->Getenv("INDIR") ; 
    else cout<<"INDIR not set, use default: "<<kInDir<<endl;	
    
    if(gSystem->Getenv("PATTERN"))   
      kPattern = gSystem->Getenv("PATTERN") ; 
    else  cout<<"PATTERN not set, use default: "<<kPattern<<endl;
    
    if(gSystem->Getenv("NFILES"))
      kFile = atoi(gSystem->Getenv("NFILES")) ;
    else cout<<"NFILES not set, use default: "<<kFile<<endl;
    
    //Check if env variables are set and are correct
    if ( kInDir  && kFile) {
      printf("Get %d files from directory %s\n",kFile,kInDir);
      if ( ! gSystem->cd(kInDir) ) {//check if ESDs directory exist
	printf("%s does not exist\n", kInDir) ;
	return ;
      }

      cout<<"INDIR   : "<<kInDir<<endl;
      cout<<"NFILES  : "<<kFile<<endl;
      cout<<"PATTERN : " <<kPattern<<endl;

      TString datafile="";
      if(kInputData == "ESD") datafile = "AliESDs.root" ;
      else if(kInputData == "AOD") datafile = "AliAOD.root" ;
      
      //Loop on ESD files, add them to chain
      Int_t event =0;
      Int_t skipped=0 ; 
      char file[120] ;
      
      for (event = 0 ; event < kFile ; event++) {
	sprintf(file, "%s/%s%d/%s", kInDir,kPattern,event,datafile.Data()) ; 
	TFile * fESD = 0 ; 
	//Check if file exists and add it, if not skip it
	if ( fESD = TFile::Open(file)) {
	  if ( fESD->Get(kTreeName) ) { 
	    printf("++++ Adding %s\n", file) ;
	    chain->AddFile(file);
	  }
	}
	else { 
	  printf("---- Skipping %s\n", file) ;
	  skipped++ ;
	}
      }
      printf("number of entries # %lld, skipped %d\n", chain->GetEntries(), skipped*100) ; 	
    }
    else {
      TString input = "AliESDs.root" ;
      cout<<">>>>>> No list added, take a single file <<<<<<<<< "<<input<<endl;
      chain->AddFile(input);
    }
    
  }// local files analysis
  
  //------------------------------
  //GRID xml files
  //-----------------------------
  else if(mode == mGRID){
    //Get colection file. It is specified by the environmental
    //variable XML

    if(gSystem->Getenv("XML") )
      kXML = gSystem->Getenv("XML");
    else
      sprintf(kXML, "collection.xml") ; 
    
    if (!TFile::Open(kXML)) {
      printf("No collection file with name -- %s -- was found\n",kXML);
      return ;
    }
    else cout<<"XML file "<<kXML<<endl;

    //Load necessary libraries and connect to the GRID
    gSystem->Load("libNetx") ; 
    gSystem->Load("libRAliEn"); 
    TGrid::Connect("alien://") ;

    //Feed Grid with collection file
    TGridCollection * collection = (TGridCollection*) TAlienCollection::Open(kXML);
    if (! collection) {
      AliError(Form("%s not found", kXML)) ; 
      return kFALSE ; 
    }
    TGridResult* result = collection->GetGridResult("",0 ,0);
   
    // Makes the ESD chain 
    printf("*** Getting the Chain       ***\n");
    for (Int_t index = 0; index < result->GetEntries(); index++) {
      TString alienURL = result->GetKey(index, "turl") ; 
      cout << "================== " << alienURL << endl ; 
      chain->Add(alienURL) ; 
    }
  }// xml analysis
  
  gSystem->ChangeDirectory(ocwd.Data());
}

