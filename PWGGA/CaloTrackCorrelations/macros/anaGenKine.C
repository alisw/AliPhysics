/* $Id:  $ */
//--------------------------------------------------
// Example macro to do analysis with the 
// analysis classes in CaloTrackCorrelations
// Can be executed with Root and AliRoot
//
// Pay attention to the options and definitions
// set in the lines below
//
//  Author : Gustavo Conesa Balbastre (INFN-LNF)
//
//-------------------------------------------------
enum anaModes {mLocal=0, mPROOF=1, mPlugin=2, mGRID=3};
//mLocal    = 0: Analyze locally files in your computer
//mPROOF    = 1: Analyze files on GRID with Plugin
//mPlugin   = 2: Analyze files on GRID with Plugin
//mGRID     = 3: Analyze files on GRID, jobs launched from aliensh

//---------------------------------------------------------------------------
// Settings to read locally several files, only for "mLocal" mode
// The different values are default, they can be set with environmental 
// variables: INDIR, PATTERN, NFILES, respectivelly

char * kInDir   = "/user/data/files/"; 
char * kPattern = ""; // Data are in files kInDir/kPattern+i 
Int_t  kFile    = 6; 

//---------------------------------------------------------------------------
// Dataset for proof analysis, mode=mPROOF
// char * kDataset = "/alice/vernet/PbPb_LHC10h_ESD";

char *  kDatasetPROOF     = "/alice/vernet/LHC11b_149646";
Int_t   kDatasetNMaxFiles = 20;
TString ccin2p3UserName   = "arbor" ;
TString alienUserName     = "narbor" ;

//---------------------------------------------------------------------------
// Collection file for grid analysis

char * kXML = "collection.xml";

//---------------------------------------------------------------------------
//Scale histograms from file. Change to kTRUE when xsection file exists
//Put name of file containing xsection 
//Put number of events per ESD file
//This is an specific case for normalization of Pythia files.
const char * kXSFileName = "pyxsec.root";

//---------------------------------------------------------------------------

//Set some default values, but used values are set in the code!

Bool_t  kMC        = kFALSE; //With real data kMC = kFALSE
TString kInputData = "ESD"; //ESD, AOD, MC, deltaAOD
Int_t   kYear      = 2011;
TString kCollision = "pp";
Bool_t  outAOD     = kFALSE; //Some tasks doesnt need it.
TString kTreeName;
TString kPass      = "";
char    kTrigger[1024];
Int_t   kRun       = 0;

//________________________
void anaGenKine(Int_t mode=mGRID)
{
  // Main
  
  //--------------------------------------------------------------------
  // Load analysis libraries
  // Look at the method below, 
  // change whatever you need for your analysis case
  // ------------------------------------------------------------------
  
  LoadLibraries(mode) ;
  //gSystem->ListLibraries();
  
  //-------------------------------------------------------------------------------------------------
  //Create chain from ESD and from cross sections files, look below for options.
  //-------------------------------------------------------------------------------------------------
  
  // Set kInputData and kTreeName looking to the kINDIR
  
  //CheckInputData(mode);
  
  // Check global analysis settings  
  
  //CheckEnvironmentVariables();
  
  //printf("*********************************************\n");
  //printf("*** Input data < %s >, pass %s, tree < %s >, MC?  < %d > ***\n",kInputData.Data(),kPass.Data(),kTreeName.Data(),kMC);
  
  //printf("*********************************************\n");
  kTreeName  = "TE";
  kInputData = "MC";
  TChain * chain   = new TChain(kTreeName) ;
  TChain * chainxs = new TChain("Xsection") ;
  
  chain  ->AddFile("galice.root");
  chainxs->AddFile("pyxsec.root");

  //CreateChain(mode, chain, chainxs); 
  
  Double_t scale = -1;
  printf("===== kMC %d, chainxs %p\n",kMC,chainxs);
  if(chainxs && chainxs->GetEntries() > 0)
  {
    Int_t nfiles = chainxs->GetEntries();
    
    //Get the cross section
    Double_t xsection = 0; 
    Float_t ntrials   = 0;
    
    GetAverageXsection(chainxs, xsection, ntrials);
    
    Int_t    nEventsPerFile = chain->GetEntries() / nfiles;
    
    Double_t trials = ntrials / nEventsPerFile ;
    
    scale = xsection/trials;
    
    printf("Get Cross section : nfiles  %d, nevents %d, nevents per file %d \n",nfiles, chain->GetEntries(),nEventsPerFile);     
    printf("                    ntrials %d, trials %2.2f, xs %2.2e, scale factor %2.2e\n", ntrials,trials,xsection,scale);
    
  } 
  
  printf("*********************************************\n");
  printf("number of entries # %lld, skipped %d\n", chain->GetEntries()) ; 	
  printf("*********************************************\n");
  
  if(!chain)
  { 
    printf("STOP, no chain available\n"); 
    return;
  }
  
  AliLog::SetGlobalLogLevel(AliLog::kError);//Minimum prints on screen
  
  //------------------------------------------
  //  Alien handler part
  //------------------------------------------
  AliAnalysisGrid *alienHandler=0x0;
  if(mode==mPlugin)
  {
    // Create and configure the alien handler plugin
    gROOT->LoadMacro("CreateAlienHandler.C");
    alienHandler = CreateAlienHandler();
    if (!alienHandler) return;
  }  
  
  //--------------------------------------
  // Make the analysis manager
  //-------------------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");
  //AliAnalysisManager::SetUseProgressBar(kTRUE);
  //mgr->SetSkipTerminate(kTRUE);
  //mgr->SetNSysInfo(1);
  
//  if(mode==mPlugin)
//  {
//    // Connect plugin to the analysis manager
//    mgr->SetGridHandler(alienHandler);
//  }
  
  // MC handler
  if((kMC || kInputData == "MC") && !kInputData.Contains("AOD"))
  {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(kFALSE);//Do not search TrackRef file
    mgr->SetMCtruthEventHandler(mcHandler);
    if( kInputData == "MC") 
    {
      cout<<"MC INPUT EVENT HANDLER"<<endl;
      mgr->SetInputEventHandler(NULL);
    }
  }
  
  
//  // AOD output handler
//  if(kInputData!="deltaAOD" && outAOD)
//  {
//    cout<<"Init output handler"<<endl;
//    AliAODHandler* aodoutHandler   = new AliAODHandler();
//    aodoutHandler->SetOutputFileName("aod.root");
//    ////aodoutHandler->SetCreateNonStandardAOD();
//    mgr->SetOutputEventHandler(aodoutHandler);
//  }
//  
//  //input
//  
//  if(kInputData == "ESD")
//  {
//    // ESD handler
//    AliESDInputHandler *esdHandler = new AliESDInputHandler();
//    esdHandler->SetReadFriends(kFALSE);
//    mgr->SetInputEventHandler(esdHandler);
//    cout<<"ESD handler "<<mgr->GetInputEventHandler()<<endl;
//  }
//  else if(kInputData.Contains("AOD"))
//  {
//    // AOD handler
//    AliAODInputHandler *aodHandler = new AliAODInputHandler();
//    mgr->SetInputEventHandler(aodHandler);
//    if(kInputData == "deltaAOD") aodHandler->AddFriend("deltaAODCaloTrackCorr.root");
//    cout<<"AOD handler "<<mgr->GetInputEventHandler()<<endl;
//  }
  
  //mgr->RegisterExternalFile("deltaAODCaloTrackCorr.root");
  mgr->SetDebugLevel(1); // For debugging, do not uncomment if you want no messages.
  
  TString outputFile = AliAnalysisManager::GetCommonFileName(); 
   
  gROOT->LoadMacro("AddTaskGenKine.C");  
  
 
  AliAnalysisTaskCaloTrackCorrelation *ana   = AddTaskGenKine(outputFile, scale);
    
  
  //-----------------------
  // Run the analysis
  //-----------------------    
  mgr->InitAnalysis();
  mgr->PrintStatus();
  
  if      (mode == mPlugin) mgr->StartAnalysis("grid");
  else if (mode == mPROOF ) mgr->StartAnalysis("proof",chain);
  else                      mgr->StartAnalysis("local",chain);
  
  cout <<" Analysis ended sucessfully "<< endl ;
  
}

//_____________________________
void  LoadLibraries(Int_t mode)
{
  
  if (mode == mPROOF) {
    //TProof::Mgr("ccalpmaster")->SetROOTVersion("ALICE_v5-27-06b");
    gROOT->LoadMacro("/afs/in2p3.fr/group/alice/laf/EnableAliRootForLAF.C");
    TProof* proof = EnableAliRootForLAF("ccaplmaster",nPROOFWorkers.Data(),ccin2p3UserName.Data(),alienUserName.Data(),"",kFALSE,kTRUE,kTRUE,"OADB:ANALYSIS:ANALYSISalice:AOD:ESD:CORRFW:STEERBase:EMCALUtils:PHOSUtils:PWGCaloTrackCorrBase:PWGGACaloTrackCorrelations");
    
    //  TProof* proof = TProof::Open("ccaplmaster",Form("workers=%s",nPROOFWorkers.Data()));
    
    //     //proof->ClearPackages();
    //     proof->UploadPackage("STEERBase");
    //     proof->UploadPackage("ESD");
    //     proof->UploadPackage("AOD");
    //     proof->UploadPackage("ANALYSIS");
    //     proof->UploadPackage("OADB");
    //     proof->UploadPackage("ANALYSISalice");
    //     proof->UploadPackage("CORRFW");
    //     //proof->UploadPackage("JETAN");
    //     proof->UploadPackage("PHOSUtils");
    //     proof->UploadPackage("EMCALUtils");
    //     proof->UploadPackage("PWGCaloTrackCorrBase");
    //     proof->UploadPackage("PWGGACaloTrackCorrelations");

    //     proof->EnablePackage("STEERBase");
    //     proof->EnablePackage("ESD");
    //     proof->EnablePackage("AOD");
    //     proof->EnablePackage("ANALYSIS");
    //     proof->EnablePackage("OADB");
    //     proof->EnablePackage("ANALYSISalice");
    //     proof->EnablePackage("CORRFW");
    //     //proof->EnablePackage("JETAN");
    //     proof->EnablePackage("PHOSUtils");
    //     proof->EnablePackage("EMCALUtils");
    //     proof->EnablePackage("PWGCaloTrackCorrBase");
    //     proof->EnablePackage("PWGGACaloTrackCorrelations");

    return;
  }  
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libMatrix");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit"); // Root + libraries to if reclusterization is done
  
  gSystem->Load("libSTEERBase");
  gSystem->Load("libGui"); // Root + libraries to if reclusterization is done
  gSystem->Load("libCDB"); // Root + libraries to if reclusterization is done
  gSystem->Load("libESD"); // Root + libraries to if reclusterization is done
  gSystem->Load("libAOD");
  gSystem->Load("libRAWDatabase"); // Root + libraries to if reclusterization is done
  gSystem->Load("libProof"); 
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTEER"); // Root + libraries to if reclusterization is done
  
  gSystem->Load("libRAWDatarec"); // Root + libraries to if reclusterization is done
  gSystem->Load("libRAWDatasim"); // Root + libraries to if reclusterization is done
  gSystem->Load("libVZERObase");  // Root + libraries to if reclusterization is done
  gSystem->Load("libVZEROrec");   // Root + libraries to if reclusterization is done
  
  gSystem->Load("libEMCALUtils");
  //SetupPar("EMCALUtils");
  gSystem->Load("libEMCALraw");  // Root + libraries to if reclusterization is done
  gSystem->Load("libEMCALbase"); // Root + libraries to if reclusterization is done
  gSystem->Load("libEMCALsim");  // Root + libraries to if reclusterization is done
  gSystem->Load("libEMCALrec");  // Root + libraries to if reclusterization is done
  //SetupPar("EMCALraw");
  //SetupPar("EMCALbase");
  //SetupPar("EMCALsim");
  //SetupPar("EMCALrec");
  
  gSystem->Load("libANALYSISalice");
  //gSystem->Load("libTender"); 
  //gSystem->Load("libTenderSupplies");
  
  gSystem->Load("libPHOSUtils");
  gSystem->Load("libEMCALUtils");
  gSystem->Load("libPWGCaloTrackCorrBase");
  gSystem->Load("libPWGGACaloTrackCorrelations");
  //SetupPar("PWGCaloTrackCorrBase");
  //SetupPar("PWGGACaloTrackCorrelations");
  
 
  //gSystem->Load("libJETAN");
  //gSystem->Load("FASTJETAN");
  //gSystem->Load("PWGJE");

  //gSystem->Load("libCORRFW");
  //gSystem->Load("libPWGGAGammaConv"); 
  //SetupPar("PWGGAGammaConv"); 
  
  // needed for plugin?
  gSystem->AddIncludePath("-I$ALICE_ROOT");
  gSystem->AddIncludePath("-I./");     
  
}

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

//______________________________________
void CheckInputData(const anaModes mode)
{
  //Sets input data and tree
  
  TString ocwd = gSystem->WorkingDirectory();
  
  //---------------------------------------
  //Local files analysis
  //---------------------------------------
  if(mode == mLocal){    
    //If you want to add several ESD files sitting in a common directory INDIR
    //Specify as environmental variables the directory (INDIR), the number of files 
    //to analyze (NFILES) and the pattern name of the directories with files (PATTERN)
    
    if(gSystem->Getenv("INDIR"))  
      kInDir = gSystem->Getenv("INDIR") ; 
    else cout<<"INDIR not set, use default: "<<kInDir<<endl;	
    
    TString sindir(kInDir);
    if     (sindir.Contains("pass1")) kPass = "pass1";
    else if(sindir.Contains("pass2")) kPass = "pass2";
    else if(sindir.Contains("pass3")) kPass = "pass3";
    
    if(gSystem->Getenv("PATTERN"))   
      kPattern = gSystem->Getenv("PATTERN") ; 
    else  cout<<"PATTERN not set, use default: "<<kPattern<<endl;
    
    cout<<"INDIR   : "<<kInDir<<endl;
    cout<<"NFILES  : "<<kFile<<endl;
    
    char fileE[120] ;   
    char fileA[120] ;   
    char fileG[120] ;
    char fileEm[120] ;   
    for (Int_t event = 0 ; event < kFile ; event++) {
      sprintf(fileE,  "%s/%s%d/AliESDs.root",    kInDir,kPattern,event) ; 
      sprintf(fileA,  "%s/%s%d/AliAOD.root",     kInDir,kPattern,event) ; 
      sprintf(fileG,  "%s/%s%d/galice.root",     kInDir,kPattern,event) ; 
      sprintf(fileEm, "%s/%s%d/embededAOD.root", kInDir,kPattern,event) ; 
      
      TFile * fESD = TFile::Open(fileE) ; 
      TFile * fAOD = TFile::Open(fileA) ; 
      
      //Check if file exists and add it, if not skip it
      if (fESD) 
      {
        kTreeName  = "esdTree";
        kInputData = "ESD";
        TFile * fG = TFile::Open(fileG);
        if(fG) { kMC = kTRUE; fG->Close();}
        else     kMC = kFALSE;
        
        // Get run number
        TTree* esdTree = (TTree*)fESD->Get("esdTree");
        AliESDEvent* esd = new AliESDEvent();
        esd->ReadFromTree(esdTree);
        esdTree->GetEvent(0);
        kRun = esd->GetRunNumber();
        
        return;
      }
      else if(fAOD)
      {
        kTreeName  = "aodTree";
        kInputData = "AOD";
        if(((TTree*) fAOD->Get("aodTree"))->GetBranch("mcparticles")) kMC=kTRUE;
        else kMC = kFALSE;
        
        // Get run number
        TTree* aodTree = (TTree*)fAOD->Get("aodTree");
        AliAODEvent* aod = new AliAODEvent();
        aod->ReadFromTree(aodTree);
        aodTree->GetEvent(0);
        kRun = aod->GetRunNumber();
        return;
      }
      else if(TFile::Open(fileEm))
      {
        kTreeName  = "aodTree";
        kInputData = "AOD";
        kMC        = kTRUE;
        
        return;
      }
      else if(TFile::Open(fileG))
      {
        kTreeName  = "TE";
        kInputData = "MC";
        kMC        = kTRUE;
        return;
      }
    }
    
    if(fESD) fESD->Close();
    if(fAOD) fAOD->Close();
    
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
    
    for (Int_t index = 0; index < result->GetEntries(); index++) {
      TString alienURL = result->GetKey(index, "turl") ; 
      cout << "================== " << alienURL << endl ; 
      
      if     (alienURL.Contains("pass1")) kPass = "pass1";
      else if(alienURL.Contains("pass2")) kPass = "pass2";
      else if(alienURL.Contains("pass3")) kPass = "pass3";
      
      kRun = AliAnalysisManager::GetRunFromAlienPath(alienURL.Data());
      printf("Run number from alien path = %d\n",kRun);
      
      TFile * fAOD = 0 ; 
      //Check if file exists and add it, if not skip it
      if (alienURL.Contains("AliESDs.root"))  
      {
        kTreeName  = "esdTree";
        kInputData = "ESD";
        alienURL.ReplaceAll("AliESDs.root","galice.root");
        if(TFile::Open(alienURL)) kMC=kTRUE;
        else kMC = kFALSE;
        return;
      }
      else if(alienURL.Contains("AliAOD.root"))
      {
        kTreeName  = "aodTree";
        kInputData = "AOD";
        fAOD = TFile::Open(alienURL);
        if(((TTree*) fAOD->Get("aodTree"))->GetBranch("mcparticles")) kMC=kTRUE;
        else kMC = kFALSE;
        return;
      }
      else if(alienURL.Contains("embededAOD.root"))
      {
        kTreeName  = "aodTree";
        kInputData = "AOD";
        kMC=kTRUE;
        return;
      }
      else if(alienURL.Contains("galice.root"))
      {
        kTreeName  = "TE";
        kInputData = "MC";
        kMC=kTRUE;
        return;
      } 
    }
  }// xml analysis
  //------------------------------
  //PROOF files
  //-----------------------------
  else if(mode == mPROOF){
    
    TFileCollection* coll  = gProof->GetDataSet(kDatasetPROOF)->GetStagedSubset();
    
    TIter iter(coll->GetList());
    
    TFileInfo* fileInfo = 0;
    while ((fileInfo = dynamic_cast<TFileInfo*> (iter())))
    {
      if (fileInfo->GetFirstUrl()) {
        TString ProofURL = fileInfo->GetFirstUrl()->GetUrl();
        cout << "================== " << ProofURL << endl ; 
        
        if     (ProofURL.Contains("pass1")) kPass = "pass1";
        else if(ProofURL.Contains("pass2")) kPass = "pass2";
        else if(ProofURL.Contains("pass3")) kPass = "pass3";
        
        kRun = AliAnalysisManager::GetRunFromAlienPath(ProofURL.Data());
        printf("Run number from alien path = %d\n",kRun);
        
        TFile * fAOD = 0 ; 
        //Check if file exists and add it, if not skip it
        if (ProofURL.Contains("AliESDs.root"))  
        {
          kTreeName  = "esdTree";
          kInputData = "ESD";
          alienURL.ReplaceAll("AliESDs.root","galice.root");
          if(TFile::Open(ProofURL)) kMC=kTRUE;
          else kMC = kFALSE;
          
          return;
        }
        else if(ProofURL.Contains("AliAOD.root"))
        {
          kTreeName  = "aodTree";
          kInputData = "AOD";
          fAOD = TFile::Open(ProofURL);
          if(((TTree*) fAOD->Get("aodTree"))->GetBranch("mcparticles")) kMC=kTRUE;
          else kMC = kFALSE;
          return;
        }
        else if(ProofURL.Contains("embededAOD.root"))
        {
          kTreeName  = "aodTree";
          kInputData = "AOD";
          kMC=kTRUE;
          return;
        }
        else if(ProofURL.Contains("galice.root"))
        {
          kTreeName  = "TE";
          kInputData = "MC";
          kMC=kTRUE;
          return;
        } 
      }
    }
  }// proof analysis
  
  gSystem->ChangeDirectory(ocwd.Data());
  
}

//_____________________________________________________________________
void CreateChain(const anaModes mode, TChain * chain, TChain * chainxs)
{
  //Fills chain with data
  TString ocwd = gSystem->WorkingDirectory();
  
  //---------------------------------------
  // Local files analysis
  //---------------------------------------
  if(mode == mLocal){    
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
      
      //if(gSystem->Getenv("XSFILE"))  
      //kXSFileName = gSystem->Getenv("XSFILE") ; 
      //else cout<<" XS file name not set, use default: "<<kXSFileName<<endl;	
      char * kGener = gSystem->Getenv("GENER");
      if(kGener) {
        cout<<"GENER "<<kGener<<endl;
        if     (!strcmp(kGener,"PYTHIA")) kXSFileName = "pyxsec.root";
        else if(!strcmp(kGener,"HERWIG")) kXSFileName = "hexsec.root";
        else cout<<" UNKNOWN GENER, use default: "<<kXSFileName<<endl;
      }
      else cout<<" GENER not set, use default xs file name: "<<kXSFileName<<endl;
      
      cout<<"INDIR   : "<<kInDir     <<endl;
      cout<<"NFILES  : "<<kFile      <<endl;
      cout<<"PATTERN : "<<kPattern   <<endl;
      cout<<"XSFILE  : "<<kXSFileName<<endl;
      
      TString datafile="";
      if     (kInputData == "ESD")        datafile = "AliESDs.root" ;
      else if(kInputData.Contains("AOD")) datafile = "AliAOD.root"  ;
      else if(kInputData == "MC")         datafile = "galice.root"  ;
      
      //Loop on ESD/AOD/MC files, add them to chain
      Int_t event =0;
      Int_t skipped=0 ; 
      char file[120] ;
      char filexs[120] ;
      
      for (event = 0 ; event < kFile ; event++) {
        sprintf(file,   "%s/%s%d/%s", kInDir,kPattern,event,datafile.Data()) ; 
        sprintf(filexs, "%s/%s%d/%s", kInDir,kPattern,event,kXSFileName) ; 
        TFile * fData = 0 ; 
        //Check if file exists and add it, if not skip it
        if ( fData = TFile::Open(file)) {
          if ( fData->Get(kTreeName) ) { 
            printf("++++ Adding %s\n", file) ;
            chain->AddFile(file);
            chainxs->Add(filexs) ; 
          }
        }
        else { 
          printf("---- Skipping %s\n", file) ;
          skipped++ ;
        }
      }
    }
    else {
      TString input = "AliESDs.root" ;
      cout<<">>>>>> No list added, take a single file <<<<<<<<< "<<input<<endl;
      chain->AddFile(input);
    }
    
  }// local files analysis
  
  //------------------------------
  // GRID xml files
  //------------------------------
  else if(mode == mGRID){
    //Get colection file. It is specified by the environmental
    //variable XML
    
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
      alienURL.ReplaceAll("AliESDs.root",kXSFileName);
      chainxs->Add(alienURL) ; 
    }
  }// xml analysis
  
  //------------------------------
  // PROOF
  //------------------------------
  else if (mode == mPROOF) {
    
    TFileCollection* ds= gProof->GetDataSet(kDatasetPROOF)->GetStagedSubset();
    
    gROOT->LoadMacro("/afs/in2p3.fr/group/alice/laf/dataset_management/CreateChainFromDataSet.C");
    chain = CreateChainFromDataSet(ds, kTreeName , kDatasetNMaxFiles);
    printf("chain has %d entries\n",chain->GetEntries());
  }
  
  gSystem->ChangeDirectory(ocwd.Data());
  
}

//______________________________
void CheckEnvironmentVariables()
{
  
  sprintf(kTrigger,"");
  
  Bool_t bRecalibrate = kFALSE;
  Bool_t bBadChannel = kFALSE;
  
  for (int i=0; i< gApplication->Argc();i++){
    
#ifdef VERBOSEARGS
    
    printf("Arg  %d:  %s\n",i,gApplication->Argv(i));
    
#endif
    
    TString sRun = "";
    
    if (!(strcmp(gApplication->Argv(i),"--trigger")))
      sprintf(trigger,gApplication->Argv(i+1));
    
    if (!(strcmp(gApplication->Argv(i),"--recalibrate")))
      bRecalibrate = atoi(gApplication->Argv(i+1));
    
    if (!(strcmp(gApplication->Argv(i),"--badchannel")))
      bBadChannel = atoi(gApplication->Argv(i+1));
    
    if (!(strcmp(gApplication->Argv(i),"--run"))){
      sRun = gApplication->Argv(i+1);
      if(sRun.Contains("LHC10")) {
        kYear = 2010;
      }
      else {
        if(kRun <=0){
          kRun = atoi(gApplication->Argv(i+1));
        }
        else printf("** Run number already set  to %d, do not set to %d\n",kRun,atoi(gApplication->Argv(i+1)));
      }//numeric run
    }//--run available
    
  }// args loop
  
  if(!sRun.Contains("LHC10")){
    if     ( kRun < 140000) 
    {
      kYear = 2010;
      if( kRun >= 136851 ) kCollision = "PbPb";
    }
    else if( kRun < 170600)
    {
      kYear = 2011;
      if( kRun >= 166500 ) kCollision = "PbPb";
    }
    else 
    {
      kYear = 2012;

    }
  }
  
  if(kMC) sprintf(kTrigger,"");
  
  printf("*********************************************\n");
  //printf("*** Settings trigger %s, recalibrate %d, remove bad channels %d, year %d, collision %s, run %d ***\n",
  //       kTrigger,bRecalibrate,bBadChannel, kYear,kCollision.Data(), kRun);
  printf("*** Settings year %d, collision %s, run %d ***\n",kYear,kCollision.Data(), kRun);
  printf("*********************************************\n");
  
}


//_________________________________________________________________
void GetAverageXsection(TTree * tree, Double_t & xs, Float_t & ntr)
{
  // Read the PYTHIA statistics from the file pyxsec.root created by
  // the function WriteXsection():
  // integrated cross section (xsection) and
  // the  number of Pyevent() calls (ntrials)
  // and calculate the weight per one event xsection/ntrials
  // The spectrum calculated by a user should be
  // multiplied by this weight, something like this:
  // TH1F *userSpectrum ... // book and fill the spectrum
  // userSpectrum->Scale(weight)
  //
  // Yuri Kharlov 19 June 2007
  // Gustavo Conesa 15 April 2008
  Double_t xsection = 0;
  UInt_t    ntrials = 0;
  xs = 0;
  ntr = 0;
  
  Int_t nfiles =  tree->GetEntries()  ;
  if (tree && nfiles > 0) {
    tree->SetBranchAddress("xsection",&xsection);
    tree->SetBranchAddress("ntrials" ,&ntrials );
    for(Int_t i = 0; i < nfiles; i++){
      tree->GetEntry(i);
      xs  += xsection ;
      ntr += ntrials ;
      cout << "xsection " <<xsection<<" ntrials "<<ntrials<<endl; 
    }
    
    xs =   xs /  nfiles;
    ntr =  ntr / nfiles;
    cout << "-----------------------------------------------------------------"<<endl;
    cout << "Average of "<< nfiles<<" files: xsection " <<xs<<" ntrials "<<ntr<<endl; 
    cout << "-----------------------------------------------------------------"<<endl;
  } 
  else cout << " >>>> Empty tree !!!! <<<<< "<<endl;
  
}



