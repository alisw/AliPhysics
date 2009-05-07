/* $Id:  $ */

//--------------------------------------------------
// Example macro to do analysis with the 
// analysis classes in PWG4PartCorr
// and JETAN at the same time
// Can be executed with Root and AliRoot
//
// Pay attention to the options and definitions
// set in the lines below
//
//  Author : Gustavo Conesa Balbastre (INFN-LNF)
//
//-------------------------------------------------
enum anaModes {mLocal, mLocalCAF,mPROOF,mGRID};
//mLocal: Analyze locally files in your computer
//mLocalCAF: Analyze locally CAF files
//mPROOF: Analyze CAF files with PROOF

//---------------------------------------------------------------------------
//Settings to read locally several files, only for "mLocal" mode
//The different values are default, they can be set with environmental 
//variables: INDIR, PATTERN, NEVENT, respectivelly
char * kInDir = "/user/data/files/"; 
char * kPattern = ""; // Data are in files kInDir/kPattern+i 
Int_t kEvent = 1; // Number of files
//---------------------------------------------------------------------------
//Collection file for grid analysis
char * kXML = "collection.xml";
//---------------------------------------------------------------------------
//Scale histograms from file. Change to kTRUE when xsection file exists
//Put name of file containing xsection 
//Put number of events per ESD file
//This is an specific case for normalization of Pythia files.
const Bool_t kGetXSectionFromFileAndScale = kFALSE ;
const char * kXSFileName = "pyxsec.root";
const Int_t kNumberOfEventsPerFile = 100; 
//---------------------------------------------------------------------------

const Bool_t kMC = kTRUE; //With real data kMC = kFALSE
const TString kInputData = "ESD";//ESD, AOD, MC
TString kTreeName = "esdTree";

void anaPartCorrJetAn(Int_t mode=mLocal, TString configName = "ConfigAnalysisGammaJetFinderCorrelation")
{
  // Main

  //--------------------------------------------------------------------
  // Load analysis libraries
  // Look at the method below, 
  // change whatever you need for your analysis case
  // ------------------------------------------------------------------
  LoadLibraries(mode) ;
  
  //-------------------------------------------------------------------------------------------------
  //Create chain from ESD and from cross sections files, look below for options.
  //------------------------------------------------------------------------------------------------- 
  if(kInputData == "ESD") kTreeName = "esdTree" ;
  else if(kInputData == "AOD") kTreeName = "aodTree" ;
  else if (kInputData == "MC") kTreeName = "TE" ;
  else {
    cout<<"Wrong  data type "<<kInputData<<endl;
    break;
  }
	
  TChain *chain   = new TChain(kTreeName) ;
  TChain * chainxs = new TChain("Xsection") ;
  CreateChain(mode, chain, chainxs);  

  if(chain){
    AliLog::SetGlobalLogLevel(AliLog::kError);//Minimum prints on screen
    
    //--------------------------------------
    // Make the analysis manager
    //-------------------------------------
    AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");
    // MC handler
    if(kMC || kInputData == "MC"){
      AliMCEventHandler* mcHandler = new AliMCEventHandler();
      mcHandler->SetReadTR(kFALSE);//Do not search TrackRef file
      mgr->SetMCtruthEventHandler(mcHandler); 
	  if( kInputData == "MC") mgr->SetInputEventHandler(NULL);
    }

    // AOD output handler
    AliAODHandler* aodoutHandler   = new AliAODHandler();
    aodoutHandler->SetOutputFileName("aod.root");
    //aodoutHandler->SetCreateNonStandardAOD();
    mgr->SetOutputEventHandler(aodoutHandler);
    
    //input
    if(kInputData == "ESD"){
      // ESD handler
      AliESDInputHandler *esdHandler = new AliESDInputHandler();
      mgr->SetInputEventHandler(esdHandler);
    }
    if(kInputData == "AOD"){
      // AOD handler
      AliAODInputHandler *aodHandler = new AliAODInputHandler();
      mgr->SetInputEventHandler(aodHandler);
    }

   //mgr->SetDebugLevel(-1); // For debugging, uncomment for lots of messages

    //-------------------------------------------------------------------------
    //Define task, put here any other task that you want to use.
    //-------------------------------------------------------------------------
	
    AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis");
    //jetana->SetDebugLevel(-1);
    mgr->AddTask(jetana);
	
    AliAnalysisTaskParticleCorrelation * taskpwg4 = new AliAnalysisTaskParticleCorrelation ("Particle");
    taskpwg4->SetConfigFileName(configName); //Default name is ConfigAnalysis
	
    mgr->AddTask(taskpwg4);
    
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("corrhistos", TList::Class(),
							      AliAnalysisManager::kOutputContainer, "histos.root");
    
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("jethistos", TList::Class(),
							      AliAnalysisManager::kOutputContainer, "histos.root");
	
    mgr->ConnectInput  (taskpwg4,     0, cinput1);
    mgr->ConnectOutput (taskpwg4,     0, coutput1 );
    mgr->ConnectOutput (taskpwg4,     1, coutput2 );
	
    mgr->ConnectInput  (jetana,     0, cinput1);
    mgr->ConnectOutput (jetana,     0, coutput1 );
    mgr->ConnectOutput (jetana,     1, coutput3 );

    //------------------------  
    //Scaling task
    //-----------------------
    Int_t nfiles = chainxs->GetEntries();
    //cout<<"Get? "<<kGetXSectionFromFileAndScale<<" nfiles "<<nfiles<<endl;
    if(kGetXSectionFromFileAndScale && nfiles > 0){
      //cout<<"Init AnaScale"<<endl;
      //Get the cross section
      Double_t xsection=0; 
      Float_t ntrials = 0;
      GetAverageXsection(chainxs, xsection, ntrials);
      
      AliAnaScale * scale = new AliAnaScale("scale") ;
      scale->Set(xsection/ntrials/kNumberOfEventsPerFile/nfiles) ;
      scale->MakeSumw2(kFALSE);
      scale->SetDebugLevel(2);
      mgr->AddTask(scale);
      
      AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("corrhistosscaled", TList::Class(),
								AliAnalysisManager::kOutputContainer, "histosscaled.root");
      mgr->ConnectInput  (scale,     0, coutput2);
      mgr->ConnectOutput (scale,     0, coutput4 );
	  
      AliAnaScale * scalejet = new AliAnaScale("scalejet") ;
      scalejet->Set(xsection/ntrials/kNumberOfEventsPerFile/nfiles) ;
      scalejet->MakeSumw2(kFALSE);
      //scalejet->SetDebugLevel(2);
      mgr->AddTask(scalejet);
      
      AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("jethistosscaled", TList::Class(),
								AliAnalysisManager::kOutputContainer, "histosscaled.root");
      mgr->ConnectInput  (scalejet,     0, coutput3);
      mgr->ConnectOutput (scalejet,     0, coutput5 );	  
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

void  LoadLibraries(const anaModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");
  
  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (mode==mLocal || mode == mLocalCAF || mode == mGRID) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
    //gSystem->Load("libSTEERBase");
    //gSystem->Load("libESD");
    //gSystem->Load("libAOD");
    //gSystem->Load("libANALYSIS");
    //gSystem->Load("libANALYSISalice");
    //gSystem->Load("libJETAN");
    //gSystem->Load("libPHOSUtils");
    //gSystem->Load("libPWG4PartCorrBase");
    //gSystem->Load("libPWG4PartCorrDep");
	
    //--------------------------------------------------------
    //If you want to use root and par files from aliroot
    //--------------------------------------------------------  
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
    SetupPar("JETAN");
    //If your analysis needs PHOS geometry uncomment following lines
//      SetupPar("PHOSUtils");
//      //Create Geometry
//    TGeoManager::Import("geometry.root") ; //need file "geometry.root" in local dir!!!!
    SetupPar("PWG4PartCorrBase");
    SetupPar("PWG4PartCorrDep");
    
  }

  //---------------------------------------------------------
  // <<<<<<<<<< PROOF mode >>>>>>>>>>>>
  //---------------------------------------------------------
  else if (mode==mPROOF) {
    //
    // Connect to proof
    // Put appropriate username here
    // TProof::Reset("proof://mgheata@lxb6046.cern.ch"); 
    TProof::Open("proof://mgheata@lxb6046.cern.ch");
    
    //    gProof->ClearPackages();
    //    gProof->ClearPackage("ESD");
    //    gProof->ClearPackage("AOD");
    //    gProof->ClearPackage("ANALYSIS"); 
    //    gProof->ClearPackage("JETAN");  
    //    gProof->ClearPackage("PWG4PartCorrDep");
    //    gProof->ClearPackage("PWG4PartCorrBase");
    
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
    // Enable JETAN
    gProof->UploadPackage("JETAN.par");
    gProof->EnablePackage("JETAN");
    // Enable PHOSUtils
    //gProof->UploadPackage("PHOSUtils.par");
    //gProof->EnablePackage("PHOSUtils");
    // Enable PartCorr analysis
    gProof->UploadPackage("PWG4PartCorrBase.par");
    gProof->EnablePackage("PWG4PartCorrBase");
    gProof->UploadPackage("PWG4PartCorrDep.par");
    gProof->EnablePackage("PWG4PartCorrDep");
    //
    gProof->ShowEnabledPackages();
  }  
  
}

void SetupPar(char* pararchivename)
{
  //Load par files, create analysis libraries
  //For testing, if par file already decompressed and modified
  //classes then do not decompress.
 
  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ; 
  TString parpar(Form("%s.par", pararchivename)) ; 
  if ( gSystem->AccessPathName(parpar.Data()) ) {
    gSystem->ChangeDirectory(gSystem->Getenv("ALICE_ROOT")) ;
    TString processline(Form(".! make %s", parpar.Data())) ; 
    gROOT->ProcessLine(processline.Data()) ;
    gSystem->ChangeDirectory(cdir) ; 
    processline = Form(".! mv $ALICE_ROOT/%s .", parpar.Data()) ;
    gROOT->ProcessLine(processline.Data()) ;
  } 
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



void CreateChain(const anaModes mode, TChain * chain, TChain * chainxs){
  //Fills chain with data
  TString ocwd = gSystem->WorkingDirectory();
  
  //-----------------------------------------------------------
  //Analysis of CAF data locally and with PROOF
  //-----------------------------------------------------------
  if(mode ==mPROOF || mode ==mLocalCAF){
    // Chain from CAF
    gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
    // The second parameter is the number of input files in the chain
    chain = CreateESDChain("ESD12001.txt", 5);  
  }
  
  //---------------------------------------
  //Local files analysis
  //---------------------------------------
  else if(mode == mLocal){    
    //If you want to add several ESD files sitting in a common directory INDIR
    //Specify as environmental variables the directory (INDIR), the number of files 
    //to analyze (NEVENT) and the pattern name of the directories with files (PATTERN)

    if(gSystem->Getenv("INDIR"))  
      kInDir = gSystem->Getenv("INDIR") ; 
    else cout<<"INDIR not set, use default: "<<kInDir<<endl;	
    
    if(gSystem->Getenv("PATTERN"))   
      kPattern = gSystem->Getenv("PATTERN") ; 
    else  cout<<"PATTERN not set, use default: "<<kPattern<<endl;
    
    if(gSystem->Getenv("NEVENT"))
      kEvent = atoi(gSystem->Getenv("NEVENT")) ;
    else cout<<"NEVENT not set, use default: "<<kEvent<<endl;
    
    //Check if env variables are set and are correct
    if ( kInDir  && kEvent) {
      printf("Get %d files from directory %s\n",kEvent,kInDir);
      if ( ! gSystem->cd(kInDir) ) {//check if ESDs directory exist
	printf("%s does not exist\n", kInDir) ;
	return ;
      }
      cout<<"INDIR : "<<kInDir<<endl;
      cout<<"NEVENT : "<<kEvent<<endl;
      cout<<"PATTERN: " <<kPattern<<endl;
		
	  TString datafile="";
      if(kInputData == "ESD") datafile = "AliESDs.root" ;
      else if(kInputData == "AOD") datafile = "aod.root" ;
      else if(kInputData == "MC")  datafile = "galice.root" ;
		
      //Loop on ESD files, add them to chain
      Int_t event =0;
      Int_t skipped=0 ; 
      char file[120] ;
      char filexs[120] ;
      
      for (event = 0 ; event < kEvent ; event++) {
	sprintf(file, "%s/%s%d/%s", kInDir,kPattern,event,datafile.Data()) ;  
	sprintf(filexs, "%s/%s%d/%s", kInDir,kPattern,event,kXSFileName) ; 
	TFile * fESD = 0 ; 
	//Check if file exists and add it, if not skip it
	if ( fESD = TFile::Open(file)) {
	  if ( fESD->Get(kTreeName) ) { 
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
      printf("number of entries # %lld, skipped %d\n", chain->GetEntries(), skipped) ; 	
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
    gSystem->Load("libNetx.so") ; 
    gSystem->Load("libRAliEn.so"); 
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
    for (Int_t index = 0; index < result->GetEntries(); index++) {
      TString alienURL = result->GetKey(index, "turl") ; 
      cout << "================== " << alienURL << endl ; 
      chain->Add(alienURL) ; 
      alienURL.ReplaceAll("AliESDs.root",kXSFileName);
      chainxs->Add(alienURL) ; 
    }
  }// xml analysis
  
  gSystem->ChangeDirectory(ocwd.Data());
}

//________________________________________________
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
    tree->SetBranchAddress("ntrials",&ntrials);
    for(Int_t i = 0; i < nfiles; i++){
      tree->GetEntry(i);
      xs += xsection ;
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



