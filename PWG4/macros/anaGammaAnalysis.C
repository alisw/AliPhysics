/* $Id$ */
/* $Log$
/* Revision 1.2  2007/12/13 09:45:45  gustavo
/* Scaling option and more comentaries added
/*
/* Revision 1.1  2007/12/07 14:13:02  gustavo
/* Example macros for execution and configuration of the analysis
/* */

//---------------------------------------------------
// Example macro to do analysis with the 
// analysis classes in PWG4Gamma
// Can be executed with Root and AliRoot
//
// Pay attention to the loading of libraries in LoadLibraries()
// in the local mode.
//
//  Author : Gustavo Conesa Balbastre (INFN-LNF)
//
//-------------------------------------------------
enum anaModes {mLocal, mLocalCAF,mPROOF,mGRID};
//mLocal: Analyze locally files in your computer
//mLocalCAF: Analyze locally CAF files
//mPROOF: Analyze CAF files with PROOF

//Settings to read locally several files
//The differnt values are default, they can be set with environmental 
//variables: OUTDIR, PATTERN, NEVENT, respectivelly
char * kInDir = "/data/"; 
char * kPattern = "Run"; // Data are in diles /data/Run0, 
// /Data/Run1 ...
Int_t kEvent = 3; // Number of files


//Scale histograms from file. Change to kTRUE when xsection file exists
//Put name of file containing xsection 
//Put number of events per ESD file
//This is an specific case for normalization of Pythia files.
const Bool_t kGetXSectionFromFileAndScale = kFALSE ;
const char * kXSFileName = "pyxsec.root";
const Int_t kNumberOfEventsPerFile = 100; 


void anaGammaAnalysis(Int_t mode=mLocal, TString configName = "ConfigGammaAnalysis")
{
  // Main
  
  //
  // Load analysis libraries
  // Look at the method below, 
  // change whatever you need for your analysis case
  // 
  LoadLibraries(mode) ;
  
  //
  //Create chain, get cross sections, look below for options.
  // 
  Double_t xsection = 0;
  Int_t ntrials = 0;
  Int_t nfiles = 0;
  TChain* chain = CreateChain(mode, xsection, ntrials, nfiles);  
  if(kGetXSectionFromFileAndScale && nfiles > 0){
    xsection /= nfiles ;
    ntrials  /= nfiles ;
    cout<<"//////////////////////////////////////////////////////////////"<<endl;
    cout<<"Average cross section: "<<xsection<<" ntrials "<<ntrials<<" files "<<nfiles<<endl;
    cout<<"//////////////////////////////////////////////////////////////"<<endl;
  }

  if(chain){
  //
  // Make the analysis manager
  //
  AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");
  // MC handler
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  // AOD output handler
  AliAODHandler* aodHandler   = new AliAODHandler();
  aodHandler->SetOutputFileName("aod.root");
  mgr->SetOutputEventHandler(aodHandler);
  // ESD handler
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdHandler);

  //mgr->SetDebugLevel(10);

  //Define task, put here any other task that you want to use.
  AliAnalysisTaskGamma * taskgamma = new AliAnalysisTaskGamma ("Gamma");
  taskgamma->SetConfigFileName(configName); //Default name is ConfigGammaAnalysis
  mgr->AddTask(taskgamma);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
                      AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
                      AliAnalysisManager::kOutputContainer, "default");
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("gammahistos", TList::Class(),
                      AliAnalysisManager::kOutputContainer, "gammahistos.root");

  mgr->ConnectInput  (taskgamma,     0, cinput1);
  mgr->ConnectOutput (taskgamma,     0, coutput1 );
  mgr->ConnectOutput (taskgamma,     1, coutput2 );

  //Scaling task
  
  if(kGetXSectionFromFileAndScale){
    AliAnaScale * scale = new AliAnaScale("scale") ;
    scale->Set(xsection/ntrials/kNumberOfEventsPerFile/nfiles) ;
    mgr->AddTask(scale);
    
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("gammahistosscaled", TList::Class(),
							      AliAnalysisManager::kOutputContainer, "gammahistosscaled.root");
    mgr->ConnectInput  (scale,     0, coutput2);
    mgr->ConnectOutput (scale,     0, coutput3 );
  }
  
  //
  // Run the analysis
  //    
  TString smode = "";
  if (mode==mLocal || mode == mLocalCAF) 
    smode = "local";
  else if (mode==mPROOF) 
    smode = "proof";
   else if (mode==mGRID) 
    smode = "grid";
  
  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis(smode.Data(),chain);
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
  
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  if (mode==mLocal || mode == mLocalCAF || mode == mGRID) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
    //gSystem->Load("libSTEERBase");
    //gSystem->Load("libESD");
    //gSystem->Load("libAOD");
    //gSystem->Load("libANALYSIS");
    //gSystem->Load("libPWG4Gamma");

    //--------------------------------------------------------
    //If you want to use your own modified classes
    //--------------------------------------------------------  
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("PWG4Gamma");
    
    //--------------------------------------------------------
    // If the modified libraries have already been compiled and 
    // you don't want to  decompress them and recompile
    //--------------------------------------------------------
    //SetupPar("STEERBase",kFALSE);
    //SetupPar("ESD",kFALSE);
    //SetupPar("AOD",kFALSE);
    //SetupPar("ANALYSIS",kFALSE);
    //SetupPar("PWG4Gamma",kFALSE);
  }
  
  // <<<<<<<<<< PROOF mode >>>>>>>>>>>>
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
    //    gProof->ClearPackage("PWG4Gamma");
    
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
    // Enable gamma jet analysis
    gProof->UploadPackage("PWG4Gamma.par");
    gProof->EnablePackage("PWG4Gamma");
    //
    gProof->ShowEnabledPackages();
  }  
  
}

void SetupPar(char* pararchivename, Bool_t decomp = kTRUE)
{

  //Load par files, create analysis libraries
  //For testing, if par file already decompressed and modified
  //classes then do not decompress.

  if (pararchivename) {
    char processline[1024];
    if(decomp){
      sprintf(processline,".! tar xvzf %s.par",pararchivename);
      gROOT->ProcessLine(processline);
    }
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


TChain * CreateChain(const anaModes mode, Double_t &xsection, Int_t &ntrials, Int_t &nfiles){
  //Creates data chain
  TChain *chain;

  TString ocwd = gSystem->WorkingDirectory();
  //Analysis of CAF data locally and with PROOF
  if(mode ==mPROOF || mode ==mLocalCAF){
    // Chain from CAF
    gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
    // The second parameter is the number of input files in the chain
    chain = CreateESDChain("ESD12001.txt", 5);  
  }

  //Local files analysis
  else if(mode == mLocal){
    chain = new TChain("esdTree") ;
    TString input = "AliESDs.root" ;
    
    //If you want to add several ESD files sitting in a common directory OUTDIR
    //Specify as environmental variables the directory (OUTDIR), the number of files 
    //to analyze (NEVENT) and the pattern name of the directories with files (PATTERN)
    if(gSystem->Getenv("OUTDIR"))  
      kInDir = gSystem->Getenv("OUTDIR") ; 
    else cout<<"OUTDIR not set, use default: "<<kInDir<<endl;	

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
      cout<<"OUTDIR : "<<kInDir<<endl;
      cout<<"NEVENT : "<<kEvent<<endl;
      cout<<"PATTERN: " <<kPattern<<endl;
      
      //Loop on ESD files, add them to chain
      Int_t event =0;
      Int_t skipped=0 ; 
      char file[120] ;
      Double_t * rv = new Double_t[2] ;
      rv[0] = rv[1] = 0.0 ;

      for (event = 0 ; event < kEvent ; event++) {
        sprintf(file, "%s/%s%d/AliESDs.root", kInDir,kPattern,event) ; 
	TFile * fESD = 0 ; 
	//Check if file exists and add it, if not skip it
	if ( fESD = TFile::Open(file)) {
	  //Get cross section if file exists
	  ReadXsection(kInDir, kPattern, event, rv) ;
	  if ( fESD->Get("esdTree") ) { 
	    printf("++++ Adding %s\n", file) ;
	    chain->AddFile(file);
	    nfiles++;
	    if(rv) {
	      cout << "xsection" << rv[0] << "; ntrials " <<rv[1]<< endl ;	      
	      xsection += rv[0] ;
	      ntrials += rv[1] ;
	    }
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
      cout<<">>>>>> No list added, take a single file <<<<<<<<< "<<input<<endl;
      chain->AddFile(input);
    }

  }// local files analysis

  //GRID xml files
  else if(mode == mGRID){
    
    //Get colection file
    TString input = "collection.xml" ;
    char kXML  [1000];
    if(gSystem->Getenv("XML") )
      kXML = gSystem->Getenv("XML");
    else
        sprintf(kXML, "collection.xml") ; 

    cout<<"XML file "<<kXML<<endl;
    
    if (!TFile::Open(kXML)) {
      printf("No collection file with name -- %s -- was found\n",kXML);
      return ;
    }

#ifdef WITHALIEN
    TGridCollection * collection =  (TGridCollection*)gROOT->ProcessLine(Form("TAlienCollection::Open(\"%s\", 0)", kXML));
    if (! collection) {
      AliError(Form("%s not found", kXML)) ; 
      return kFALSE ; 
    }
    
    TGridResult* result = collection->GetGridResult("",0 ,0);
    TList* analysisfilelist = result->GetFileInfoList();
    
    // Makes the ESD chain 
    printf("*** Getting the Chain       ***\n");
    chain->AddFileInfoList(analysisfilelist);
    
#endif

  }// xml analysis
  
  gSystem->ChangeDirectory(ocwd.Data());
  
  return chain ;
}

//________________________________________________
void ReadXsection(const char * inDir, const char * pattern, const Int_t event, Double_t * rv)
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

  Double_t xsection;
  UInt_t    ntrials;

  char cfile[100] ;
  sprintf(cfile, "%s/%s%d/%s", inDir, pattern,event,kXSFileName) ;
  TFile *file = new TFile(cfile,"readonly");

  if ( ! file ) {
    AliInfo(Form("Cross section file not found %s", cfile)) ;
    rv=0x0 ;
  }
  else{
    TTree *tree = file->Get("Xsection");
    if (tree) {
      tree->SetBranchAddress("xsection",&xsection);
      tree->SetBranchAddress("ntrials",&ntrials);
      tree->GetEntry(0);
      rv[0] = xsection ;
      rv[1] = ntrials ;
    } 
    else
      rv = 0x0 ;
  }
}
