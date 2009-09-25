/* $Id:  $ */
//--------------------------------------------------
// Example macro to do multi-platform JETAN/FASTJET analysis
// Can be executed with Root and AliRoot
//
// Configured by options and definitions set in the lines below and
// additional external configuration files and environement variables.
//
// Author: K. Read
//
//-------------------------------------------------
enum anaModes {mLocal, mLocalCAF, mPROOF, mGRID, mPLUGIN};
//mLocal    = 0: Analyze locally files in your computer
//mLocalCAF = 1: Analyze locally CAF files
//mPROOF    = 2: Analyze CAF files with PROOF
//mGRID     = 3: Analyze files on GRID
//mPLUGIN   = 4: Analyze files on GRID with AliEn plugin

//---------------------------------------------------------------------------
//Settings to read locally several files, only for "mLocal" mode
//The different values are default, they can be set with environmental 
//variables: INDIR, PATTERN, NFILES, respectively
//char * kInDir = "/afs/cern.ch/user/k/kread/public/data/"; 
char * kInDir = "/user/data/files";
char * kPattern = ""; // Data are in files kInDir/kPattern+i
Int_t kFile = 1; // Number of files
//---------------------------------------------------------------------------
//Collection file for grid analysis
char * kXML = "collection.xml";
//---------------------------------------------------------------------------
//Data directory for PROOF analysis
char * kmydataset = "/COMMON/COMMON/LHC09a4_run8101X";
//char * kmydataset = "/PWG4/mcosenti/LHC08d10_ppElectronB_Jets#esdTree";
//---------------------------------------------------------------------------
//Scale histograms from file. Change to kTRUE when xsection file exists
//Put name of file containing xsection 
//Put number of events per ESD file
//This is an specific case for normalization of Pythia files.
const Bool_t kGetXSectionFromFileAndScale = kTRUE ;
const char * kXSFileName = "pyxsec.root";
const Int_t kNumberOfEventsPerFile = 200; 
//---------------------------------------------------------------------------

const Bool_t kMC = kTRUE; //With real data kMC = kFALSE
const TString kInputData = "ESD";//ESD, AOD, MC
TString kTreeName = "esdTree";
//const   Bool_t kMergeAODs = kTRUE; //uncomment for AOD merging
const Bool_t kMergeAODs = kFALSE; //uncomment for no AOD merging
Int_t sevent = 0;
	
Int_t mode = mLocal;
char sconfig1[1024] = "ConfigPWG4AODtoAOD";        //"ConfigAnalysis";
char sconfig2[1024] = "ConfigJetAnalysisFastJet.C";//"ConfigAnalysis";
char sconfig3[1024] = "ConfigAnalysisElectron";    //"ConfigAnalysis";

void anaJete()
{
  // Main

  //Process environmental variables from command line:
  ProcessEnvironment();	
  printf("Final    Variables: mode %d, config1 %s, config2 %s, config3 %s, sevent %d\n", mode,sconfig1,sconfig2,sconfig3,sevent);

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

  TChain * chain   = new TChain(kTreeName) ;
  TChain * chainxs = new TChain("Xsection") ;

  if (mode==mLocal || mode==mLocalCAF || mode == mGRID) {
    CreateChain(mode, chain, chainxs);  
  }

  if( chain || mode==mPROOF || mode==mPLUGIN ){
    AliLog::SetGlobalLogLevel(AliLog::kError);//Minimum prints on screen
    
    //--------------------------------------
    // Make the analysis manager
    //-------------------------------------
    AliAnalysisManager *mgr  = new AliAnalysisManager("Jet Manager", "Jet Manager");

    if( mode == mPLUGIN ){
      // Create and configure the alien handler plugin
      if (!AliAnalysisGrid::CreateToken()) return NULL;
      AliAnalysisAlien *plugin = new AliAnalysisAlien();
      plugin->SetRunMode("submit");
      //Uncomment the following 3 lines to permit auto xml creation
      //plugin->SetGridDataDir("/alice/sim/PDC_08b/LHC08d10/"); //dummy
      //plugin->SetDataPattern("AliESDs.root"); //dummy
      //plugin->AddRunNumber(30010); //dummy
      plugin->AddDataFile("mycollect.xml");
      plugin->SetGridWorkingDir("work3");
      plugin->SetAdditionalLibs("anaJet.C ConfigJetAnalysisFastJet.C ConfigAnalysisElectron.C ANALYSIS.par ANALYSISalice.par AOD.par ESD.par STEERBase.par JETAN.par FASTJETAN.par");
      plugin->SetJDLName("anaJet.jdl");
      plugin->SetExecutable("anaJet.sh");
      plugin->SetOutputFiles("histos.root");
      AliAnalysisGrid *alienHandler = plugin;
      if (!alienHandler) return;

      // Connect plug-in to the analysis manager
      mgr->SetGridHandler(alienHandler);
    }

    // MC handler
    if(kMC || kInputData == "MC"){
      AliMCEventHandler* mcHandler = new AliMCEventHandler();
      mcHandler->SetReadTR(kFALSE);//Do not search TrackRef file
      mgr->SetMCtruthEventHandler(mcHandler);
      if( kInputData == "MC") mgr->SetInputEventHandler(NULL);
    }

    // AOD output handler
    AliAODHandler* aodoutHandler   = new AliAODHandler();
    aodoutHandler->SetOutputFileName("aodoutput.root");
    if(kMergeAODs)aodoutHandler->SetCreateNonStandardAOD();
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
      if(kMergeAODs){
        char path[1024];
        sprintf(path,"AliAOD.root");
        if(gSystem->Getenv("SIMPATH"))
            sprintf(path,"%s/AliAOD.root",gSystem->Getenv("SIMPATH"));
        cout<<"Config: Second input file: "<<path<<endl;
        aodHandler->SetMergeEvents(kTRUE);
        aodHandler->AddFriend(path);
        cout<<"Config: Starting event for second input file: "<<sevent<<endl;
        aodHandler->SetMergeOffset(sevent);
      }
    }

    mgr->SetDebugLevel(3); // For debugging, do not uncomment if you want no messages.


    const Bool_t kDoESDFilter = kTRUE; //need this for JETAN with input ESDs
  //const Bool_t kDoESDFilter = kFALSE;
    if(kInputData == "ESD" && kDoESDFilter){
      printf("Applying ESD filter cuts appropriate for jet analysis\n");
      //
      // Set of cuts
      // 
      // standard
      AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
      //esdTrackCutsL->SetMinNClustersTPC(50);
      //esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
      //esdTrackCutsL->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
      //esdTrackCutsL->SetRequireTPCRefit(kTRUE);
      //esdTrackCutsL->SetMinNsigmaToVertex(3); //keep commented out
      //esdTrackCutsL->SetRequireSigmaToVertex(kTRUE); //keep commented out
      //esdTrackCutsL->SetAcceptKinkDaughters(kFALSE);
      //
      // hard
      AliESDtrackCuts* esdTrackCutsH = new AliESDtrackCuts("AliESDtrackCuts", "Hard");
      esdTrackCutsH->SetMinNClustersTPC(100);
      esdTrackCutsH->SetMaxChi2PerClusterTPC(2.0);
      esdTrackCutsH->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
      esdTrackCutsH->SetRequireTPCRefit(kTRUE);
      //esdTrackCutsH->SetMinNsigmaToVertex(2);
      //esdTrackCutsH->SetRequireSigmaToVertex(kTRUE);
      esdTrackCutsH->SetAcceptKinkDaughters(kFALSE);
      //
      AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
      trackFilter->AddCuts(esdTrackCutsL);
      //  trackFilter->AddCuts(esdTrackCutsH);
      //
      AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
      esdfilter->SetTrackFilter(trackFilter);
      esdfilter->SetDebugLevel(10);
      mgr->AddTask(esdfilter);
    }


    //-------------------------------------------------------------------------
    //Define task, put here any other task that you want to use.
    //-------------------------------------------------------------------------

    //
    // Jet analysis
    //
    AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis",chain);
    jetana->SetDebugLevel(2);
    jetana->SetConfigFile(sconfig2);  //Default ConfigJetAnalysisFastJet
    //Uncommenting the following line produces too many AddAtAndExpand warnings for now.
    if(kMergeAODs)jetana->ReadAODFromOutput();  //Uncomment when AOD merging
    mgr->AddTask(jetana);

    //
    // electron analysis
    //
    AliAnalysisTaskParticleCorrelation * taskpwg4 = new AliAnalysisTaskParticleCorrelation ("Particle");
    taskpwg4->SetConfigFileName("ConfigAnalysisElectron"); //Default name is ConfigAnalysisElectron
    mgr->AddTask(taskpwg4);
    
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TList::Class(),
							      AliAnalysisManager::kOutputContainer, "histos.root");


    if(kInputData == "ESD" && kDoESDFilter){
      mgr->ConnectInput  (esdfilter,  0, cinput1  );
      mgr->ConnectOutput (esdfilter,  0, coutput1 );
    }

    mgr->ConnectInput  (jetana,    0, cinput1  );
    mgr->ConnectOutput (jetana,    0, coutput1 );
    mgr->ConnectOutput (jetana,    1, coutput2 );
  
    mgr->ConnectInput  (taskpwg4, 0, cinput1  );
    mgr->ConnectOutput (taskpwg4, 0, coutput1 );
    mgr->ConnectOutput (taskpwg4, 1, coutput2 );
 

    //------------------------  
    //Scaling task
    //-----------------------
    Int_t nfiles = chainxs->GetEntries();
    Int_t nevents = chain->GetEntries();
    cout<<"Get? "<<kGetXSectionFromFileAndScale<<" nfiles "<<nfiles<<" nevents "<<nevents<<endl;
    if(kGetXSectionFromFileAndScale && nfiles > 0){
      //cout<<"Init AnaScale"<<endl;
      //Get the cross section
      Double_t xsection=0; 
      Float_t ntrials = 0;
      GetAverageXsection(chainxs, xsection, ntrials);
      
      AliAnaScale * scale = new AliAnaScale("scale") ;

     cout<<"Scale factor "<<xsection/ntrials/kNumberOfEventsPerFile/nfiles<<"  "<<xsection/ntrials/nevents<<endl;
      scale->Set(xsection/ntrials/nevents) ;
      scale->MakeSumw2(kTRUE);//If you want histograms with error bars set to kTRUE
      scale->SetDebugLevel(2);
      mgr->AddTask(scale);
      
      AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("histosscaled", TList::Class(),
								AliAnalysisManager::kOutputContainer, "histosscaled.root");
      mgr->ConnectInput  (scale,     0, coutput2);
      mgr->ConnectOutput (scale,     0, coutput3 );
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
    else if (mode==mPLUGIN) 
      smode = "grid";
    
    //mgr->ResetAnalysis();
    mgr->InitAnalysis();
    mgr->PrintStatus();
    if (mode==mPROOF)
      mgr->StartAnalysis(smode.Data(),kmydataset,1500,0);
    else if (mode==mPLUGIN)
      mgr->StartAnalysis(smode.Data());
    else
      mgr->StartAnalysis(smode.Data(),chain,20);

    cout <<" Analysis ended sucessfully "<< endl ;

  }
  else cout << "Chain was not produced ! "<<endl;
  
}

void  LoadLibraries(const anaModes mode) {
  
  
  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (mode==mLocal || mode == mLocalCAF || mode == mGRID || mode == mPLUGIN) {
    bool usepar = true;

    //--------------------------------------
    // Load the needed libraries most of them already loaded by aliroot
    //--------------------------------------
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libXMLIO.so");
    gSystem->Load("libCGAL");
    gSystem->Load("libfastjet");
    gSystem->Load("libSISConePlugin");

    if(usepar){
      //--------------------------------------------------------
      //If you want to use root and par files from aliroot
      //--------------------------------------------------------  
      SetupPar("STEERBase");
      SetupPar("ESD");
      SetupPar("AOD");
      SetupPar("ANALYSIS");
      SetupPar("ANALYSISalice");
      cerr<<"Now Loading JETAN"<<endl;
      SetupPar("JETAN");
      cerr<<"Done Loading JETAN"<<endl;
      cerr<<"Now Loading FASTJETAN"<<endl;
      SetupPar("FASTJETAN");
      cerr<<"Done Loading FASTJETAN"<<endl;
      SetupPar("PWG4PartCorrBase");
      SetupPar("PWG4PartCorrDep");
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
      gSystem->Load("libJETAN");
      gSystem->Load("libFASTJETAN");
      gSystem->Load("libPWG4PartCorrBase");
      gSystem->Load("libPWG4PartCorrDep");
    }

	
  }

  //---------------------------------------------------------
  // <<<<<<<<<< PROOF mode >>>>>>>>>>>>
  //---------------------------------------------------------
  else if (mode==mPROOF) {
    //
    // Connect to proof
    // Put appropriate username here
    // char* myproofname = "alicecaf";
    char* myproofname = "kread@localhost";

    //TProof::Reset("proof://kread@lxb6046.cern.ch");
    //TProof::Reset("proof://myproofname);
    //TProof::Reset("myproofname",kTRUE);
    gEnv->SetValue("XSec.GSI.DelegProxy","2");   
    //TProof::Mgr(myproofname)->ShowROOTVersions();
    //TProof::Mgr(myproofname)->SetROOTVersion("v5-23-04");
    TProof::Open(myproofname);

    // gProof->ClearPackages();
    // gProof->SetLogLevel(5);
    // gProof->ClearPackage("STEERBase");
    // gProof->ClearPackage("ESD");
    // gProof->ClearPackage("AOD");
    // gProof->ClearPackage("ANALYSIS");
    // gProof->ClearPackage("ANALYSISalice");
    // gProof->ClearPackage("JETAN");
    // gProof->ClearPackage("FASTJETAN");
    // gProof->ShowEnabledPackages();

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
    // Enable the Analysis Package
    gProof->UploadPackage("ANALYSISalice.par");
    gProof->EnablePackage("ANALYSISalice");
    // Enable JETAN analysis
    gProof->UploadPackage("JETAN.par");
    gProof->EnablePackage("JETAN");
    // Enable FASTJETAN analysis
    gProof->UploadPackage("FASTJETAN.par");
    gProof->EnablePackage("FASTJETAN");

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
  //Analysis of CAF data locally
  //-----------------------------------------------------------
  if(mode == mLocalCAF){
    // Read the input list of files and add them to the chain
    TString line;
    ifstream in;
    in.open("ESDlist.txt");
    while (in.good()) {
      in >> line;
      if (line.Length() == 0) continue;
      // cout << " line = " << line << endl;
      chain->Add(line);
    }
  }
  
  //---------------------------------------
  //Local files analysis
  //---------------------------------------
  else if(mode == mLocal){
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
	if(!strcmp(kGener,"PYTHIA")) kXSFileName = "pyxsec.root";
	else if(!strcmp(kGener,"HERWIG")) kXSFileName = "hexsec.root";
	else cout<<" UNKNOWN GENER, use default: "<<kXSFileName<<endl;
      }

      else cout<<" GENER not set, use default xs file name: "<<kXSFileName<<endl;

      cout<<"INDIR : "<<kInDir<<endl;
      cout<<"NFILES : "<<kFile<<endl;
      cout<<"PATTERN: " <<kPattern<<endl;
      cout<<"XSFILE  : "<<kXSFileName<<endl;
      
      TString datafile="";
      if(kInputData == "ESD") datafile = "AliESDs.root" ;
      else if(kInputData == "AOD") datafile = "AliAOD.root" ;
      else if(kInputData == "MC")  datafile = "galice.root" ;
      
      //Loop on ESD files, add them to chain
      Int_t event =0;
      Int_t skipped=0 ; 
      char file[120] ;
      char filexs[120] ;
      
      for (event = 0 ; event < kFile ; event++) {
	sprintf(file, "%s/%s%d/%s", kInDir,kPattern,event,datafile.Data()) ; 
	sprintf(filexs, "%s/%s%d/%s", kInDir,kPattern,event,kXSFileName) ; 
	TFile * fESD = 0 ; 
	//Check if file exists and add it, if not skip it
	if ( fESD = TFile::Open(file)) {
	  if ( fESD->Get(kTreeName) ) { 
	    printf("++++ Adding %s\n", file) ;
	    chain->AddFile(file);
	    if(kGetXSectionFromFileAndScale)chainxs->Add(filexs) ; 
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

void ProcessEnvironment(){

	if (gSystem->Getenv("MODE"))
		mode = atoi(gSystem->Getenv("MODE"));

	if (gSystem->Getenv("CONFIG1"))
		sprintf(sconfig1,gSystem->Getenv("CONFIG1"));

	if (gSystem->Getenv("CONFIG2"))
		sprintf(sconfig2,gSystem->Getenv("CONFIG2"));

	if (gSystem->Getenv("CONFIG3"))
		sprintf(sconfig3,gSystem->Getenv("CONFIG3"));

	if (gSystem->Getenv("SEVENT"))
		sevent = atoi (gSystem->Getenv("SEVENT"));
	
	printf("PROCESS: Variables: mode %d, config1 %s, config2 %s, config3 %s, sevent %d\n", mode,sconfig1,sconfig2,sconfig3,sevent);

	
}
