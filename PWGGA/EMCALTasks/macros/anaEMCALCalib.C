// $Id$
//--------------------------------------------------
// Example macro to do EMCAL calibration analysis
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
//variables: INDIR, PATTERN, NFILES, respectivelly
char * kInDir = "/user/data/files/"; 
char * kPattern = ""; // Data are in files kInDir/kPattern+i 
Int_t kFile = 1; // Number of files
//---------------------------------------------------------------------------
//Collection file for grid analysis
char * kXML = "collection.xml";
//---------------------------------------------------------------------------

const TString kInputData = "ESD"; //ESD, AOD, MC
TString kTreeName = "esdTree";

void anaEMCALCalib(Int_t mode=mLocal)
{
  // Main
  char cmd[200] ; 
  sprintf(cmd, ".! rm -rf aod.root pi0calib.root") ; 
  gROOT->ProcessLine(cmd) ; 
  //--------------------------------------------------------------------
  // Load analysis libraries
  // Look at the method below, 
  // change whatever you need for your analysis case
  // ------------------------------------------------------------------
  LoadLibraries(mode) ;
  //gSystem->Unload("libPWGGAEMCALTasks.so");
  //Try to set the new library
  //gSystem->Load("./PWGGAEMCALTasks/libPWGGAEMCALTasks.so");
  //gSystem->ListLibraries();
  
  //-------------------------------------------------------------------------------------------------
  //Create chain from ESD and from cross sections files, look below for options.
  //-------------------------------------------------------------------------------------------------
  if(kInputData == "ESD") kTreeName = "esdTree" ;
  else if(kInputData == "AOD") kTreeName = "aodTree" ;
  else {
    cout<<"Wrong  data type "<<kInputData<<endl;
    break;
  }
  
  TChain *chain       = new TChain(kTreeName) ;
  CreateChain(mode, chain);  
  
  if(chain){
    AliLog::SetGlobalLogLevel(AliLog::kError);//Minimum prints on screen
    
    //--------------------------------------
    // Make the analysis manager
    //-------------------------------------
    AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");
    
    //input
    if(kInputData == "ESD")
    {
      // ESD handler
      AliESDInputHandler *esdHandler = new AliESDInputHandler();
      mgr->SetInputEventHandler(esdHandler);
      esdHandler->SetReadFriends(kFALSE);
      cout<<"ESD handler "<<mgr->GetInputEventHandler()<<endl;
    }
    if(kInputData == "AOD")
    {
      // AOD handler
      AliAODInputHandler *aodHandler = new AliAODInputHandler();
      mgr->SetInputEventHandler(aodHandler);
      cout<<"AOD handler "<<mgr->GetInputEventHandler()<<endl;
      
    }
    
    // mgr->SetDebugLevel(1);
    
    //-------------------------------------------------------------------------
    //Define task, put here any other task that you want to use.
    //-------------------------------------------------------------------------
    
/*    
    if(kInputData == "ESD"){
      
      gROOT->LoadMacro("AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
      
    }
*/    
    
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
    
    AliAnalysisTaskEMCALPi0CalibSelection * pi0calib = new AliAnalysisTaskEMCALPi0CalibSelection ("EMCALPi0Calibration");
    if(kInputData == "ESD") pi0calib->SelectCollisionCandidates(); 
    //pi0calib->SetDebugLevel(10); 
    //pi0calib->UseFilteredEventAsInput();
    pi0calib->SetClusterMinEnergy(0.5);
    pi0calib->SetClusterMaxEnergy(10.);
    pi0calib->SetAsymmetryCut(0.5);
    pi0calib->SetClusterMinNCells(0);
    pi0calib->SetNCellsGroup(0);
    pi0calib->SwitchOnSameSM();
    //pi0calib->SwitchOnOldAODs();
    
/*
    TGeoHMatrix *matrix[4];
    
    double rotationMatrix[4][9] = {-0.014585, -0.999892, -0.002031, 0.999892, -0.014589,  0.001950, -0.001979, -0.002003,  0.999996,
      -0.014585,  0.999892,  0.002031, 0.999892,  0.014589, -0.001950, -0.001979,  0.002003, -0.999996,
      -0.345861, -0.938280, -0.003412, 0.938281, -0.345869,  0.001950, -0.003010, -0.002527,  0.999992,
      -0.345861,  0.938280,  0.003412, 0.938281,  0.345869, -0.001950, -0.003010,  0.002527, -0.999992};
    
    double translationMatrix[4][3] = {0.367264,    446.508738,  175.97185+0.3,
      1.078181,    445.826258, -174.026758+0.3,
      -153.843916, 418.304256,  175.956905+0.8,
      -152.649580, 417.621779, -174.040392+0.8};
    for(int j=0; j<4; j++)
    {
      matrix[j] = new TGeoHMatrix();
      matrix[j]->SetRotation(rotationMatrix[j]);
      matrix[j]->SetTranslation(translationMatrix[j]);
      matrix[j]->Print();
      pi0calib->SetGeometryMatrixInSM(matrix[j],j);
    }
    
    
    pi0calib->SwitchOnLoadOwnGeometryMatrices();
*/
    
    pi0calib->SwitchOnClusterCorrection();
    AliEMCALRecoUtils * reco = pi0calib->GetEMCALRecoUtils();
    reco->SetParticleType(AliEMCALRecoUtils::kPhoton);
    reco->SetW0(4.5);
    
    //---------------------
    // Geometry alignment
    //---------------------

    //reco->SetPositionAlgorithm(AliEMCALRecoUtils::kUnchanged);
    reco->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);
    
    //reco->SetMisalTransShift(0,1.134);   reco->SetMisalTransShift(1,8.2); reco->SetMisalTransShift(2,1.197);
    //reco->SetMisalTransShift(3,-3.093);  reco->SetMisalTransShift(4,6.82);reco->SetMisalTransShift(5,1.635);
    
    //reco->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerIndex);
    //reco->SetMisalTransShift(0,1.08);   reco->SetMisalTransShift(1,8.35); reco->SetMisalTransShift(2,1.12);
    //reco->SetMisalRotShift(3,-8.05);    reco->SetMisalRotShift(4,8.05);  
    //reco->SetMisalTransShift(3,-0.42);  reco->SetMisalTransShift(5,1.55);
    
    //---------------------
    // Non linearity
    //---------------------
    
    reco->SetNonLinearityFunction(AliEMCALRecoUtils::kNoCorrection);
    //reco->SetNonLinearityFunction(AliEMCALRecoUtils::kPi0GammaGamma);
    //reco->SetNonLinearityParam(0,1.04);     reco->SetNonLinearityParam(1,-0.1445);
    //reco->SetNonLinearityParam(2,1.046);    
    
    //     reco->SetNonLinearityFunction(AliEMCALRecoUtils::kPi0GammaConversion);
    //     reco->SetNonLinearityParam(0,1.033);     reco->SetNonLinearityParam(1,0.0566186);
    //     reco->SetNonLinearityParam(2,0.982133);    
    
    
    //      reco->SetNonLinearityFunction(AliEMCALRecoUtils::kPi0MC);
    //      reco->SetNonLinearityParam(0,1.001);     reco->SetNonLinearityParam(1,-0.01264);
    //      reco->SetNonLinearityParam(2,-0.03632);    
    //      reco->SetNonLinearityParam(3,0.1798);     reco->SetNonLinearityParam(4,-0.522);
    
    
    //---------------------
    // Recalibration
    //---------------------

    // 
    //     reco->SwitchOnRecalibration();
    //      TFile * f = new TFile("RecalibrationFactors.root","read");
    //      TH2F * h0 = (TH2F*)f->Get("EMCALRecalFactors_SM0")->Clone();
    //      TH2F * h1 = (TH2F*)f->Get("EMCALRecalFactors_SM1")->Clone();
    //      TH2F * h2 = (TH2F*)f->Get("EMCALRecalFactors_SM2")->Clone();
    //      TH2F * h3 = (TH2F*)f->Get("EMCALRecalFactors_SM3")->Clone();
    
    //      reco->SetEMCALChannelRecalibrationFactors(0,h0);
    //      reco->SetEMCALChannelRecalibrationFactors(1,h1);
    //      reco->SetEMCALChannelRecalibrationFactors(2,h2);
    //      reco->SetEMCALChannelRecalibrationFactors(3,h3);
    
    //---------------------
    // Bad channels removal
    //---------------------

/*    
    reco->SwitchOnBadChannelsRemoval();
    
    // SM0
    reco->SetEMCALChannelStatus(0,3,13);  reco->SetEMCALChannelStatus(0,44,1); reco->SetEMCALChannelStatus(0,3,13); 
    reco->SetEMCALChannelStatus(0,20,7);  reco->SetEMCALChannelStatus(0,38,2);   
    // SM1
    reco->SetEMCALChannelStatus(1,4,7);   reco->SetEMCALChannelStatus(1,4,13);  reco->SetEMCALChannelStatus(1,9,20); 
    reco->SetEMCALChannelStatus(1,14,15); reco->SetEMCALChannelStatus(1,23,16); reco->SetEMCALChannelStatus(1,32,23); 
    reco->SetEMCALChannelStatus(1,37,5);  reco->SetEMCALChannelStatus(1,40,1);  reco->SetEMCALChannelStatus(1,40,2);
    reco->SetEMCALChannelStatus(1,40,5);  reco->SetEMCALChannelStatus(1,41,0);  reco->SetEMCALChannelStatus(1,41,1);
    reco->SetEMCALChannelStatus(1,41,2);  reco->SetEMCALChannelStatus(1,41,4);
    // SM2        
    reco->SetEMCALChannelStatus(2,14,15); reco->SetEMCALChannelStatus(2,18,16); reco->SetEMCALChannelStatus(2,18,17); 
    reco->SetEMCALChannelStatus(2,18,18); reco->SetEMCALChannelStatus(2,18,20); reco->SetEMCALChannelStatus(2,18,21); 
    reco->SetEMCALChannelStatus(2,18,23); reco->SetEMCALChannelStatus(2,19,16); reco->SetEMCALChannelStatus(2,19,17); 
    reco->SetEMCALChannelStatus(2,19,19); reco->SetEMCALChannelStatus(2,19,20); reco->SetEMCALChannelStatus(2,19,21); 
    reco->SetEMCALChannelStatus(2,19,22);
    //SM3
    reco->SetEMCALChannelStatus(3,4,7);
*/
    
    reco->SetNumberOfCellsFromEMCALBorder(1);
    
    //reco->Print("");
    
    pi0calib->PrintInfo();
    mgr->AddTask(pi0calib);
    
    AliAnalysisDataContainer *coutput2 = 
    mgr->CreateContainer("pi0calib", TList::Class(), AliAnalysisManager::kOutputContainer, "pi0calib.root");
    
    AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer("Cuts", TList::Class(), 
                                                               AliAnalysisManager::kOutputContainer, "pi0calib.root");
    
    mgr->ConnectInput  (pi0calib,     0, cinput1);
    mgr->ConnectOutput (pi0calib, 1, coutput2 );
    mgr->ConnectOutput (pi0calib, 2, cout_cuts);
    
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
    
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libANALYSISalice.so");
    //TGeoManager::Import("geometry.root") ; //need file "geometry.root" in local dir!!!!
    gSystem->Load("libPWGGAEMCALTasks.so");
    //SetupPar("PWGGAEMCALTasks");
    
    //--------------------------------------------------------
    //If you want to use root and par files from aliroot
    //--------------------------------------------------------  
    /*     
	   SetupPar("STEERBase");
	   SetupPar("ESD");
	   SetupPar("AOD");
	   SetupPar("ANALYSIS");
	   SetupPar("ANALYSISalice");  
	   SetupPar("PHOSUtils");
	   SetupPar("EMCALUtils");
	   //Create Geometry
	   TGeoManager::Import("geometry.root") ; //need file "geometry.root" in local dir!!!!
	   SetupPar("PWGGAEMCALTasks");
     */
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



void CreateChain(const anaModes mode, TChain * chain){
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
      else if(kInputData == "AOD") datafile = "aod.root" ;
      else if(kInputData == "MC")  datafile = "galice.root" ;
      
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
    }
  }// xml analysis
  
  gSystem->ChangeDirectory(ocwd.Data());
}

