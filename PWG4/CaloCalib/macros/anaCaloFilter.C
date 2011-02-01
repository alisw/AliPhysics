/* $Id:  $ */
//--------------------------------------------------
// Example macro to do Calorimeters filtering
// copy ESDs into AODs
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
Bool_t kUsePhysSel = kTRUE;

void anaCaloFilter(Int_t mode=mLocal)
{
  // Main
  char cmd[200] ; 
  sprintf(cmd, ".! rm -rf AliAOD.root") ; 
  gROOT->ProcessLine(cmd) ; 
  //--------------------------------------------------------------------
  // Load analysis libraries
  // Look at the method below, 
  // change whatever you need for your analysis case
  // ------------------------------------------------------------------
  LoadLibraries(mode) ;
  //gSystem->ListLibraries(); 
  //gSystem->Unload("libPWG4CaloCalib.so"); 
  //gSystem->Unload("libEMCALUtils.so");
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

    // AOD output handler
    AliAODHandler* aodoutHandler   = new AliAODHandler();
    aodoutHandler->SetOutputFileName("AliAOD.root");
    ////aodoutHandler->SetCreateNonStandardAOD();
    mgr->SetOutputEventHandler(aodoutHandler);
    
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
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
    
    // ESD physics selection task
    if(kInputData == "ESD" && kUsePhysSel){
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
    }
    
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
    
    AliAnalysisTaskCaloFilter * filter = new AliAnalysisTaskCaloFilter("EMCALFilter");
    if(kUsePhysSel)filter->SelectCollisionCandidates(); 
    filter->SetCaloFilter(AliAnalysisTaskCaloFilter::kBoth); //kPHOS or kBoth
    filter->SwitchOnClusterCorrection();
    //filter->SetDebugLevel(10);
    AliEMCALRecoUtils * reco = filter->GetEMCALRecoUtils();
    reco->SetParticleType(AliEMCALRecoUtils::kPhoton);
    reco->SetW0(4.5);
        
    reco->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);

    TGeoHMatrix *matrix[4];
    
    double rotationMatrix[4][9] = {-0.014587, -0.999892, -0.002031, 0.999892, -0.014591,  0.001979, -0.002009, -0.002002,  0.999996,
				 -0.014587,  0.999892,  0.002031, 0.999892,  0.014591, -0.001979, -0.002009,  0.002002, -0.999996,
				 -0.345864, -0.938278, -0.003412, 0.938276, -0.345874,  0.003010, -0.004004, -0.002161,  0.999990,
				 -0.345861,  0.938280,  0.003412, 0.938276,  0.345874, -0.003010, -0.004004,  0.002161, -0.999990};
    
    double translationMatrix[4][3] = {0.351659,    447.576446,  176.269742,
				      1.062577,    446.893974, -173.728870,
				      -154.213287, 419.306156,  176.753692,
				      -153.018950, 418.623681, -173.243605};
    for(int j=0; j<4; j++)
      {
	matrix[j] = new TGeoHMatrix();
	matrix[j]->SetRotation(rotationMatrix[j]);
	matrix[j]->SetTranslation(translationMatrix[j]);
	matrix[j]->Print();
	filter->SetEMCALGeometryMatrixInSM(matrix[j],j);
      }
    
    
    filter->SwitchOnLoadOwnEMCALGeometryMatrices();
    
    reco->SetNonLinearityFunction(AliEMCALRecoUtils::kNoCorrection);

    //Time dependent corrections    
    //Recover file from alien  /alice/cern.ch/user/g/gconesab/TimeDepCorrectionDB
    reco->SwitchOnTimeDepCorrection();
    char cmd[200] ;
    sprintf(cmd, ".!tar xvfz CorrectionFiles.tgz") ;
    gROOT->ProcessLine(cmd) ;

    //Recalibration factors
    //Recover the file from alien  /alice/cern.ch/user/g/gconesab/RecalDB
    reco->SwitchOnRecalibration();
    TFile * f = new TFile("RecalibrationFactors.root","read");
    TH2F * h0 = (TH2F*)f->Get("EMCALRecalFactors_SM0");
    TH2F * h1 = (TH2F*)f->Get("EMCALRecalFactors_SM1");
    TH2F * h2 = (TH2F*)f->Get("EMCALRecalFactors_SM2");
    TH2F * h3 = (TH2F*)f->Get("EMCALRecalFactors_SM3");
    
    reco->SetEMCALChannelRecalibrationFactors(0,h0);
    reco->SetEMCALChannelRecalibrationFactors(1,h1);
    reco->SetEMCALChannelRecalibrationFactors(2,h2);
    reco->SetEMCALChannelRecalibrationFactors(3,h3);
    
    //Bad channels
    //Recover the file from alien  /alice/cern.ch/user/g/gconesab/BadChannelsDB
    reco->SwitchOnBadChannelsRemoval();
    reco->SwitchOnDistToBadChannelRecalculation();
    TFile * fbad = new TFile("BadChannels.root","read");
    TH2I * hbad0 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod0");
    TH2I * hbad1 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod1");
    TH2I * hbad2 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod2");
    TH2I * hbad3 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod3");
    reco->SetEMCALChannelStatusMap(0,hbad0);
    reco->SetEMCALChannelStatusMap(1,hbad1);
    reco->SetEMCALChannelStatusMap(2,hbad2);
    reco->SetEMCALChannelStatusMap(3,hbad3);

    //reco->Print("");
    filter->PrintInfo(); 
    mgr->AddTask(filter);
        
    //AliAnalysisDataContainer *cout_cuts2 = mgr->CreateContainer("Cuts", TList::Class(), 
    //					       AliAnalysisManager::kOutputContainer, "pi0calib.root");
    
    mgr->ConnectInput  (filter,  0, cinput1);
    mgr->ConnectOutput (filter, 0, coutput1 );
    //mgr->ConnectOutput (filter, 2, cout_cuts2);

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
  
  sprintf(cmd, ".! rm -rf CorrectionFiles") ;
  
}

void  LoadLibraries(const anaModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");
  gSystem->Load("libMatrix.so");
  gSystem->Load("libPhysics.so");  

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
    gSystem->Load("libPHOSUtils.so"); 
    //gSystem->Load("libEMCALUtils.so");  
    SetupPar("EMCALUtils"); 
    SetupPar("PWG4CaloCalib"); 
    
    //TGeoManager::Import("geometry.root") ; //need file "geometry.root" in local dir!!!!
    //gSystem->Load("libPWG4CaloCalib.so");
   
    /*
      //     gSystem->Load("libPWG4omega3pi.so");
      //     gSystem->Load("libCORRFW.so");
      //     gSystem->Load("libPWG3base.so");
      //     gSystem->Load("libPWG3muon.so");
 */
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
	   SetupPar("PWG4CaloCalib");
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
      else if(kInputData == "AOD") datafile = "AliAOD.root" ;
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

