enum analysisModes {mLocal,mLocalPAR,mPROOF,mGrid,mGridPAR};
enum analysisTypes {mESD,mAOD,mMC,mMCESD};

//
class AliAnalysisGrid;
class AliAnalysisTaskBF;
class AliBalance;

//________________________________________________________________________//
void runBalanceFunctionInpp(Int_t mode = mLocal, 
			    Int_t type = mAOD,
			    Bool_t DATA = kFALSE) {
  // Time:
  TStopwatch timer;
  timer.Start();
  
  //Check analysis mode
  if((mode < 0) || (mode > 4)) {
    Printf("Analysis mode not recognized!");
    Printf("You can select out of 0: local, 1: local with par files, 2: proof, 3: grid, 4: grid with par files");
    return;
  }
  
  //Check analysis type
  if((type < 0) || (type > 3)) {
    Printf("Analysis type not recognized!");
    Printf("You can select out of 0: ESD, 1: AOD, 2: MC (stack), 3: MC (from the ESD)");
    return;
  }
  
  // Load needed libraries:
  LoadLibraries(mode);
  
 // Create and configure the AliEn plug-in:
  if(mode == mGrid || mode == mGridPAR) {
    gROOT->LoadMacro("CreateAlienHandler.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandler(runListFileName);  
    if (!alienHandler) return;
  }
  // Chains:   
  if(mode == mLocal || mode == mLocalPAR) {
    TChain* chain = 0x0;
    if((type == mESD)||(type == mMCESD))  
      chain = new TChain("esdTree");
    else if(type == mAOD)
      chain = new TChain("aodTree");
    else if(type == mMC)
      chain = new TChain("TE");

    TString filename;
    for(Int_t i = 1; i < 20; i++) {
      filename = "/data/alice2/pchrist/pp/LHC10c/0.9TeV/Data/";
      filename += "Set"; filename += i; 
      if((type == mESD)||(type == mMCESD))  
	filename += "/AliESDs.root";
      else if(type == mAOD)
	filename += "/AliAOD.root";
      else if(type == mMC)
     	filename += "/galice.root";

      chain->Add(filename.Data());
    }
  }
  //Proof
  if(mode == mPROOF) {
    gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_ROOT")));
  }
  
  // analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager("balanceFunctionManager");
  if(mode == mGrid || mode == mGridPAR)
    mgr->SetGridHandler(alienHandler);
    
  // input handler (ESD or AOD)
  AliVEventHandler* inputH = NULL;
  if((type == mESD)||(type == mMCESD))  
    inputH = new AliESDInputHandler();
  else if(type == mAOD)
    inputH = new AliAODInputHandler();
  mgr->SetInputEventHandler(inputH);
    
  // mc event handler
  if((type == mMC) || (type == mMCESD)) {
    AliMCEventHandler* mchandler = new AliMCEventHandler();
    // Not reading track references
    mchandler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mchandler);
  }   

  if(type != mAOD){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
  }

  //Add the BF task (all centralities)
  gROOT->LoadMacro("AddTaskBalanceFunctionInpp.C"); 
  AddTaskBalanceFunctionInpp();

  // enable debug printouts
  mgr->SetDebugLevel(2);
  mgr->SetUseProgressBar(1,100);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  // start analysis
  if(mode == mLocal || mode == mLocalPAR) 
    mgr->StartAnalysis("local",chain);
  else if(mode == mPROOF) 
    mgr->StartAnalysis("proof",dataDir,nRuns,offset);
  else if(mode == mGrid || mode == mGridPAR) 
    mgr->StartAnalysis("grid");

  // Print real and CPU time used for analysis:  
  timer.Stop();
  timer.Print();
}

//=============================================================//
void LoadLibraries(const analysisModes mode) {  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libCore.so");        
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");

  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (mode==mLocal || mode==mGrid || mode == mGridPAR) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libPWG2ebye.so");
    // Use AliRoot includes to compile our task
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
  }
  
  else if (mode == mLocalPAR) {
    //--------------------------------------------------------
    //If you want to use root and par files from aliroot
    //--------------------------------------------------------  
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
    SetupPar("PWG2ebye");
}
  
  //---------------------------------------------------------
  // <<<<<<<<<< PROOF mode >>>>>>>>>>>>
  //---------------------------------------------------------
  else if (mode==mPROOF) {
    // Connect to proof
    printf("*** Connect to PROOF ***\n");
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    // Put appropriate username here
    TProof::Open("alice-caf.cern.ch");
    //TProof::Open("skaf.saske.sk");
    //TProof::Open("prf000-iep-grid.saske.sk");

    gProof->EnablePackage("VO_ALICE@AliRoot::v4-21-12-AN");
  }  
  
} // end of void LoadLibraries(const anaModes mode)

//======================================================================//
void SetupPar(char* pararchivename) {
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
    processline = Form(".! mv /tmp/%s .", parpar.Data()) ;
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

} // end of void SetupPar(char* pararchivename) 
