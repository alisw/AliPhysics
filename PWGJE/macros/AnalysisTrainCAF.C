//______________________________________________________________________________
void AnalysisTrainCAF(Int_t nEvents = 10000, Int_t nOffset = 0, char *ds = "/PWG4/kleinb/LHC09a1_test500")
{
  // Example of running analysis train in CAF. To run in debug mode:
  //  - export ROOTSYS=debug  on your local client
  //  - un-comment gProof->ClearPackages()
  //  - un-comment lines with debugging info


  Bool_t debug         = kTRUE;
  Bool_t useMC         = kTRUE;
  Bool_t readTR        = kFALSE;
  Bool_t bPROOF        = kFALSE;
  Bool_t bLOCALPAR     = kFALSE;  // flag that swtiches on loading of local par files insead of loading libs, needed for grid and local testing

    
  Int_t iAODanalysis   = 1;
  Int_t iAODhandler    = 1;
  Int_t iESDfilter     = 1;  // Only active if iAODanalysis=0
  Int_t iJETAN         = 1;
  Int_t iJETANESD      = 0;
  Int_t iJETANMC       = 0;
  Int_t iJETANMC2       = 0;
  Int_t iFASTJET     = 1;
  Int_t iDIJETAN       = 0;
  Int_t iPWG4SPECTRUM  = 0;
  Int_t iPWG4JFSYSTEMATICS  = 0;
  Int_t iPWG4JETCORRECTION  = 0;
  Int_t iPWG4THREEJETS  = 0;
  Int_t iPWG4UE        = 0;
  Int_t iPWG4PID        = 0;

  if (iAODanalysis) {
    useMC = kFALSE;
    readTR = kFALSE;
    iESDfilter = 0;
  }    
  if (iJETAN) iESDfilter=1;
  if (iESDfilter) iAODhandler=1;
  
  // Dataset from CAF
  TString dataset(ds);
  TChain *chain = 0;
  // CKB quick hack for local analysis
  gROOT->LoadMacro("CreateESDChain.C");
  TChain *chain = CreateChain("aodTree",ds,1);
  //  TChain *chain = CreateChain("esdTree",ds,100);
  //  chain = new TChain("aodTree");
  //  chain->Add("/Users/kleinb/bigdisk/1/LHC09a3/001/AliAOD.root");

 
  printf("==================================================================\n");
  printf("===========    RUNNING ANALYSIS TRAIN IN CAF MODE    =============\n");
  printf("==================================================================\n");
  if (iAODanalysis) printf("=  AOD analysis on dataset: %s\n", dataset.Data());
  else              printf("=  ESD analysis on dataset: %s\n", dataset.Data());
  if (iESDfilter)   printf("=  ESD filter                                                     =\n");
  if (iJETAN)       printf("=  Jet analysis from AOD                                         =\n");
  if (iJETANESD)    printf("=  Jet analysis from ESD                                         =\n");
  if (iJETANMC)     printf("=  Jet analysis from Kinematics                                  =\n");
  if (iJETANMC2)     printf("=  Jet analysis 2 from Kinematics                               =\n");
  if (iFASTJET)     printf("=  Loading FastJet                               =\n");
  if (iDIJETAN)     printf("=  DiJet analysis                                                 =\n");
  if (iPWG4SPECTRUM)printf("=  PWG4 Jet spectrum analysis                                    =\n");
  if (iPWG4JFSYSTEMATICS)printf("=  PWG4 Jet Finder systematics                                   =\n");
  if (iPWG4JETCORRECTION)printf("=  PWG4 Jet Correction                                   =\n");
  if (iPWG4THREEJETS)printf("=  PWG4 Three Jets                                   =\n");

  if (iPWG4UE)      printf("=  PWG4 UE                                                        =\n");
  printf("==================================================================\n");
  if (useMC) printf(":: use MC    TRUE\n");
  else       printf(":: use MC    FALSE\n");
  if (readTR) printf(":: read TR   TRUE\n");
  else        printf(":: read TR   FALSE\n");
  if (debug) printf(":: debugging TRUE\n");
  else       printf(":: debugging FALSE\n");
    
  // Load common libraries
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  if(iFASTJET){
    gSystem->Load("libCGAL");
    gSystem->Load("libfastjet");
    gSystem->Load("libsiscone");
    gSystem->Load("libSISConePlugin");
  }


  // Reset user processes if CAF if not responding anymore
  // TProof::Reset("alicecaf"); 
  // One may enable a different ROOT version on CAF
  
  //  const char* proofNode = "localhost";
  const char* proofNode = "alicecaf";
  
  // Connect to proof
  if(bPROOF){
    TProof::Mgr(proofNode)->ShowROOTVersions();
    // TProof::Mgr(proofNode)->SetROOTVersion("v5-21-01-alice_dbg");
    TProof::Open(proofNode); 
    
    // Clear packages if changing ROOT version on CAF or local
     gProof->ClearPackages();
    // Enable proof debugging if needed
    //    gProof->SetLogLevel(5);
    // To debug the train in PROOF mode, type in a root session:
    // root[0] TProof::Mgr("lxb6064")->GetSessionLogs()->Display("*",0,10000);
    // Common packages
    // --- Enable the STEERBase Package
    gProof->UploadPackage("STEERBase.par");
    gProof->EnablePackage("STEERBase");	   
    // --- Enable the ESD Package	   
    gProof->UploadPackage("ESD.par");	   
    gProof->EnablePackage("ESD");	   
    // --- Enable the AOD Package	   
    gProof->UploadPackage("AOD.par");	   
    gProof->EnablePackage("AOD");	   
    // --- Enable the ANALYSIS Package	   
    gProof->UploadPackage("ANALYSIS.par"); 
    gProof->EnablePackage("ANALYSIS");	   
    // --- Enable the ANALYSISalice Package
    gProof->UploadPackage("ANALYSISalice.par");
    gProof->EnablePackage("ANALYSISalice");
    
      
    // --- Enable the JETAN Package
    if (iJETAN||iJETANESD||iJETANMC||iJETANMC2) {
      gProof->UploadPackage("JETAN.par");
      gProof->EnablePackage("JETAN");
      if(iFASTJET){
	gProof->UploadPackage("FASTJETAN.par");
	gProof->EnablePackage("FASTJETAN"); 
      }
    }   
    // --- Enable particle correlation analysis
    if (iPWG4UE||iPWG4SPECTRUM||iPWG4JFSYSTEMATICS||iPWG4JETCORRECTION||iPWG4THREEJETS) {
      gProof->UploadPackage("JETAN.par");
      gProof->EnablePackage("JETAN");
      gProof->UploadPackage("PWG4JetTasks.par");
      gProof->EnablePackage("PWG4JetTasks");
    }   
    
  }
  else{

    //  
    // We are local or on grid
    // access remote files in lcoal case as well so open alien connection

    /*
    printf("*** Connect to AliEn ***\n");
    TGrid::Connect("alien://");

    chain = CreateChainFromCollection("wn.xml","esdTree",2); 
    */

    if(bLOCALPAR){
      SetupPar("STEERBase");
      SetupPar("ESD");	 
      SetupPar("AOD");	   
      SetupPar("ANALYSIS"); 
      SetupPar("ANALYSISalice");
      if (iJETAN||iJETANESD||iJETANMC||iJETANMC2){
	SetupPar("JETAN");	   
	if(iFASTJET)	SetupPar("FASTJETAN");	   
      }
      if (iPWG4UE||iPWG4SPECTRUM||iPWG4JFSYSTEMATICS){
	SetupPar("JETAN");	   
	SetupPar("PWG4JetTasks");
      }
    }
    else{
      Printf("Loading Local libs");
      gSystem->Load("libSTEERBase");
      gSystem->Load("libESD");
      gSystem->Load("libAOD");
      gSystem->Load("libANALYSIS");
      gSystem->Load("libANALYSISalice");  
      // --- Enable the JETAN Package
      if (iJETAN||iJETANESD||iJETANMC||iJETANMC2){
	gSystem->Load("libJETAN");
	if(iFASTJET)gSystem->Load("libFASTJETAN");
      }
      // --- Enable particle correlation analysis
      if (iPWG4UE||iPWG4SPECTRUM||iPWG4JFSYSTEMATICS||iPWG4THREEJETS){
	gSystem->Load("libJETAN");
	gSystem->Load("libPWG4JetTasks"); 
      }
    }

  }


    // Make the analysis manager
    AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "A test setup for the analysis train");
    if (iAODanalysis) {
    // AOD input handler
       AliAODInputHandler *aodH = new AliAODInputHandler();
       mgr->SetInputEventHandler(aodH);
    } else {   
    // ESD input handler
       AliESDInputHandler *esdHandler = new AliESDInputHandler();
       mgr->SetInputEventHandler(esdHandler);
//       esdHandler->SetInactiveBranches("FMD CaloCluster");
    }
    // Monte Carlo handler
    if (useMC && !iAODanalysis) {
       AliMCEventHandler* mcHandler = new AliMCEventHandler();
       mgr->SetMCtruthEventHandler(mcHandler);
       mcHandler->SetReadTR(readTR); 
    }   
    // Top container for input 
    AliAnalysisDataContainer *cinput = 0;
    
    cinput = mgr->GetCommonInputContainer();

    // This container is managed by the AOD handler
    AliAnalysisDataContainer *cout_aod = 0;
    if (iAODhandler) {
       // AOD output handler
       AliAODHandler* aodHandler   = new AliAODHandler();
       //      aodHandler->SetFillAOD(kFALSE);
       mgr->SetOutputEventHandler(aodHandler);       
       aodHandler->SetOutputFileName(Form("AliAODs_pwg4_%07d-%07d.root",nOffset,nOffset+nEvents));
       cout_aod = mgr->GetCommonOutputContainer();
       cout_aod->SetSpecialOutput();
    }   

    // Debugging if needed
    if (debug) mgr->SetDebugLevel(10);
    //    AliLog::EnableDebug(kTRUE);
    AliLog::SetGlobalLogLevel(1);


    if (iESDfilter && !iAODanalysis) {
      gSystem->Load("libCORRFW");
      gSystem->Load("libPWGmuon");

      gROOT->LoadMacro(Form("%s/ANALYSIS/macros/AddTaskESDFilter.C",gSystem->ExpandPathName("${ALICE_ROOT}")));
      //      gROOT->LoadMacro("AddTaskESDfilter.C");
      AliAnalysisTaskESDfilter *esdfilter = AddTaskESDFilter();
      Printf("esdFilter %p",esdfilter); 
    }   
    // Jet analysis from the AOD
    if (iJETAN) {
      gROOT->LoadMacro("AddTaskJets.C");
      //      AliAnalysisTaskJets *jetanaAOD  = AddTaskJets("AOD","UA1",0.4);
      //      AliAnalysisTaskJets *jetanaAOD  = AddTaskJets("AOD","UA1",0.4);
      //      jetanaAOD->SetNonStdBranch("jetsAOD_UA1");    
      AliAnalysisTaskJets *jetanaAOD  = AddTaskJets();
      Int_t i = AddTaskJetsDelta();
    }   
    // JETANALYSIS from the ESD
    if (iJETANESD && !iAODanalysis) {
      gROOT->LoadMacro("AddTaskJets.C");
      AliAnalysisTaskJets *jetanaESD = AddTaskJets("ESD","UA1");
      jetanaESD->SetDebugLevel(0);
      jetanaESD->SetNonStdBranch("jetsESD");    
    }   
    // Jet analysisMC
    if (iJETANMC ){ 
      gROOT->LoadMacro("AddTaskJets.C");
      //      AliAnalysisTaskJets *jetanaMC =  AddTaskJets("AODMC","UA1",0.4);
      AliAnalysisTaskJets *jetanaMC =  AddTaskJets("AODMC","UA1",0.4);
      jetanaMC->SetDebugLevel(0);
      jetanaMC->SetNonStdBranch("jetsMC_UA1");
    }   
    if (iJETANMC2 ){ 
      gROOT->LoadMacro("AddTaskJets.C");
      //      AliAnalysisTaskJets *jetanaMC2 = AddTaskJets("AODMC2","UA1",0.4);
      AliAnalysisTaskJets *jetanaMC2 = AddTaskJets("AODMC2","UA1",0.4);
      jetanaMC2->SetDebugLevel(0);
      jetanaMC2->SetNonStdBranch("jetsMC2_UA1");
    }   
    // Dijet analysis
    if(iDIJETAN){
      gROOT->LoadMacro("AddTaskDiJets.C");
      AliAnalysisTaskDiJets *dijetana  = AddTaskDiJets();
    }
    if (iPWG4SPECTRUM) {
      gROOT->LoadMacro("AddTaskJetSpectrum2.C");
      AliAnalysisTaskJetSpectrum2* pwg4spec = AddTaskJetSpectrum2();
      pwg4spec->SetAODInput(kTRUE);
      pwg4spec->SetBranchRec("jets");
      pwg4spec->SetAnalysisType(0);
      pwg4spec->SetDebugLevel(0);
    }   
    if (iPWG4JFSYSTEMATICS) {
      gROOT->LoadMacro("AddTaskJFSystematics.C");
      AliAnalysisTaskJFSystematics* pwg4jfs = AddTaskJFSystematics("jetsMC","jets");
      pwg4jfs->SetAODInput(kTRUE);
      pwg4jfs->SetDebugLevel(0);
    }   
    if (iPWG4JETCORRECTION) {
      gROOT->LoadMacro("AddTaskJetCorrections.C");
      AliAnalysisTaskJetCorrections* pwg4jc = AddTaskJetCorrections();
      pwg4jc->SetDebugLevel(11);
    }   
    if (iPWG4THREEJETS) {
      gROOT->LoadMacro("AddTaskThreeJets.C");
      AliAnalysisTaskThreeJets* pwg4jjj = AddTaskThreeJets();
      pwg4jjj->SetDebugLevel(11);
    }   
    if (iPWG4UE) {
      gROOT->LoadMacro("AddTaskUE.C");
      AliAnalysisTaskUE* ueana = AddTaskUE();
    }   
    if(iPWG4PID){
      gROOT->LoadMacro("AddTaskPWG4PidDetEx.C");
      AliAnalysisTaskPWG4PidDetEx *taskPid = AddTaskPWG4PidDetEx();
      taskPid->SetDebugLevel(0);
  }
    // Run the analysis
    //    
    if (mgr->InitAnalysis()) {
      mgr->PrintStatus();
      Printf("Chain with %d entries",chain->GetEntries()); 
       if(bPROOF)mgr->StartAnalysis("proof",dataset.Data(), nEvents,nOffset);
       else mgr->StartAnalysis("local",chain,nEvents);
    }   
}
TChain *CreateChainFromCollection(const char* xmlfile, const char *treeName="esdTree",Int_t nFiles = 0)
{
// Create a chain from an alien collection.                                                                          
   TGridCollection * myCollection  = gGrid->OpenCollection(xmlfile);

   if (!myCollection) {
      ::Error("CreateChainSingle", "Cannot create an AliEn collection from %s", xmlfile) ;
     return NULL ;
   }

  TChain* chain = new TChain(treeName);
  myCollection->Reset() ;
  Int_t iCount = 0;
  while ( myCollection->Next() ){
    if(nFiles!=0)iCount++;
    if(iCount > nFiles)break;
    chain->Add(myCollection->GetTURL("")) ;
    Printf("Adding %s",myCollection->GetTURL(""));
  }
  chain->ls();
  return chain;
}


void SetupPar(char* pararchivename)
{
  //Load par files, create analysis libraries                                                         
  //For testing, if par file already decompressed and modified                                        
  //classes then do not decompress.                                                                   

  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ;
  TString parpar(Form("%s.par", pararchivename)) ;
  /*
  if ( gSystem->AccessPathName(parpar.Data()) ) {
    gSystem->ChangeDirectory(gSystem->Getenv("ALICE_ROOT")) ;
    TString processline(Form(".! make %s", parpar.Data())) ;
    gROOT->ProcessLine(processline.Data()) ;
    gSystem->ChangeDirectory(cdir) ;
    processline = Form(".! mv /tmp/%s .", parpar.Data()) ;
    gROOT->ProcessLine(processline.Data()) ;
  }
  */

  if (!gSystem->AccessPathName(pararchivename) ) {
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
