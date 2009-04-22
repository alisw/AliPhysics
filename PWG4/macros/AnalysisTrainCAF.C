//______________________________________________________________________________
void AnalysisTrainCAF(Int_t nEvents = 10000, Int_t nOffset = 0, char *ds = "/PWG4/kleinb/LHC08r_jetjet50")
{
  // Example of running analysis train in CAF. To run in debug mode:
  //  - export ROOTSYS=debug  on your local client
  //  - un-comment gProof->ClearPackages()
  //  - un-comment lines with debugging info


  Bool_t debug         = kTRUE;
  Bool_t useMC         = kTRUE;
  Bool_t readTR        = kFALSE;
  Bool_t bPROOF        = kFALSE;
  Bool_t bOLD          = kFALSE;  // a flag to be compatible with the older AF, to be removed ASA grid and proof are updated


    
  Int_t iAODanalysis   = 0;
  Int_t iAODhandler    = 1;
  Int_t iESDfilter     = 1;  // Only active if iAODanalysis=0
  Int_t iJETAN         = 1;
  Int_t iJETANESD      = 1;
  Int_t iJETANMC       = 1;
  Int_t iDIJETAN       = 1;
  Int_t iPWG4SPECTRUM  = 1;
  Int_t iPWG4UE        = 1;
  Int_t iPWG4PID        = 1;

  if (iAODanalysis) {
    useMC = kFALSE;
    readTR = kFALSE;
    iESDfilter = 0;
    iMUONfilter = 0;
  }    
  if (iJETAN) iESDfilter=1;
  if (iESDfilter) iAODhandler=1;
  
  // Dataset from CAF
  TString dataset(ds);
  TChain *chain = 0;
  // CKB quick hack for local analysis
  gROOT->LoadMacro("CreateESDChain.C");
  TChain *chain = CreateESDChain("tmp.txt",1000);

 
  printf("==================================================================\n");
  printf("===========    RUNNING ANALYSIS TRAIN IN CAF MODE    =============\n");
  printf("==================================================================\n");
  if (iAODanalysis) printf("=  AOD analysis on dataset: %s\n", dataset.Data());
  else              printf("=  ESD analysis on dataset: %s\n", dataset.Data());
  if (iESDfilter)   printf("=  ESD filter                                                    =\n");
  if (iJETAN)       printf("=  Jet analysis from AOD                                         =\n");
  if (iJETANESD)    printf("=  Jet analysis from ESD                                         =\n");
  if (iJETANMC)     printf("=  Jet analysis from Kinematics                                  =\n");
  if (iDIJETAN)     printf("=  DiJet analysis                                                =\n");
  if (iPWG4SPECTRUM)printf("=  PWG4 Jet spectrum analysis                                    =\n");
  if (iPWG4UE)      printf("=  PWG4 UE                                                       =\n");
  printf("==================================================================\n");
  if (useMC) printf(":: use MC    TRUE\n");
  else       printf(":: use MC    FALSE\n");
  if (readTR) printf(":: read TR   TRUE\n");
  else        printf(":: read TR   FALSE\n");
  if (debug) printf(":: debugging TRUE\n");
  else       printf(":: debugging FALSE\n");
    
  // Load common libraries
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  

  // Reset user processes if CAF if not responding anymore
  // TProof::Reset("alicecaf"); 
  // One may enable a different ROOT version on CAF
  
  const char* proofNode = "localhost";
  
  // Connect to proof
  if(bPROOF){
    TProof::Mgr(proofNode)->ShowROOTVersions();
    // TProof::Mgr(proofNode)->SetROOTVersion("v5-21-01-alice_dbg");
    TProof::Open(proofNode); 
    
    // Clear packages if changing ROOT version on CAF or local
    // gProof->ClearPackages();
    // Enable proof debugging if needed
    //    gProof->SetLogLevel(5);
    // To debug the train in PROOF mode, type in a root session:
    // root[0] TProof::Mgr("lxb6064")->GetSessionLogs()->Display("*",0,10000);
    // Common packages
    // --- Enable the STEERBase Package
    gProof->UploadPackage("${ALICE_ROOT}/STEERBase.par");
    gProof->EnablePackage("STEERBase");
    // --- Enable the ESD Package
    gProof->UploadPackage("${ALICE_ROOT}/ESD.par");
    gProof->EnablePackage("ESD");
    // --- Enable the AOD Package
    gProof->UploadPackage("${ALICE_ROOT}/AOD.par");
    gProof->EnablePackage("AOD");
    // --- Enable the ANALYSIS Package
    gProof->UploadPackage("${ALICE_ROOT}/ANALYSIS.par");
    gProof->EnablePackage("ANALYSIS");
    // --- Enable the ANALYSISalice Package
    gProof->UploadPackage("${ALICE_ROOT}/ANALYSISalice.par");
    gProof->EnablePackage("ANALYSISalice");
    
      
    // --- Enable the JETAN Package
    if (iJETAN||iJETANESD||iJETANMC) {
      gProof->UploadPackage("${ALICE_ROOT}/JETAN.par");
      gProof->EnablePackage("JETAN");
    }   
    // --- Enable particle correlation analysis
    if (iPWG4UE||iPWG4SPECTRUM) {
      gProof->UploadPackage("${ALICE_ROOT}/PWG4JetTasks.par");
      gProof->EnablePackage("PWG4JetTasks");
    }   
    
  }
  else{

    //  
    // We are local or on grid
    // access remote files in lcoal case as well so open alien conection

    /*
    printf("*** Connect to AliEn ***\n");
    TGrid::Connect("alien://");

    chain = CreateChainFromCollection("wn.xml","esdTree",2); 
    */

    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    
    
    // --- Enable the JETAN Package
    if (iJETAN||iJETANESD||iJETANMC) gSystem->Load("libJETAN");
    // --- Enable particle correlation analysis
    if (iPWG4UE||iPWG4SPECTRUM)gSystem->Load("libPWG4JetTasks"); 
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
    
    if(bOLD){
      cinput = mgr->CreateContainer("cInput",TChain::Class(),AliAnalysisManager::kInputContainer);
    } 
    else { 
     cinput = mgr->GetCommonInputContainer();
    }
    // This container is managed by the AOD handler
    AliAnalysisDataContainer *cout_aod = 0;
    if (iAODhandler) {
       // AOD output handler
       AliAODHandler* aodHandler   = new AliAODHandler();
       aodHandler->SetFillAOD(kFALSE);
       mgr->SetOutputEventHandler(aodHandler);       
       aodHandler->SetOutputFileName(Form("AliAODs_pwg4_%07d-%07d.root",nOffset,nOffset+nEvents));
       if(bOLD){
	 cout_aod = mgr->CreateContainer("cAOD", TTree::Class(),AliAnalysisManager::kOutputContainer, "default");
       }
       else{
	 cout_aod = mgr->GetCommonOutputContainer();
       }
       cout_aod->SetSpecialOutput();
    }   

    // Debugging if needed
    if (debug) mgr->SetDebugLevel(0);
    //    AliLog::EnableDebug(kTRUE);
    AliLog::SetGlobalLogLevel(0);


    if (iESDfilter && !iAODanalysis) {
      gROOT->LoadMacro("AddTaskESDfilter.C");
      AliAnalysisTaskESDfilter *esdfilter = 0;
      if(bOLD)esdfilter = AddTaskESDfilter(mgr,cinput,cout_aod);
      else esdfilter = AddTaskESDfilter();
    }   
    // Jet analysis from the AOD
    if (iJETAN) {
      gROOT->LoadMacro("AddTaskJets.C");
      AliAnalysisTaskJets *jetanaAOD = 0;
      if(bOLD)jetanaAOD = AddTaskJets("AOD","UA1",mgr,cinput);
      else AliAnalysisTaskJets *jetanaAOD = AddTaskJets("AOD","UA1");
    }   
    // JETANALYSIS from the ESD
    if (iJETANESD && !iAODanalysis) {
      gROOT->LoadMacro("AddTaskJets.C");
      AliAnalysisTaskJets *jetanaESD = 0;
      if(bOLD)jetanaESD = AddTaskJets("ESD","UA1",mgr,cinput);
      else jetanaESD = AddTaskJets("ESD","UA1");
      jetanaESD->SetDebugLevel(10);
      jetanaESD->SetNonStdBranch("jetsESD");    
    }   
    // Jet analysisMC
    if (iJETANMC && useMC) {
      gROOT->LoadMacro("AddTaskJets.C");
      AliAnalysisTaskJets *jetanaMC = 0;
      if(bOLD)jetanaMC = AddTaskJets("MC","UA1",mgr,cinput);
      else jetanaMC = AddTaskJets("MC","UA1");
      jetanaMC->SetDebugLevel(10);
      jetanaMC->SetNonStdBranch("jetsMC");
    }   
    // Dijet analysis
    if(iDIJETAN){
      gROOT->LoadMacro("AddTaskDiJets.C");
      AliAnalysisTaskDiJets *dijetana = 0;
      if(bOLD) dijetana = AddTaskDiJets(mgr,cinput);
      else dijetana = AddTaskDiJets();
    }
    if (iPWG4SPECTRUM) {
      gROOT->LoadMacro("AddTaskJetSpectrum.C");
      AliAnalysisTaskJetSpectrum* pwg4spec = 0;
      if(bOLD)pwg4spec = AddTaskJetSpectrum(mgr,cinput);
      else pwg4spec = AddTaskJetSpectrum();
    }   
    if (iPWG4UE) {
      gROOT->LoadMacro("AddTaskUE.C");
      AliAnalysisTaskUE* ueana = 0;
      if(bOLD)ueana = AddTaskUE(mgr,cinput);
      else ueana = AddTaskUE();
    }   
    if(iPWG4PID){
      gROOT->LoadMacro("AddTaskPWG4PidDetEx.C");
      AliAnalysisTaskPWG4PidDetEx *taskPid = 0;
      if(bOLD)taskPid = AddTaskPWG4PidDetEx(mgr,cinput);
      else taskPid = AddTaskPWG4PidDetEx();
      taskPid->SetDebugLevel(10);
  }
    // Run the analysis
    //    
    if (mgr->InitAnalysis()) {
       mgr->PrintStatus();
       if(bPROOF)mgr->StartAnalysis("proof",dataset.Data(), nEvents,nOffset);
       else mgr->StartAnalysis("local",chain);
    }   
}
TChain *CreateChainFromCollection(const char* xmlfile, const char *treeName="esdTree",Int_t nFiles = 0)
{
// Create a chain from an alien collection.                                                                          
   TAlienCollection * myCollection  = TAlienCollection::Open(xmlfile);

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
