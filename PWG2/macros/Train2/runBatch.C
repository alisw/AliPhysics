void runBatch() {
  TStopwatch timer;
  timer.Start();

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");
  gSystem->Load("libPhysics.so");

  int runFemto = 1;
  int runSpectraProtons = 1;
  int runSpectraV0 = 1;
  int runFlow = 1;
  int runResonances = 1;
  int runEvChar = 1;
  int runKink = 1;
  int runUnicor = 1;
  int runFMDanalysis = 1;

  //____________________________________________________//
  //_____________Setting up STEERBase.par_______________//
  //____________________________________________________//
  //  setupPar("STEERBase");
  gSystem->Load("libSTEERBase.so");

  //____________________________________________________//
  //_____________Setting up ESD.par_____________________//
  //____________________________________________________//
  //  setupPar("ESD");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");

  //____________________________________________________//
  //_____________Setting up AOD.par_____________________//
  //____________________________________________________//
  //  setupPar("AOD");
  gSystem->Load("libAOD.so");

  //_________________________________________________________//
  //_____________Setting up ANALYSIS.par_____________________//
  //_________________________________________________________//
  //  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS.so");

  //_________________________________________________________//
  //_____________Setting up ANALYSISalice.par________________//
  //_________________________________________________________//
  //  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice.so");

  //____________________________________________________//
  //_____________Setting up CORRFW library______________//
  //____________________________________________________//
  //  setupPar("CORRFW");
  gSystem->Load("libCORRFW.so");
  
  //____________________________________________________//
  //_____________Setting up PWG2AOD.par_________________//
  //____________________________________________________//
  setupPar("PWG2AOD");
  gSystem->Load("libPWG2AOD.so");
  
  if (runFemto) {
    //____________________________________________________//
    //_____________Setting up PWG2femtoscopy.par__________//
    //____________________________________________________//
    setupPar("PWG2femtoscopy");
    gSystem->Load("libPWG2femtoscopy.so");
    
    //____________________________________________________//
    //_____________Setting up PWG2femtoscopyUser.par______//
    //____________________________________________________//
    setupPar("PWG2femtoscopyUser");
    gSystem->Load("libPWG2femtoscopyUser.so");
  }
  
  if (runSpectraProtons) {
    //____________________________________________________//
    //_____________Setting up PWG2spectra library ________//
    //____________________________________________________//
    gSystem->Load("libPWG2spectra.so");
  }

  if (runSpectraV0) {
    //____________________________________________________//
    //_____________Setting up PWG2spectra library ________//
    //____________________________________________________//
    gSystem->Load("libPWG2spectra.so");
  }

  if (runFlow) {
    //____________________________________________________//
    //_____________Setting up PWG2flowCommon.par__________//
    //____________________________________________________//
    setupPar("PWG2flowCommon");
    gSystem->Load("libPWG2flowCommon.so");

    //____________________________________________________//
    //_____________Setting up PWG2flowTasks.par___________//
    //____________________________________________________//
    setupPar("PWG2flowTasks");
    gSystem->Load("libPWG2flowTasks.so");
  }

  if (runResonances) {
    //____________________________________________________//
    //_____________Setting up PWG2resonances.par__________//
    //____________________________________________________//
    setupPar("PWG2resonances");
    gSystem->Load("libPWG2resonances.so");
  }
  
  if (runEvChar) {
    //____________________________________________________//
    //_____________Setting up PWG2evchar library__________//
    //____________________________________________________//
    setupPar("PWG2evchar");
    gSystem->Load("libPWG2evchar.so");
  }
  
  if (runKink) {
    //____________________________________________________//
    //_____________Setting up PWG2evchar library__________//
    //____________________________________________________//
    setupPar("PWG2kink");
    gSystem->Load("libPWG2kink.so");
  }
  
  if (runUnicor) {
    //____________________________________________________//
    //_____________Setting up PWG2evchar library__________//
    //____________________________________________________//
    setupPar("PWG2unicor");
    gSystem->Load("libPWG2unicor.so");
  }
  
  if (runFMDanalysis) {
    //____________________________________________________//
    //_____________Setting up PWG2evchar library__________//
    //____________________________________________________//
    setupPar("FMDanalysis");
    gSystem->Load("libFMDanalysis.so");
  }
  
  //ANALYSIS PART
  const char *collectionfile="wn.xml";
  //____________________________________________//
  //Usage of event tags
  AliTagAnalysis *analysis = new AliTagAnalysis();
  TChain *chain = 0x0;
  chain = analysis->GetChainFromCollection(collectionfile,"esdTree");
  
  //  const char *collectionfile="../../LHC09a4/esd.LHC09a4.81305.mini.list";
  //  gROOT->LoadMacro("CreateESDChain.C");
  //  chain = CreateESDChain(collectionfile);
  
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliESDInputHandler* esdH = new AliESDInputHandler;
  //  esdH->SetInactiveBranches("FMD CaloCluster");
  mgr->SetInputEventHandler(esdH);  

  AliMCEventHandler *mcH = new AliMCEventHandler;
  mgr->SetMCtruthEventHandler(mcH);
  if (runFemto) {
    //____________________________________________//
    // 1st task - FEMTOSCOPY
    
    gROOT->LoadMacro("AddTaskFemto.C");
    AliAnalysisTaskFemto *taskfemto = AddTaskFemto();
  }

  if (runSpectraProtons) {
    //____________________________________________//
    // 2nd task - SPECTRA protons
    
    gROOT->LoadMacro("AddTaskProtons.C");
    AliAnalysisTaskProtons *taskprotons = AddTaskProtons();
  }

  if (runFlow) {
    //____________________________________________//
    // 3rd task - FLOW

    // Flow analysis method can be:(set to kTRUE or kFALSE)
    Bool_t SP     = kTRUE;
    Bool_t LYZ1   = kTRUE;
    Bool_t LYZ2   = kFALSE;
    Bool_t LYZEP  = kFALSE;
    Bool_t GFC    = kFALSE;
    Bool_t QC     = kTRUE;
    Bool_t FQD    = kTRUE;
    Bool_t MCEP   = kTRUE;
    
    Bool_t METHODS[] = {SP,LYZ1,LYZ2,LYZEP,GFC,QC,FQD,MCEP};
    
    // Analysis type can be ESD, AOD, MC, ESDMC0, ESDMC1
    const TString type = "ESD";
    
    // Boolean to fill/not fill the QA histograms
    Bool_t QA = kTRUE;   
    
    // Boolean to use/not use weights for the Q vector
    Bool_t WEIGHTS[] = {kFALSE,kFALSE,kFALSE}; //Phi, v'(pt), v'(eta)

    gROOT->LoadMacro("AddTaskFlow.C");
    AliAnalysisTaskFlowEvent* taskFE = AddTaskFlow(type,METHODS,QA,WEIGHTS);
  }

  if (runResonances) {
    //____________________________________________//
    // 4th task - RESONANCES
    
    int useMC = 1;

    gROOT->LoadMacro("AddAnalysisTaskRsn.C");
    //    AliAnalysisTaskFemto *taskfemto = AddTaskFemto();
    AddAnalysisTaskRsn(AliLog::kInfo, "rsn.root", useMC);
  }

  if (runEvChar) {
    //____________________________________________//
    // 5th task - EVENT CHARACTERIZARION
    
    gROOT->LoadMacro("AddTaskSPDdNdEta.C");
    AliAnalysisTaskSPDdNdEta *taskspddndeta = AddTaskSPDdNdEta();
  }

  if (runSpectraV0) {
    //____________________________________________//
    // 6th, 7th, 8th tasks - SPECTRA V0

    // cascades
    gROOT->LoadMacro("AddTaskCheckCascade.C");
    AliAnalysisTaskCheckCascade *taskcheckcascade = AddTaskCheckCascade(0);      

    // v0's
    gROOT->LoadMacro("AddTaskCheckV0.C");
    AliAnalysisTaskCheckV0 *taskcheckV0 = AddTaskCheckV0();

    // strangeness
    gROOT->LoadMacro("AddTaskStrange.C");
    AliAnalysisTaskStrange *taskstrange = AddTaskStrange();
  }

  if (runKink) {
    //____________________________________________//
    // 9th, 10th, 11th tasks - KINK
    gROOT->LoadMacro("AddTaskKink.C");
    AliAnalysisKinkESDMC *taskkink = AddTaskKink();

    gROOT->LoadMacro("AddTaskKinkResonance.C");
    AliAnalysisTaskKinkResonance *taskkinkres = AddTaskKinkResonance();

    gROOT->LoadMacro("AddTaskKinkResonanceLikeSign.C");
    AliResonanceKinkLikeSign *taskkinklikesign = AddTaskKinkResonanceLikeSign();
  }

  if (runUnicor) {
    //____________________________________________//
    // 12th task - UNICOR
    gROOT->LoadMacro("AddTaskUnicor.C");
    AliAnalysisTaskUnicor *taskunicor = AddTaskUnicor();
  }

  if (runFMDanalysis) {
    //____________________________________________//
    // 13th task - FMD
    gROOT->LoadMacro("AddTaskFMD.C");
    AliFMDAnalysisTaskSE *taskfmd = AddTaskFMD();
  }

  //____________________________________________//
  // Running the train

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}

Int_t setupPar(const char* pararchivename) {
  ///////////////////
  // Setup PAR File//
  ///////////////////
  if (pararchivename) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
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
    
    gSystem->ChangeDirectory("../");
  }

  return 1;
}
