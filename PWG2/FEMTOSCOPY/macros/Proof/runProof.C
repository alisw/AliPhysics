void runProof(const char *chainlistfile, int nEvents=0, int offset=0, int domc=0, int prooflite=0) {
  TStopwatch timer;
  timer.Start();
  
  printf("*** Open PROOF ***");
  if (prooflite) 
    TProof::Open("");
  else
    TProof::Open("alicecaf");

  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");

  gProof->UploadPackage("STEERBase.par");
  gProof->EnablePackage("STEERBase");
  gProof->UploadPackage("ESD.par");
  gProof->EnablePackage("ESD");
  gProof->UploadPackage("AOD.par");
  gProof->EnablePackage("AOD");
  gProof->UploadPackage("ANALYSIS.par");
  gProof->EnablePackage("ANALYSIS");
  gProof->UploadPackage("ANALYSISalice.par");
  gProof->EnablePackage("ANALYSISalice");
  gProof->UploadPackage("PWG2AOD.par");
  gProof->EnablePackage("PWG2AOD");
  gProof->UploadPackage("PWG2femtoscopy.par");
  gProof->EnablePackage("PWG2femtoscopy");
  gProof->UploadPackage("PWG2femtoscopyUser.par");
  gProof->EnablePackage("PWG2femtoscopyUser");
    
  gSystem->SetIncludePath("-I$ROOTSYS/include -I./STEERBase/ -I./ESD/ -I./AOD/ -I./ANALYSIS/ -I./ANALYSISalice/ -I./PWG2AOD/AOD -I./PWG2femtoscopy/FEMTOSCOPY/AliFemto -I./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser");
  gProof->Exec(".L AddTaskFemto.C",kTRUE);
  gROOT->LoadMacro("AddTaskFemto.C");
  cout << "Loaded AddTaskFemto macro "<< endl;

  gProof->ShowEnabledPackages();

  //ANALYSIS PART
  TChain *chain = 0x0;
  if (prooflite) {
    gROOT->LoadMacro("CreateESDChain.C");
    chain = CreateESDChain(chainlistfile,10000);
  }

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliESDInputHandler* esdH = new AliESDInputHandler;

  esdH->SetInactiveBranches("FMD CaloCluster");
  mgr->SetInputEventHandler(esdH);  

  if (domc) {
    AliMCEventHandler *mcH = new AliMCEventHandler;
    mgr->SetMCtruthEventHandler(mcH);
  }
  //____________________________________________//
  // 1st Pt task

  AliAnalysisTaskFemto *taskfemto = AddTaskFemto();

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if (prooflite) 
    mgr->StartAnalysis("proof",chain,2000000,0);
  else 
    mgr->StartAnalysis("proof",chainlistfile,nEvents,offset);

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
