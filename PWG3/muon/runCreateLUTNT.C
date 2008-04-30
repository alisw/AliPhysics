void runCreateLUTNT() {

  gSystem->Load("libVMC.so");
  gSystem->Load("libMinuit.so");

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
 
  gSystem->Load("libProofPlayer.so");

  //____________________________________________________//
  //_____________Setting up STEERBase.par_____________________//
  //____________________________________________________//
  setupPar("STEERBase");

  //____________________________________________________//
  //_____________Setting up ESD.par_____________________//
  //____________________________________________________//
  setupPar("ESD");

  //____________________________________________________//
  //_____________Setting up AOD.par_____________________//
  //____________________________________________________//
  setupPar("AOD");

  //_____________________________________________________________//
  //_____________Setting up ANALYSIS.par_________________________//
  //_____________________________________________________________//
  setupPar("ANALYSIS");
  
  gROOT->LoadMacro("AliAnalysisTaskLUT.cxx+");
  gROOT->LoadMacro("analysisTaskLUTNT.C");

  TStopwatch timer;
  timer.Start();

  analysisTaskLUTNT();

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
        Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
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
