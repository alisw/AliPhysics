void runAnalysis() {
  TStopwatch timer;
  timer.Start();

  //____________________________________________________//
  //_____________Setting up ESD.par_____________________//
  //____________________________________________________//
  const char* pararchivename1 = "ESD";
  //////////////////////////////////////////
  // Libraries required to load
  //////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // Setup PAR File
  if (pararchivename1) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename1);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchivename1);

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
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");

  //_____________________________________________________________//
  //_____________Setting up ANALYSIS.par_____________________//
  //_____________________________________________________________//
  const char* pararchivename2 = "ANALYSIS";
  //////////////////////////////////////////
  // Libraries required to load
  //////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // Setup PAR File
  if (pararchivename2) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename2);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchivename2);

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
  gSystem->Load("libANALYSIS.so");

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://"); 

  gROOT->LoadMacro("AliAnalysisTaskPt.cxx+");  
  gROOT->LoadMacro("demoInteractive.C");
  demoInteractive();

  timer.Stop();
  timer.Print();
}
