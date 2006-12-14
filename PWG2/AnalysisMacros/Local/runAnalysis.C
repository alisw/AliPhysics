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
  //_____________Setting up ANALYSIS_NEW.par_____________________//
  //_____________________________________________________________//
  const char* pararchivename2 = "ANALYSIS_NEW";
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
  gSystem->Load("libANALYSIS_NEW.so");

  gROOT->LoadMacro("AliAnalysisTaskPt.cxx+");  
  gROOT->LoadMacro("demoLocal.C");
  demoLocal();

  timer.Stop();
  timer.Print();
}
