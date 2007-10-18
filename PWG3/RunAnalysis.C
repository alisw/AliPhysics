//
// Example to generate the PWG3base library from a par file
// par file is created in ALICE_ROOT directory via make PGW3base.par
// To run the generation of AOD you also need the PWG0base.par file
// Gines Martinez, Nantes oct 2007
//
// Copy the par file to your working directory and execute root
// .L $ALICE_ROOT/PWG3/RunAnalysis.C
// setupPar("PWG3base")
// setupPar("PWG0base") 
// Now you can run you analysis macro
// .L RunAODGeneration.C
// 
   
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
