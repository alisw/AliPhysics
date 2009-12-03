void LoadLibraries(Bool_t useParFiles=kFALSE) {

  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so"); 
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG3base.so");
  gSystem->Load("libPWG3vertexingHF.so");
  gSystem->Load("libPWG3muon.so");
 

  if(useParFiles) {
    setupPar("STEERBase");
    setupPar("ESD");
    setupPar("AOD");
    setupPar("ANALYSIS");
    setupPar("ANALYSISalice");
    setupPar("CORRFW");  
    setupPar("PWG3base");
    setupPar("PWG3vertexingHF");
    setupPar("PWG3muon");
  }

  return;
}
//------------------------------------------------------------------------
Int_t setupPar(const char* pararchivename) {
  ///////////////////
  // Setup PAR File//
  ///////////////////

  if (pararchivename) {
    char processline[1024];
    TString base = gSystem->BaseName(pararchivename);
    TString dir  = gSystem->DirName(pararchivename);
    TString ocwd = gSystem->WorkingDirectory();
    // Move to dir where the par files are and unpack 
    gSystem->ChangeDirectory(dir.Data());
    sprintf(processline,".! tar xvzf %s.par",base.Data());
    gROOT->ProcessLine(processline);
    // Move to par folder                           
    gSystem->ChangeDirectory(base.Data());

    // check for BUILD.sh and execute                
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");

      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
        Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
	gSystem->ChangeDirectory(ocwd.Data());
        return -1;
      }
    }
    // check for SETUP.C and execute                
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      // If dir not empty, set the full include path 
      if (dir.Length()) {
	sprintf(processline, ".include %s", pararchivename);
	gROOT->ProcessLine(processline);
      }
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
    gSystem->ChangeDirectory(ocwd.Data());
  }
  return 1;
}
