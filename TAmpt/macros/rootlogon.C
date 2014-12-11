{
  gSystem->Load("libTree");
  gSystem->Load("libGui");
  gSystem->Load("libVMC");
  gSystem->Load("libMinuit");
  gSystem->Load("libMinuit2");
  gSystem->Load("libProof");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libSTEER");
  gSystem->Load("libEVGEN");
  gSystem->Load("libFASTSIM");
  gSystem->Load("libampt");
  gSystem->Load("libTAmpt");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/PYTHIA6");
  if (gSystem->Getenv("TMPDIR")) 
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  if (!gAlice)
    new AliRun("gAlice","The ALICE Off-line Simulation Framework");
  if (1) {
    delete gRandom;
    gRandom = new TRandom3(0);
  }
}

