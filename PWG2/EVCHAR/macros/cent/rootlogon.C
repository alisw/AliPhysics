// $Id$
{
  gSystem->Load("libTree.so");
  gSystem->Load("libGui.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so"); 
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCDB.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libMinuit2.so");
  gSystem->Load("libProof.so");
  gSystem->Load("libRAWDatabase.so");
  gSystem->Load("libSTEER.so");
  gSystem->Load("libEVGEN.so");
  gSystem->Load("libFASTSIM.so");
  if (1) {
    gSystem->Load("libhijing.so");
    gSystem->Load("libTHijing.so");
  } else {
    gSystem->Load("libampt.so");
    gSystem->Load("libTAmpt.so");
  }
  gSystem->Load("libEGPythia6.so");
  gSystem->Load("libpythia6.so");
  gSystem->Load("libAliPythia6.so");
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT");
  if (gSystem->Getenv("TMPDIR")) 
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  if (!gAlice)
    new AliRun("gAlice","The ALICE Off-line Simulation Framework");
  if (1) {
    delete gRandom;
    gRandom = new TRandom(10);
  }

  gROOT->LoadMacro("createHijingGlauberTestTree.C+g");
}
