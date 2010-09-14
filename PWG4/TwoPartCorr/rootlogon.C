{
  gSystem->Load("libTree.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so"); 
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  if (gSystem->Getenv("ALICE_ROOT_ORIG")) 
    gSystem->AddIncludePath("-I$ALICE_ROOT_ORIG/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  if (gSystem->Getenv("TMPDIR")) 
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  if (1) {
      delete gRandom;
      gRandom = new TRandom3(0);
   }
  gROOT->LoadMacro("MyTreeClasses.C+g");
}



