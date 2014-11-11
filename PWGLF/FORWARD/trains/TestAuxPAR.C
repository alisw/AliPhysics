void
TestAuxPAR()
{
  gROOT->Macro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/scripts/LoadLibs.C");
  gSystem->Load("libProof");
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include");

  gROOT->LoadMacro("Railway.C++");
  gROOT->LoadMacro("ParUtilities.C++");

  TList files;
  files.Add(new TObjString("AAFRailway.C"));
  files.Add(new TObjString("GridRailway.C"));
  files.Add(new TObjString("analysis2/trains/../ForwardAODConfig.C"));

  ParUtilities::MakeAuxFilePAR(files, "test", true);

  TProof::Open("lite://");
  gProof->UploadPackage("test.par");
  gProof->EnablePackage("test.par");
}

  
