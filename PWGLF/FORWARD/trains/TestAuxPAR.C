void
TestAuxPAR()
{
  gROOT->Macro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/scripts/LoadLibs.C");
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include");
  gSystem->AddIncludePath("-I${ALICE_PHYSICS}/include");

  gROOT->LoadMacro("Railway.C++");
  gROOT->LoadMacro("ParUtilities.C++");

  TList files;
  files.Add(new TObjString("LocalRailway.C"));
  files.Add(new TObjString("GridRailway.C"));
  files.Add(new TObjString("analysis2/trains/../ForwardAODConfig.C"));

  ParUtilities::MakeAuxFilePAR(files, "test", true);

  ParUtilities::Load("test");
}
