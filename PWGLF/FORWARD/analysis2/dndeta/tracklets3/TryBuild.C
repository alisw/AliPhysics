void
TryBuild()
{
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include -I${ALICE_PHYSICS}/include");
  gROOT->LoadMacro("AliAODTracklet.C+g");
  gROOT->LoadMacro("AliTrackletAODUtils.C+g");
  gROOT->LoadMacro("AliTrackletAODTask.C+g");
  gROOT->LoadMacro("AliTrackletAODdNdeta.C+g");
}
