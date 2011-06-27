{
  gSystem->Load("libGui");
  gROOT->LoadMacro("$ALICE_ROOT.trunk/FMD/scripts/Compile.C");
  Compile("$ALICE_ROOT.trunk/FMD/scripts/DrawCalib.C","g");
  DrawCalib();
}
