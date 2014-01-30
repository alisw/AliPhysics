void runcheck(){
  gROOT->Macro("loadlibsrec.C");
  gROOT->Macro("${ALICE_ROOT}/STEER/CheckESD.C");
}
