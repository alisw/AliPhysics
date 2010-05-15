
int runAll(const char *data="./", const char *grp="./", bool runPhos = true, bool runEmcal = true, bool runTM = true)
{
  gSystem->Load("libAliHLTCalo");
  gROOT->LoadMacro("rec_hlt_calo.C");
  gROOT->LoadMacro("read_HLT_ESDs.C");
  rec_hlt_calo("./", "./", runPhos, runEmcal, runTM);
  read_HLT_ESDs();
}
