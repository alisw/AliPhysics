
int runPhos(const char *data="./", const char *grp="./")
{
  gSystem->Load("libAliHLTCalo");
  gROOT->LoadMacro("rec_hlt_calo_phos.C");
  gROOT->LoadMacro("read_HLT_ESDs.C");
  rec_hlt_calo_phos();
  read_HLT_ESDs();
}
