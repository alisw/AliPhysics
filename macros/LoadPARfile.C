void LoadPARfile(TString parfile, Bool_t clean=kFALSE) {
  parfile.ReplaceAll(".par", "");
  const char *pardir = "unpacked_pars";
  gSystem->mkdir(pardir);
  TString buf;
  if (clean) {
    buf.Form("rm -rf %s/%s", pardir, parfile.Data());
    gSystem->Exec(buf.Data());
  }
  buf.Form("tar -C %s -xf %s.par",
           pardir, parfile.Data(), pardir, parfile.Data());
  gSystem->Exec(buf.Data());
  buf.Form("%s/%s", pardir, parfile.Data());
  TString owd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(buf.Data());
  buf.Form("test -x PROOF-INF/BUILD.sh && PROOF-INF/BUILD.sh",
           pardir, parfile.Data());
  gSystem->Exec(buf.Data());
  buf.Form("PROOF-INF/SETUP.C", pardir, parfile.Data());
  if (!gSystem->AccessPathName(buf.Data())) gROOT->Macro(buf.Data());
  gSystem->ChangeDirectory(owd.Data());
}
