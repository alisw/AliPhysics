void TestLoadFastJet() {
  TString libs = "libCGAL libfastjet libsiscone libsiscone_spherical "
                 "libfastjetplugins libfastjettools libfastjetcontribfragile "
                 "libPWGJEEMCALJetTasks";
  TString tok;
  Ssiz_t from;
  while (libs.Tokenize(tok, from)) {
    if (!tok.BeginsWith("lib")) tok = "lib" + tok;
    cout << "Loading " << tok << endl;
    if (gSystem->Load(tok.Data())) gSystem->Exit(1);
  }
  cout << "All libraries loaded" << endl;
  gSystem->Exit(0);
}
