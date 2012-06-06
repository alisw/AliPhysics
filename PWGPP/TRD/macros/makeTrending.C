void makeTrending(const Char_t *fl, const Char_t *db = "$ALICE_ROOT/PWGPP/TRD/data/TRD.Trend.root", Bool_t relative=kFALSE)
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGPP.so");
  gSystem->Load("libPWGmuon.so");

  AliTRDtrendingManager *tm = AliTRDtrendingManager::Instance();
  tm->Load(db);
  tm->SetRelativeMeanSigma(relative);
  tm->MakeTrends(fl);
  return;
}

