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

  TH1 *h1(NULL); TObjArray *trend = new TObjArray(100);
  h1 = tm->MakeTrends(fl, trend);

  TFile *fOut = TFile::Open("Trend_TRDgraph.root", "RECREATE");
  h1->SetDirectory(fOut); h1->Write();
  trend->Write();
  fOut->Close();
  return;
}

