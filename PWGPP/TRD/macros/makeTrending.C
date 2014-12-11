void makeTrending(const Char_t *fl, Bool_t relative=kFALSE, const Char_t *db = "$ALICE_ROOT/PWGPP/TRD/data/TRD.Trend.root")
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTender.so");
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

