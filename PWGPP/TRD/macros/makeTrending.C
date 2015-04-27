void makeTrending(const Char_t *fl, Bool_t relative=kFALSE, const Char_t *db = "$ALICE_PHYSICS/PWGPP/TRD/data/TRD.Trend.root")
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGPP");
  gSystem->Load("libPWGmuon");

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

