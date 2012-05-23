void makeTrending(const Char_t *fl)
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGPP.so");
  gSystem->Load("libPWGmuon.so");

  AliTRDtrendingManager *tm = AliTRDtrendingManager::Instance();
  tm->Load();
  tm->MakeTrends(fl);
  return /*task*/;
}

