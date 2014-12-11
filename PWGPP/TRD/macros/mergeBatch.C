//______________________________________________________
void mergeBatch(const Char_t *mark, const Char_t *files, const Int_t nfiles=20, const Int_t first=0, Bool_t kSVN=kTRUE, Bool_t kCLEAR=kFALSE)
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");
  gSystem->Load("libPWGmuon");

  Int_t ntry(0);
  while(AliTRDpwgppHelper::MergeBatch(mark, files, nfiles, first, kSVN, kCLEAR) && ntry<5) ntry++;
}
