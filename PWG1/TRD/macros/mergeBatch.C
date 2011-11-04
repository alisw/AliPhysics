//______________________________________________________
void mergeBatch(const Char_t *mark, const Char_t *files, const Int_t nfiles=20, const Int_t first=0, Bool_t kSVN=kTRUE, Bool_t kCLEAR=kFALSE)
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libPWG1.so");

  Int_t ntry(0);
  while(AliTRDpwg1Helper::MergeBatch(mark, files, nfiles, first, kSVN, kCLEAR) && ntry<5) ntry++;
}
