// -*- C++ -*-

struct ScanData {
  Int_t fillNumber;
  Int_t minNumberOfTracks;
  const char* vtxFileName;
  const char* sepFileName;
  const char* scanName;

  const char* scanName1;
  Int_t scanType1;
  Double_t t1,t2;
  Double_t offset1;

  const char* scanName2;
  Int_t scanType2;
  Double_t t3,t4;
  Double_t offset2;
} ;

void CheckCopyFile(TString fn) {
  TFile *f = TFile::Open(fn);
  if (f) {
    f->Close();
    return;
  }
  if (!gGrid)
    TGrid::Connect("alien://");
  gSystem->Unlink(fn);
  TObjArray *a = fn.Tokenize("/");
  TString url = TString::Format("alien:///alice/cern.ch/user/d/dcaffarr/%s", a->At(a->GetEntries()-1)->GetName());
  delete a;
  Printf("url=%s", url.Data());
  TFile::Cp(url, TString("file://")+fn);
}
