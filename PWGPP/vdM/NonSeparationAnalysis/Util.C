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

TH1* DrawFrame(TVirtualPad* c,
	       Double_t yMin, Double_t yMax,
	       const char* title,
	       const char* yTitle,
	       Bool_t logY=kFALSE,
	       Double_t xMin=-0.65,
	       Double_t xMax= 0.65,
	       TString label="",
	       const char* xTitle="Separation [mm]",
	       Double_t bottomMargin=0.15) {
  c->SetLeftMargin(0.17);
  c->SetBottomMargin(bottomMargin);
  TH1* hf = c->DrawFrame(xMin, yMin, xMax, yMax);
  hf->SetTitle(title);
  hf->GetXaxis()->SetTitle(xTitle);
  hf->GetYaxis()->SetTitle(yTitle);
  hf->SetLabelSize(0.06, "X");
  hf->SetLabelSize(0.06, "Y");
  hf->SetTitleSize(0.08, "XY");
  hf->SetTitleOffset(0.9, "X");
  hf->SetTitleOffset(1.1, "Y");
  if (logY)
    c->SetLogy();

  if (label != "") {
    TLatex *tl = new TLatex;
    tl->SetNDC();
    tl->SetTextSize(0.05);
    tl->SetTextAlign(23);
    tl->DrawLatex(0.55, 0.88, label);
  }
  return hf;
}
