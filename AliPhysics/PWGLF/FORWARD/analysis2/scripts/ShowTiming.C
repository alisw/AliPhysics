void showParts(TProfile* p, const char* title)
{
  if (!p) return;
  Printf(title);
  
  Int_t    n   = p->GetNbinsX();
  Double_t sum = 0;
  for (Int_t i = 1; i < n; i++) sum += p->GetBinContent(i);

  Double_t total = p->GetBinContent(n);
  Printf("sum %f total %f diff %f", sum, total, total-sum);

  for (Int_t i = 1; i < n; i++) { 
    Double_t cc = p->GetBinContent(i); 
    Printf("  //  %-24s: fraction of sum %4.1f%%   of total %4.1f%%", 
	   p->GetXaxis()->GetBinLabel(i), 100*cc/sum, 100*cc/total); 
  }
}
  
void ShowTiming(const char* filename)
{
  TFile *_file0 = TFile::Open(filename);
  new TBrowser;
  TList* forward;
  gDirectory->GetObject("Forward",forward);
  TList* forwardResults;
  gDirectory->GetObject("ForwardResults",forwardResults);
  TProfile* t   = (TProfile*)forward->FindObject("timing");
  TList*    dc  = (TList*)forward->FindObject("fmdDensityCalculator");
  TProfile* dct = (TProfile*)dc->FindObject("timing");
  TProfile* tr  = (TProfile*)forwardResults->FindObject("timing");
  
  showParts(t, "Full task");
  showParts(dct, "Density calculator");

  if (t) {
    TCanvas* c1 = new TCanvas("c1", "c1");
    t->Draw("s text30");
  }

  if (dct) {
    TCanvas* c2 = new TCanvas("c2", "c2");
    dct->Draw("s text30");
  }

  if (tr) {
    TCanvas* c3 = new TCanvas("c3", "c3");
    tr->Draw("s text30");
  }

}
