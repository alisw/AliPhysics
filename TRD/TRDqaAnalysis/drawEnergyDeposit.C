
void drawEnergyDeposit(const char *filename) {
  
  gROOT->SetStyle("Plain");
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadLeftMargin(0.07);
  gStyle->SetPadBottomMargin(0.07);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
  
  TFile *file = new TFile(filename, "READ");
  if (file->IsZombie()) return;
  
  TCanvas *c = new TCanvas();
  c->Divide(5, 4, 0.001, 0.001);
  
  for(int i=0; i<5; i++) {
    TH1 *hist = (TH1*)file->Get(Form("probPos%d", i));
    c->cd(i+1);
    gPad->SetLogy();
    hist->Draw();
  }
  
 for(int i=0; i<5; i++) {
    TH1 *hist = (TH1*)file->Get(Form("ptSigPos%d", i));
    c->cd(i+6);
    gPad->SetLogx();
    hist->Draw("col");
  }

 for(int i=0; i<5; i++) {
    TH1 *hist = (TH1*)file->Get(Form("probNeg%d", i));
    c->cd(i+11);
    gPad->SetLogy();
    hist->Draw();
  }

 for(int i=0; i<5; i++) {
    TH1 *hist = (TH1*)file->Get(Form("ptSigNeg%d", i));
    c->cd(i+16);
    gPad->SetLogx();
    hist->Draw("col");
  }
}
