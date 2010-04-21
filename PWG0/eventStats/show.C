void show(const char* fileName = "event_stat.root")
{
  if (TString(fileName).BeginsWith("alien://"))
    TGrid::Connect("alien://");

  TFile::Open(fileName);
  
  c = new TCanvas("c", "c", 1200, 600);
  hist = (TH1*) gFile->Get("physics_selection/fHistStatistics");
  hist->SetStats(0);
  hist->Draw("TEXT");
  hist->SetTitle("");
  hist->SetMarkerSize(1.5);
  c->SetLeftMargin(0.20);
  c->SetBottomMargin(0.2);
  c->SaveAs("stat.png");
  
  c = new TCanvas("c2", "c2", 800, 600);
  c->SetLeftMargin(0.25);
  hist = (TH1*) gFile->Get("physics_selection/fHistBunchCrossing");
  hist->SetStats(0);
  hist->SetTitle("");
  hist->Draw("COL");
  c->SaveAs("bc.png");

  // BG
  gStyle->SetStatX(0.87);
  gStyle->SetStatY(0.93);

  TH2F * hBB  = (TH2F*) gDirectory->Get("physics_selection/background_identification/hCvsTCuts_CINT1B-ABCE-NOPF-ALL");
  TH2F * hBE  = (TH2F*) gDirectory->Get("physics_selection/background_identification/hCvsTCuts_CINT1C-ABCE-NOPF-ALL");
  TH2F * hBEA = (TH2F*) gDirectory->Get("physics_selection/background_identification/hCvsTCuts_CINT1A-ABCE-NOPF-ALL");

  if(!hBB || !hBE || !hBEA) {
    printf("WARNING: no BG histos\n");
  }
  else
  {
    hBE->Add(hBEA);
    hBE->SetTitle("CINT1A-ABCE-NOPF-ALL + CINT1C-ABCE-NOPF-ALL");
    c = new TCanvas("cbg","cbg",800,400);
    c->Divide(2,1);
    c->cd(1);
    hBB->Draw("colz");
    c->cd(2);
    hBE->Draw("colz");
    c->Print("bg.png");
  }
}
