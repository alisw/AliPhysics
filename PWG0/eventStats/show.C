void show(const char* fileName = "event_stat.root")
{
  if (TString(fileName).BeginsWith("alien://"))
    TGrid::Connect("alien://");

  TFile::Open(fileName);
  
  c = new TCanvas;
  hist = (TH1*) gFile->Get("physics_selection/fHistStatistics");
  hist->SetStats(0);
  hist->Draw("TEXT");
  c->SetLeftMargin(0.25);
  c->SetBottomMargin(0.2);
  c->SaveAs("stat.png");
  
  c = new TCanvas;
  c->SetLeftMargin(0.25);
  hist = (TH1*) gFile->Get("physics_selection/fHistBunchCrossing");
  hist->SetStats(0);
  hist->Draw("TEXT");
  c->SaveAs("bc.png");
}
