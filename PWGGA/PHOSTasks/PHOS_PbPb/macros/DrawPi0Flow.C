DrawPi0Flow()
{
  gStyle->SetOptStat(0);
  TFile * f = new TFile("Pi0Flow_000167920.root") ;
  TList *histoList = (TList*)f->Get("PHOSPi0Flow");

  TH1F *hev = (TH1F*)histoList->FindObject("hTotSelEvents") ;
  TH2F *hCenPHOS  = (TH2F*)histoList->FindObject("hCenPHOS") ;
  TH2F *hPi0All_cen0 = (TH2F*)histoList->FindObject("hPi0All_cen0");
  TH2F *hPi0Allcore_cen0 = (TH2F*)histoList->FindObject("hPi0Allcore_cen0");

  TCanvas *c1 = new TCanvas("c1","Event selection");
  hev->SetXTitle("Event selection");
  hev->SetYTitle("N_{events}");
  hev->Draw();
  c1->Print("PHOS_EvSel.eps");

  TCanvas *c2 = new TCanvas("c2","PHOS multiplicity vs centrality");
  hCenPHOS->SetXTitle("centrality (%)");
  hCenPHOS->SetYTitle("Number of PHOS clusters");
  hCenPHOS->Draw("colz");
  c2->Print("PHOS_MultCentrality.eps");

  TCanvas *c3 = new TCanvas("c3","gg mass vs pt, PID=All");
  c3->SetLogz();
  hPi0All_cen0->SetTitle("PID: All");
  hPi0All_cen0->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
  hPi0All_cen0->SetYTitle("p_{T} (GeV/c)");
  hPi0All_cen0->Draw("colz");
  c3->Print("PHOS_MggAll.eps");

  TCanvas *c4 = new TCanvas("c4","gg mass vs pt, PID=Allcore");
  c4->SetLogz();
  hPi0Allcore_cen0->SetTitle("PID: Allcore");
  hPi0Allcore_cen0->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})"); 
  hPi0Allcore_cen0->SetYTitle("p_{T} (GeV/c)");
  hPi0Allcore_cen0->Draw("colz");
  c4->Print("PHOS_MggAllcore.eps");

  //  hPi0Allcore_cen0->ProjectionY()->Draw();

}
