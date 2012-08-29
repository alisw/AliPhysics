DrawPi0Flow(const TString filename = "Pi0Flow_000167920.root")
{
  gStyle->SetOptStat(0);
  TFile * f = new TFile(filename) ;
  TList *histoList = (TList*)f->Get("PHOSPi0Flow");

  TH1F *hev        = (TH1F*)histoList->FindObject("hTotSelEvents") ;
  TH2F* hZvertex   = (TH2F*)histoList->FindObject("hZvertex");
  TH2F* hCenTrack  = (TH2F*)histoList->FindObject("hCenTrack");
  TH2F *hCenPHOS   = (TH2F*)histoList->FindObject("hCenPHOS") ;
  TH2F* phiRPV0A   = (TH2F*)histoList->FindObject("phiRPV0A");
  TH2F* phiRPV0C   = (TH2F*)histoList->FindObject("phiRPV0C");
  TH2F* hCellNXZM1 = (TH2F*)histoList->FindObject("hCellNXZM1");
  TH2F* hCellNXZM2 = (TH2F*)histoList->FindObject("hCellNXZM2");
  TH2F* hCellNXZM3 = (TH2F*)histoList->FindObject("hCellNXZM3");
  TH2F *hPi0All_cen0     = (TH2F*)histoList->FindObject("hPi0All_cen0");
  TH2F *hPi0Allcore_cen0 = (TH2F*)histoList->FindObject("hPi0Allcore_cen0");

  //-----------------------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1","Event selection");
  hev->GetXaxis()->SetBinLabel(1,"Total");
  hev->GetXaxis()->SetBinLabel(2,"Total");
  hev->GetXaxis()->SetBinLabel(3,"Goog vertex");
  hev->GetXaxis()->SetBinLabel(4,"No pileup");
  hev->GetXaxis()->SetBinLabel(5,"Good centrality");
  hev->GetXaxis()->SetBinLabel(6,"Good RP");
  hev->GetXaxis()->SetBinLabel(7,"MC");
  hev->GetXaxis()->SetBinLabel(8,"QA filled");
  hev->GetXaxis()->SetBinLabel(9,"Clusters filled");
  hev->GetXaxis()->SetBinLabel(10,"Pi0 filled");
  hev->GetXaxis()->SetBinLabel(11,"Mixed filled");
  hev->GetXaxis()->SetBinLabel(12,"Lists updated");
  hev->SetYTitle("N_{events}");
  hev->Draw();
  c1->Print("PHOS_EvSel.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c2 = new TCanvas("c2","Track multiplicity vs centrality");
  hCenTrack->SetXTitle("centrality (%)");
  hCenTrack->SetYTitle("Number of tracks");
  hCenTrack->Draw("colz");
  c2->Print("TrackMultCentrality.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c3 = new TCanvas("c3","PHOS multiplicity vs centrality");
  hCenPHOS->SetXTitle("centrality (%)");
  hCenPHOS->SetYTitle("Number of PHOS clusters");
  hCenPHOS->Draw("colz");
  c3->Print("PHOSMultCentrality.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c4 = new TCanvas("c4","RP phi");
  phiRPV0A1 = phiRPV0A->ProjectionX();
  phiRPV0C1 = phiRPV0C->ProjectionX();
  phiRPV0A1->SetLineColor(kRed);
  phiRPV0C1->SetLineColor(kBlue);
  phiRPV0A1->SetXTitle("RP #phi");
  phiRPV0C1->SetXTitle("RP #phi");
  phiRPV0A1->Draw();
  phiRPV0C1->Draw("same");
  leg = new TLegend(0.7,0.8,0.89,0.89);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->AddEntry(phiRPV0A1,"V0A","l");
  leg->AddEntry(phiRPV0C1,"V0C","l");
  leg->Draw();
  c4->Print("V0RPphi.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c5 = new TCanvas("c5","XZ in M1");
  c5->SetLogz();
  hCellNXZM1->SetTitle("Cell occupancy in module 4");
  hCellNXZM1->SetXTitle("X (cells)");
  hCellNXZM1->SetYTitle("X (cells)");
  hCellNXZM1->Draw("colz");
  c5->Print("PHOS_XYM4.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c6 = new TCanvas("c6","XZ in M2");
  c6->SetLogz();
  hCellNXZM2->SetTitle("Cell occupancy in module 3");
  hCellNXZM2->SetXTitle("X (cells)");
  hCellNXZM2->SetYTitle("X (cells)");
  hCellNXZM2->Draw("colz");
  c6->Print("PHOS_XYM3.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c7 = new TCanvas("c7","XZ in M3");
  c7->SetLogz();
  hCellNXZM3->SetTitle("Cell occupancy in module 2");
  hCellNXZM3->SetXTitle("X (cells)");
  hCellNXZM3->SetYTitle("X (cells)");
  hCellNXZM3->Draw("colz");
  c7->Print("PHOS_XYM2.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c8 = new TCanvas("c8","gg mass vs pt, PID=All");
  c8->SetLogz();
  hPi0All_cen0->SetTitle("PID: All");
  hPi0All_cen0->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
  hPi0All_cen0->SetYTitle("p_{T} (GeV/c)");
  hPi0All_cen0->Draw("colz");
  c3->Print("PHOS_MggAll.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c9 = new TCanvas("c9","gg mass vs pt, PID=Allcore");
  c9->SetLogz();
  hPi0Allcore_cen0->SetTitle("PID: Allcore");
  hPi0Allcore_cen0->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})"); 
  hPi0Allcore_cen0->SetYTitle("p_{T} (GeV/c)");
  hPi0Allcore_cen0->Draw("colz");
  c9->Print("PHOS_MggAllcore.eps");

}
