{
  //=========Macro generated from canvas: c3/cr
  //=========  (Sat Oct 24 14:06:41 2009) by ROOT version5.23/01
  TCanvas *c3 = new TCanvas("c3", "cr",0,0,1000,600);
  gStyle->SetOptStat(0);
  c3->Range(0,0,1,1);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetBorderSize(2);
  c3->SetFrameBorderMode(0);

  TGraphBentErrors *grbe = new TGraphBentErrors(4);
  grbe->SetName("Graph");
  grbe->SetTitle("Graph");
  grbe->SetFillColor(1);
  grbe->SetMarkerColor(2);
  grbe->SetMarkerStyle(21);
  grbe->SetPoint(0,10,1163);
  grbe->SetPointError(0,0,0,319.8685,242.6435,0,0,0,0);
  grbe->SetPoint(1,20,1473);
  grbe->SetPointError(1,0,0,355.2512,271.8543,0,0,0,0);
  grbe->SetPoint(2,40,2172);
  grbe->SetPointError(2,0,0,490.1327,376.0158,0,0,0,0);
  grbe->SetPoint(3,80,3329);
  grbe->SetPointError(3,0,0,800.9482,426.2069,0,0,0,0);

  TH1 *Graph1 = new TH1F("Graph1","",100,3,87);
  Graph1->SetMinimum(10);
  Graph1->SetMaximum(5000);
  Graph1->SetDirectory(0);
  Graph1->SetStats(0);
  Graph1->GetXaxis()->SetTitle("p_{T}^{track} (GeV/c)");
  //  Graph1->GetXaxis()->SetTitleOffset(1.5);
  Graph1->GetYaxis()->SetTitle("Rejection Power");
  //Graph1->GetYaxis()->SetTitleOffset(1.5);
  grbe->SetHistogram(Graph1);

  grbe->Draw("apl");

  grbe = new TGraphBentErrors(4);
  grbe->SetName("Grapha");
  grbe->SetTitle("Grapha");
  grbe->SetFillColor(1);
  grbe->SetMarkerColor(4);
  grbe->SetMarkerStyle(21);
  grbe->SetPoint(0,10,214);
  grbe->SetPointError(0,0,0,23.04028,20.36199,0,0,0,0);
  grbe->SetPoint(1,20,600);
  grbe->SetPointError(1,0,0,86.54205,72.02112,0,0,0,0);
  grbe->SetPoint(2,40,651);
  grbe->SetPointError(2,0,0,77.49317,70.38893,0,0,0,0);
  grbe->SetPoint(3,80,1705);
  grbe->SetPointError(3,0,0,249.8343,221.2924,0,0,0,0);

  TH1 *Graph2 = new TH1F("Graph2","Rejection Power",100,3,87);
  Graph2->SetMinimum(17.42646);
  Graph2->SetMaximum(2099.826);
  Graph2->SetDirectory(0);
  Graph2->SetStats(0);
  Graph2->GetXaxis()->SetTitle("p_{T}^{track} (GeV/c)");
  Graph2->GetYaxis()->SetTitle("Rejection Power");
  grbe->SetHistogram(Graph2);

  grbe->Draw("pl");

  grbe = new TGraphBentErrors(2);
  grbe->SetName("Graphb");
  grbe->SetTitle("Graphb");
  grbe->SetFillColor(1);
  grbe->SetMarkerColor(6);
  grbe->SetMarkerStyle(29);
  grbe->SetMarkerSize(2);
  grbe->SetPoint(0,40,686);
  grbe->SetPointError(0,0,0,0,0,0,0,0,0);
  grbe->SetPoint(1,80,1210);
  grbe->SetPointError(1,0,0,0,0,0,0,0,0);

  TH1 *Graph3 = new TH1F("Graph3","Rejection Power",100,36,84);
  Graph3->SetMinimum(633.6);
  Graph3->SetMaximum(1262.4);
  Graph3->SetDirectory(0);
  Graph3->SetStats(0);
  Graph3->GetXaxis()->SetTitle("p_{T}^{track} (GeV/c)");
  Graph3->GetYaxis()->SetTitle("Rejection Power");
  grbe->SetHistogram(Graph3);

  grbe->Draw("pl");

  TLegend *leg = new TLegend(0.15,0.7,0.7,0.9,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  TLegendEntry *entry=leg->AddEntry("Graph","80% e^{-} Efficiency - Simulation","p");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(21);
  entry->SetMarkerSize(1);
  entry=leg->AddEntry("Grapha","90% e^{-} Efficiency - Simulation","p");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(21);
  entry->SetMarkerSize(1);
  entry=leg->AddEntry("Graphb","90% e^{-} Efficiency - Test Beam","p");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(21);
  entry->SetMarkerSize(1);
  leg->Draw();
  //  c3->Modified();
  //c3->cd();
  //c3->SetSelected(c3);
}
