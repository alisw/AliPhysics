

{
  //=========Macro generated from canvas: c7/7 efficiency per pt
  //=========  (Sat Oct 24 13:36:23 2009) by ROOT version5.14/00
  TCanvas *c7 = new TCanvas("c7", "7 efficiency per pt",158,233,600,400);
  c7->Range(-5,-0.08875,55,0.89875);
  c7->SetFillColor(0);
  c7->SetBorderMode(0);
  c7->SetBorderSize(2);
  c7->SetFrameBorderMode(0);
  c7->SetFrameBorderMode(0);

  TGraphErrors *gre = new TGraphErrors(5);
  gre->SetName("Graph");
  gre->SetTitle("PID efficiency");
  gre->SetFillColor(1);
  gre->SetMarkerStyle(20);
  gre->SetPoint(0,5,0.140749);
  gre->SetPointError(0,0,0);
  gre->SetPoint(1,15,0.559721);
  gre->SetPointError(1,0,0);
  gre->SetPoint(2,25,0.465536);
  gre->SetPointError(2,0,0);
  gre->SetPoint(3,35,0.418012);
  gre->SetPointError(3,0,0);
  gre->SetPoint(4,45,0.341103);
  gre->SetPointError(4,0,0);

  TH1 *Graph1 = new TH1F("Graph1","PID efficiency",100,1,49);
  Graph1->SetMinimum(0.01);
  Graph1->SetMaximum(0.8);
  Graph1->SetDirectory(0);
  Graph1->SetStats(0);
  Graph1->GetXaxis()->SetTitle("pT (GeV/c)");
  Graph1->GetYaxis()->SetTitle("efficiency");
  gre->SetHistogram(Graph1);

  gre->Draw("apl");

  gre = new TGraphErrors(5);
  gre->SetName("Graph2");
  gre->SetTitle("Graph2");
  gre->SetFillColor(1);
  gre->SetMarkerColor(2);
  gre->SetMarkerStyle(20);
  gre->SetPoint(0,5,0.136597);
  gre->SetPointError(0,0,0);
  gre->SetPoint(1,15,0.525284);
  gre->SetPointError(1,0,0);
  gre->SetPoint(2,25,0.425577);
  gre->SetPointError(2,0,0);
  gre->SetPoint(3,35,0.397009);
  gre->SetPointError(3,0,0);
  gre->SetPoint(4,45,0.316323);
  gre->SetPointError(4,0,0);
  gre->Draw("pl");

  TLegend *leg = new TLegend(0.5,0.6,0.85,0.85,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetFillStyle(1001);
  TLegendEntry *entry=leg->AddEntry("Graph","Identified Electrons (EMCAL)","p");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1);
  entry->SetTextAlign(12);
  entry->SetTextColor(1);
  entry=leg->AddEntry("Graph2","Identified Electrons (EMCAL+TRD+TPC)","p");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(2);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1);
  entry->SetTextAlign(12);
  entry->SetTextColor(1);
  leg->Draw();
  TLatex *   tex = new TLatex(0.4,0.91,"LHC09b4, 5.5Tev, b-jet");
  tex->SetNDC();
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.4,0.85,"oct20");
  tex->SetNDC();
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.15,0.85,"version = 99");
  tex->SetNDC();
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.15,0.75,"nEvents=406580");
  tex->SetNDC();
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  TPaveText *pt = new TPaveText(0.01,0.940161,0.217919,0.995,"blNDC");
  pt->SetName("title");
  pt->SetBorderSize(1);
  text = pt->AddText("PID efficiency");
  pt->Draw();
  c7->Modified();
  c7->cd();
  c7->SetSelected(c7);
}
