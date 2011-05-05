void
dndeta_final(Double_t max=6)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleFont(132, "xyz");
  gStyle->SetTitleSize(0.1, "xyz");
  gStyle->SetTitleOffset(0.4, "y");
  gStyle->SetTitleOffset(0.8, "x");
  gStyle->SetLabelFont(132, "xyz");
  gStyle->SetLabelSize(0.08, "xyz");
  gStyle->SetNdivisions(212, "x");
  gStyle->SetNdivisions(208, "y");
  gStyle->SetTextFont(132);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  // gStyle->SetFillColor(0);
  // gStyle->SetFillStyle(0);
  
  TCanvas* c = new TCanvas("c", "c", 900, 900);
  c->SetFillColor(0);
  c->SetFillStyle(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.02);
  c->SetBottomMargin(0.15);
  c->Divide(1,3,0,0);
  
  // --- INEL --------------------------------------------------------
  TVirtualPad* p = c->cd(1);
  p->SetGridx();
  p->SetRightMargin(.01);
  THStack* inel     = new THStack("inel", "INEL");
  TLatex*  inelT    = new TLatex(1-p->GetRightMargin()-.01, 
				 1-p->GetTopMargin()-.01, 
				 "INEL");
  inelT->SetNDC();
  inelT->SetTextAlign(33);
  inelT->SetTextSize(0.12);
  TLegend* inelL    = new TLegend(.3, .02, .8, .4);
  inelL->SetBorderSize(0);
  inelL->SetNColumns(2);
  inelL->SetFillColor(0);
  inelL->SetFillStyle(0);
  TLegendEntry* e = inelL->AddEntry("d1", "Forward", "lp");
  e->SetMarkerColor(kRed+2);
  e->SetMarkerStyle(29);
  e = inelL->AddEntry("d2", "Central", "lp");
  e->SetMarkerColor(kMagenta+2);    
  e->SetMarkerStyle(29);
  e = inelL->AddEntry("d3", "Data", "lp");
  e->SetMarkerStyle(29);
  e = inelL->AddEntry("d4", "Mirrored data", "lp");
  e->SetMarkerStyle(30);
  e = inelL->AddEntry("d5", "Systematic error", "f");
  e->SetFillColor(kGray); 
  e->SetLineColor(kGray);
  e->SetLineWidth(0);
  e->SetFillStyle(3001);
  
  gROOT->LoadMacro("export_pp_0900GeV_INEL_m10p10cm_000100000ev.C");
  export_pp_0900GeV_INEL_m10p10cm_000100000ev(inel, inelL, 20);
  export_pp_0900GeV_INEL_m10p10cm_000100000ev(inel, inelL, 21);
  export_pp_0900GeV_INEL_m10p10cm_000100000ev(inel, inelL, 22);
  inel->Draw("nostack e1");
  inel->GetHistogram()->SetYTitle("#frac{1}{N}#frac{dN_{ch}}{d#eta}");
  inel->GetHistogram()->SetXTitle("#eta");
  inel->GetHistogram()->GetYaxis()->SetDecimals();
  inelL->Draw();
  inelT->Draw();

  // --- INEL>0 ------------------------------------------------------
  p = c->cd(2);
  p->SetGridx();
  p->SetRightMargin(.01);
  THStack* inelgt0     = new THStack("inelgt0", "INEL>0");
  TLatex*  inelgt0T    = new TLatex(1-p->GetRightMargin()-.01, 
				    1-p->GetTopMargin()-.01, 
				    "INEL>0");
  inelgt0T->SetNDC();
  inelgt0T->SetTextAlign(33);
  inelgt0T->SetTextSize(0.12);
  gROOT->LoadMacro("export_pp_0900GeV_INEL_m10p10cm_000100000ev.C");
  export_pp_0900GeV_INEL_m10p10cm_000100000ev(inelgt0, 0, 20);
  export_pp_0900GeV_INEL_m10p10cm_000100000ev(inelgt0, 0, 21);
  export_pp_0900GeV_INEL_m10p10cm_000100000ev(inelgt0, 0, 22);
  inelgt0->Draw("nostack e1");
  inelgt0->GetHistogram()->SetXTitle("#eta");
  inelgt0->GetHistogram()->GetYaxis()->SetDecimals();
  inelgt0T->Draw();

  // --- NSD ---------------------------------------------------------
  p = c->cd(3);
  p->SetGridx();
  p->SetRightMargin(.01);
  THStack* nsd     = new THStack("nsd", "NSD");
  TLatex*  nsdT    = new TLatex(1-p->GetRightMargin()-.01, 
				1-p->GetTopMargin()-.01, 
				"NSD");
  nsdT->SetNDC();
  nsdT->SetTextAlign(33);
  nsdT->SetTextSize(0.12);
  gROOT->LoadMacro("export_pp_0900GeV_NSD_m10p10cm_000100000ev.C");
  export_pp_0900GeV_NSD_m10p10cm_000100000ev(nsd, 0, 20);
  export_pp_0900GeV_NSD_m10p10cm_000100000ev(nsd, 0, 21);
  export_pp_0900GeV_NSD_m10p10cm_000100000ev(nsd, 0, 22);
  nsd->Draw("nostack e1");
  nsd->GetHistogram()->SetXTitle("#eta");
  nsd->GetHistogram()->GetYaxis()->SetDecimals();
  nsdT->Draw();

  c->cd();
  c->SaveAs("dndeta_final.png");
}

  
