void
DrawDeltas()
{
  TFile* file0   = TFile::Open("dt_middle_none/trdt.root");
  TFile* file1   = TFile::Open("mc_middle_PtPidStr/trmc.root");
  TFile* file2   = TFile::Open("mc_middle_none/trmc.root");
  TList* l0      = (TList*)file0->Get("clist");
  TH1*   data    = (TH1*)  l0->FindObject("b0_TrData_WDist");
  TList* l1      = (TList*)file1->Get("clist");
  TH1*   reweigh = (TH1*)  l1->FindObject("b0_TrData_WDist");
  TList* l2      = (TList*)file2->Get("clist");
  TH1*   none    = (TH1*)  l2->FindObject("b0_TrData_WDist");

  data->Scale(1./data->Integral());
  data->SetMarkerColor(kRed+2);
  data->SetLineColor(kRed+2);
  data->SetMarkerStyle(20);
  data->SetTitle("Real data");
  
  reweigh->Scale(1./reweigh->Integral());
  reweigh->SetMarkerColor(kBlue+2);
  reweigh->SetLineColor(kBlue+2);
  reweigh->SetMarkerStyle(25);
  reweigh->SetTitle("Simulated data, reweighed");

  none->Scale(1./none->Integral());
  none->SetLineColor(kGreen+2);
  none->SetMarkerColor(kGreen+2);
  none->SetMarkerStyle(24);
  none->SetTitle("Simulated data");
    
  THStack* stack = new THStack("stack","");
  stack->Add(data);
  stack->Add(none);
  stack->Add(reweigh);
    
  TH1* rnone     = none   ->Clone("rnone");    rnone   ->Divide(data);
  TH1* rreweigh  = reweigh->Clone("rreweigh"); rreweigh->Divide(data);
  rnone   ->SetTitle("Simulated/Real data");
  rreweigh->SetTitle("Simulated (reweighed)/Real data");
  THStack* ratios = new THStack("ratios","");
  ratios->Add(rnone);
  ratios->Add(rreweigh);

  TCanvas* c= new TCanvas("c","c",900,1000);
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.01);
  c->Divide(1,2,0,0);

  c->cd(1); stack ->Draw("nostack");
  c->cd(2); ratios->Draw("nostack");
  c->GetPad(1)->SetLogy();   
  // c->GetPad(1)->SetLogx();  
  // c->GetPad(2)->SetLogx();  
  c->GetPad(1)->SetRightMargin(0.01);
  c->GetPad(2)->SetRightMargin(0.01);

  stack->GetHistogram()->SetXTitle("\\Delta");
  stack->GetHistogram()->SetYTitle("P(\\Delta) \\hbox{ Norm. to }"
				   "\\int \\mathrm{d}\\Delta");

  ratios->GetHistogram()->SetXTitle("\\Delta");
  ratios->GetHistogram()->SetYTitle("Ratio to real data");
  
  TLegend* l = c->GetPad(1)->BuildLegend(.6,.6,.98,.98);
  l->SetFillStyle(0);
  l->SetBorderSize(0);

  l = c->GetPad(2)->BuildLegend(.6,.12,.98,.6);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  TLine* ll = new TLine(data->GetXaxis()->GetXmin(), 1,
			data->GetXaxis()->GetXmax(), 1);
  ll->SetLineColor(data->GetLineColor());
  ll->SetLineStyle(7);
  ll->SetLineWidth(2);
  ll->Draw();

  c->SaveAs("deltas.png");
}
