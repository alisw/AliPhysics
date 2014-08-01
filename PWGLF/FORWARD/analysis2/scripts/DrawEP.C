void
DrawEP()
{
  TFile *file = TFile::Open("AnalysisResults.root");
  TList* l    = static_cast<TList*>(gDirectory->Get("Forward"));
  TList* ep   = static_cast<TList*>(l->FindObject("fmdEventPlaneFinder"));
  ep->ls();

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  TCanvas* c = new TCanvas("epFMD", "From FMD", 1000,1000);
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.01);
  c->Divide(1,2);
    
  TVirtualPad* p = c->cd(1);
  p->Divide(2,1);

  TVirtualPad* q = p->cd(1);
  q->SetTopMargin(0.01);
  q->SetRightMargin(0.01);
  
  THStack* epS    = new THStack("psiR", "From FMD");
  TH1*     epAll = 0;
  TH1*     epA   = 0;
  TH1*     epC   = 0;
  epS->Add(epAll = static_cast<TH1*>(ep->FindObject("epFMD")));
  epS->Add(epA   = static_cast<TH1*>(ep->FindObject("epFMDA")));
  epS->Add(epC   = static_cast<TH1*>(ep->FindObject("epFMDC")));
  epAll->SetMarkerStyle(20);
  epAll->SetMarkerColor(epAll->GetLineColor());
  epA->SetMarkerStyle(21);
  epA->SetMarkerColor(epA->GetLineColor());
  epC->SetMarkerStyle(22);
  epC->SetMarkerColor(epC->GetLineColor());
  epS->Draw("nostack");
  epS->GetHistogram()->SetXTitle("#Psi_{R} [radians]");
  
  TLegend* leg = new TLegend(0.11, 0.11, .98, .4);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetNColumns(3);
  leg->AddEntry(epAll, "Full FMD",        "p");
  leg->AddEntry(epA,   "A-side (#eta>0)", "p");
  leg->AddEntry(epC,   "C-side (#eta<0)", "p");
  leg->Draw();

  q = p->cd(2);
  q->SetTopMargin(0.01);
  q->SetRightMargin(0.01);
  q->SetLogy();

  THStack* diff = new THStack("diff", "Different to others");
  TH1* dSelf  = 0;
  TH1* dTPC   = 0;
  TH1* dVZERO = 0;
  diff->Add(dSelf  = static_cast<TH1*>(ep->FindObject("diffFMDAC")));
  diff->Add(dVZERO = static_cast<TH1*>(ep->FindObject("diffFMDVZERO")));
  diff->Add(dTPC   = static_cast<TH1*>(ep->FindObject("diffFMDTPC")));
  dSelf->SetMarkerStyle(20);
  dSelf->SetMarkerColor(dSelf->GetLineColor());
  dTPC->SetMarkerStyle(21);
  dTPC->SetMarkerColor(dTPC->GetLineColor());
  dVZERO->SetMarkerStyle(22);
  dVZERO->SetMarkerColor(dVZERO->GetLineColor());
  diff->Draw("nostack hist");
  diff->GetHistogram()->SetXTitle("#delta#Psi_{R} [radians]");

  leg = new TLegend(0.6, 0.6, .98, .98);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  // leg->SetNColumns(3);
  leg->AddEntry(epAll, "#Psi_{R,A}-#Psi_{R,C}",       "p");
  leg->AddEntry(epA,   "#Psi_{R,FMD}-#Psi_{R,TPC}",   "p");
  leg->AddEntry(epC,   "#Psi_{R,FMD}-#Psi_{R,VZERO}", "p");
  leg->Draw();
  
  p = c->cd(2);
  p->Divide(3,1);
  
  q = p->cd(1); q->SetTopMargin(0); q->SetRightMargin(0.13);
  ep->FindObject("corrFMDAC")->Draw("colz");

  q = p->cd(2); q->SetTopMargin(0); q->SetRightMargin(0.13);
  ep->FindObject("corrFMDTPC")->Draw("colz");

  q = p->cd(3); q->SetTopMargin(0); q->SetRightMargin(0.13);
  ep->FindObject("corrFMDVZERO")->Draw("colz");
  
  c->Print("ep_from_fmd.png");
}

