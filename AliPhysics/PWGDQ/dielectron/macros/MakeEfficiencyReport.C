void SetupStyle();

void MakeEfficiencyReport(const char* outputFile="JpsiEffReport.pdf",
                          const char* histos="jpsi_HistosSE.root",
                          const char* cf="jpsi_CF.root")
{
  //
  // Make a pdf file with the efficiency report
  //

  SetupStyle();
  
  AliDielectronCFdraw d(cf);
  d.SetRangeUser("PairType",1,1);
  d.SetRangeUser("Y",-.89,.9,"0");
  
  
  TFile f("jpsi_HistosSE.root");
  
  AliDielectronHistos h;
  TIter nextHists((TList*)f.Get("Dielectron_Histos"));
  
  TPaveText pt(.02,.6,.98,.8);
  TText *t1=pt.AddText("");
  TText *t2=pt.AddText("");
  
  TCanvas *c1=new TCanvas;
  
  TPDF p(outputFile);

  //
  // Efficiency plots
  //
  t1->SetTitle("Efficiency plots");
  t2->SetTitle("Pair");
  pt.Draw();
  c1->Update();
  
//   d.DrawEfficiency("Pt","1;2;4;10;16;8",0);
  d.DrawEfficiency("Pt","4;16;22",0);
  c1->Update();
  

  //
  // Efficiency as a function of the impact parameter
  //
//   c1->Clear();
  Int_t refStep=16;
  const Int_t impSteps=5;
  Double_t impMax[impSteps]={.05,.04,.03,.02,.01};
  TString title=d.GetCFContainer()->GetStepTitle(refStep);
//   d.DrawEfficiency("Pt:Y","1;2;4;10;16;8",0);
  d.DrawEfficiency("Pt",Form("%d",refStep),0,"sameleg2");
//   c1->Update();
  
  for (Int_t i=0; i<impSteps;++i){
    d.SetRangeUser("Leg1_ImpactParXY",-1*impMax[i],impMax[i]);
    d.SetRangeUser("Leg2_ImpactParXY",-1*impMax[i],impMax[i]);
    d.GetCFContainer()->SetStepTitle(refStep, (title+Form(" |leg dXY|<%.2f",impMax[i])).Data());
    d.DrawEfficiency("Pt",Form("%d",refStep),0,"same+leg2");
  }

  d.UnsetRangeUser("Leg1_ImpactParXY");
  d.UnsetRangeUser("Leg2_ImpactParXY");
  d.GetCFContainer()->SetStepTitle(refStep, title.Data());
  ((TH1*)gPad->GetListOfPrimitives()->FindObject("eff_16/00"))->SetMaximum(1.);

  c1->Update();


  c1->Clear();
  d.SetRangeUser("Y",-.89,.9);
//   d.DrawEfficiency("Pt:Y","1;2;4;10;16;8",0,"colz2");
  d.DrawEfficiency("Pt:Y","4;16;22",0,"colz2");
  c1->Update();
  c1->Clear();

  //
  // Inv mass plots
  //
//   t1->SetTitle("Invariant Mass plots");
//   t2->SetTitle("");
//   pt.Draw();
//   c1->Update();
  
  
  //
  // Make QA info
  //
  
  t1->SetTitle("QA summary plots for");
  THashList *list=0x0;
  while ( (list=(THashList*)nextHists()) ){
    h.SetHistogramList(*list);
    t2->SetTitle(list->GetName());
    pt.Draw();
    c1->Update();
    h.Draw();
    c1->Clear();
  }
  p.Close();
  delete c1;
}

void SetupStyle()
{
  const Int_t NCont=255;
  
  TStyle *st = new TStyle("mystyle","mystyle");
  gROOT->GetStyle("Plain")->Copy((*st));
  st->SetTitleX(0.1);
  st->SetTitleW(0.8);
  st->SetTitleH(0.08);
  st->SetStatX(.9);
  st->SetStatY(.9);
  st->SetNumberContours(NCont);
  st->SetPalette(1,0);
  st->SetOptStat("erm");
  st->SetOptFit(0);
  st->SetGridColor(kGray+1);
  st->SetPadGridX(kTRUE);
  st->SetPadGridY(kTRUE);
  st->SetPadTickX(kTRUE);
  st->SetPadTickY(kTRUE);
  st->cd();
  
  const Int_t NRGBs = 5;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  
}