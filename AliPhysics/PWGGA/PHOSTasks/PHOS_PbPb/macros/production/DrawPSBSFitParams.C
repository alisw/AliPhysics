#include "TCanvas.h"



namespace RawProduction {
  class Output;
}


void DrawAbs(const RawProduction::Output& output, const char* trigger, int cent, TStringToken& names) 
{
  names.NextToken();

  const char* pid = "All";
  char canvName[256] = Form("PSBS_Abs_%s_c%03i_%s_%s", trigger, cent, pid, names.Data());
  TCanvas* canv = new TCanvas(canvName, canvName);  
  TLegend* leg = new TLegend(0.6,0.8,0.95,0.95);

  TH1* hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, pid, names.Data()));
  if(names.Contains("mr")) hist->SetTitle(Form("Peak Position, %s, %s, %s", trigger, RawProduction::GetCentString(cent), pid));
  if(names.Contains("sr")) hist->SetTitle(Form("Peak Width, %s, %s, %s", trigger, RawProduction::GetCentString(cent), pid));
  if(names.Contains("mr")) hist->GetYaxis()->SetRangeUser(0.12, 0.15);
  if(names.Contains("mr")) hist->GetYaxis()->SetTitle("Peak #mu");
  if(names.Contains("sr")) hist->GetYaxis()->SetRangeUser(0., 0.012);
  if(names.Contains("sr")) hist->GetYaxis()->SetTitle("Peak #sigma");
  hist->GetXaxis()->SetTitle("p_{T}");
  //Printf(hist->GetTitle());
  hist->SetMarkerStyle(21);
  //hist->SetMarkerSize(1.5);
  hist->SetMarkerColor(kBlack);
  hist->SetLineColor(kBlack);
  hist->Draw();
  leg->AddEntry(hist, "Pol1, Ratio", "lep");
  
  int marker = 21;
  Color_t color[3] = {kRed, kBlue, kMagenta};
  while( names.NextToken() ) {
    hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "All", names.Data()));
    //Printf(hist->GetName());
    hist->SetMarkerStyle(++marker);
    hist->SetMarkerColor(color[marker-22]);
    hist->SetLineColor(color[marker-22]);
    hist->Draw("same");
    char legName[256] = "";
    if(names.Contains("1")) sprintf(legName, "Pol1");
    if(names.Contains("2")) sprintf(legName, "Pol2");
    if( marker <23 ) sprintf(legName, "%s, Ratio", legName);
    leg->AddEntry(hist, legName, "lep");
  }
  
  hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, pid, names.Data()));
  hist->Draw("same");

  leg->Draw();
  
  canv->SaveAs(Form("imgs/%s.png", canvName));
  canv->SaveAs(Form("imgs/%s.pdf", canvName));
}

void DrawRatios(const RawProduction::Output& output, const char* trigger, int cent)
{
  TStringToken graphs("mr1;sr1;mr1r;sr1r;mr2;sr2;mr2r;sr2r", ";");
  while( graphs.NextToken() ) {
    const char* graph = graphs.Data();

    // ratio
    TCanvas* canv = new TCanvas(Form("%s_c%03i_%s_ratio", trigger, cent, graph), Form("%s_c%03i_%s_ratio", trigger, cent, graph));
    TLegend* leg = new TLegend(0.6,0.8,0.95,0.95);

    hAll = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "All", graph));


    TH1* hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "CPV", graph));
    hist = (TH1*) hist->Clone(Form("dAll_%s", hist->GetName()));
    if(graphs.Contains("mr1")) hist->SetTitle( "Peak Position, Pol1" );
    if(graphs.Contains("mr2")) hist->SetTitle( "Peak Position, Pol2" );
    if(graphs.Contains("mr1r")) hist->SetTitle( "Peak Position, Pol1, Ratio" );
    if(graphs.Contains("mr2r")) hist->SetTitle( "Peak Position, Pol2, Ratio" );
    if(graphs.Contains("sr1")) hist->SetTitle( "Peak Width, Pol1" );
    if(graphs.Contains("sr2")) hist->SetTitle( "Peak Width, Pol2" );
    if(graphs.Contains("sr1r")) hist->SetTitle( "Peak Width, Pol1, Ratio" );
    if(graphs.Contains("sr2r")) hist->SetTitle( "Peak Width, Pol2, Ratio" );
    hist->SetTitle(Form("%s, %s, centrality: %s", hist->GetTitle(), trigger, RawProduction::GetCentString(cent)));
    if(graphs.Contains("mr")) hist->GetYaxis()->SetRangeUser(0.9, 1.05);
    if(graphs.Contains("mr")) hist->GetYaxis()->SetTitle("Peak #mu");
    if(graphs.Contains("sr")) hist->GetYaxis()->SetRangeUser(0.4, 1.4);
    if(graphs.Contains("sr")) hist->GetYaxis()->SetTitle("Peak #sigma");
    hist->GetXaxis()->SetTitle("p_{T}");
    hist->Divide(hAll);
    hist->SetMarkerStyle(22);
    hist->SetMarkerColor(kRed);
    hist->SetLineColor(kRed);
    hist->Draw();
    leg->AddEntry(hist, "CPV/All", "lep");

    hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "Disp", graph));
    hist = (TH1*) hist->Clone(Form("dAll_%s", hist->GetName()));
    hist->Divide(hAll);
    hist->SetMarkerStyle(23);
    hist->SetMarkerColor(kBlue);
    hist->SetLineColor(kBlue);
    hist->Draw("same");
    leg->AddEntry(hist, "Disp/All", "lep");

    hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "Both", graph));
    hist = (TH1*) hist->Clone(Form("dAll_%s", hist->GetName()));
    hist->Divide(hAll);
    hist->SetMarkerStyle(21);
    hist->SetMarkerColor(kMagenta);
    hist->SetLineColor(kMagenta);
    hist->Draw("same");
    leg->AddEntry(hist, "Both/All", "lep");

    
    leg->Draw();

    canv->SaveAs(Form("imgs/PSBSFits_ratio_%s_c%03i_%s.png", trigger, cent, graph));
    canv->SaveAs(Form("imgs/PSBSFits_ratio_%s_c%03i_%s.pdf", trigger, cent, graph));
    



    // core, ratio
    canv = new TCanvas(Form("%s_c%03i_%s_ratio_core", trigger, cent, graph), Form("%s_c%03i_%s_ratio_core", trigger, cent, graph));
    leg = new TLegend(0.6,0.8,0.95,0.95);

    hAll = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "All", graph));

    hAllcore = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "Allcore", graph));

    hist = hAllcore;
    hist = (TH1*) hist->Clone(Form("dAll_%s", hist->GetName()));
    hist->Divide(hAll);
    if(graphs.Contains("r1")) hist->SetTitle( Form("%s, Pol1", hist->GetTitle()) );
    if(graphs.Contains("r2")) hist->SetTitle( Form("%s, Pol2", hist->GetTitle()) );
    hist->SetTitle(Form("%s, %s, centrality: %s", hist->GetTitle(), trigger, RawProduction::GetCentString(cent)));
    if(graphs.Contains("mr")) hist->GetYaxis()->SetRangeUser(0.9, 1.05);
    if(graphs.Contains("mr")) hist->GetYaxis()->SetTitle("Peak #mu");
    if(graphs.Contains("sr")) hist->GetYaxis()->SetRangeUser(0.4, 1.4);
    if(graphs.Contains("sr")) hist->GetYaxis()->SetTitle("Peak #sigma");
    hist->GetXaxis()->SetTitle("p_{T}");
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kBlack);
    hist->SetLineColor(kBlack);
    hist->Draw();
    leg->AddEntry(hist, "Allcore/All", "lep");

    hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "CPVcore", graph));
    hist = (TH1*) hist->Clone(Form("dAllcore_%s", hist->GetName()));
    hist->Divide(hAllcore);
    hist->SetMarkerStyle(22);
    hist->SetMarkerColor(kRed);
    hist->SetLineColor(kRed);
    hist->Draw("same");
    leg->AddEntry(hist, "CPVcore/Allcore", "lep");

    hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "Dispcore", graph));
    hist = (TH1*) hist->Clone(Form("dAllcore_%s", hist->GetName()));
    hist->Divide(hAllcore);
    hist->SetMarkerStyle(23);
    hist->SetMarkerColor(kBlue);
    hist->SetLineColor(kBlue);
    hist->Draw("same");
    leg->AddEntry(hist, "Dispcore/Allcore", "lep");

    hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "Bothcore", graph));
    hist = (TH1*) hist->Clone(Form("dAllcore_%s", hist->GetName()));
    hist->Divide(hAllcore);
    hist->SetMarkerStyle(21);
    hist->SetMarkerColor(kMagenta);
    hist->SetLineColor(kMagenta);
    hist->Draw("same");
    leg->AddEntry(hist, "Bothcore/Allcore", "lep");

    
    leg->Draw();

    canv->SaveAs(Form("imgs/PSBSFits_ratio_core_%s_c%03i_%s.png", trigger, cent, graph));
    canv->SaveAs(Form("imgs/PSBSFits_ratio_core_%s_c%03i_%s.pdf", trigger, cent, graph));
    
  }
}


void DrawPSBSFitMethodeParams(const RawProduction::Output& output, const char* trigger, int cent)
{
  TStringToken mrst("mr1r;mr2r;mr1;mr2", ";");
  DrawAbs(output, trigger, cent, mrst);

  TStringToken srst("sr1r;sr2r;sr1;sr2", ";");
  DrawAbs(output, trigger, cent, srst);
  
    
  DrawRatios(output, trigger, cent);
}

void DrawPSBSFitParams()
{
  gROOT->LoadMacro("MakeRawProduction.C+g");
  RawProduction::Output output;
  gStyle->SetOptStat(0);
  
  
  DrawPSBSFitMethodeParams(output, "kMB", -10);
  DrawPSBSFitMethodeParams(output, "kPHOSPb", -10);

  DrawPSBSFitMethodeParams(output, "kCentral", -1);
  DrawPSBSFitMethodeParams(output, "kMB", -1);
  DrawPSBSFitMethodeParams(output, "kPHOSPb", -1);

  DrawPSBSFitMethodeParams(output, "kSemiCentral", -11);
  DrawPSBSFitMethodeParams(output, "kMB", -11);
  DrawPSBSFitMethodeParams(output, "kPHOSPb", -11);

  DrawPSBSFitMethodeParams(output, "kPHOSPb", -6);
  DrawPSBSFitMethodeParams(output, "kMB", -6);
}
