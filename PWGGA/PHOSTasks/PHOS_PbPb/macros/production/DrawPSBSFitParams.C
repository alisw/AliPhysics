#include "TCanvas.h"


const char* pids[12] = {"Bothcore", "Allcore", "Allwou", "Disp2", "Dispwou","Dispcore", "CPV2", "CPVcore", "Both", "Disp", "CPV", "All"};
Int_t markers[12] =  {  25,          24,       24,        32,     30,         32,         26,        26,     21,     23,     22,    20};
Float_t markerSizes[12] =  {  1.2,   0.8,      1.2,     0.8,      1.2,        1.6,        0.8,       1.2,     1.,     1,       1,     1.  };
Int_t colors[12] =    {kMagenta-4, kGray, kGray+1, kBlue-9, kBlue-7, kBlue-4,   kRed-7,   kRed-4, kMagenta+1, kBlue, kRed, kBlack};
TH1* hists[12] = {0};
int from =0;

namespace RawProduction {
  class Output;
}


void DrawPSBSFitMethodeParams(const RawProduction::Output& output, const char* trigger, int cent)
{
  TStringToken graphs("mr1;sr1;mr1r;sr1r;mr2;sr2;mr2r;sr2r", ";");
  while( graphs.NextToken() ) {
    const char* graph = graphs.Data();
    
    TCanvas* canv = new TCanvas(Form("%s_c%03i_%s", trigger, cent, graph), Form("%s_c%03i_%s", trigger, cent, graph));
    TH1* hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, pids[from], graph));
    hists[from] = hist;
    
    if(graphs.Contains("mr")) hist->GetYaxis()->SetRangeUser(0.11, 0.15);
    if(graphs.Contains("sr")) hist->GetYaxis()->SetRangeUser(0, 0.012);
    
    if(graphs.Contains("r1")) hist->SetTitle( Form("%s, Pol1", hist->GetTitle()) );
    if(graphs.Contains("r2")) hist->SetTitle( Form("%s, Pol1", hist->GetTitle()) );
    hist->SetMarkerStyle(markers[from]);
    hist->SetMarkerColor(colors[from]);
    hist->SetLineColor(colors[from]);
    hist->Draw();
    
    
    for(int index=from+1; index<12; index++  ) {
      hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, pids[index], graph));
      hists[index] = hist;
      hist->SetMarkerStyle(markers[index]);
      hist->SetMarkerSize(markerSizes[index]);
      hist->SetMarkerColor(colors[index]);
      hist->SetLineColor(colors[index]);
      hist->Draw("same");
      
    }
    
    TLegend* leg = new TLegend(0.8,0.6,0.95,0.95);
    for(int index=11; index>=from; index--  )
      leg->AddEntry(hists[index], pids[index], "lep");
    leg->Draw();
    canv->SaveAs(Form("imgs/PSBSFits_%s_c%03i_%s.png", trigger, cent, graph));
    canv->SaveAs(Form("imgs/PSBSFits_%s_c%03i_%s.pdf", trigger, cent, graph));



    // // core
    // canv = new TCanvas(Form("%s_c%03i_%s_core", trigger, cent, graph), Form("%s_c%03i_%s_core", trigger, cent, graph));
    // leg = new TLegend(0.8,0.6,0.95,0.95);

    // hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "Allcore", graph));
    // if(graphs.Contains("mr")) hist->GetYaxis()->SetRangeUser(0.11, 0.15);
    // if(graphs.Contains("sr")) hist->GetYaxis()->SetRangeUser(0, 0.012);
    // if(graphs.Contains("r1")) hist->SetTitle( Form("%s, Pol1", hist->GetTitle()) );
    // if(graphs.Contains("r2")) hist->SetTitle( Form("%s, Pol2", hist->GetTitle()) );
    // hist->SetMarkerStyle(20);
    // hist->SetMarkerColor(kBlack);
    // hist->SetLineColor(kBlack);
    // hist->Draw();
    // leg->AddEntry(hist, "Allcore", "lep");

    // hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "CPVcore", graph));
    // hist->SetMarkerStyle(22);
    // //hist->SetMarkerSize();
    // hist->SetMarkerColor(kRed);
    // hist->SetLineColor(kRed);
    // hist->Draw("same");
    // leg->AddEntry(hist, "CPVcore", "lep");

    // hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "Dispcore", graph));
    // hist->SetMarkerStyle(23);
    // //hist->SetMarkerSize();
    // hist->SetMarkerColor(kBlue);
    // hist->SetLineColor(kBlue);
    // hist->Draw("same");
    // leg->AddEntry(hist, "Dispcore", "lep");

    // hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "Bothcore", graph));
    // hist->SetMarkerStyle(21);
    // //hist->SetMarkerSize();
    // hist->SetMarkerColor(kMagenta);
    // hist->SetLineColor(kMagenta);
    // hist->Draw("same");
    // leg->AddEntry(hist, "Bothcore", "lep");

    
    // leg->Draw();

    // canv->SaveAs(Form("imgs/PSBSFits_core_%s_c%03i_%s.png", trigger, cent, graph));
    // canv->SaveAs(Form("imgs/PSBSFits_core_%s_c%03i_%s.pdf", trigger, cent, graph));


    // core, ratio
    canv = new TCanvas(Form("%s_c%03i_%s_ratio_core", trigger, cent, graph), Form("%s_c%03i_%s_ratio_core", trigger, cent, graph));
    leg = new TLegend(0.6,0.8,0.95,0.95);

    hAll = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "All", graph));
    if(graphs.Contains("r1")) hist->SetTitle( Form("%s, Pol1", hist->GetTitle()) );
    if(graphs.Contains("r2")) hist->SetTitle( Form("%s, Pol2", hist->GetTitle()) );

    hAllcore = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "Allcore", graph));

    hist = hAllcore;
    hist = (TH1*) hist->Clone(Form("dAll_%s", hist->GetName()));
    hist->Divide(hAll);
    if(graphs.Contains("mr")) hist->GetYaxis()->SetRangeUser(0.9, 1.05);
    if(graphs.Contains("sr")) hist->GetYaxis()->SetRangeUser(0.4, 1.4);
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

void DrawPSBSFitParams()
{
  gROOT->LoadMacro("MakeRawProduction.C+g");
  RawProduction::Output output;
  gStyle->SetOptStat(0);
  
  
  DrawPSBSFitMethodeParams(output, "kCentral", -1);
  DrawPSBSFitMethodeParams(output, "kMB", -10);
  DrawPSBSFitMethodeParams(output, "kPHOSPb", -10);
  
}
