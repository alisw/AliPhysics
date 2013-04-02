#include "TCanvas.h"
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TLegend.h>
#include <TString.h>

namespace RawProduction {
  class Output;
}

TF1* GetEfficency(const char* trigger, int cent, const char* pid)
{
  if( -1 == cent )
    cent = 0;
  char fileName[256]; sprintf(fileName, "eff_pi0_%s.root", trigger);
  char funcName[256]; sprintf(funcName, "eff_pi0_%s_cen%i", pid, cent);
  
  TFile* file = TFile::Open(fileName, "READ");
  TF1* func = (TF1*)file->Get(funcName);
  return func;
}

// void DrawAbs(const RawProduction::Output& output, const char* trigger, int cent, TStringToken& names) 
// {
//   names.NextToken();
// 
//   const char* pid = "All";
//   char canvName[256] = Form("PSBS_Abs_%s_c%03i_%s_%s", trigger, cent, pid, names.Data());
//   TCanvas* canv = new TCanvas(canvName, canvName);  
//   TLegend* leg = new TLegend(0.6,0.8,0.95,0.95);
// 
//   TH1* hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, pid, names.Data()));
//   if(names.Contains("mr")) hist->SetTitle(Form("Peak Position, %s, %s, %s", trigger, RawProduction::GetCentString(cent), pid));
//   if(names.Contains("sr")) hist->SetTitle(Form("Peak Width, %s, %s, %s", trigger, RawProduction::GetCentString(cent), pid));
//   if(names.Contains("mr")) hist->GetYaxis()->SetRangeUser(0.12, 0.15);
//   if(names.Contains("mr")) hist->GetYaxis()->SetTitle("Peak #mu");
//   if(names.Contains("sr")) hist->GetYaxis()->SetRangeUser(0., 0.012);
//   if(names.Contains("sr")) hist->GetYaxis()->SetTitle("Peak #sigma");
//   hist->GetXaxis()->SetTitle("p_{T}");
//   //Printf(hist->GetTitle());
//   hist->SetMarkerStyle(21);
//   //hist->SetMarkerSize(1.5);
//   hist->SetMarkerColor(kBlack);
//   hist->SetLineColor(kBlack);
//   hist->Draw();
//   leg->AddEntry(hist, "Pol1, Ratio", "lep");
//   
//   int marker = 21;
//   Color_t color[3] = {kRed, kBlue, kMagenta};
//   while( names.NextToken() ) {
//     hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "All", names.Data()));
//     //Printf(hist->GetName());
//     hist->SetMarkerStyle(++marker);
//     hist->SetMarkerColor(color[marker-22]);
//     hist->SetLineColor(color[marker-22]);
//     hist->Draw("same");
//     char legName[256] = "";
//     if(names.Contains("1")) sprintf(legName, "Pol1");
//     if(names.Contains("2")) sprintf(legName, "Pol2");
//     if( marker <23 ) sprintf(legName, "%s, Ratio", legName);
//     leg->AddEntry(hist, legName, "lep");
//   }
//   
//   hist = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, pid, names.Data()));
//   hist->Draw("same");
// 
//   leg->Draw();
//   
//   canv->SaveAs(Form("imgs/%s.png", canvName));
//   canv->SaveAs(Form("imgs/%s.pdf", canvName));
// }

void DrawRatios(const RawProduction::Output& output, const char* trigger, int cent)
{
  const int nProdTypes = 2;
  enum ProdType { RAW=0, EffCor=1 };
  
  for(int prodType=0; prodType < nProdTypes; ++prodType) {  
    TStringToken graphs("yr1 yr1int yr2 yr2int", " ");
    while ( graphs.NextToken() ) {
      char graphName[32]; sprintf(graphName, "%s", graphs.Data());

      // All
      TH1* hAll = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "All", graphName));
      char newName[256] = Form("E*%s_All", graphName);
      hAll = (TH1*) hAll->Clone(newName);
      hAll->SetTitle(Form("%s, %s, %s, %s", graphName, trigger, RawProduction::GetCentString(cent), graphName));
      hAll->GetXaxis()->SetTitle("IM_{#gamma #gamma}");
      hAll->SetLineColor(kBlack);
      hAll->SetMarkerColor(kBlack);
      hAll->SetMarkerStyle(20);
      hAll->GetXaxis()->SetLabelFont(63); //font in pixels
      hAll->GetXaxis()->SetLabelSize(16); //in pixels
      hAll->GetYaxis()->SetLabelFont(63); //font in pixels
      hAll->GetYaxis()->SetLabelSize(16); //in pixels

      // Allcore
      TH1* hAllcore = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "Allcore", graphName));
      sprintf(newName, "E*%s_Allcore", graphName);
      hAllcore = (TH1*) hAllcore->Clone(newName);
      hAllcore->SetLineColor(kCyan);
      hAllcore->SetMarkerColor(kCyan);
      hAllcore->SetMarkerStyle(25);
  
      // CPV
      TH1* hCPV = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "CPV", graphName));
      sprintf(newName, "E*%s_CPV", graphName);
      hCPV = (TH1*) hCPV->Clone(newName);
      hCPV->SetLineColor(kBlue);
      hCPV->SetMarkerColor(kBlue);
      hCPV->SetMarkerStyle(22);

      // CPVcore
      TH1* hCPVcore = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "CPVcore", graphName));
      sprintf(newName, "E*%s_CPVcore", graphName);
      hCPVcore = (TH1*) hCPVcore->Clone(newName);
      hCPVcore->SetLineColor(kGreen);
      hCPVcore->SetMarkerColor(kGreen);
      hCPVcore->SetMarkerStyle(26);

      // Both
      TH1* hBoth = output.GetHistogram(Form("%s/c%03i/%s/%s", trigger, cent, "Both", graphName));
      sprintf(newName, "E*%s_Both", graphName);
      hBoth = (TH1*) hBoth->Clone(newName);
      hBoth->SetLineColor(kRed);
      hBoth->SetMarkerColor(kRed);
      hBoth->SetMarkerStyle(33);

      // Efficiancy correction
      if( EffCor == prodType ) {
	hAll->GetYaxis()->SetTitle("#frac{d^{2}N_{#pi^{0}}}{p_{T}dp_{T}dy N_{ev}}");
	hAll->Divide(GetEfficency(trigger, cent, "All"));
	hAllcore->Divide(GetEfficency(trigger, cent, "Allcore"));
	hCPV->Divide(GetEfficency(trigger, cent, "CPV"));
	hCPVcore->Divide(GetEfficency(trigger, cent, "CPVcore"));
	hBoth->Divide(GetEfficency(trigger, cent, "Both"));
      }
      
      // hAll / hAllcore
      TH1* hAll_hAll = hAll->Clone(Form("%s_div_hAll", hAll->GetName()));
      hAll_hAll->Divide(hAll);
      // hAllcore / hAll
      TH1* hAllcore_hAll = hAllcore->Clone(Form("%s_div_hAll", hAllcore->GetName()));
      hAllcore_hAll->Divide(hAll);
      // hCPV / hAll
      TH1* hCPV_hAll = hCPV->Clone(Form("%s_div_hAll", hCPV->GetName()));
      hCPV_hAll->Divide(hAll);
      // hAllcore / hCPVcore
      TH1* hCPVcore_hAllcore = hCPVcore->Clone(Form("%s_div_Allcore", hCPVcore->GetName()));
      hCPVcore_hAllcore->Divide(hAllcore);
      // hBoth / hAll
      TH1* hBoth_hAll = hBoth->Clone(Form("%s_div_hAll", hBoth->GetName()));
      hBoth_hAll->Divide(hAll);

  
      TCanvas* canv = new TCanvas(Form("%s_c%03i_%s_ratio", trigger, cent, graphName), Form("%s_c%03i_%s_ratio", trigger, cent, graphName));
  
      canv->cd();
      TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
      pad1->SetBottomMargin(0);
      pad1->SetLogy();
      pad1->Draw();
      pad1->cd();
      hAll->GetYaxis()->SetTitleSize(0.06);
      hAll->GetYaxis()->SetTitleOffset(0.6);
      hAll->GetYaxis()->SetRangeUser(2.0e-7, 90.);
      hAll->DrawCopy();
      hAllcore->DrawCopy("same");
      hCPV->DrawCopy("same");
      hCPVcore->DrawCopy("same");
      hBoth->DrawCopy("same");
      TLegend* leg1 = new TLegend(0.7,0.6,0.85,0.88);
      leg1->AddEntry(hAll, "All", "lep");
      leg1->AddEntry(hAllcore, "Allcore", "lep");
      leg1->AddEntry(hCPV, "CPV", "lep");
      leg1->AddEntry(hCPVcore, "CPVcore", "lep");
      leg1->AddEntry(hBoth, "Both", "lep");
      leg1->Draw();
  
      canv->cd();
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
      pad2->SetTopMargin(0);
      pad2->Draw();
      pad2->cd();
      pad2->SetGridy();
      hAll_hAll->GetYaxis()->SetRangeUser(0.2, 1.4);
      hAll_hAll->GetYaxis()->SetTitle("");
      hAll_hAll->SetTitle("");
      hAll_hAll->DrawCopy("AXIS");
      hAll_hAll->DrawCopy("AXIGsame");
      hAllcore_hAll->DrawCopy("same");
      hCPV_hAll->DrawCopy("same");
      hCPVcore_hAllcore->DrawCopy("same");
      hBoth_hAll->DrawCopy("same");
      TLegend* leg2 = new TLegend(0.7,0.63,0.85,0.98);
      //leg2->AddEntry(hAll_hAll, "All/All", "lep");
      leg2->AddEntry(hAllcore_hAll, "Allcore/All", "lep");
      leg2->AddEntry(hCPV_hAll, "CPV/All", "lep");
      leg2->AddEntry(hCPVcore_hAllcore, "CPVcore/Allcore", "lep");
      leg2->AddEntry(hBoth_hAll, "Both/All", "lep");
      leg2->Draw();
  
      if(RAW == prodType) {
	canv->SaveAs(Form("imgs/prod_ratios_%s_%i_%s_RAW.pdf", trigger, cent, graphName));      
	canv->SaveAs(Form("imgs/prod_ratios_%s_%i_%s_RAW.png", trigger, cent, graphName));      
      } else if (EffCor == prodType ) {
	canv->SaveAs(Form("imgs/prod_ratios_%s_%i_%s.pdf", trigger, cent, graphName));      
	canv->SaveAs(Form("imgs/prod_ratios_%s_%i_%s.png", trigger, cent, graphName));      
      }
      delete canv;
    }
  }
    
}


void DrawProductions(const RawProduction::Output& output, const char* trigger, int cent)
{
  // TStringToken mrst("mr1r;mr2r;mr1;mr2", ";");
  // DrawAbs(output, trigger, cent, mrst);
    
  DrawRatios(output, trigger, cent);
}

void DrawProduction()
{
  gROOT->LoadMacro("MakeRawProduction.C+g");
  RawProduction::Output output;
  gStyle->SetOptStat(0);
  
  
//   DrawProductions(output, "kMB", -10);
//   DrawProductions(output, "kPHOSPb", -10);

  DrawProductions(output, "kCentral", -1);
//   DrawProductions(output, "kMB", -1);
//   DrawProductions(output, "kPHOSPb", -1);
// 
//   DrawProductions(output, "kSemiCentral", -11);
//   DrawProductions(output, "kMB", -11);
//   DrawProductions(output, "kPHOSPb", -11);
// 
//   DrawProductions(output, "kPHOSPb", -6);
//   DrawProductions(output, "kMB", -6);
}
