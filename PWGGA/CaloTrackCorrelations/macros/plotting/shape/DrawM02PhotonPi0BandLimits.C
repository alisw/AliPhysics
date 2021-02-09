/// \file DrawM02PhotonPi0BandLimits.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Plot Shower shape vs E 2D histogram with limits
///
/// Macro to plot cluster 2D histograms
/// shower shape long axis vs energy and the pi0 and photon selection bands.
/// Several performance plots produced with this or similar macro.
///

//---------------------------------------------------------
// Set includes and declare methods for compilation

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TFile.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TObject.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLine.h>
#include <TF1.h>
#include <TMath.h>
#include <TLatex.h>
#include <TPaveText.h>

#endif

/// Main method
/// \param filepath  : string with input file path
/// \param filename  : string with input file name
/// \param dirname   : string with TDirectoyFileName file name
/// \param listname  : string with TList file name
/// \param histoname : string 2D histogram name, it should be M02 vs E or pT
///
void DrawM02PhotonPi0BandLimits
(
 TString filepath  = "data/", 
 TString filename  = "LHC11cd_EMC7",  
 TString dirname   = "",  
 TString listname  = "",  
 TString histoname = "AnaPhoton_hLam0E"  
)
{
  gStyle->SetOptTitle(0);
  gStyle->SetTitleOffset(2.0,"Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  gStyle->SetPadRightMargin(0.17);
  gStyle->SetPadLeftMargin(0.08);
  gStyle->SetPadTopMargin(0.02);
  //gStyle->SetPadBottomMargin(0.15);
  
  //
  // Open file and get histogram
  //
  TFile * file = TFile::Open(Form("%s/%s.root",filepath.Data(),filename.Data()));
  if ( !file ) { printf("No file %s\n",filename.Data()) ; return ; }

  TH2F  * hM02 = 0;
  
  if ( dirname!="" ) 
  {
    TDirectoryFile * dir = (TDirectoryFile *) file->Get(dirname);
    if ( !dir  ) { printf("No directory %s\n", dirname .Data()) ; return ; }
    
    TList * list = (TList*) dir->Get(listname);
    if ( !list ) { printf("No list %s\n"    ,  listname.Data()) ; return ; }

    hM02 = (TH2F*) list->FindObject(histoname);
  }
  else if ( listname!="" ) 
  {
    TList * list = (TList*) file->Get(listname);
    if ( !list ) { printf("No list %s\n"    ,  listname.Data()) ; return ; }
    
    hM02 = (TH2F*) list->FindObject(histoname);
  }
  else
    hM02 = (TH2F*) file->Get(histoname);
  
  //
  // Normalize 2D to bin energy
  //
  Int_t nbinsy = hM02->GetNbinsY();
  for(Int_t j = 1; j <= hM02->GetNbinsX(); j++)
  {
    TH1D* temp1 = hM02->ProjectionY(Form("Bin%d",j),j,j);
  
    Float_t scale1 = temp1 -> Integral(-1,-1);
  
    for(Int_t i = 1; i <= nbinsy; i++)
    {
      //printf("NLM2: i %d, j %d;  content %f / scale %f = %f\n",i,j,hNLM2->GetBinContent(j,i),scale2,hNLM2->GetBinContent(j,i)/scale2);
      if(scale1>0)
      {
        hM02->SetBinContent(j,i, hM02->GetBinContent(j,i)/scale1);
        hM02->SetBinError  (j,i, hM02->GetBinError  (j,i)/scale1);
      }
      else
      {
        hM02->SetBinContent(j,i, 0);
        hM02->SetBinError  (j,i, 0);
      }
    }
  }
  
  //
  // Selection band min/max lines
  // Based on AN145: https://alice-notes.web.cern.ch/node/145
  //
  TF1 *lM02MinNLM1 = new TF1("M02MinNLM1","exp(2.135-0.245*x)",6,13.6);
  lM02MinNLM1->SetLineColor(6);
  lM02MinNLM1->SetLineWidth(3);
  
  TF1 *lM02MinNLM2 = new TF1("M02MinNLM2","exp(2.135-0.245*x)",6,13.6);
  lM02MinNLM2->SetLineColor(6);
  lM02MinNLM2->SetLineWidth(3);

  TF1 *lM02MaxNLM1 = new TF1("M02MaxNLM1","exp(0.0662-0.0201*x)-0.0955+0.00186*x[0]+9.91/x[0]",6,100);
  lM02MaxNLM1->SetLineColor(6);
  lM02MaxNLM1->SetLineWidth(3);
  
  TF1 *lM02MaxNLM2 = new TF1("M02MaxNLM2","exp(0.353-0.0264*x)-0.524+0.00559*x[0]+21.9/x[0]",6,100);
  lM02MaxNLM2->SetLineColor(6);
  lM02MaxNLM2->SetLineWidth(3);

  TLine *l03 = new TLine(13.5,0.3,40,0.3);
  l03->SetLineColor(6);
  l03->SetLineWidth(3);
  l03->SetLineStyle(1);

  TLine *lGammaMax = new TLine(5,0.27,40,0.27);
  lGammaMax->SetLineColor(kRed);
  lGammaMax->SetLineWidth(3);
  lGammaMax->SetLineStyle(2);
  
  TLine *lGammaMin = new TLine(5,0.1,40,0.1);
  lGammaMin->SetLineColor(kRed);
  lGammaMin->SetLineWidth(3);
  lGammaMin->SetLineStyle(2);

  //
  // Plot
  //
  TCanvas * c = new TCanvas("cM02","M02",500,500);
  
  gPad->SetLogz();
  
  //hM02->Rebin2D(1,1);
  hM02->SetAxisRange(5,40,"X");
  hM02->SetAxisRange(0,2.4,"y");
  hM02->SetZTitle("#it{E} bin normalized to integral");
  hM02->SetXTitle("#it{E} (GeV)");
  hM02->SetYTitle("#sigma^{2}_{long}");
  hM02->SetTitleOffset(1.5,"Z");
  hM02->SetTitleOffset(1.05,"Y");

  hM02->Draw("colz");

  lM02MinNLM1->Draw("same");
  //lM02MinNLM2->Draw("same");
  //lM02MaxNLM1->Draw("same");
  lM02MaxNLM2->Draw("same");
  l03->Draw("same");
  lGammaMin->Draw("same");
  lGammaMax->Draw("same");

  TLegend *l2 = new TLegend(0.5,0.65,0.77,0.75);
  l2->SetTextSize(0.035);
  l2->AddEntry(lM02MaxNLM1,"#pi^{0} band limits","L");
  l2->AddEntry(lGammaMin,"#gamma band limits","L");
  l2->SetLineColor(0);

  l2->Draw();

  TPaveText *pt = new TPaveText(0.45,0.77,0.82,0.97,"brNDC");
  pt->SetTextSize(0.035);
  pt->AddText("pp, #sqrt{#it{s}} = 7 TeV");
  pt->AddText("EMCal-L0 trigger");
  //pt->AddLine(.0,.5,1.,.5);
//  pt->AddText("ALICE Performance");
//  pt->AddText("08.09.2016");
  pt->SetFillColor(kWhite);
  pt->SetLineColor(1);
  pt->Draw("same");
  
  c->Print("M02_BandLimits.eps");
  //c->Print(Form("%s/M02_BandLimits_%s.eps",filepath.Data(),filename.Data()));
}
