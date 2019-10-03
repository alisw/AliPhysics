#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"

#include <Riostream.h>

//
// Multiply two histos of variable bin width
//
TH1 * MultiplyHistos(TH1 *h1, TH1 *h2) {
  
  Int_t nbins1 = h1->GetNbinsX();
  Int_t nbins2 = h2->GetNbinsX();
  Double_t binwidth1 = h1->GetBinWidth(1);
  Double_t binwidth2 = h2->GetBinWidth(1);

  if ( nbins1!=nbins2 || binwidth1!=binwidth2 ) {
    cout << "Histos do not have the same binning, do not seem compatible" << endl;
    return NULL;
  }

  //
  // Get the bins & limits
  Double_t *limits = new Double_t[nbins1+1];
  Double_t xlow=0., binwidth=0.;
  for (Int_t i=1; i<=nbins1; i++) {
    binwidth = h1->GetBinWidth(i);
    xlow = h1->GetBinLowEdge(i);
    limits[i-1] = xlow;
  }
  limits[nbins1] = xlow + binwidth;

  TH1D *hMultiply = new TH1D("hMultiply","hMultiply",nbins1,limits);

  Double_t value=0., err=0.;
  for (Int_t ibin=1; ibin<=nbins1; ibin++) {
    value = h1->GetBinContent(ibin) * h2->GetBinContent(ibin);
    err = value * TMath::Sqrt(  (h1->GetBinError(ibin)/h1->GetBinContent(ibin)) * (h1->GetBinError(ibin)/h1->GetBinContent(ibin))  +
				(h2->GetBinError(ibin)/h2->GetBinContent(ibin)) * (h2->GetBinError(ibin)/h2->GetBinContent(ibin))   );
    hMultiply->SetBinContent(ibin,value);
    hMultiply->SetBinError(ibin,err);
  }
  
  return (TH1*)hMultiply;
}

//
// Main function
// 
void ComputeEfficiencyInputFromTwoSteps (const char* recolhc10d3filename="Distributions.root",
					 const char* recolhc10d3histoname="RECPIDpt",
					 const char* simuAcclhc10d3filename="Distributions.root",
					 const char* simuAcclhc10d3histoname="MCAccpt",
					 const char* simuLimAcclhc10d4filename="",
					 const char* simuLimAcclhc10d4histoname="MCLimAccpt",
					 const char* simuAcclhc10d4filename="",
					 const char* simuAcclhc10d4histoname="MCAccpt",
					 const char* outfilename="ComputeEfficiencyInputFromTwoSteps.root") 
{

  TFile *recolhc10d3file = new TFile(recolhc10d3filename,"read");
  TH1D *hrecolhc10d3 = (TH1D*)recolhc10d3file->Get(recolhc10d3histoname);

  TFile *simuAcclhc10d3file = new TFile(simuAcclhc10d3filename,"read");
  TH1D *hsimuAcclhc10d3 = (TH1D*)simuAcclhc10d3file->Get(simuAcclhc10d3histoname);

  TFile *simuLimAcclhc10d4file = new TFile(simuLimAcclhc10d4filename,"read");
  TH1D *hsimuLimAcclhc10d4 = (TH1D*)simuLimAcclhc10d4file->Get(simuLimAcclhc10d4histoname);

  TFile *simuAcclhc10d4file = new TFile(simuAcclhc10d4filename,"read");
  TH1D *hsimuAcclhc10d4 = (TH1D*)simuAcclhc10d4file->Get(simuAcclhc10d4histoname);


  TFile *out = new TFile(outfilename,"recreate");
  TH1D *hRecoPIDCorr = (TH1D*)MultiplyHistos(hrecolhc10d3,hsimuAcclhc10d4);
  hRecoPIDCorr->SetNameTitle("hRecoPIDCorr","hRecoPIDCorr");
  TH1D *hSimuCorr = (TH1D*)MultiplyHistos(hsimuAcclhc10d3,hsimuLimAcclhc10d4);
  hSimuCorr->SetNameTitle("hSimuCorr","hSimuCorr");

  out->cd(); 
  hRecoPIDCorr->Write();
  hSimuCorr->Write();
  out->Close();

}
