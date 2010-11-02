/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*   Origin: Alberica Toia, CERN, Alberica.Toia@cern.ch                   */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  class to determine centrality percentiles from 2D distributions          // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TH1D.h>
#include <TString.h>
#include <TFile.h>
#include <TMath.h>
#include <TROOT.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <vector>
#include <TMap.h>
#include <map>
#include "AliCentralityByFunction.h"


ClassImp(AliCentralityByFunction)  
 
//______________________________________________________________________________

Double_t fitf_pol2(Double_t* x, Double_t* par) {
  Double_t  fitValue  = 
    (par[0]+
     par[1]*x[0]+
     par[2]*x[0]*x[0]);
  return fitValue;
}
Double_t fitf_pol3(Double_t* x, Double_t* par) {
  Double_t  fitValue  = 
    (par[0]+
     par[1]*x[0]+
     par[2]*x[0]*x[0]+
     par[3]*x[0]*x[0]*x[0]);
  return fitValue;
}
Double_t fitf_pol4(Double_t* x, Double_t* par) {
  Double_t  fitValue  = 
    (par[0]+
     par[1]*x[0]+
     par[2]*x[0]*x[0]+
     par[3]*x[0]*x[0]*x[0]+
     par[4]*x[0]*x[0]*x[0]*x[0]);
  return fitValue;
}
Double_t fitf_pol6(Double_t* x, Double_t* par) {
  Double_t  fitValue  = 
    (par[0]+
     par[1]*x[0]+
     par[2]*x[0]*x[0]+
     par[3]*x[0]*x[0]*x[0]+
     par[4]*x[0]*x[0]*x[0]*x[0]+
     par[5]*x[0]*x[0]*x[0]*x[0]*x[0]+
     par[6]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
  return fitValue;
}

AliCentralityByFunction::AliCentralityByFunction() {
  // standard constructor
  fitter["fitf_pol2"] = new TF1("fitf_pol2",fitf_pol2,0,1,3);
  fitter["fitf_pol2"]->SetLineColor(kRed); 

  fitter["fitf_pol3"] = new TF1("fitf_pol3",fitf_pol3,0,1,4);
  fitter["fitf_pol3"]->SetLineColor(kRed); 

  fitter["fitf_pol4"] = new TF1("fitf_pol4",fitf_pol4,0,1,5);
  fitter["fitf_pol4"]->SetLineColor(kRed); 

  fitter["fitf_pol6"] = new TF1("fitf_pol6",fitf_pol6,0,1,6);
  fitter["fitf_pol6"]->SetLineColor(kRed); 
}

AliCentralityByFunction::~AliCentralityByFunction() {
  // destructor
}

void AliCentralityByFunction::AddHisto(TString name) {
  histnames.push_back(name);
}

void AliCentralityByFunction::SetPercentileFile(TString filename) {
  outrootfilename = filename;
}

void AliCentralityByFunction::SetPercentileCrossSection(Float_t xsec) {
  percentXsec = xsec;
}

void AliCentralityByFunction::SetFitFunction(TString distribution, TString func, Double_t xmin, Double_t xmax) {
  fitfunc[distribution] = func;
  fitter[fitfunc[distribution]]->SetRange(xmin,xmax);
}

void AliCentralityByFunction::MakePercentiles(TString infilename) {
  TH1D *hpercentile;
  
  // open inrootfile, outrootfile
  inrootfile  = new TFile(infilename);
  outrootfile = new TFile(outrootfilename,"RECREATE");
  
  // loop over all distribution names  
  vector<TString>::const_iterator hni;
  for(hni=histnames.begin(); hni!=histnames.end(); hni++) {
    hpercentile = FitHisto(*hni);
    outrootfile->cd();
    hpercentile->Write();
  }
  // close inrootfile, outrootfile
  inrootfile->Close();
  outrootfile->Close();
  
}

 TH1D *AliCentralityByFunction::FitHisto(TString hdistributionName) {
  TH2D *hdist  = (TH2D*) (inrootfile->Get(hdistributionName)); 
  TProfile *profile =hdist->ProfileX(); 
  //  fitter[fitfunc[hdistributionName]]->SetRange(0,profile->GetBinCenter(profile->GetNbinsX()));
  profile->Fit(fitter[fitfunc[hdistributionName]], "RNM");

  outrootfile->cd();
  profile->Write();
  fitter[fitfunc[hdistributionName]]->Write(hdistributionName.Append("_fit"));

  return MakePercentHisto(hdist);
}

TH1D * AliCentralityByFunction::MakePercentHisto(TH2D *histo) {
  TH2D *htemp = (TH2D*)histo->Clone("htemp");
  TString hdistributionName = htemp->GetTitle();

  const int NUM_DIVISION=500;
  
  TH1D *hpercent  = new TH1D("","",NUM_DIVISION,0,htemp->GetXaxis()->GetBinCenter(htemp->GetNbinsX()));
  hpercent->Reset();
  
  double Xpoint, Ypoint, Xcomp, Ycomp, count;
  
  double Xmax = htemp->GetXaxis()->GetBinCenter(htemp->GetNbinsX());
  double Xmin = 0;
  
  double delta = (Xmax - Xmin) / NUM_DIVISION;
  double slopePerp;   // slope of the perpendicular line
  double slopeFunction;

  cout << "Start Percentile Histo for distribution " << hdistributionName << " with fit function " << fitfunc[hdistributionName] << endl;

  // move perpendicular line along fitting curve
  // and count number of points on the right
  for (int i = 0; i <= NUM_DIVISION; i++) {
    Xcomp = Xmin + i * delta;
    Ycomp = fitter[fitfunc[hdistributionName]]->Eval(Xcomp);
    count = 0.0;
    slopeFunction = fitter[fitfunc[hdistributionName]]->Derivative(Xcomp,NULL,0.0001);
    slopePerp = -1.0 /slopeFunction;
    //
    // equation of perpendicular line
    // (y-Ycomp)/(x-Xcomp) = slopePerp
    //
    for (int ibiny = 1; ibiny < htemp->GetNbinsY(); ibiny++) {
      Ypoint = htemp->GetYaxis()->GetBinCenter(ibiny);

      double XonLine =(Ypoint - Ycomp)/slopePerp + Xcomp;
      //      double YonLine = slopePerp*(XonLine - Xcomp) +Ycomp;
      
      for (int ibinx = 1; ibinx < htemp->GetNbinsX(); ibinx++) {
	Xpoint = htemp->GetXaxis()->GetBinCenter(ibinx);
	//
	// (XonLine,Ypoint) lies on the perpendicular
	//
	 if ((Xpoint > XonLine && Xpoint > Xcomp - 0.5)||
	     (Xpoint > Xcomp + 0.2)) {
		  
	  count  += htemp->GetBinContent(ibinx,ibiny);
	}
      }
    }
    count = count/htemp->Integral() * 100.0;   // change to percentage
    hpercent->SetBinContent(i, count);
  }

  hpercent->SetName(hdistributionName.Append("_percentile"));

  return hpercent;
}



