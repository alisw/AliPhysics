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
AliCentralityByFunction::AliCentralityByFunction() :
  finrootfile(0),
  foutrootfilename(0),
  foutrootfile(0),
  fhistnames(),
  fpercentXsec(0),
  fitfunc(),
  fitter()
{
  // standard constructor
  fitter["fitf_pol2"] = new TF1("pol2",0,1,3);
  fitter["fitf_pol2"]->SetLineColor(kRed); 

  fitter["fitf_pol3"] = new TF1("pol3",0,1,4);
  fitter["fitf_pol3"]->SetLineColor(kRed); 

  fitter["fitf_pol4"] = new TF1("pol4",0,1,5);
  fitter["fitf_pol4"]->SetLineColor(kRed); 

  fitter["fitf_pol6"] = new TF1("pol6",0,1,6);
  fitter["fitf_pol6"]->SetLineColor(kRed); 
}

void AliCentralityByFunction::SetFitFunction(TString distribution, TString func, Double_t xmin, Double_t xmax) 
{
  // Set fit function
  fitfunc[distribution] = func;
  fitter[fitfunc[distribution]]->SetRange(xmin,xmax);
}

void AliCentralityByFunction::MakePercentiles(TString infilename) 
{
  // Make percentile bins

  TH1D *hpercentile;
  
  // open inrootfile, outrootfile
  finrootfile  = new TFile(infilename);
  foutrootfile = new TFile(foutrootfilename,"RECREATE");
  
  // loop over all distribution names  
  std::vector<TString>::const_iterator hni;
  for(hni=fhistnames.begin(); hni!=fhistnames.end(); hni++) {
    hpercentile = FitHisto(*hni);
    foutrootfile->cd();
    hpercentile->Write();
  }
  // close inrootfile, outrootfile
  finrootfile->Close();
  foutrootfile->Close();
}

TH1D *AliCentralityByFunction::FitHisto(TString hdistributionName) 
{
  // Fit histogram
  TH2D *hdist  = (TH2D*) (finrootfile->Get(hdistributionName)); 
  TProfile *profile =hdist->ProfileX(); 
  //  fitter[fitfunc[hdistributionName]]->SetRange(0,profile->GetBinCenter(profile->GetNbinsX()));
  profile->Fit(fitter[fitfunc[hdistributionName]], "RNM");

  foutrootfile->cd();
  profile->Write();
  fitter[fitfunc[hdistributionName]]->Write(hdistributionName.Append("_fit"));

  return MakePercentHisto(hdist);
}

TH1D * AliCentralityByFunction::MakePercentHisto(TH2D *histo) 
{
  TH2D *htemp = (TH2D*)histo->Clone("htemp");
  TString hdistributionName = htemp->GetTitle();

  const int num_DIVISION=500;
  
  TH1D *hpercent  = new TH1D("","",num_DIVISION,0,htemp->GetXaxis()->GetBinCenter(htemp->GetNbinsX()));
  hpercent->Reset();
  
  double xpoint, ypoint, xcomp, ycomp, count;
  
  double xmax = htemp->GetXaxis()->GetBinCenter(htemp->GetNbinsX());
  double xmin = 0;
  
  double delta = (xmax - xmin) / num_DIVISION;
  double slopePerp;   // slope of the perpendicular line
  double slopeFunction;

  std::cout << "Start Percentile Histo for distribution " << hdistributionName << " with fit function " << fitfunc[hdistributionName] << std::endl;

  // move perpendicular line along fitting curve
  // and count number of points on the right
  for (int i = 0; i <= num_DIVISION; i++) {
    xcomp = xmin + i * delta;
    ycomp = fitter[fitfunc[hdistributionName]]->Eval(xcomp);
    count = 0.0;
    slopeFunction = fitter[fitfunc[hdistributionName]]->Derivative(xcomp,NULL,0.0001);
    slopePerp = -1.0 /slopeFunction;
    //
    // equation of perpendicular line
    // (y-Ycomp)/(x-Xcomp) = slopePerp
    //
    for (int ibiny = 1; ibiny < htemp->GetNbinsY(); ibiny++) {
      ypoint = htemp->GetYaxis()->GetBinCenter(ibiny);

      double xonLine =(ypoint - ycomp)/slopePerp + xcomp;
      //      double YonLine = slopePerp*(XonLine - Xcomp) +Ycomp;
      
      for (int ibinx = 1; ibinx < htemp->GetNbinsX(); ibinx++) {
	xpoint = htemp->GetXaxis()->GetBinCenter(ibinx);
	//
	// (XonLine,Ypoint) lies on the perpendicular
	//
	 if ((xpoint > xonLine && xpoint > xcomp - 0.5)||
	     (xpoint > xcomp + 0.2)) {
		  
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
