// ----------------------------------------------------------------------
//                     AliBWTools
// 
// This class provides some tools which can be useful in the analsis
// of spectra, to fit or transform histograms. See the comments of the
// individual methods for details
//
// Author: M. Floris (CERN)
// ----------------------------------------------------------------------

#include "AliBWTools.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include "TMinuit.h"
#include "AliLog.h"
#include <iostream>

using namespace std;


ClassImp(AliBWTools)

TF1 * AliBWTools::fdNdptForETCalc = 0;

AliBWTools::AliBWTools() {
  // ctor
}

AliBWTools::~AliBWTools(){
  // dtor
}

TH1 * AliBWTools::GetdNdmtFromdNdpt(const TH1 * hpt, Double_t mass) {
  // convert the x axis from pt to mt. Assumes you have 1/pt dNdpt in the histo you start with

  Int_t nbins = hpt->GetNbinsX();
  Float_t * xbins = new Float_t[nbins+1];
  for(Int_t ibins = 0; ibins <= nbins; ibins++){
    xbins[ibins] = TMath::Sqrt(hpt->GetBinLowEdge(ibins+1)*hpt->GetBinLowEdge(ibins+1) +
    			       mass *mass) - mass;
    // // xbins[ibins] = TMath::Sqrt(hpt->GetBinLowEdge(ibins+1)*hpt->GetBinLowEdge(ibins+1) +
    // // 			       mass *mass);
    //    cout << ibins << " "<< xbins[ibins]  << endl;

  }

  TH1D * hmt = new TH1D(TString(hpt->GetName())+"_mt",
		      TString(hpt->GetName())+"_mt",
		      nbins, xbins);
  for(Int_t ibins = 1; ibins <= nbins; ibins++){
    hmt->SetBinContent(ibins, hpt->GetBinContent(ibins));
    hmt->SetBinError(ibins,   hpt->GetBinError(ibins));

  }

  hmt->SetXTitle("m_{t} - m_{0} (GeV/c^{2})");
  hmt->SetYTitle("1/m_{t} dN/dm_{t} (a.u.)");
  hmt->SetMarkerStyle(hpt->GetMarkerStyle());
  hmt->SetMarkerColor(hpt->GetMarkerColor());
  hmt->SetLineColor(hpt->GetLineColor());

  return hmt;

}

TH1 * AliBWTools::GetdNdptFromdNdmt(const TH1 * hmt, Double_t mass) {
  // convert the x axis from mt to pt. Assumes you have 1/mt dNdmt in the histo you start with

  Int_t nbins = hmt->GetNbinsX();
  Float_t * xbins = new Float_t[nbins+1];
  for(Int_t ibins = 0; ibins <= nbins; ibins++){
    xbins[ibins] = TMath::Sqrt((hmt->GetBinLowEdge(ibins+1)+mass)*(hmt->GetBinLowEdge(ibins+1)+mass) -
    			       mass *mass);
    xbins[ibins] = Float_t(TMath::Nint(xbins[ibins]*100))/100;
    // // xbins[ibins] = TMath::Sqrt(hmt->GetBinLowEdge(ibins+1)*hmt->GetBinLowEdge(ibins+1) +
    // // 			       mass *mass);
    cout << ibins << " "<< xbins[ibins]  << endl;

  }

  TH1D * hptL = new TH1D(TString(hmt->GetName())+"_pt",
		      TString(hmt->GetName())+"_pt",
		      nbins, xbins);
  for(Int_t ibins = 1; ibins <= nbins; ibins++){
    hptL->SetBinContent(ibins, hmt->GetBinContent(ibins));
    hptL->SetBinError(ibins,   hmt->GetBinError(ibins));

  }

  hptL->SetXTitle("p_{t} (GeV/c)");
  hptL->SetYTitle("1/p_{t} dN/dp_{t} (a.u.)");
  hptL->SetMarkerStyle(hmt->GetMarkerStyle());
  hptL->SetMarkerColor(hmt->GetMarkerColor());
  hptL->SetLineColor(hmt->GetLineColor());

  return hptL;

}


TH1 * AliBWTools::GetdNdPtFromOneOverPt(const TH1 * h1Pt) {

  // convert an histo from 1/pt dNdpt to dNdpt, using the pt at the center of the bin


  TH1 * hPt = (TH1 *) h1Pt->Clone((TString(h1Pt->GetName()) + "_inv").Data());
  hPt->Reset();

  Int_t nbinx = hPt->GetNbinsX();

  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){

    Double_t cont = h1Pt->GetBinContent(ibinx);
    Double_t err  = h1Pt->GetBinError(ibinx);
    
    Double_t pt   = h1Pt->GetBinCenter(ibinx);
    
    if(pt > 0) {
      cont *= pt;
      err  *= pt;
    } else {
      cont = 0;
      err  = 0;
    }

    hPt->SetBinContent(ibinx, cont);
    hPt->SetBinError(ibinx, err);
    
  }

  hPt->SetXTitle("p_{t} (GeV)");
  hPt->SetYTitle("dN/dp_{t} (GeV^{-2})");

  return hPt;    

}




TH1 * AliBWTools::GetOneOverPtdNdPt(const TH1 * hPt) {

  // convert an histo from dNdpt to 1/pt dNdpt, using the pt at the center of the bin

  TH1 * h1Pt = (TH1 *) hPt->Clone((TString(hPt->GetName()) + "_inv").Data());
  h1Pt->Reset();

  Int_t nbinx = h1Pt->GetNbinsX();

  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){

    Double_t cont = hPt->GetBinContent(ibinx);
    Double_t err  = hPt->GetBinError(ibinx);
    
    Double_t pt   = hPt->GetBinCenter(ibinx);
    
    if(pt > 0) {
      cont /= pt;
      err  /= pt;
    } else {
      cont = 0;
      err  = 0;
    }

    h1Pt->SetBinContent(ibinx, cont);
    h1Pt->SetBinError(ibinx, err);
    
  }

  h1Pt->SetXTitle("p_{t} (GeV)");
  h1Pt->SetYTitle("1/p_{t} dN/dp_{t} (GeV^{-2})");

  return h1Pt;    

}


TGraphErrors * AliBWTools::GetGraphFromHisto(const TH1F * h, Bool_t binWidth) {
  // Convert a histo to a graph
  // if binWidth is true ex is set to the bin width of the histos, otherwise it is set to zero
  Int_t nbin = h->GetNbinsX();

  TGraphErrors * g = new TGraphErrors();
  Int_t ipoint = 0;
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    Double_t xerr = binWidth ? h->GetBinWidth(ibin)/2 : 0;
    if (h->GetBinContent(ibin)) {
      g->SetPoint     (ipoint,   h->GetBinCenter(ibin), h->GetBinContent(ibin));
      g->SetPointError(ipoint,   xerr,  h->GetBinError(ibin));
      ipoint++;
    }
  }
  
  g->SetMarkerStyle(h->GetMarkerStyle());
  g->SetMarkerColor(h->GetMarkerColor());
  g->SetLineColor(h->GetLineColor());
  g->SetLineStyle(h->GetLineStyle());
  g->SetLineWidth(h->GetLineWidth());

  g->SetTitle(h->GetTitle());
  g->SetName(TString("g_")+h->GetName());

  return g;

}

TH1F * AliBWTools::GetHistoFromGraph(const TGraphErrors * g, const TH1F* hTemplate) {

  // convert a graph to histo with the binning of hTemplate.
  // Warning: the template should be chosen
  // properly: if you have two graph points in the same histo bin this
  // won't work!

  TH1F * h = (TH1F*) hTemplate->Clone(TString("h_")+g->GetName());
  h->Reset();
  Int_t npoint = g->GetN();
  //  g->Print();
  for(Int_t ipoint = 0; ipoint < npoint; ipoint++){
    Float_t x  = g->GetX() [ipoint];
    Float_t y  = g->GetY() [ipoint];
    Float_t ey = g->GetEY()[ipoint];
    Int_t bin = h->FindBin(x);
    //    cout << "bin: "<< bin << " -> " << x << ", "<< y <<", " << ey << endl;
    
    h->SetBinContent(bin,y);
    h->SetBinError  (bin,ey);
  }
 
  h->SetMarkerStyle(g->GetMarkerStyle());
  h->SetMarkerColor(g->GetMarkerColor());
  h->SetLineColor  (g->GetLineColor());

 
  return h;
}

TGraphErrors * AliBWTools::ConcatenateGraphs(const TGraphErrors * g1,const TGraphErrors * g2){

  // concatenates two graphs

  Int_t npoint1=g1->GetN();
  Int_t npoint2=g2->GetN();

  TGraphErrors * gClone = (TGraphErrors*) g1->Clone();

//   for(Int_t ipoint = 0; ipoint < npoint1; ipoint++){
//     gClone->SetPointError(ipoint,0,g1->GetEY()[ipoint]);

//   }
  for(Int_t ipoint = 0; ipoint < npoint2; ipoint++){
    gClone->SetPoint(ipoint+npoint1,g2->GetX()[ipoint],g2->GetY()[ipoint]);
    gClone->SetPointError(ipoint+npoint1,g2->GetEX()[ipoint],g2->GetEY()[ipoint]);
    //    gClone->SetPointError(ipoint+npoint1,0,g2->GetEY()[ipoint]);
  }

  gClone->GetHistogram()->GetXaxis()->SetTimeDisplay(1);
  gClone->SetTitle(TString(gClone->GetTitle())+" + "+g2->GetTitle());
  gClone->SetName(TString(gClone->GetName())+"_"+g2->GetName());

  return gClone;
}


TH1F * AliBWTools::Combine3HistosWithErrors(const TH1 * h1,  const TH1 * h2,  const TH1* h3, 
					    TH1 * he1,  TH1 * he2,  TH1 * he3, 
					    const TH1* htemplate, Int_t statFrom, 
					    Float_t renorm1, Float_t renorm2, Float_t renorm3,
					    TH1 ** hSyst, Bool_t errorFromBinContent) {

  // Combines 3 histos (h1,h2,h3), weighting by the errors provided in
  // he1,he2,he3, supposed to be the independent systematic errors.
  // he1,he2,he3 are also assumed to have the same binning as h1,h2,h3
  // The combined histo must fit the template provided (no check is performed on this)
  // The histogram are supposed to come from the same (nearly) sample
  // of tracks. statFrom tells the stat error of which of the 3 is
  // going to be assigned to the combined
  // Optionally, it is possible to rescale any of the histograms.
  // if hSyst is give, the histo is filled with combined syst error vs pt
  // if errorFromBinContent is true, the weights are taken from the he* content rather than errors
  TH1F * hcomb = (TH1F*) htemplate->Clone(TString("hComb_")+h1->GetName()+"_"+h2->GetName()+h3->GetName());

  // TODO: I should have used an array for h*local...

  // Clone histos locally to rescale them
  TH1F * h1local = (TH1F*) h1->Clone();
  TH1F * h2local = (TH1F*) h2->Clone();
  TH1F * h3local = (TH1F*) h3->Clone();
  h1local->Scale(renorm1);
  h2local->Scale(renorm2);
  h3local->Scale(renorm3);

  const TH1 * hStatError = 0;
  if (statFrom == 0)      hStatError = h1; 
  else if (statFrom == 1) hStatError = h2; 
  else if (statFrom == 2) hStatError = h3; 
  else {
    Printf("AliBWTools::Combine3HistosWithErrors: wrong value of the statFrom parameter");
    return NULL;
  }
  Printf("AliBWTools::Combine3HistosWithErrors: improve error on combined");
  // Loop over all bins and take weighted mean of all points
  Int_t nBinComb = hcomb->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nBinComb; ibin++){
    Int_t ibin1 = h1local->FindBin(hcomb->GetBinCenter(ibin));
    Int_t ibin2 = h2local->FindBin(hcomb->GetBinCenter(ibin));
    Int_t ibin3 = h3local->FindBin(hcomb->GetBinCenter(ibin));
    Int_t ibinError = -1; // bin used to get stat error

    if (statFrom == 0)      ibinError = ibin1; 
    else if (statFrom == 1) ibinError = ibin2; 
    else if (statFrom == 2) ibinError = ibin3; 
    else Printf("AliBWTools::Combine3HistosWithErrors: wrong value of the statFrom parameter");


    Double_t y[3];
    Double_t ye[3];
    y[0]  = h1local->GetBinContent(ibin1);
    y[1]  = h2local->GetBinContent(ibin2);
    y[2]  = h3local->GetBinContent(ibin3);
    if (errorFromBinContent) {
      ye[0] = he1->GetBinContent(he1->FindBin(hcomb->GetBinCenter(ibin)));
      ye[1] = he2->GetBinContent(he2->FindBin(hcomb->GetBinCenter(ibin)));
      ye[2] = he3->GetBinContent(he3->FindBin(hcomb->GetBinCenter(ibin)));
    } else {
      ye[0] = he1->GetBinError(ibin1);
      ye[1] = he2->GetBinError(ibin2);
      ye[2] = he3->GetBinError(ibin3);
    }
    // Set error to 0 if content is 0 (means it was not filled)
    if(!h1local->GetBinContent(ibin1)) ye[0] = 0;
    if(!h2local->GetBinContent(ibin2)) ye[1] = 0;
    if(!h3local->GetBinContent(ibin3)) ye[2] = 0;
    
    // Compute weighted mean
    //    cout << "Bins:  "<< ibin1 << " " << ibin2 << " " << ibin3 << endl;    
    Double_t mean, err;
    WeightedMean(3,y,ye,mean,err);


    // Fill combined
    hcomb->SetBinContent(ibin,mean);
    Double_t statError = 0;
    if (hStatError->GetBinContent(ibinError)) {
      //      cout << "def" << endl;
      statError = hStatError->GetBinError(ibinError)/hStatError->GetBinContent(ibinError)*hcomb->GetBinContent(ibin);
    }
    else if (h1local->GetBinContent(ibin1)) {
      //      cout << "1" << endl;
      statError = h1local->GetBinError(ibin1)/h1local->GetBinContent(ibin1)*hcomb->GetBinContent(ibin);
    }
    else if (h2local->GetBinContent(ibin2)) {
      //      cout << "2" << endl;
      statError = h2local->GetBinError(ibin2)/h2local->GetBinContent(ibin2)*hcomb->GetBinContent(ibin);
    }
    else if (h3local->GetBinContent(ibin3)) {
      //      cout << "3" << endl;
      statError = h3local->GetBinError(ibin3)/h3local->GetBinContent(ibin3)*hcomb->GetBinContent(ibin);
    }
    hcomb->SetBinError  (ibin,statError);
    if(hSyst) (*hSyst)->SetBinContent(ibin,err);
    //    cout << "BIN " << ibin << " " << mean << " " << statError << endl;

  }

  hcomb->SetMarkerStyle(hStatError->GetMarkerStyle());
  hcomb->SetMarkerColor(hStatError->GetMarkerColor());
  hcomb->SetLineColor  (hStatError->GetLineColor());

  hcomb->SetXTitle(hStatError->GetXaxis()->GetTitle());
  hcomb->SetYTitle(hStatError->GetYaxis()->GetTitle());

  delete h1local;
  delete h2local;
  delete h3local;

  return hcomb;
}


void AliBWTools::GetMeanDataAndExtrapolation(const TH1 * hData, TF1 * fExtrapolation, Double_t &mean, Double_t &error, Float_t min, Float_t max){
  // Computes the mean of the combined data and extrapolation in a
  // given range, use data where they are available and the function
  // where data are not available
  // ERROR from DATA ONLY is returned in this version! 
  //
  Printf("AliBWTools::GetMeanDataAndExtrapolation: WARNING from data only");
  Float_t minData    = GetLowestNotEmptyBinEdge (hData);
  Int_t minDataBin   = GetLowestNotEmptyBin     (hData);
  Float_t maxData    = GetHighestNotEmptyBinEdge(hData);
  Int_t maxDataBin   = GetHighestNotEmptyBin    (hData);
  Double_t integral  = 0; 
  mean      = 0;
  error     = 0; 
  if (min < minData) {
    // Compute integral exploiting root function to calculate moments, "unnormalizing" them
    mean     += fExtrapolation->Mean(min,minData)*fExtrapolation->Integral(min,minData);
    integral += fExtrapolation->Integral(min,minData);
    cout << " Low "<< mean << " " << integral << endl;
    
  }
  
  if (max > maxData) {
    // Compute integral exploiting root function to calculate moments, "unnormalizing" them
    mean     += fExtrapolation->Mean(maxData,max)*fExtrapolation->Integral(maxData,max);
    integral += fExtrapolation->Integral(maxData,max);
    cout << " Hi "<< mean << " " << integral << endl;
  } 
  Float_t err2 = 0;
  
  for(Int_t ibin = minDataBin; ibin <= maxDataBin; ibin++){
    if(hData->GetBinCenter(ibin) < min) continue;
    if(hData->GetBinCenter(ibin) > max) continue;
    mean     = mean + (hData->GetBinCenter(ibin) *  hData->GetBinWidth(ibin)* hData->GetBinContent(ibin));
    err2     = err2 + TMath::Power(hData->GetBinError(ibin) * hData->GetBinCenter(ibin) *  hData->GetBinWidth(ibin),2);
    integral = integral + hData->GetBinContent(ibin) * hData->GetBinWidth(ibin);
  }
  cout << " Data "<< mean << " " << integral << endl;
  
  mean = mean/integral;
  error = TMath::Sqrt(err2)/integral;


}

TH1F * AliBWTools::CombineHistos(const TH1 * h1, TH1 * h2, const TH1* htemplate, Float_t renorm1){
  // Combine two histos. This assumes the histos have the same binning
  // in the overlapping region. It computes the arithmetic mean in the
  // overlapping region and assigns as an error the relative error
  // h2. TO BE IMPROVED

  TH1F * hcomb = (TH1F*) htemplate->Clone(TString(h1->GetName())+"_"+h2->GetName());

  TH1F * h1local = (TH1F*) h1->Clone();
  h1local->Scale(renorm1);
  
  Int_t nBinComb = hcomb->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nBinComb; ibin++){
    Int_t ibin1 = h1local->FindBin(hcomb->GetBinCenter(ibin));
    Int_t ibin2 = h2->FindBin(hcomb->GetBinCenter(ibin));
    
      if (!h1local->GetBinContent(ibin1) && !h2->GetBinContent(ibin2) ) {
	// None has data: go to next bin
	hcomb->SetBinContent(ibin,0);
	hcomb->SetBinError  (ibin,0);	
      } else if(h1local->GetBinContent(ibin1) && !h2->GetBinContent(ibin2)) {
	// take data from h1local:
	hcomb->SetBinContent(ibin,h1local->GetBinContent(ibin1));
	hcomb->SetBinError  (ibin,h1local->GetBinError(ibin1));
      } else if(!h1local->GetBinContent(ibin1) && h2->GetBinContent(ibin2)) {
	// take data from h2:
	hcomb->SetBinContent(ibin,h2->GetBinContent(ibin2));
	hcomb->SetBinError  (ibin,h2->GetBinError(ibin2));
      }  else {
	hcomb->SetBinContent(ibin,(h1local->GetBinContent(ibin1) +h2->GetBinContent(ibin2))/2);
	//	hcomb->SetBinError  (ibin,h1local->GetBinError(ibin1)/h1local->GetBinContent(ibin1)*hcomb->GetBinContent(ibin));
	hcomb->SetBinError  (ibin,h2->GetBinError(ibin2)/h2->GetBinContent(ibin2)*hcomb->GetBinContent(ibin));
      }


  }
  

  hcomb->SetMarkerStyle(h1local->GetMarkerStyle());
  hcomb->SetMarkerColor(h1local->GetMarkerColor());
  hcomb->SetLineColor  (h1local->GetLineColor());

  hcomb->SetXTitle(h1local->GetXaxis()->GetTitle());
  hcomb->SetYTitle(h1local->GetYaxis()->GetTitle());
  delete h1local;
  return hcomb;  

}


void AliBWTools::GetFromHistoGraphDifferentX(const TH1F * h, TF1 * f, TGraphErrors ** gBarycentre, TGraphErrors ** gXlw) {

  // Computes the baycentre in each bin and the xlw as defined in NIMA
  // 355 - 541 using f. Returs 2 graphs with the same y content of h
  // but with a different x (barycentre and xlw)

  Int_t nbin = h->GetNbinsX();
  
  (*gBarycentre) = new TGraphErrors();
  (*gXlw)        = new TGraphErrors();

  Int_t ipoint = 0;
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    Float_t min = h->GetBinLowEdge(ibin);
    Float_t max = h->GetBinLowEdge(ibin+1);
    Double_t xerr = 0;
    Double_t xbar = f->Mean(min,max);
    // compute xLW
    Double_t temp = 1./(max-min) * f->Integral(min,max);
    Double_t epsilon   = 0.000000001;
    Double_t increment = 0.0000000001;
    Double_t xLW = min;

    while ((f->Eval(xLW)- temp) > epsilon) {
      xLW += increment;
      if(xLW > max) {
	Printf("Cannot find xLW");
	break;
      }
    }
      
    if (h->GetBinContent(ibin)!=0) {
      (*gBarycentre)->SetPoint     (ipoint,   xbar, h->GetBinContent(ibin));
      (*gBarycentre)->SetPointError(ipoint,   xerr, h->GetBinError(ibin));
      (*gXlw)       ->SetPoint     (ipoint,   xLW,  h->GetBinContent(ibin));
      (*gXlw)       ->SetPointError(ipoint,   xerr, h->GetBinError(ibin));
      ipoint++;
    }
  }
  
  (*gBarycentre)->SetMarkerStyle(h->GetMarkerStyle());
  (*gBarycentre)->SetMarkerColor(h->GetMarkerColor());
  (*gBarycentre)->SetLineColor(h->GetLineColor());

  (*gBarycentre)->SetTitle(h->GetTitle());
  (*gBarycentre)->SetName(TString("g_")+h->GetName());

  (*gXlw)->SetMarkerStyle(h->GetMarkerStyle());
  (*gXlw)->SetMarkerColor(h->GetMarkerColor());
  (*gXlw)->SetLineColor(h->GetLineColor());
  (*gXlw)->SetTitle(h->GetTitle());
  (*gXlw)->SetName(TString("g_")+h->GetName());


}


Float_t AliBWTools::GetMean(TH1F * h, Float_t min, Float_t max, Float_t * error) {

  // Get the mean of histo in a range; root is not reliable in sub
  // ranges with variable binning.  
  Int_t minBin = h->FindBin(min);
  Int_t maxBin = h->FindBin(max-0.00001);

  Float_t mean = 0 ;
  Float_t integral = 0;
  Float_t err2 = 0;
  for(Int_t ibin = minBin; ibin <= maxBin; ibin++){
    mean     = mean + (h->GetBinCenter(ibin) *  h->GetBinWidth(ibin)* h->GetBinContent(ibin));
    err2     = err2 + TMath::Power(h->GetBinError(ibin) * h->GetBinCenter(ibin) *  h->GetBinWidth(ibin),2);
    integral = integral + h->GetBinContent(ibin) * h->GetBinWidth(ibin);
  }
  
  Float_t value = mean/integral;  
  if (error) (*error) = TMath::Sqrt(err2);
  return value;


}

void AliBWTools::GetMean(TF1 * func, Float_t &mean, Float_t &error, Float_t min, Float_t max, Int_t normPar) {

  // Get the mean of function in a range; If normPar is >= 0, it means
  // that the function is defined such that par[normPar] is its
  // integral.  In this case the error on meanpt can be calculated
  // correctly. Otherwise, the function is normalized in get moment,
  // but the error is not computed correctly.

  return GetMoment("fMean", TString("x*")+func->GetExpFormula(), func, mean, error, min, max, normPar);

}

void AliBWTools::GetMeanSquare(TF1 * func, Float_t &mean, Float_t &error, Float_t min, Float_t max, Int_t normPar) {

  // Get the mean2 of function in a range;  If normPar is >= 0, it means
  // that the function is defined such that par[normPar] is its
  // integral.  In this case the error on meanpt can be calculated
  // correctly. Otherwise, the function is normalized in get moment,
  // but the error is not computed correctly.

  return GetMoment("fMean2", TString("x*x*")+func->GetExpFormula(), func, mean, error, min, max, normPar);


}

void AliBWTools::GetMoment(TString name, TString var, TF1 * func, Float_t &mean, Float_t &error, Float_t min, Float_t max, Int_t normPar) {

  // returns the integral of a function derived from func by prefixing some operation.
  // the derived function MUST have the same parameter in the same order
  // Used as a base method for mean and mean2
  //  If normPar is >= 0, it means that the function is defined such
  // that par[normPar] is its integral.  In this case the error on
  // meanpt can be calculated correctly. Otherwise, the function is
  // normalized using its numerical integral, but the error is not computed
  // correctly. 

  // TODO:
  // - improve to propagate error also in the case you need the
  //   numerical integrals (will be rather slow)
  // - this version assumes that func is defined using a
  //   TFormula. Generalize to the case of a C++ function

  if (normPar<0) Printf("AliBWTools::GetMoment: Warning: If normalization is required, the error may bot be correct");
  if (!strcmp(func->GetExpFormula(),"")) Printf("AliBWTools::GetMoment: Warning: Empty formula in the base function");
  Int_t npar = func->GetNpar();

  // Definition changes according to the value of normPar
  TF1 * f = normPar < 0 ? 
    new TF1(name, var) :	              // not normalized
    new TF1(name, var+Form("/[%d]",normPar)); // normalized with par normPar

  // integr is used to normalize if no parameter is provided
  Double_t integr  = normPar < 0 ? func->Integral(min,max) : 1;
  
  // The parameter of the function used to compute the mean should be
  // the same as the parent function: fixed if needed and they should
  // also have the same errors.

  //  cout << "npar :" << npar << endl;
  
  for(Int_t ipar = 0; ipar < npar; ipar++){
    Double_t parmin, parmax;
    Double_t value = func->GetParameter(ipar);
    f->SetParameter(ipar,value);
    func->GetParLimits(ipar, parmin, parmax);
    if ( parmin == parmax )   {
      //      if ( parmin || (parmin == 1 && !value) ) { // not sure we I check parmin == 1 here. 
      if ( parmin || (TMath::Abs(parmin-1)<0.000001 && !value) ) { // not sure we I check parmin == 1 here. Changed like this because of coding conventions. Does it still work? FIXME
	f->FixParameter(ipar,func->GetParameter(ipar));
	//	cout << "Fixing " << ipar << "("<<value<<","<<parmin<<","<<parmax<<")"<<endl;
      }       
      else {
	f->SetParError (ipar,func->GetParError(ipar) );
	//	cout << "Setting Err" << ipar << "("<<func->GetParError(ipar)<<")"<<endl;      
      }
    }
    else {
      f->SetParError (ipar,func->GetParError(ipar) );
      //      cout << "Setting Err" << ipar << "("<<func->GetParError(ipar)<<")"<<endl;      
    }
  }

//   f->Print();
//   cout << "----" << endl;
//   func->Print();

  mean  = normPar < 0 ? f->Integral     (min,max)/integr : f->Integral     (min,max);
  error = normPar < 0 ? f->IntegralError(min,max)/integr : f->IntegralError(min,max);
//   cout << "Mean " << mean <<"+-"<< error<< endl;
//   cout << "Integral Error "  << error << endl;
  
}



Bool_t AliBWTools::Fit (TH1 * h1, TF1* func, Float_t min, Float_t max) { 
  
  // Fits h1 with func, preforms several checks on the quality of the
  // fit and tries to improve it. If the fit is not good enough, it
  // returs false.

  Double_t amin; Double_t edm; Double_t errdef; Int_t nvpar; Int_t nparx;
  TVirtualFitter *fitter;
  cout << "--- Fitting : " << h1->GetName() << " ["<< h1->GetTitle() <<"] ---"<< endl;

  h1-> Fit(func,"IME0","",min,max);      
  Int_t fitResult = h1-> Fit(func,"IME0","",min,max);      
//   h1-> Fit(func,"0","",min,max);      
//   Int_t fitResult = h1-> Fit(func,"0IE","",min,max);      
  

// From TH1:
// The fitStatus is 0 if the fit is OK (i.e no error occurred).  The
// value of the fit status code is negative in case of an error not
// connected with the minimization procedure, for example when a wrong
// function is used.  Otherwise the return value is the one returned
// from the minimization procedure.  When TMinuit (default case) or
// Minuit2 are used as minimizer the status returned is : fitStatus =
// migradResult + 10*minosResult + 100*hesseResult +
// 1000*improveResult.  TMinuit will return 0 (for migrad, minos,
// hesse or improve) in case of success and 4 in case of error (see
// the documentation of TMinuit::mnexcm). So for example, for an error
// only in Minos but not in Migrad a fitStatus of 40 will be returned.
// Minuit2 will return also 0 in case of success and different values
// in migrad minos or hesse depending on the error. See in this case
// the documentation of Minuit2Minimizer::Minimize for the
// migradResult, Minuit2Minimizer::GetMinosError for the minosResult
// and Minuit2Minimizer::Hesse for the hesseResult.  If other
// minimizers are used see their specific documentation for the status
// code returned.  For example in the case of Fumili, for the status
// returned see TFumili::Minimize.
 

  if( gMinuit->fLimset ) {
    Printf("ERROR: AliBWTools: Parameters at limits");
    return kFALSE;
  } 


  ///
  fitter = TVirtualFitter::GetFitter();   
  Int_t  fitStat = fitter->GetStats(amin, edm, errdef, nvpar, nparx);  

  if( ( (fitStat < 3  && gMinuit->fCstatu != "UNCHANGED ")|| (edm > 1e6) || (fitResult !=0 && fitResult < 4000) ) && 
      TString(gMinuit->fCstatu) != "SUCCESSFUL"  &&
      TString(gMinuit->fCstatu) != "CONVERGED "  ) {
    if(fitStat < 3 && gMinuit->fCstatu != "UNCHANGED ") {
      Printf("WARNING: AliBWTools: Cannot properly compute errors");
      if (fitStat == 0) Printf(" not calculated at all");
      if (fitStat == 1) Printf(" approximation only, not accurate");
      if (fitStat == 2) Printf(" full matrix, but forced positive-definite");
    }

    if (edm > 1e6) {

      Printf("WARNING: AliBWTools: Huge EDM  (%f): Fit probably failed!", edm);
    }
    if (fitResult != 0) {
      Printf("WARNING: AliBWTools: Fit Result (%d)", fitResult);
    }
      
    Printf ("AliBWTools: Trying Again with Strategy = 2");

    gMinuit->Command("SET STRATEGY 2"); // more effort
    fitResult = 0;
    fitResult = h1-> Fit(func,"0","",min,max);      
    fitResult = h1-> Fit(func,"IME0","",min,max);      
    fitResult = h1-> Fit(func,"IME0","",min,max);      
      
    fitter = TVirtualFitter::GetFitter();   
  
    fitStat = fitter->GetStats(amin, edm, errdef, nvpar, nparx);  

    if(fitStat < 3 && gMinuit->fCstatu != "UNCHANGED ") {
      Printf("ERROR: AliBWTools: Cannot properly compute errors");
      if (fitStat == 0) Printf(" not calculated at all");
      if (fitStat == 1) Printf(" approximation only, not accurate");
      if (fitStat == 2) Printf(" full matrix, but forced positive-definite");
      cout << "[" <<gMinuit->fCstatu<<"]" << endl;
      return kFALSE;
    }

    if (edm > 1e6) {
      Printf("ERROR: AliBWTools: Huge EDM  (%f): Fit probably failed!", edm);

      return kFALSE;
    }
    if (fitResult != 0) {
      Printf("ERROR: AliBWTools: Fit Result (%d)", fitResult);

      return kFALSE;
    }
      
    gMinuit->Command("SET STRATEGY 1"); // back to normal value

  }

  cout << "---- FIT OK ----" << endl;
  
  return kTRUE;
  
}

Int_t AliBWTools::GetLowestNotEmptyBin(const TH1*h) {

  // Return the index of the lowest non empty bin in the histo h

  Int_t nbin = h->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    if(h->GetBinContent(ibin)>0) return ibin;
  }
  
  return -1;

}

Int_t AliBWTools::GetHighestNotEmptyBin(const TH1*h) {

  // Return the index of the highest non empty bin in the histo h

  Int_t nbin = h->GetNbinsX();
  for(Int_t ibin = nbin; ibin > 0; ibin--){
    if(h->GetBinContent(ibin)>0) return ibin;
  }
  
  return -1;

}

void AliBWTools::GetResiduals(const TGraphErrors * gdata, const TF1 * func, TH1F ** hres, TGraphErrors ** gres) {

  // Returns a graph of residuals vs point and the res/err distribution

  Int_t npoint = gdata->GetN();

  (*gres) =new TGraphErrors();
  (*hres) = new TH1F(TString("hres_")+gdata->GetName()+"-"+func->GetName(),
                  TString("hres_")+gdata->GetName()+"-"+func->GetName(),
                  20,-5,5);


  for(Int_t ipoint = 0; ipoint < npoint; ipoint++){
    Float_t x   = gdata->GetX()[ipoint];
    Float_t res = (gdata->GetY()[ipoint] - func->Eval(x))/func->Eval(x);
    Float_t err = gdata->GetEY()[ipoint]/func->Eval(x);
    (*hres)->Fill(res/err);
    (*gres)->SetPoint(ipoint, x, res/err);
    //    (*gres)->SetPointError(ipoint, 0, err);
    
  }
  
  (*gres)->SetMarkerStyle(gdata->GetMarkerStyle());
  (*gres)->SetMarkerColor(gdata->GetMarkerColor());
  (*gres)->SetLineColor  (gdata->GetLineColor());
  (*gres)->GetHistogram()->GetYaxis()->SetTitle("(data-function)/function");
  (*hres)->SetMarkerStyle(gdata->GetMarkerStyle());
  (*hres)->SetMarkerColor(gdata->GetMarkerColor());
  (*hres)->SetLineColor  (gdata->GetLineColor());



}

void AliBWTools::GetResiduals(const TH1F* hdata, const TF1 * func, TH1F ** hres, TH1F ** hresVsBin) {

  // Returns an histo of residuals bin by bin and the res/err distribution

  if (!func) {
    Printf("AliBWTools::GetResiduals: No function provided");
    return;
  }
  if (!hdata) {
    Printf("AliBWTools::GetResiduals: No data provided");
    return;
  }

  (*hresVsBin) = (TH1F*) hdata->Clone(TString("hres_")+hdata->GetName());
  (*hresVsBin)->Reset();
  (*hres) = new TH1F(TString("hres_")+hdata->GetName()+"-"+func->GetName(),
		     TString("hres_")+hdata->GetName()+"-"+func->GetName(),
		     20,-5,5);

  Int_t nbin = hdata->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    if(!hdata->GetBinContent(ibin)) continue;
    Float_t res = (hdata->GetBinContent(ibin) - func->Eval(hdata->GetBinCenter(ibin)) ) / 
      func->Eval(hdata->GetBinCenter(ibin));
    Float_t err = hdata->GetBinError  (ibin) /  func->Eval(hdata->GetBinCenter(ibin));
    (*hresVsBin)->SetBinContent(ibin,res);
    (*hresVsBin)->SetBinError  (ibin,err);
    (*hres)->Fill(res/err);
    
  }
  
  (*hresVsBin)->SetMarkerStyle(hdata->GetMarkerStyle());
  (*hresVsBin)->SetMarkerColor(hdata->GetMarkerColor());
  (*hresVsBin)->SetLineColor  (hdata->GetLineColor()  );
  (*hresVsBin)->GetYaxis()->SetTitle("(data-function)/function");
  (*hres)->SetMarkerStyle(hdata->GetMarkerStyle());
  (*hres)->SetMarkerColor(hdata->GetMarkerColor());
  (*hres)->SetLineColor  (hdata->GetLineColor()  );

}

void AliBWTools::GetYield(TH1* h,  TF1 * f, Double_t &yield, Double_t &yieldError, Float_t min, Float_t max,
			  Double_t *partialYields, Double_t *partialYieldsErrors){

  // Returns the yield extracted from the data in the histo where
  // there are points and from the fit to extrapolate, in the given
  // range.

  // Partial yields are also returned if the corresponding pointers are non null

  Int_t bin1 = h->FindBin(min);
  Int_t bin2 = h->FindBin(max);
  Float_t bin1Edge = GetLowestNotEmptyBinEdge (h);
  Float_t bin2Edge = GetHighestNotEmptyBinEdge(h);

  Double_t integralFromHistoError ;
  Double_t integralFromHisto = DoIntegral(h,bin1,bin2,-1,-1,-1,-1,integralFromHistoError,"width",1);
  
  Double_t integralBelow      = min < bin1Edge ? f->Integral(min,bin1Edge) : 0;
  Double_t integralBelowError = min < bin1Edge ? f->IntegralError(min,bin1Edge) : 0;
  Double_t integralAbove      = max > bin2Edge ? f->Integral(bin2Edge,max) : 0;
  Double_t integralAboveError = max > bin2Edge ? f->IntegralError(bin2Edge,max) : 0;

//   cout << "GetYield INFO" << endl;
//   cout << " " << bin1Edge << " " << bin2Edge << endl;  
//   cout << " " << integralFromHisto      << " " << integralBelow      << " " << integralAbove      << endl;
//   cout << " " << integralFromHistoError << " " << integralBelowError << " " << integralAboveError << endl;
  
  if(partialYields) {
    partialYields[0] = integralFromHisto;
    partialYields[1] = integralBelow;
    partialYields[2] = integralAbove;
  }
  if(partialYieldsErrors) {
    partialYieldsErrors[0] = integralFromHistoError;
    partialYieldsErrors[1] = integralBelowError;
    partialYieldsErrors[2] = integralAboveError;
  }
  yield      = integralFromHisto+integralBelow+integralAbove;
  yieldError = TMath::Sqrt(integralFromHistoError*integralFromHistoError+
			   integralBelowError*integralBelowError+
			   integralAboveError*integralAboveError);

}

TGraphErrors * AliBWTools::DivideGraphByFunc(const TGraphErrors * g, const TF1 * f, Bool_t invert){ 

  // Divides g/f. If invert == true => f/g

  TGraphErrors * gRatio = new TGraphErrors();
  Int_t npoint = g->GetN();
  for(Int_t ipoint = 0; ipoint < npoint; ipoint++){
    Double_t x = g->GetX()[ipoint];
    Double_t ratio  = invert ? f->Eval(x)/g->GetY()[ipoint] :g->GetY()[ipoint]/f->Eval(x);
    gRatio->SetPoint     (ipoint, x, ratio);
    if(f->Eval(x) && strcmp(g->ClassName(),"TGraphAsymmErrors")) gRatio->SetPointError(ipoint, 0, g->GetEY()[ipoint]/f->Eval(x));
    //    cout << x << " " << g->GetY()[ipoint] << " " << f->Eval(x) << endl;
    
  }
  gRatio->SetMarkerStyle(20);
  //gRatio->Print();
  return gRatio;

}

TGraphErrors * AliBWTools::Divide2Graphs(const TGraphErrors * g1, const TGraphErrors * g2){ 

  // Divides g1/g2, looks for point with very close centers
  Int_t ipoint=0;
  TGraphErrors * gRatio = new TGraphErrors();
  Int_t npoint1 = g1->GetN();
  Int_t npoint2 = g2->GetN();
  for(Int_t ipoint1 = 0; ipoint1 < npoint1; ipoint1++){
    Double_t x1 = g1->GetX()[ipoint1];
    for(Int_t ipoint2 = 0; ipoint2 < npoint2; ipoint2++){
      Double_t x2 = g2->GetX()[ipoint2];
      if((TMath::Abs(x1-x2)/(x1+x2)*2)<0.01) {
	Double_t ratio   = g2->GetY()[ipoint2]  ? g1->GetY()[ipoint1]/g2->GetY()[ipoint2] : 0;
	Double_t eratio  = g2->GetY()[ipoint2]  ? 
	  TMath::Sqrt(g1->GetEY()[ipoint1]*g1->GetEY()[ipoint1]/g1->GetY()[ipoint1]/g1->GetY()[ipoint1] + 
		      g2->GetEY()[ipoint2]/g2->GetY()[ipoint2]/g2->GetY()[ipoint2] ) * ratio
	  : 0;
	gRatio->SetPoint     (ipoint, x1, ratio);
	gRatio->SetPointError(ipoint, 0, eratio);
	ipoint++;
	cout << ipoint << " [" << x1 << "] " <<  g1->GetY()[ipoint1] << "/" << g2->GetY()[ipoint2] << " = " << ratio <<"+-"<<eratio<< endl;
	
    //    cout << x << " " << g->GetY()[ipoint] << " " << f->Eval(x) << endl;
      }
    
    }
  }
  gRatio->SetMarkerStyle(20);
  //gRatio->Print();
  return gRatio;

}

TGraphErrors * AliBWTools::DivideGraphByHisto(const TGraphErrors * g, TH1 * h, Bool_t invert){ 

  // Divides g/h. If invert == true => h/g

  Bool_t skipError = kFALSE;
  if(!strcmp(g->ClassName(),"TGraph")) skipError = kTRUE;
  if(!strcmp(g->ClassName(),"TGraphAsymmErrors")) skipError = kTRUE;
  if(skipError) {
    Printf("WARNING: Skipping graph errors");
  }
  TGraphErrors * gRatio = new TGraphErrors();
  Int_t npoint = g->GetN();
  for(Int_t ipoint = 0; ipoint < npoint; ipoint++){
    Double_t xj  = g->GetX()[ipoint];
    Double_t yj  = g->GetY()[ipoint];
    Double_t yje = skipError ? 0 : g->GetEY()[ipoint];

    Int_t binData = h->FindBin(xj);
    Double_t yd   = h->GetBinContent(binData);
    Double_t yde  = h->GetBinError(binData);
    Double_t xd   = h->GetBinCenter(binData);
    

     
    if (!yd) continue;
    
    if (TMath::Abs(xd-xj)/TMath::Abs(xd) > 0.01){
      Printf( "WARNING: bin center (%f)  and x graph (%f) are more than 1 %% away, skipping",xd,xj );
      continue;
      
    }

    Double_t ratio = invert ? yd/yj : yj/yd;

    gRatio->SetPoint(ipoint, xj, ratio);
    gRatio->SetPointError(ipoint, 0, TMath::Sqrt(yde*yde/yd/yd + yje*yje/yj/yj)*ratio);
    //    gRatio->SetPointError(ipoint, 0, yje/yj * ratio);
  }

  return gRatio;


}

TH1F * AliBWTools::DivideHistoByFunc(TH1F * h, TF1 * f, Bool_t invert){ 

  // Divides h/f. If invert == true => f/g
  // Performs the integral of f on the bin range to perform the ratio
  // Returns a histo with the same binnig as h

  // Prepare histo for ratio
  TH1F * hRatio = (TH1F*) h->Clone(TString("hRatio_")+h->GetName()+"_"+f->GetName());
  hRatio->Reset();
  // Set y title
  if(!invert) hRatio->SetYTitle(TString(h->GetName())+"/"+f->GetName());
  else        hRatio->SetYTitle(TString(f->GetName())+"/"+h->GetName());

  // Loop over all bins
  Int_t nbin = hRatio->GetNbinsX();

  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    Double_t yhisto = h->GetBinContent(ibin);
    Double_t yerror = h->GetBinError(ibin);
    Double_t xmin   = h->GetBinLowEdge(ibin);
    Double_t xmax   = h->GetBinLowEdge(ibin+1);
    Double_t yfunc  = f->Integral(xmin,xmax)/(xmax-xmin);
    Double_t ratio = invert ? yfunc/yhisto : yhisto/yfunc ;
    Double_t error = yerror/yfunc  ;
    hRatio->SetBinContent(ibin,ratio);
    hRatio->SetBinError  (ibin,error);
  }

  return hRatio;

}

void AliBWTools::WeightedMean(Int_t npoints, const Double_t *x, const Double_t *xerr, Double_t &mean, Double_t &meanerr){

  // Performs the weighted mean of npoints numbers in x with errors in xerr

  mean = 0;
  meanerr = 0;

  Double_t sumweight = 0;

  for (Int_t ipoint = 0; ipoint < npoints; ipoint++){
    
    Double_t xerr2 = xerr[ipoint]*xerr[ipoint];
    if(xerr2>0){
      //      cout << "xe2 " << xerr2 << endl;
      Double_t weight = 1. / xerr2;
      sumweight += weight;
      mean += weight * x[ipoint];
    }// else cout << " Skipping " << ipoint << endl;
    
  }


  if(sumweight){
    mean /= sumweight;
    meanerr = TMath::Sqrt(1./ sumweight);
  }
  else {
    //    cout << " No sumweight" << endl;
    mean = 0;
    meanerr = 0;
  }

  
}

TH1 * AliBWTools::GetRelativeError(TH1 * h){
  // Returns an histogram with the same binning as h, filled with the relative error bin by bin
  TH1 * hrel = (TH1*) h->Clone(TString(h->GetName())+"_rel");
  hrel->Reset();
  Int_t nbin = hrel->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    hrel->SetBinContent(ibin,h->GetBinError(ibin)/h->GetBinContent(ibin));
    hrel->SetBinError(ibin,0);
  }
  
  return hrel;
}


void AliBWTools::GetValueAndError(TH1 * hdest, const TH1 * hvalue, const TH1 * herror, Bool_t isPercentError) {
  
  // Put into source, bin-by-bin, the values from hvalue and the
  // errors from content from herror. 
  // Used mainly to combine histos of systemati errors with their spectra
  // Set isPercentError to kTRUE if the error is given in % 

  if(hdest == NULL){ 
    Printf("AliBWTools::GetValueAndError Errror: hdest is null");
    return;
  }


  Int_t nbin  = hdest->GetNbinsX();
  Int_t nBinSourceVal = hvalue->GetNbinsX();
  Int_t nBinSourceErr = herror->GetNbinsX();
  
  for(Int_t iBinDest = 1; iBinDest <= nbin; iBinDest++){
    Float_t lowPtDest=hdest->GetBinLowEdge(iBinDest);
    Float_t binWidDest=hdest->GetBinWidth(iBinDest);
    // Loop over Source bins and find overlapping bins to Dest
    // First value then error
    // Value
    Bool_t foundValue = kFALSE;
    for(Int_t iBinSourceVal=1; iBinSourceVal<=nBinSourceVal; iBinSourceVal++){
      Float_t lowPtSource=  hvalue->GetBinLowEdge(iBinSourceVal) ;
      Float_t binWidSource= hvalue->GetBinWidth(iBinSourceVal);
      if(TMath::Abs(lowPtDest-lowPtSource)<0.001 && TMath::Abs(binWidSource-binWidDest)<0.001){
	Double_t content = hvalue->GetBinContent(iBinSourceVal);
	hdest->SetBinContent(iBinDest, content);
	foundValue = kTRUE;
	break;
      }
    }
    // if (!foundValue){
    //   Printf("AliBWTools::GetValueAndError: Error: cannot find matching value source bin for destination %d",iBinDest);
    // }

    // Error
    Bool_t foundError = kFALSE;
    for(Int_t iBinSourceErr=1; iBinSourceErr<=nBinSourceErr; iBinSourceErr++){
      Float_t lowPtSource=  herror->GetBinLowEdge(iBinSourceErr) ;
      Float_t binWidSource= herror->GetBinWidth(iBinSourceErr);
      if(TMath::Abs(lowPtDest-lowPtSource)<0.001 && TMath::Abs(binWidSource-binWidDest)<0.001){
	Double_t error = herror->GetBinContent(iBinSourceErr);
	//	cout << "-> " << iBinDest << " " << error << " " << hdest->GetBinContent(iBinDest) << endl;
	
	hdest->SetBinError(iBinDest, isPercentError ? error * hdest->GetBinContent(iBinDest) : error);
	foundError=kTRUE;
	break;
      }      
    }
    // if (!foundError ){
    //   Printf("AliBWTools::GetValueAndError: Error: cannot find matching error source bin for destination %d",iBinDest);
    // }
  }
  

}

void AliBWTools::AddHisto(TH1 * hdest, const TH1* hsource, Bool_t getMirrorBins) {

  // Adds hsource to hdest bin by bin, even if they have a different
  // binning If getMirrorBins is true, it takes the negative bins
  // (Needed because sometimes the TPC uses the positive axis for
  // negative particles and the possitive axis for positive
  // particles).


  if (hdest == NULL) {
    Printf("Error: hdest is NULL\n");
    return;
  } 
  if (hsource == NULL) {
    Printf("Error: hsource is NULL\n");
    return;
  } 

  Int_t nBinSource = hsource->GetNbinsX();
  Int_t nBinDest = hdest->GetNbinsX();

  // Loop over destination bins, 
  for(Int_t iBinDest=1; iBinDest<=nBinDest; iBinDest++){
    Float_t lowPtDest=hdest->GetBinLowEdge(iBinDest);
    Float_t binWidDest=hdest->GetBinWidth(iBinDest);
    // Loop over Source bins and find overlapping bins to Dest
    Bool_t found = kFALSE;
    for(Int_t iBinSource=1; iBinSource<=nBinSource; iBinSource++){      
      Float_t lowPtSource= getMirrorBins ? -hsource->GetBinLowEdge(iBinSource)+hsource->GetBinWidth(iBinSource) : hsource->GetBinLowEdge(iBinSource) ;
      Float_t binWidSource= hsource->GetBinWidth(iBinSource)  ;
      if(TMath::Abs(lowPtDest-lowPtSource)<0.001 && TMath::Abs(binWidSource-binWidDest)<0.001){
	Float_t dest=hdest->GetBinContent(iBinDest);
	Float_t source=hsource->GetBinContent(iBinSource);
	Float_t edest=hdest->GetBinError(iBinDest);
	Float_t esource=hsource->GetBinError(iBinSource);
	Double_t cont=dest+source;
	Double_t econt=TMath::Sqrt(edest*edest+esource*esource);
	hdest->SetBinContent(iBinDest,cont);
	hdest->SetBinError  (iBinDest,econt);
	found = kTRUE;
	
	break;
      }
    }
    // if (!found){
    //   Printf("Error: cannot find matching source bin for destination %d",iBinDest);
    // }
  }


}

void AliBWTools::GetHistoCombinedErrors(TH1 * hdest, const TH1 * h1) {

  // Combine the errors of hdest with the errors of h1, summing in
  // quadrature. Results are put in hdest. Histograms are assumed to
  // have the same binning

  Int_t nbin = hdest->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    Double_t e1 = hdest->GetBinError(ibin);
    Double_t e2 = h1->GetBinError(ibin);
    hdest->SetBinError(ibin, TMath::Sqrt(e1*e1+e2*e2));
  }
  
  
}

TH1F * AliBWTools::DivideHistosDifferentBins(const TH1F* h1, const TH1F* h2) {
  // Divides 2 histos even if they have a different binning. Finds
  // overlapping bins and divides them
  
  // 1. clone histo
  TH1F * hRatio = new TH1F(*h1);
  Int_t nBinsH1=h1->GetNbinsX();
  Int_t nBinsH2=h2->GetNbinsX();
  // Loop over H1 bins, 
  for(Int_t iBin=1; iBin<=nBinsH1; iBin++){
    hRatio->SetBinContent(iBin,0.);
    hRatio->SetBinContent(iBin,0.);
    Float_t lowPtH1=h1->GetBinLowEdge(iBin);
    Float_t binWidH1=h1->GetBinWidth(iBin);
    // Loop over H2 bins and find overlapping bins to H1
    for(Int_t jBin=1; jBin<=nBinsH2; jBin++){
      Float_t lowPtH2=h2->GetBinLowEdge(jBin);
      Float_t binWidH2=h2->GetBinWidth(jBin);
      if(TMath::Abs(lowPtH1-lowPtH2)<0.001 && TMath::Abs(binWidH2-binWidH1)<0.001){
	Float_t numer=h1->GetBinContent(iBin);
	Float_t denom=h2->GetBinContent(jBin);
	Float_t enumer=h1->GetBinError(iBin);
	Float_t edenom=h2->GetBinError(jBin);
	Double_t ratio=0.;
	Double_t eratio=0.;
	if(numer>0. && denom>0.){
	  ratio=numer/denom;
	  eratio=ratio*TMath::Sqrt((enumer/numer)*(enumer/numer)+(edenom/denom)*(edenom/denom));
	}
	hRatio->SetBinContent(iBin,ratio);
	hRatio->SetBinError(iBin,eratio);
	break;
      }
    }
  }
  return hRatio;
}

Double_t AliBWTools::DoIntegral(TH1* h, Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Double_t & error ,
				Option_t *option, Bool_t doError) 
{
   // function to compute integral and optionally the error  between the limits
   // specified by the bin number values working for all histograms (1D, 2D and 3D)
  // copied from TH! to fix a bug still present in 5-27-06b
   Int_t nbinsx = h->GetNbinsX();
   if (binx1 < 0) binx1 = 0;
   if (binx2 > nbinsx+1 || binx2 < binx1) binx2 = nbinsx+1;
   if (h->GetDimension() > 1) {
      Int_t nbinsy = h->GetNbinsY();
      if (biny1 < 0) biny1 = 0;
      if (biny2 > nbinsy+1 || biny2 < biny1) biny2 = nbinsy+1;
   } else {
      biny1 = 0; biny2 = 0;
   }
   if (h->GetDimension() > 2) {
      Int_t nbinsz = h->GetNbinsZ();
      if (binz1 < 0) binz1 = 0;
      if (binz2 > nbinsz+1 || binz2 < binz1) binz2 = nbinsz+1;
   } else {
      binz1 = 0; binz2 = 0;
   }

   //   - Loop on bins in specified range
   TString opt = option;
   opt.ToLower();
   Bool_t width   = kFALSE;
   if (opt.Contains("width")) width = kTRUE;


   Double_t dx = 1.;
   Double_t dy = 1.;
   Double_t dz = 1.;
   Double_t integral = 0;
   Double_t igerr2 = 0;
   for (Int_t binx = binx1; binx <= binx2; ++binx) {
     if (width) dx = h->GetXaxis()->GetBinWidth(binx);
     for (Int_t biny = biny1; biny <= biny2; ++biny) {
       if (width) dy = h->GetYaxis()->GetBinWidth(biny);
       for (Int_t binz = binz1; binz <= binz2; ++binz) {
	 if (width) dz = h->GetZaxis()->GetBinWidth(binz);
	 Int_t bin  = h->GetBin(binx, biny, binz);
	 if (width) integral += h->GetBinContent(bin)*dx*dy*dz;
	 else       integral += h->GetBinContent(bin);
	 if (doError) {
	   if (width)  igerr2 += h->GetBinError(bin)*h->GetBinError(bin)*dx*dy*dz*dx*dy*dz;
	   else        igerr2 += h->GetBinError(bin)*h->GetBinError(bin);
	 }
	 //	 cout << h->GetBinContent(bin) << " " <<  h->GetBinError(bin) << " " << dx*dy*dz << " "  << integral << " +- " << igerr2 << endl;
	 
       }
     }
   }
   
   if (doError) error = TMath::Sqrt(igerr2);
   return integral;
}

Double_t AliBWTools::dMtdptFunction(Double_t *x, Double_t *p) {

  // Computes the dmt/dptdeta function using the dN/dpt function
  // This is a protected function used internally by GetdMtdy to integrate dN/dpt function using mt as a weight
  // The mass of the particle is given as p[0]
  Double_t pt   = x[0];
  Double_t mass = p[0]; 
  Double_t mt = TMath::Sqrt(pt*pt + mass*mass);
  Double_t jacobian = pt/mt;
  if(!fdNdptForETCalc){
    Printf("AliBWTools::dMtdptFunction: ERROR: fdNdptForETCalc not defined");
    return 0;
  }
  Double_t dNdpt = fdNdptForETCalc->Eval(pt);
  return dNdpt*mt*jacobian; // FIXME: do I have to normalize somehow?
  
}

Double_t AliBWTools::GetdMtdEta(TH1 *hData, TF1 * fExtrapolation, Double_t mass) {
  // Computes dMtdEta integrating dN/dptdy with the proper weights and jacobian.
  Printf("WARNING ALIBWTOOLS::GetdMtdEta: ONLY USING FUNCTION FOR THE TIME BEING");

  // Assign the fiunction used internally by dMtdptFunction
  fdNdptForETCalc = fExtrapolation;
  // Create the function to be integrated
  TF1 * funcdMtdPt = new TF1 ("funcdMtdPt", dMtdptFunction, 0.0, 20, 1);
  funcdMtdPt->SetParameter(0,mass);
  // Integrate it
  Double_t dMtdEta = funcdMtdPt->Integral(0,100);
  // Clean up
  fdNdptForETCalc=0;
  delete funcdMtdPt;
  //return 
  return dMtdEta;

}
