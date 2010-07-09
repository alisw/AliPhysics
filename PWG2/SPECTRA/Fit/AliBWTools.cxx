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

TF1 * AliBWTools::fFuncForNormalized = 0;

ClassImp(AliBWTools)

AliBWTools::AliBWTools() {
  // ctor
}

AliBWTools::~AliBWTools(){
  // dtor
}

TH1 * AliBWTools::GetdNdmtFromdNdpt(TH1 * hpt, Double_t mass) {
  // convert the x axis from pt to mt. Assumes you have 1/pt dNdpt in the histo you start with

  Int_t nbins = hpt->GetNbinsX();
  Float_t * xbins = new Float_t[nbins+1];
  for(Int_t ibins = 0; ibins <= nbins; ibins++){
    xbins[ibins] = TMath::Sqrt(hpt->GetBinLowEdge(ibins+1)*hpt->GetBinLowEdge(ibins+1) +
			       mass *mass) - mass;
//     xbins[ibins] = TMath::Sqrt(hpt->GetBinLowEdge(ibins+1)*hpt->GetBinLowEdge(ibins+1) +
// 			       mass *mass);
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
  
  return hmt;

}

TH1 * AliBWTools::GetdNdPtFromOneOverPt(TH1 * h1Pt) {

  // convert an histo from 1/pt dNdpt to dNdpt, using the pt at the center of the bin


  TH1 * hPt = (TH1 *) h1Pt->Clone((TString(h1Pt->GetName()) + "_inv").Data());
  hPt->Reset();

  Int_t nbinx = hPt->GetNbinsX();

  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){

    Double_t cont = h1Pt->GetBinContent(ibinx);
    Double_t err  = h1Pt->GetBinError(ibinx);
    
    Double_t pt   = h1Pt->GetBinCenter(ibinx);
    
    if(pt != 0) {
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




TH1 * AliBWTools::GetOneOverPtdNdPt(TH1 * hPt) {

  // convert an histo from dNdpt to 1/pt dNdpt, using the pt at the center of the bin

  TH1 * h1Pt = (TH1 *) hPt->Clone((TString(hPt->GetName()) + "_inv").Data());
  h1Pt->Reset();

  Int_t nbinx = h1Pt->GetNbinsX();

  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){

    Double_t cont = hPt->GetBinContent(ibinx);
    Double_t err  = hPt->GetBinError(ibinx);
    
    Double_t pt   = hPt->GetBinCenter(ibinx);
    
    if(pt != 0) {
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


TGraphErrors * AliBWTools::GetGraphFromHisto(TH1F * h, Bool_t binWidth) {
  // Convert a histo to a graph
  // if binWidth is true ex is set to the bin width of the histos, otherwise it is set to zero
  Int_t nbin = h->GetNbinsX();

  TGraphErrors * g = new TGraphErrors();
  Int_t ipoint = 0;
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    Double_t xerr = binWidth ? h->GetBinWidth(ibin)/2 : 0;
    if (h->GetBinContent(ibin)!=0) {
      g->SetPoint     (ipoint,   h->GetBinCenter(ibin), h->GetBinContent(ibin));
      g->SetPointError(ipoint,   xerr,  h->GetBinError(ibin));
      ipoint++;
    }
  }
  
  g->SetMarkerStyle(h->GetMarkerStyle());
  g->SetMarkerColor(h->GetMarkerColor());
  g->SetLineColor(h->GetLineColor());

  g->SetTitle(h->GetTitle());
  g->SetName(TString("g_")+h->GetName());

  return g;

}

TH1F * AliBWTools::GetHistoFromGraph(TGraphErrors * g, TH1F* hTemplate) {

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

TGraphErrors * AliBWTools::ConcatenateGraphs(TGraphErrors * g1,TGraphErrors * g2){

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


TH1F * AliBWTools::Combine3HistosWithErrors(TH1 * h1,  TH1 * h2,  TH1* h3, 
					    TH1 * he1, TH1 * he2, TH1 * he3, 
					    TH1* htemplate, Int_t statFrom, 
					    Float_t renorm1, Float_t renorm2, Float_t renorm3) {

  // Combines 3 histos (h1,h2,h3), weighting by the errors provided in
  // he1,he2,he3, supposed to be the independent systematic errors.
  // he1,he2,he3 are also assumed to have the same binning as h1,h2,h3
  // The combined histo must fit the template provided (no check is performed on this)
  // The histogram are supposed to come from the same (nearly) sample
  // of tracks. statFrom tells the stat error of which of the 3 is
  // going to be assigned to the combined
  // Optionally, it is possible to rescale any of the histograms.

  TH1F * hcomb = (TH1F*) htemplate->Clone(TString("hComb_")+h1->GetName()+"_"+h2->GetName()+h3->GetName());

  // TODO: I should have used an array for h*local...

  // Clone histos locally to rescale them
  TH1F * h1local = (TH1F*) h1->Clone();
  TH1F * h2local = (TH1F*) h2->Clone();
  TH1F * h3local = (TH1F*) h3->Clone();
  h1local->Scale(renorm1);
  h2local->Scale(renorm2);
  h3local->Scale(renorm3);

  TH1 * hStatError = 0;
  if (statFrom == 0)      hStatError = h1; 
  else if (statFrom == 1) hStatError = h2; 
  else if (statFrom == 2) hStatError = h3; 
  else Printf("AliBWTools::Combine3HistosWithErrors: wrong value of the statFrom parameter");
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
    ye[0] = he1->GetBinError(ibin1);
    ye[1] = he2->GetBinError(ibin2);
    ye[2] = he3->GetBinError(ibin3);
 
    // Set error to 0 if content is 0 (means it was not filled)
    if(h1local->GetBinContent(ibin1)==0) ye[0] = 0;
    if(h2local->GetBinContent(ibin2)==0) ye[1] = 0;
    if(h3local->GetBinContent(ibin3)==0) ye[2] = 0;
    
    // Compute weighted mean
    Double_t mean, err;
    WeightedMean(3,y,ye,mean,err);


    // Fill combined
    // TODO: return error from weighted mean somehow...
    hcomb->SetBinContent(ibin,mean);
    Double_t statError = 0;
    if (hStatError->GetBinContent(ibinError) != 0) {
      //      cout << "def" << endl;
      statError = hStatError->GetBinError(ibinError)/hStatError->GetBinContent(ibinError)*hcomb->GetBinContent(ibin);
    }
    else if (h1local->GetBinContent(ibin1) != 0) {
      //      cout << "1" << endl;
      statError = h1local->GetBinError(ibin1)/h1local->GetBinContent(ibin1)*hcomb->GetBinContent(ibin);
    }
    else if (h2local->GetBinContent(ibin2) != 0) {
      //      cout << "2" << endl;
      statError = h2local->GetBinError(ibin2)/h2local->GetBinContent(ibin2)*hcomb->GetBinContent(ibin);
    }
    else if (h3local->GetBinContent(ibin3) != 0) {
      //      cout << "3" << endl;
      statError = h3local->GetBinError(ibin3)/h3local->GetBinContent(ibin3)*hcomb->GetBinContent(ibin);
    }
    hcomb->SetBinError  (ibin,statError);

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


TH1F * AliBWTools::CombineHistos(TH1 * h1, TH1 * h2, TH1* htemplate, Float_t renorm1){

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
    
      if (h1local->GetBinContent(ibin1) == 0 && h2->GetBinContent(ibin2) == 0) {
	// None has data: go to next bin
	hcomb->SetBinContent(ibin,0);
	hcomb->SetBinError  (ibin,0);	
      } else if(h1local->GetBinContent(ibin1) != 0 && h2->GetBinContent(ibin2) == 0) {
	// take data from h1local:
	hcomb->SetBinContent(ibin,h1local->GetBinContent(ibin1));
	hcomb->SetBinError  (ibin,h1local->GetBinError(ibin1));
      } else if(h1local->GetBinContent(ibin1) == 0 && h2->GetBinContent(ibin2) != 0) {
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


void AliBWTools::GetFromHistoGraphDifferentX(TH1F * h, TF1 * f, TGraphErrors ** gBarycentre, TGraphErrors ** gXlw) {

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
    // compute x_LW
    Double_t temp = 1./(max-min) * f->Integral(min,max);
    Double_t epsilon   = 0.000000001;
    Double_t increment = 0.0000000001;
    Double_t x_lw = min;

    while ((f->Eval(x_lw)- temp) > epsilon) {
      x_lw += increment;
      if(x_lw > max) {
	Printf("Cannot find x_lw");
	break;
      }
    }
      
    if (h->GetBinContent(ibin)!=0) {
      (*gBarycentre)->SetPoint     (ipoint,   xbar, h->GetBinContent(ibin));
      (*gBarycentre)->SetPointError(ipoint,   xerr, h->GetBinError(ibin));
      (*gXlw)       ->SetPoint     (ipoint,   x_lw, h->GetBinContent(ibin));
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


Float_t AliBWTools::GetMean(TH1F * h, Float_t min, Float_t max) {

  // Get the mean of histo in a range; root is not reliable in sub ranges with variable binning.

  Int_t minBin = h->FindBin(min);
  Int_t maxBin = h->FindBin(max);

  Float_t mean = 0 ;
  Float_t integral = 0;
  for(Int_t ibin = minBin; ibin < maxBin; ibin++){
    mean     = mean + (h->GetBinCenter(ibin) *  h->GetBinWidth(ibin)* h->GetBinContent(ibin));
    integral = integral + h->GetBinContent(ibin) * h->GetBinWidth(ibin);
  }
  
  return mean/integral;


}

void AliBWTools::GetMean(TF1 * func, Float_t &mean, Float_t &error, Float_t min, Float_t max) {

  // Get the mean of function in a range; 

  return GetMoment("fMean", "x*", func, mean, error, min,max);

}

void AliBWTools::GetMeanSquare(TF1 * func, Float_t &mean, Float_t &error, Float_t min, Float_t max) {

  // Get the mean2 of function in a range; 

  return GetMoment("fMean2", "x*x*", func, mean, error, min,max);


}

void AliBWTools::GetMoment(TString name, TString var, TF1 * func, Float_t &mean, Float_t &error, Float_t min, Float_t max) {

  // returns the integral of a function derived from func by prefixing some operation.
  // Used as a base method for mean and mean2
  Printf("AliBWTools::GetMoment: Error on <pt> is not correct!!! It is overestimated, fix required");
  Int_t npar = func->GetNpar();

  TF1 * f = new TF1(name, var+func->GetName());	// FIXME
//   fFuncForNormalized = func;// TMP: computing mean pt
//   TF1 * f = new TF1(name,GetNormalizedFunc,0,10,npar);// FIXME
//   for(Int_t ipar = 0; ipar < npar; ipar++){ // FIXME
//     f->SetParameter(ipar,func->GetParameter(ipar)); // FIXME
//   }
  
  
  // The parameter of the function used to compute the mean should be
  // the same as the parent function: fixed if needed and they should
  // also have the same errors.

  //  cout << "npar :" << npar << endl;
  
  for(Int_t ipar = 0; ipar < npar; ipar++){
    Double_t parmin, parmax;
    Double_t value = func->GetParameter(ipar);
    func->GetParLimits(ipar, parmin, parmax);
    if ( parmin == parmax )   {
      if ( parmin != 0 || (parmin == 1 && value == 0) ) {
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
  //  f->Print();
//   mean  = f->Integral     (min,max)/func->Integral(min,max);
//   error = f->IntegralError(min,max)/func->Integral(min,max);
  mean  = f->Integral     (min,max);
  error = f->IntegralError(min,max);
//   cout << "Mean " << mean <<"+-"<< error<< endl;
//   cout << "Integral Error "  << error << endl;
  
}

Double_t AliBWTools::GetNormalizedFunc(double * x, double* p){

  // Static function used to provide normalized pointer to a function

  Int_t npar = fFuncForNormalized->GetNpar();
  for(Int_t ipar = 0; ipar < npar; ipar++){ // FIXME
    fFuncForNormalized->SetParameter(ipar,p[ipar]); // FIXME
  }

  return x[0]*fFuncForNormalized->Eval(x[0])/fFuncForNormalized->Integral(0,100);
  
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

Int_t AliBWTools::GetLowestNotEmptyBin(TH1*h) {

  // Return the index of the lowest non empty bin in the histo h

  Int_t nbin = h->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    if(h->GetBinContent(ibin)>0) return ibin;
  }
  
  return -1;

}

Int_t AliBWTools::GetHighestNotEmptyBin(TH1*h) {

  // Return the index of the highest non empty bin in the histo h

  Int_t nbin = h->GetNbinsX();
  for(Int_t ibin = nbin; ibin > 0; ibin--){
    if(h->GetBinContent(ibin)>0) return ibin;
  }
  
  return -1;

}

void AliBWTools::GetResiduals(TGraphErrors * gdata, TF1 * func, TH1F ** hres, TGraphErrors ** gres) {

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

void AliBWTools::GetResiduals(TH1F* hdata, TF1 * func, TH1F ** hres, TH1F ** hresVsBin) {

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
    if(hdata->GetBinContent(ibin)==0) continue;
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

void AliBWTools::GetYield(TH1* h, TF1 * f, Double_t &yield, Double_t &yieldError, Float_t min, Float_t max,
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
  Double_t integralFromHisto = h->IntegralAndError(bin1,bin2,integralFromHistoError,"width");
  
  Double_t integralBelow      = min < bin1Edge ? f->Integral(min,bin1Edge) : 0;
  Double_t integralBelowError = min < bin1Edge ? f->IntegralError(min,bin1Edge) : 0;
  Double_t integralAbove      = max > bin2Edge ? f->Integral(bin2Edge,max) : 0;
  Double_t integralAboveError = max > bin2Edge ? f->Integral(bin2Edge,max) : 0;

  cout << "GetYield INFO" << endl;
  cout << " " << bin1Edge << " " << bin2Edge << endl;  
  cout << " " << integralFromHisto      << " " << integralBelow      << " " << integralAbove      << endl;
  cout << " " << integralFromHistoError << " " << integralBelowError << " " << integralAboveError << endl;
  
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

TGraphErrors * AliBWTools::DivideGraphByFunc(TGraphErrors * g, TF1 * f, Bool_t invert){ 

  // Divides g/f. If invert == true => f/g

  TGraphErrors * gRatio = new TGraphErrors();
  Int_t npoint = g->GetN();
  for(Int_t ipoint = 0; ipoint < npoint; ipoint++){
    Double_t x = g->GetX()[ipoint];
    Double_t ratio  = invert ? f->Eval(x)/g->GetY()[ipoint] :g->GetY()[ipoint]/f->Eval(x);
    gRatio->SetPoint     (ipoint, x, ratio);
    gRatio->SetPointError(ipoint, 0, g->GetEY()[ipoint]/f->Eval(x));
    //    cout << x << " " << g->GetY()[ipoint] << " " << f->Eval(x) << endl;
    
  }
  gRatio->SetMarkerStyle(20);
  //gRatio->Print();
  return gRatio;

}

TGraphErrors * AliBWTools::DivideGraphByHisto(TGraphErrors * g, TH1 * h, Bool_t invert){ 

  // Divides g/h. If invert == true => h/g


  TGraphErrors * gRatio = new TGraphErrors();
  Int_t npoint = g->GetN();
  for(Int_t ipoint = 0; ipoint < npoint; ipoint++){
    Double_t xj  = g->GetX()[ipoint];
    Double_t yj  = g->GetY()[ipoint];
    Double_t yje = g->GetEY()[ipoint];

    Int_t binData = h->FindBin(xj);
    Double_t yd   = h->GetBinContent(binData);
    Double_t yde  = h->GetBinError(binData);
    Double_t xd   = h->GetBinCenter(binData);
    
    //    cout << TMath::Abs((xd-xj)/xd) << endl;
    

     
    if (yd == 0) continue;
    //    if (binData == 28 ) cout << TMath::Abs(xd-xj)/TMath::Abs(xd) << endl;
    
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
    }
    
  }


  if(sumweight){
    mean /= sumweight;
    meanerr = TMath::Sqrt(1./ sumweight);
  }
  else {
    mean = 0;
    meanerr = 0;
  }

  
}
