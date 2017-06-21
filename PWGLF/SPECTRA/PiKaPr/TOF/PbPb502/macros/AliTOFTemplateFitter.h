#ifndef ALITOFTEMPLATEFITTER
#define ALITOFTEMPLATEFITTER
//Flags to set the modes to be used
#define USECDECONVOLUTION//Cholesky-like
#define USEFITFUNCTIONS//Fit functions
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <RooAddPdf.h>
#include <RooRealVar.h>
#include <RooChi2Var.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooGaussianTail.h>
#include <RooHistPdf.h>
#include <RooAbsPdf.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TFractionFitter.h>
#include <TObjArray.h>
#include <TOFsignal.C>
#include <AliUtilTOFParams.h>
#include <UtilMessages.h>
#include "TMath.h"
#ifdef USECDECONVOLUTION
#include "AliCDeconv.h"
#endif
#endif

using namespace RooFit;

////////////////////////////////////////////////////////////////////////////
///                                                                       //
///                                                                       //
/// Utilities for template fitting e.g. yield extraction                  //
///                                                                       //
///                                                                       //
/// Authors:                                                              //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                        //
////////////////////////////////////////////////////////////////////////////

Double_t CHI2 = -1;
TCanvas *cFit = 0x0;//TCanvas of the single fit

//Utilities
//_________________________________________________________________________________________________
Double_t ComputeChi2(const TH1 *hdata, const TH1 *hfit, const Double_t xlow, const Double_t xhigh){//Macro to compute chi2
  if(xhigh < xlow){
    Warningmsg("ComputeChi2", Form("Range for Chi2 not well defined [%f, %f]", xlow, xhigh));
    return -1;
  }
  TH1 *h1 = (TH1*)hdata->Clone("h1");
  TH1 *h2 = (TH1*)hfit->Clone("h2");
  for (Int_t bin = 1; bin < h1->GetXaxis()->FindBin(xlow); bin++) h1->SetBinContent(bin, 0);
  for (Int_t bin = h1->GetXaxis()->FindBin(xlow) + 1; bin <= h1->GetNbinsX(); bin++) h1->SetBinContent(bin, 0);
  for (Int_t bin = 1; bin < h2->GetXaxis()->FindBin(xlow); bin++) h2->SetBinContent(bin, 0);
  for (Int_t bin = h2->GetXaxis()->FindBin(xlow) + 1; bin <= h2->GetNbinsX(); bin++) h2->SetBinContent(bin, 0);
  
  const Double_t chi2 = h1->Chi2Test(h2, "UUCHI2/NDF"); 
  delete h1;
  delete h2;
  return chi2;
}

//_________________________________________________________________________________________________
Bool_t PerformFitWithTFF(TH1F* hData, TObjArray* mc, Double_t *range, Double_t *fitrange, TArrayD &fraction, TArrayD &fractionErr, TObjArray *& prediction){
  //TH1F* hData             -> Histogram of data to be fitted
  //TObjArray* mc           -> Templates for the fit
  //Double_t *range         -> Integration range over which to compute fractions
  //Double_t *fitrange      -> Fit range
  //TArrayD &fraction       -> Fit results: fractions computed in range
  //TArrayD &fractionErr    -> Fit results: errors on fractions
  //TObjArray *& prediction -> Fit results: Fitted templates and also global fit result (ntemplates+1 objects)
  
  const Int_t ntemplates = mc->GetEntries();
  Double_t min = 100;
  Int_t fitlow = 1, fithigh = hData->GetNbinsX();
  if(fitrange[0] < fitrange[1]){
    fitlow = hData->GetXaxis()->FindBin(fitrange[0]);
    fithigh = hData->GetXaxis()->FindBin(fitrange[1]);
  }
  
  Infomsgcolor("PerformFitWithTFF", Form("Data h. summary: name (%s)  nbins (%i) first bin (%.3f) last bin (%.3f) entries (%.2f) eff. (%.2f)", hData->GetName(), hData->GetNbinsX(), hData->GetXaxis()->GetBinLowEdge(1), hData->GetXaxis()->GetBinUpEdge(hData->GetNbinsX()), hData->GetEntries(), hData->GetEffectiveEntries()), cyanTxt);
  for(Int_t c = 0; c<ntemplates; c++) Infomsgcolor("PerformFitWithTFF", Form("Template h. #%i/%i summary: name (%s)  nbins (%i) first bin (%.3f) last bin (%.3f) entries (%.2f) eff. (%.2f)", c+1, ntemplates, ((TH1F*)mc->At(c))->GetName(), ((TH1F*)mc->At(c))->GetNbinsX(), ((TH1F*)mc->At(c))->GetXaxis()->GetBinLowEdge(1), ((TH1F*)mc->At(c))->GetXaxis()->GetBinUpEdge(((TH1F*)mc->At(c))->GetNbinsX()), ((TH1F*)mc->At(c))->GetEntries(), ((TH1F*)mc->At(c))->GetEffectiveEntries()), cyanTxt);
  
  TFractionFitter* fit = new TFractionFitter(hData, mc, "Q"); // initialise
  for(Int_t c = 0; c<ntemplates; c++){
    if (fraction[c] > 0) {
      Infomsg("PerformFitWithTFF", Form("Setting maximum for template %i to %f", c, fraction[c]));
      fit->Constrain(c, 0.000, fraction[c]); // constrain fraction to be between 0 and required maximum
    }
    else fit->Constrain(c, 0.000, 1.0); // constrain fraction to be between 0 and 1
  } 
  fit->GetFitter()->SetMaxIterations(100000000);
  
  Int_t StopSign = 0;//Check if the fit is feasible
  
  TH1F *htemplates[ntemplates];
  for(Int_t temp = 0; temp < ntemplates; temp++){//Check on templates
    htemplates[temp] = (TH1F*) mc->At(temp)->Clone(Form("Input template %i", temp));
    htemplates[temp]->SetTitle(Form("Input template %i", temp));
    htemplates[temp]->SetLineStyle(2);
    if(htemplates[temp]->Integral(fitlow, fithigh) < min){
      Warningmsg("PerformFitWithTFF", Form("%s does NOT have sufficient entries (%.0f/%.0f) but it has %.0f entries", htemplates[temp]->GetName(), htemplates[temp]->Integral(fitlow, fithigh), min, htemplates[temp]->GetEntries()));
      StopSign = -1;
    }
  }
  if(hData->Integral(fitlow, fithigh) < min){
    Warningmsg("PerformFitWithTFF", Form("%s does NOT have sufficient entries (%.0f/%.0f) but it has %.0f entries", hData->GetName(), hData->Integral(), min, hData->GetEntries()));
    StopSign  = -1; //Check on data
  }
  
  //fit->SetRangeX(200, 1000);                    // use only the first 15 bins in the fit
  if(fitrange[0] < fitrange[1]) fit->SetRangeX(hData->GetXaxis()->FindBin(fitrange[0]), hData->GetXaxis()->FindBin(fitrange[1]));
  Int_t status = -1;
  if(StopSign == 0){
    Infomsgcolor("PerformFitWithTFF", Form("I've checked that the fit is feasible (non empty histograms) ---> Fitting!!"), cyanTxt);
    status = fit->Fit();               // perform the fit
    if(0){//status == 4){
      delete fit;
      Warningmsg("PerformFitWithTFF", "Fitting Again due to poor fitting conditions before (fit previously failed!)");
      
      fit = new TFractionFitter(hData, mc, "Q");
      fit->GetFitter()->SetMaxIterations(100000000);
      
      for(Int_t c = 0; c<ntemplates; c++){
        for(Int_t cc = 1 ; cc <= htemplates[c]->GetNbinsX(); cc++){//Excluding empty bins from fit
          if(htemplates[c]->GetBinContent(cc) == 0 && fit->IsExcluded(cc) == kFALSE) fit->ExcludeBin(cc);
        }
      }
      status = fit->Fit();               // perform the fit AGAIN
    }
  }
  else{
    Warningmsg("PerformFitWithTFF", "Fit not performed since minimum requirements are not fulfilled:\nTemplates: ");
    for(Int_t c = 0; c<ntemplates; c++) cout<<" htemplates"<<c<<" ("<<htemplates[c]->Integral()<<")"<<flush;
    cout<<" hData ("<<hData->Integral()<<")"<<endl;
  }
  
  Infomsgcolor("PerformFitWithTFF", Form("Fit status: %i (0 means success!)", status), cyanTxt);
  
  Double_t value[ntemplates], error[ntemplates];//Fit results and errors
  Double_t yield[ntemplates], yielderror[ntemplates];//Fit results and errors
  
  //Drawing section for data and templates
  hData->SetLineColor(1);
  hData->DrawCopy("Ep");
  for(Int_t temp = 0; temp < ntemplates; temp++) htemplates[temp]->DrawCopy("same");
  
  if(status == 0){                       // check on fit status
    
    //Get the fit result
    TH1F* result = (TH1F*) fit->GetPlot();
    //Get the fitted template
    TH1F* MCPred[ntemplates];
    for(Int_t temp = 0; temp < ntemplates; temp++){
      MCPred[temp] = (TH1F*)fit->GetMCPrediction(temp);
      MCPred[temp]->SetLineColor(htemplates[temp]->GetLineColor());
    }
    
    //Drawing section after FIT
    result->SetTitle("Fit result");
    result->SetName("Fit result");
    result->SetLineColor(kMagenta);
    result->SetMarkerColor(kMagenta);
    result->DrawCopy("pesame");
    prediction->Add((TH1F *) result->Clone(Form("result")));
    
    
    //Getting results section
    Double_t factor, dataIntegral, tempIntegral, sum = 0;
    
    dataIntegral = hData->GetSumOfWeights();
    
    for(Int_t temp = 0; temp < ntemplates; temp++){
      fit->GetResult(temp, value[temp], error[temp]);
      tempIntegral = MCPred[temp]->GetSumOfWeights();
      sum += tempIntegral*value[temp];
      factor = dataIntegral/tempIntegral*value[temp];
      MCPred[temp]->SetLineStyle(2);
      MCPred[temp]->SetTitle(Form("Fitted Template %i unscaled", temp));
      MCPred[temp]->SetName(Form("Fitted Template %i unscaled", temp));
      //       MCPred[temp]->DrawCopy("same");
      MCPred[temp]->SetLineStyle(1);
      MCPred[temp]->Scale(factor);
      MCPred[temp]->SetTitle(Form("Fitted Template %i", temp));
      MCPred[temp]->SetName(Form("Fitted Template %i", temp));
      MCPred[temp]->DrawCopy("same");
      prediction->Add((TH1F *) MCPred[temp]->Clone(Form("prediction%i", temp)));
      Infomsgcolor("PerformFitWithTFF", Form("Fraction %i: %f (%f) scaling factor: %f / %f * %f  = %f", temp, value[temp], error[temp], dataIntegral, tempIntegral, value[temp], factor), yellowTxt);
      
      if(range[0] >= range[1]) Fatalmsg("PerformFitWithTFF", "Integration range badly defined");
      
      const Int_t binlow = MCPred[temp]->FindBin((1.-0.0001)*range[0]);
      const Int_t binup = MCPred[temp]->FindBin((1.-0.0001)*range[1]);
      if(temp == 0){
        Infomsgcolor("PerformFitWithTFF", Form("The sum of the integrals of the %i functions from the fit is %f while the Integral of the data is %f, the difference is %f", ntemplates, sum, dataIntegral, sum-dataIntegral), yellowTxt);
        
        Infomsgcolor("PerformFitWithTFF", Form("Now obtaining the total yield in the range [%f, %f]", MCPred[temp]->GetXaxis()->GetBinLowEdge(binlow), MCPred[temp]->GetXaxis()->GetBinUpEdge(binup)), cyanTxt);
      }
      yield[temp] = MCPred[temp]->IntegralAndError(binlow, binup, yielderror[temp]);
      Infomsgcolor("PerformFitWithTFF", Form("Template %i from Integral: %f (%f)", temp, yield[temp], yielderror[temp]), yellowTxt);
    }
    
  }
  else{
    Warningmsg("PerformFitWithTFF", Form("Fit not performed successfully so all fractions are set equal"));
    for(Int_t temp = 0; temp < ntemplates; temp++) yield[temp] = 1;
  }
  
  if(fit){
    delete fit;
    fit = 0;
  }
  
  Double_t Den = 0, eDen = 0, Num = 0, eNum = 0, Ratio = 0, RatioErr = 0;
  for(Int_t temp = 0; temp < ntemplates; temp++){
    Den += yield[temp];
    eDen += yielderror[temp]*yielderror[temp];
  }
  eDen = TMath::Sqrt(eDen);
  
  if(Den > 0){
    
    for(Int_t temp = 0; temp < ntemplates; temp++){
      Num = yield[temp];
      eNum = yielderror[temp];
      Ratio = Num/Den;
      
      if(Num>0) RatioErr = Ratio*TMath::Sqrt( eNum*eNum/(Num*Num) + eDen*eDen/(Den*Den));
      else RatioErr = 0;
      //       fraction[temp] = Ratio;
      //       fractionErr[temp] = RatioErr;
      fraction[temp] = yield[temp];
      fractionErr[temp] = yielderror[temp];
    }
  }
  
  TString f[2];
  for(Int_t temp = 0; temp < ntemplates; temp++){
    f[0] += Form(" Template %i = %f(%f)", temp, fraction[temp], fractionErr[temp]);
    f[1] += Form(" Template %i = %f(%f)", temp, value[temp], error[temp]);
  }
  
  Infomsgcolor("PerformFitWithTFF", Form("Fractions from Integral: %s", f[0].Data()), cyanTxt);
  Infomsgcolor("PerformFitWithTFF", Form("Fractions from fit     :  %s\n", f[1].Data()), cyanTxt);
  
  if(status == 0) return kTRUE;
  else return kFALSE;
}

//_________________________________________________________________________________________________
Bool_t PerformFitWithRooFit(TH1F* hData, TObjArray *mc, Double_t *range, Double_t *fitrange, TArrayD &fraction, TArrayD &fractionErr, TObjArray *& prediction, Double_t &chi2 = CHI2){
  const Int_t ntemplates = mc->GetEntries();
  Int_t templatetype[ntemplates];//Template type all set to 0 for TH1F by default but can be also 1 for TF1
  const Bool_t normtoone = kFALSE;
  const Bool_t recursive = kFALSE;
  const Bool_t chi2fit = kFALSE;
  const Double_t min = 100;
  const Bool_t showfit = kTRUE;//Flag to draw the input data before the fit in a new TCanvas
  const Int_t tempchi2 = chi2;//Value of the template on which to compute chi2
  
  //If fitrange[1] < fitrange[0] you are requiring that the fitrange is decided in the range where all templates are positive
  Double_t minx = fitrange[0], maxx = fitrange[1];
  
  Infomsg("PerformFitWithRooFit", Form("Fitting data with %i templates between %f and %f", ntemplates, minx, maxx));
  Infomsgcolor("PerformFitWithRooFit", Form("Data histogram summary: name (%s)  nbins (%i) first bin (%.3f) last bin (%.3f)", hData->GetName(), hData->GetNbinsX(), hData->GetXaxis()->GetBinLowEdge(1), hData->GetXaxis()->GetBinUpEdge(hData->GetNbinsX())), cyanTxt);
  for(Int_t c = 0; c<ntemplates; c++){
    const TString cname = mc->At(c)->ClassName();
    if(cname.Contains("TH1")){
      templatetype[c] = 0;
      TH1F * h = (TH1F*)mc->At(c);
      Infomsgcolor("PerformFitWithRooFit", Form("Template histogram #%i/%i summary: name (%s)  nbins (%i) first bin (%.3f) last bin (%.3f)", c+1, ntemplates, h->GetName(), h->GetNbinsX(), h->GetXaxis()->GetBinLowEdge(1), h->GetXaxis()->GetBinUpEdge(h->GetNbinsX())), cyanTxt);
    }
    else if(cname.Contains("TF1")){
      templatetype[c] = 1;
      TF1 * f = (TF1*)mc->At(c);
      Double_t r[2];
      f->GetRange(r[0], r[1]);
      Int_t n = f->GetNpar();
      Infomsgcolor("PerformFitWithRooFit", Form("Template function #%i/%i summary: name (%s)  first bin (%.3f) last bin (%.3f) #parameters %i", c+1, ntemplates, f->GetName(), r[0], r[1], n), cyanTxt);
    }
    else Fatalmsg("PerformFitWithRooFit", Form("Input class %s for template %i is not recognized!", cname.Data(), c));
  }
  
  Int_t StopSign = 0;//Check if the fit is feasible
  
  TVirtualPad *currpad = gPad;
  if(showfit){//TCanvas for debug purposes
    TLine linemin(minx, 0, minx, hData->GetMaximum());
    TLine linemax(maxx, 0, maxx, hData->GetMaximum());
    linemin.SetLineColor(kRed);
    linemax.SetLineColor(kRed);
    
    if(cFit == 0x0){
      cFit = new TCanvas("TemplateRooFit", "Fit with roofit");
    }
    cFit->Clear();
    cFit->Divide(2);
    
    cFit->cd(1);
    gPad->SetLogy();
    hData->DrawCopy()->GetXaxis()->SetRangeUser(minx, maxx);
    linemin.DrawClone();
    linemax.DrawClone();
    
    cFit->cd(2);
    gPad->SetLogy();
    if(hData->GetEffectiveEntries() > 0) hData->DrawNormalized()->GetXaxis()->SetRangeUser(minx, maxx);
    else hData->DrawCopy()->GetXaxis()->SetRangeUser(minx, maxx);
    linemin.DrawClone();
    linemax.DrawClone();
  }
  
  //Define xvariable
  RooRealVar x("x", Form("%s", hData->GetXaxis()->GetTitle()), range[0], range[1]);//Variable of the x axis
  
  TH1F *htemplates[ntemplates];
  RooAbsPdf *ftemplates[ntemplates];
  const Int_t npar = 3;//number of parameters
  const TString parname[npar] = {"Mean", "Sigma", "Tail"};
  RooRealVar fpar[ntemplates][npar];
  const Bool_t fitsigma = kFALSE;
  
  RooRealVar mean("mean", "mean", 0, fitsigma ? -3 : -50, fitsigma ? +3 : 50., "");
  RooRealVar sigma("sigma", "sigma", fitsigma ? 1 : 80, fitsigma ? 0.5 : 70, fitsigma ? 1.5 : 120, "");
  RooRealVar tail("tail", "tail", 1.25, 0.5, 1.5, "");
  
  for(Int_t temp = 0; temp < ntemplates; temp++){//Check on templates
    htemplates[temp] = 0x0;
    ftemplates[temp] = 0x0;
    
    if(templatetype[temp] == 0){//TH1
      htemplates[temp] = (TH1F*) mc->At(temp)->Clone(Form("IT_%s%i", ((TH1F*) mc->At(temp))->GetName(), temp));
      htemplates[temp]->SetTitle(Form("Input template %i", temp));
      htemplates[temp]->SetLineStyle(2);
      if(htemplates[temp]->Integral() < min){
        Warningmsg("PerformFitWithRooFit", Form("%s does NOT have sufficient entries (%.0f/%.0f) but it has %.0f entries", htemplates[temp]->GetName(), htemplates[temp]->Integral(), min, htemplates[temp]->GetEntries()));
        StopSign = -1;
      }
      if(showfit){
        cFit->cd(1);
        htemplates[temp]->DrawCopy("same");
        cFit->cd(2);
        htemplates[temp]->DrawNormalized("same");
        cFit->Flush();
        cFit->Update();
        gSystem->ProcessEvents();
      }
    }
    else{//RooAbsPdf
      TF1 * f = (TF1*) mc->At(temp);
      if(npar != f->GetNpar()) Fatalmsg("PerformFitWithRooFit", "Different number of parameters");
      
      Double_t p[npar], pLow[npar], pUp[npar];
      
      TString parstring = "";
      for(Int_t par = 0; par < npar; par++){
        p[par] = f->GetParameter(par);
        f->GetParLimits(par, pLow[par], pUp[par]);
        parstring += Form("%f [%f, %f] ", p[par], pLow[par], pUp[par]);
        fpar[temp][par] = RooRealVar(parname[par], parname[par], p[par], pLow[par], pUp[par], "");
        
      }
      Infomsg("PerformFitWithRooFit", Form("Defining new fit functions %s with parameters %s", f->GetName(), parstring.Data()));
      ftemplates[temp] = new RooGaussianTail("signal", "bkg1", x, mean, sigma, tail);
      
      if(showfit){
        cFit->cd(1);
        f->DrawCopy("same");
        cFit->cd(2);
        f->DrawCopy("same");
      }
      
    }
  }
  if(hData->Integral() < min){//Check on data
    Warningmsg("PerformFitWithRooFit", Form("%s does NOT have sufficient entries (%.0f/%.0f) but it has %.0f entries", hData->GetName(), hData->Integral(), min, hData->GetEntries()));
    StopSign  = -1; //Check on data
  }
  currpad->cd();
  
  if(StopSign < 0){
    for(Int_t temp = 0; temp < ntemplates; temp++) if(ftemplates[temp]) delete ftemplates[temp];
    return kFALSE;
  }
  
  if(minx > maxx){//Setting extremes of the fit
    Int_t minbin = 1;
    Int_t maxbin = 1;
    Int_t bin = 1;
    Bool_t first = kTRUE;
    for (Int_t temp = 0; temp < ntemplates; temp++) {
      if(templatetype[temp] != 0) continue;
      
      if(first){
        minbin = htemplates[temp]->FindFirstBinAbove(1);
        maxbin = htemplates[temp]->FindLastBinAbove(1);
        first = kFALSE;
      }
      else{
        bin = htemplates[temp]->FindFirstBinAbove(1);
        if(minbin > bin) minbin = bin;
        
        bin = htemplates[temp]->FindLastBinAbove(1);
        if(maxbin < bin) maxbin = bin;
      }
    }
    minx = hData->GetXaxis()->GetBinLowEdge(minbin);
    maxx = hData->GetXaxis()->GetBinUpEdge(maxbin);
    
  }
  
  RooRealVar *frac[ntemplates];
  for (Int_t temp = 0; temp < ntemplates; temp++) {
    frac[temp] = new RooRealVar(Form("frac%i", temp), Form("frac%i", temp), fraction[temp] != 0 ? TMath::Abs(fraction[temp]) : 1., 0., ((Double_t) hData->GetEntries()));
    if(fraction[temp] < 0){
      Infomsg("PerformFitWithRooFit", Form("Requested to fix the fraction for template %i to %f", temp, -fraction[temp]));
      frac[temp]->setConstant(kTRUE);
    }
    
  }
  
  if(range[0] >= range[1]) Fatalmsg("PerformFitWithRooFit", "Range is not well defined");
  
  RooDataHist roohData(Form("roohData%s", hData->GetName()), Form("roohData%s", hData->GetTitle()), x, hData);//Data histogram
  
  RooDataHist *roohtemplates[ntemplates];
  RooHistPdf *roohPdf[ntemplates];
  
  for(Int_t temp = 0; temp < ntemplates; temp++){//Get the templates
    roohtemplates[temp] = 0x0;
    roohPdf[temp] = 0x0;
    
    switch ((templatetype[temp])) {
      case 0:
      {
        roohtemplates[temp] = new RooDataHist(Form("Roo%s", htemplates[temp]->GetName()), Form("Roo %s", htemplates[temp]->GetTitle()), x, htemplates[temp]);
        
        roohPdf[temp] = new RooHistPdf (Form("PDF%s", htemplates[temp]->GetName()), Form("PDF%s", htemplates[temp]->GetTitle()), x, *roohtemplates[temp]);
        break;
      }
      case 1:
      {
        break;
      }
      default:
      Fatalmsg("PerformFitWithRooFit", "Cannot find good histogram type");
      break;
    }
    
  }
  
  RooAddPdf *convolution = 0x0;//Prepare the convolution
  RooArgList templist("templist");//List of the templates
  RooArgList flist("flist");//List of the parameters
  if(StopSign == 0){
    
    if(normtoone){
      if(recursive){
        if(ntemplates == 1){
          convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), RooArgList(*roohPdf[0]), RooArgList(*frac[0]));
          
        }
        else if(ntemplates == 2){
          RooAbsPdf * step = roohPdf[1];
          convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), RooArgList(*roohPdf[0], *step), RooArgList(*frac[0]));
          
        }
        else if(ntemplates == 3){
          RooAbsPdf * step = new RooAddPdf(Form("step%s", hData->GetName()), Form("step%s", hData->GetName()), RooArgList(*roohPdf[1], *roohPdf[2]), RooArgList(*frac[1]));
          convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), RooArgList(*roohPdf[0], *step), RooArgList(*frac[0]));
          
        }
        else if(ntemplates == 4){
          RooAbsPdf * step_1 = new RooAddPdf(Form("step_1%s", hData->GetName()), Form("step_1%s", hData->GetName()), RooArgList(*roohPdf[2], *roohPdf[3]), RooArgList(*frac[2]));
          RooAbsPdf * step = new RooAddPdf(Form("step%s", hData->GetName()), Form("step%s", hData->GetName()), RooArgList(*roohPdf[1], *step_1), RooArgList(*frac[1]));
          convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), RooArgList(*roohPdf[0], *step), RooArgList(*frac[0]));
          
        }
        else if(ntemplates == 5){
          RooAbsPdf * step_1 = new RooAddPdf(Form("step_1%s", hData->GetName()), Form("step_1%s", hData->GetName()), RooArgList(*roohPdf[3], *roohPdf[4]), RooArgList(*frac[3]));
          RooAbsPdf * step_2 = new RooAddPdf(Form("step_2%s", hData->GetName()), Form("step_2%s", hData->GetName()), RooArgList(*roohPdf[2], *step_1), RooArgList(*frac[2]));
          RooAbsPdf * step = new RooAddPdf(Form("step%s", hData->GetName()), Form("step%s", hData->GetName()), RooArgList(*roohPdf[1], *step_2), RooArgList(*frac[1]));
          convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), RooArgList(*roohPdf[0], *step), RooArgList(*frac[0]));
          
          
        }
        else if(ntemplates == 6){
          RooAbsPdf * step_1 = new RooAddPdf(Form("step_1%s", hData->GetName()), Form("step_1%s", hData->GetName()), RooArgList(*roohPdf[4], *roohPdf[5]), RooArgList(*frac[4]));
          RooAbsPdf * step_2 = new RooAddPdf(Form("step_2%s", hData->GetName()), Form("step_2%s", hData->GetName()), RooArgList(*roohPdf[3], *step_1), RooArgList(*frac[3]));
          RooAbsPdf * step_3 = new RooAddPdf(Form("step_3%s", hData->GetName()), Form("step_3%s", hData->GetName()), RooArgList(*roohPdf[2], *step_2), RooArgList(*frac[2]));
          RooAbsPdf * step = new RooAddPdf(Form("step%s", hData->GetName()), Form("step%s", hData->GetName()), RooArgList(*roohPdf[1], *step_3), RooArgList(*frac[1]));
          convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), RooArgList(*roohPdf[0], *step), RooArgList(*frac[0]));
          
        }
        else{
          Fatalmsg("PerformFitWithRooFit", "Number of templates not yet implemented");
          convolution = 0x0;
        }
      }
      else{
        if(ntemplates == 1) convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), RooArgList(*roohPdf[0]), RooArgList(*frac[0]));
        else if(ntemplates == 2) convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), RooArgList(*roohPdf[0], *roohPdf[1]), RooArgList(*frac[0]));
        else if(ntemplates == 3) convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), RooArgList(*roohPdf[0], *roohPdf[1], *roohPdf[2]), RooArgList(*frac[0], *frac[1]));
        else if(ntemplates == 4) convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), RooArgList(*roohPdf[0], *roohPdf[1], *roohPdf[2], *roohPdf[3]), RooArgList(*frac[0], *frac[1], *frac[2]));
        else if(ntemplates == 5) convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), RooArgList(*roohPdf[0], *roohPdf[1], *roohPdf[2], *roohPdf[3], *roohPdf[4]), RooArgList(*frac[0], *frac[1], *frac[2], *frac[3]));
        else if(ntemplates == 6) convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), RooArgList(*roohPdf[0], *roohPdf[1], *roohPdf[2], *roohPdf[3], *roohPdf[4], *roohPdf[5]), RooArgList(*frac[0], *frac[1], *frac[2], *frac[3], *frac[4]));
        else{
          Fatalmsg("PerformFitWithRooFit", "Number of templates not yet implemented");
          convolution = 0x0;
        }
        
      }
    }
    else{
      for(Int_t temp = 0; temp < ntemplates; temp++){
        switch ((templatetype[temp])) {
          case 0:
          {
            templist.add(*roohPdf[temp]);
            break;
          }
          case 1:
          {
            templist.add(*ftemplates[temp]);
            break;
          }
          default:
          Fatalmsg("PerformFitWithRooFit", "Cannot find good histogram type");
          break;
        }
        flist.add(*frac[temp]);
        
      }
      convolution = new RooAddPdf(Form("convolution%s", hData->GetName()), Form("convolution%s", hData->GetName()), templist, flist);
      
    }
    Infomsg("PerformFitWithRooFit", "Fitting!!");
    convolution->Print();
    
    RooFitResult *result = 0x0;
    
    if(chi2fit) result = convolution->chi2FitTo(roohData, Save(), Binned(kTRUE),  Range(minx, maxx));
    else result = convolution->fitTo(roohData, Save(), Range(minx, maxx), Extended(kTRUE), SumW2Error(kFALSE), Verbose(kFALSE), PrintEvalErrors(10));
    
    if(normtoone){
      Double_t sum = 0;
      for(Int_t temp = 0; temp < ntemplates-1; temp++) sum += frac[temp]->getVal();
      frac[ntemplates-1]->setVal(1-sum);
    }
    Infomsg("PerformFitWithRooFit", "Printing results:");
    if(result) result->Print("V");
    //     if(!result) return kFALSE;
    
    TH1F *htemp = (TH1F*)convolution->createHistogram(Form("roofitpredictiontotal%s", hData->GetName()), x);
    htemp->SetTitle("Fitted result");
    htemp->SetLineColor(kMagenta);
    htemp->Scale(normtoone ? hData->GetEntries() : hData->GetBinWidth(hData->GetNbinsX()/2));
    prediction->Add(htemp);
    
    //     convolution->plotOn(dataFrame, LineColor(kGreen));
    for(Int_t temp = 0; temp < ntemplates; temp++){
      //       convolution->plotOn(dataFrame, Components(*roohPdf[temp]), LineColor(htemplates[temp]->GetLineColor()));
      htemp = (TH1F*)roohPdf[temp]->createHistogram(Form("roofitprediction%s%i", hData->GetName(), temp), x);
      htemp->SetTitle(Form("Fitted Template %i", temp));
      htemp->Scale((normtoone ? hData->GetEntries() : 1)*frac[temp]->getVal());
      htemp->SetLineColor(htemplates[temp]->GetLineColor());
      prediction->Add(htemp);
    }
    
  }
  
  
  //   dataFrame->Draw();
  
  hData->SetLineColor(kBlack);
  hData->DrawCopy();
  for(Int_t temp = 0; temp < prediction->GetEntries(); temp++){
    TH1F *htemp = (TH1F*)prediction->At(temp);
    //     if(temp == 0) htemp->DrawCopy();
    htemp->DrawCopy("same");
    
  }
  
  TString f;
  
  for(Int_t temp = 0; temp < ntemplates; temp++){
    htemplates[temp]->SetLineStyle(2);
    htemplates[temp]->DrawCopy("same");
    
    fraction[temp] = frac[temp]->getVal();
    fractionErr[temp] = frac[temp]->getError();
    
    f += Form(" Template %i = %f(%f)", temp, fraction[temp], fractionErr[temp]);
  }
  Infomsgcolor("PerformFitWithRooFit", Form("Fractions from fit     :  %s\n", f.Data()), cyanTxt);
  cout<<"Entries : "<<hData->GetEntries()<<endl;
  
  const Double_t chi2range[2] = {tempchi2 >= 0 ? htemplates[tempchi2]->GetXaxis()->GetBinLowEdge(htemplates[tempchi2]->FindFirstBinAbove(1)) : minx, tempchi2 >= 0 ? htemplates[tempchi2]->GetXaxis()->GetBinUpEdge(htemplates[tempchi2]->FindLastBinAbove(1)) : maxx};
  RooDataHist* subrange = static_cast<RooDataHist*>(roohData.reduce(Form("x>%f", chi2range[0])));
  subrange = static_cast<RooDataHist*>(subrange->reduce(Form("x<%f", chi2range[1])));
  //   RooDataHist* subrange = static_cast<RooDataHist*>(roohData.reduce(Range(chi2range[0], chi2range[1])));
  //   subrange->Print("v");
  RooChi2Var Chi2("Chi2", "Chi2", *convolution, *subrange, Extended(), Range(chi2range[0], chi2range[1]));
  chi2 = Chi2.getVal()/(subrange->numEntries()-ntemplates);
  cout<<chi2<<" wrt "<<ComputeChi2(hData, (TH1F*)prediction->At(0), chi2range[0], chi2range[1])<<endl;
  
  if(convolution) delete convolution;
  for(Int_t temp = 0; temp < ntemplates; temp++){
    if(htemplates[temp]) delete htemplates[temp];
    if(ftemplates[temp]) delete ftemplates[temp];
    if(roohtemplates[temp]) delete roohtemplates[temp];
    if(roohPdf[temp]) delete roohPdf[temp];
    if(frac[temp]) delete frac[temp];
  }
  
  return kTRUE;
  
}

//_________________________________________________________________________________________________
Bool_t PerformFitWithFunctions(TH1F* hData, TObjArray *func, TF1 *funcsum, Double_t *range, Double_t *fitrange, TArrayD &fraction, TArrayD &fractionErr, TObjArray *& prediction){
  #ifdef USEFITFUNCTIONS
  const Int_t ntemplates = func->GetEntries();
  TF1 *singlefun[ntemplates];
  Double_t *parameters[ntemplates];
  Int_t nparameters[ntemplates];
  const Int_t nparameterssum = funcsum->GetNpar();
  const Double_t min = 100;
  const Bool_t showfit = kTRUE;//Flag to draw the input data before the fit in a new TCanvas
  
  //If fitrange[1] < fitrange[0] you are requiring that the fitrange is decided in the range where all templates are positive
  Double_t minx = fitrange[0], maxx = fitrange[1];
  
  Infomsg("PerformFitWithFunctions", Form("Fitting data with %i templates between %f and %f", ntemplates, minx, maxx));
  Infomsgcolor("PerformFitWithFunctions", Form("Data histogram summary: name (%s)  nbins (%i) first bin (%.3f) last bin (%.3f)", hData->GetName(), hData->GetNbinsX(), hData->GetXaxis()->GetBinLowEdge(1), hData->GetXaxis()->GetBinUpEdge(hData->GetNbinsX())), cyanTxt);
  
  //Single functions
  Int_t ntot = 0;
  for(Int_t c = 0; c<ntemplates; c++){
    const TString cname = func->At(c)->ClassName();
    if(cname.Contains("TF1")){
      ntot++;
      singlefun[c] = (TF1*)func->At(c);
      Double_t r[2];
      singlefun[c]->GetRange(r[0], r[1]);
      Int_t n = singlefun[c]->GetNpar();
      nparameters[c] = n;
      parameters[c] = new Double_t[n];
      Infomsgcolor("PerformFitWithFunctions", Form("Template function #%i/%i summary: name (%s)  first bin (%.3f) last bin (%.3f) #parameters %i", c+1, ntemplates, singlefun[c]->GetName(), r[0], r[1], n), cyanTxt);
    }
    else Fatalmsg("PerformFitWithFunctions", Form("Input class %s for template %i is not recognized!", cname.Data(), c));
  }
  
  if(ntot != ntemplates) Fatalmsg("PerformFitWithFunctions", "Different number of templates than expected");
  
  //Sum functions
  if(funcsum){
    TF1 * f = funcsum;
    Double_t r[2];
    f->GetRange(r[0], r[1]);
    Int_t nsum = 0;
    for(Int_t c = 0; c<ntemplates; c++) nsum += nparameters[c];
    if(nparameterssum != nsum) Fatalmsg("PerformFitWithFunctions", "Different number of parameters");
    Infomsgcolor("PerformFitWithFunctions", Form("Template function Sum summary: name (%s)  first bin (%.3f) last bin (%.3f) #parameters %i", f->GetName(), r[0], r[1], nparameterssum), cyanTxt);
  }
  else Fatalmsg("PerformFitWithFunctions", Form("Cannot find input sum function!"));
  
  Int_t StopSign = 0;//Check if the fit is feasible
  
  
  if(hData->Integral() < min){//Check on data
    Warningmsg("PerformFitWithFunctions", Form("%s does NOT have sufficient entries (%.0f/%.0f) but it has %.0f entries", hData->GetName(), hData->Integral(), min, hData->GetEntries()));
    StopSign  = -1; //Check on data
  }
  
  if(StopSign < 0) return kFALSE;
  
  TVirtualPad *currpad = gPad;
  TCanvas *cFit = 0x0;
  if(showfit){//TCanvas for debug purposes
    TLine * linemin = new TLine(minx, 0, minx, hData->GetMaximum());
    TLine * linemax = new TLine(maxx, 0, maxx, hData->GetMaximum());
    linemin->SetLineColor(kRed);
    linemax->SetLineColor(kRed);
    
    cFit = new TCanvas("TemplateRooFit", "Fit with functions");
    cFit->Divide(2);
    cFit->cd(1);
    gPad->SetLogy();
    hData->DrawCopy()->GetXaxis()->SetRangeUser(minx, maxx);
    for(Int_t c = 0; c<ntemplates; c++) singlefun[c]->DrawCopy("same");
    funcsum->DrawCopy("same");
    linemin->Draw();
    linemax->Draw();
    
    cFit->cd(2);
    gPad->SetLogy();
    hData->DrawNormalized()->GetXaxis()->SetRangeUser(minx, maxx);
    for(Int_t c = 0; c<ntemplates; c++) singlefun[c]->DrawCopy("same");
    funcsum->DrawCopy("same");
    linemin->Draw();
    linemax->Draw();
  }
  currpad->cd();
  
  Double_t frac[ntemplates];
  Double_t fracerr[ntemplates];
  
  if(range[0] >= range[1]) Fatalmsg("PerformFitWithFunctions", "Range is not well defined");
  
  TH1F *datahisto = (TH1F*) hData->Clone(Form("tobefitted%s", hData->GetName()));
  
  const TString fitopt = "0 M S"/*"0 M L S"*/;
  //Fit with single functions
  for(Int_t c = 0; c<ntemplates; c++){
    TFitResultPtr fitresults = datahisto->Fit(singlefun[c], fitopt, "", fitrange[0], fitrange[1]);
    if(fitresults) Infomsg("PerformFitWithFunctions", "Fit correctly performed");
    else continue;
    for(Int_t i = 0; i < nparameters[c]; i++){
      Infomsg("PerformFitWithFunctions", Form("Getting Parameter %i: %s -> %f", i, singlefun[c]->GetParName(i), singlefun[c]->GetParameter(i)));
      parameters[c][i] = singlefun[c]->GetParameter(i);
    }
  }
  
  //Fit with sum function
  if(funcsum){
    //Setting start parameters from other fits
    for(Int_t c = 0; c<ntemplates; c++){
      for(Int_t i = 0; i<nparameters[c]; i++){
        funcsum->SetParameter(c*nparameters[c] + i, parameters[c][i]);
      }
    }
    for(Int_t i = 0; i<funcsum->GetNpar(); i++) cout<<"Sum function prepared with par "<<i<<" ("<<funcsum->GetParName(i)<<") "<<funcsum->GetParameter(i)<<endl;
    
    TFitResultPtr fitresults = datahisto->Fit(funcsum, fitopt, "", fitrange[0], fitrange[1]);
    if(fitresults){
      //Setting start parameters from convoluted fit
      for(Int_t c = 0; c<ntemplates; c++){
        for(Int_t i = 0; i<nparameters[c]; i++){
          singlefun[c]->SetParameter(i, funcsum->GetParameter(c*nparameters[c] + i));
          singlefun[c]->SetParError(i, funcsum->GetParError(c*nparameters[c] + i));
        }
      }
      Infomsg("PerformFitWithFunctions", "Fit correctly performed");
      
    }
    else return kFALSE;
    
    for(Int_t i = 0; i < nparameterssum; i++){
      Infomsg("PerformFitWithFunctions", Form("Getting Parameter %i: %s -> %f", i, funcsum->GetParName(i), funcsum->GetParameter(i)));
    }
    
  }
  
  TString f;
  prediction->Add(ConvertFunToHisto(funcsum, hData, Form("funcpredictiontotal%s", hData->GetName())));
  const Double_t width = datahisto->GetXaxis()->GetBinWidth(datahisto->GetXaxis()->FindBin((range[0] + range[1])/2));
  for(Int_t temp = 0; temp < ntemplates; temp++){
    prediction->Add(ConvertFunToHisto(singlefun[temp], hData, Form("funcprediction%s%i", hData->GetName(), temp)));
    fraction[temp] = singlefun[temp]->Integral(range[0], range[1])/width;
    fractionErr[temp] = singlefun[temp]->IntegralError(range[0], range[1])/width;
    
    f += Form(" Template %i = %f(%f)", temp, fraction[temp], fractionErr[temp]);
  }
  Infomsgcolor("PerformFitWithFunctions", Form("Fractions from fit     :  %s\n", f.Data()), cyanTxt);
  cout<<"Entries : "<<hData->GetEntries()<<endl;
  
  return kTRUE;
  #else 
  return kFALSE;
  #endif  
}

//_________________________________________________________________________________________________
Bool_t PerformFitWithCD(TH1F* hData, TObjArray* mc, Double_t *range, Double_t *fitrange, TArrayD &fraction, TArrayD &fractionErr, TObjArray *& prediction){
  #ifdef USECDECONVOLUTION
  const Int_t ntemplates = mc->GetEntries();
  Double_t min = 100;
  
  Infomsgcolor("PerformFitWithCD", Form("Data histogram summary: name (%s)  nbins (%i) first bin (%.3f) last bin (%.3f)", hData->GetName(), hData->GetNbinsX(), hData->GetXaxis()->GetBinLowEdge(1), hData->GetXaxis()->GetBinUpEdge(hData->GetNbinsX())), cyanTxt);
  for(Int_t c = 0; c<ntemplates; c++) Infomsgcolor("PerformFitWithCD", Form("Template histogram #%i/%i summary: name (%s)  nbins (%i) first bin (%.3f) last bin (%.3f)", c+1, ntemplates, ((TH1F*)mc->At(c))->GetName(), ((TH1F*)mc->At(c))->GetNbinsX(), ((TH1F*)mc->At(c))->GetXaxis()->GetBinLowEdge(1), ((TH1F*)mc->At(c))->GetXaxis()->GetBinUpEdge(((TH1F*)mc->At(c))->GetNbinsX())), cyanTxt);
  
  AliCDeconv* fit = new AliCDeconv(); // initialise
  fit->SetData(hData);
  fit->SetTemplates(mc);
  
  Int_t StopSign = 0;//Check if the fit is feasible
  
  TH1F *htemplates[ntemplates];
  for(Int_t temp = 0; temp < ntemplates; temp++){//Check on templates
    htemplates[temp] = (TH1F*) mc->At(temp)->Clone(Form("Input template %i", temp));
    htemplates[temp]->SetTitle(Form("Input template %i", temp));
    htemplates[temp]->SetLineStyle(2);
    if(htemplates[temp]->Integral() < min){
      Warningmsg("PerformFitWithCD", Form("%s does NOT have sufficient entries (%.0f/%.0f) but it has %.0f entries", htemplates[temp]->GetName(), htemplates[temp]->Integral(), min, htemplates[temp]->GetEntries()));
      StopSign = -1;
    }
  }
  if(hData->Integral() < min){
    Warningmsg("PerformFitWithCD", Form("%s does NOT have sufficient entries (%.0f/%.0f) but it has %.0f entries", hData->GetName(), hData->Integral(), min, hData->GetEntries()));
    StopSign  = -1; //Check on data
  }
  
  //fit->SetRangeX(200, 1000);                    // use only the first 15 bins in the fit
  //   if(fitrange[0] < fitrange[1]) fit->SetRangeX(hData->GetXaxis()->FindBin(fitrange[0]), hData->GetXaxis()->FindBin(fitrange[1]));
  Int_t status = -1;
  if(StopSign == 0){
    Infomsgcolor("PerformFitWithCD", Form("I've checked that the fit is feasible (non empty histograms) ---> Fitting!!"), cyanTxt);
    status = 0;
    fit->Run(kFALSE);               // perform the fit
  }
  else{
    Warningmsg("PerformFitWithCD", "Fit not performed since minimum requirements are not fulfilled:\nTemplates: ");
    for(Int_t c = 0; c<ntemplates; c++) cout<<" htemplates"<<c<<" ("<<htemplates[c]->Integral()<<")"<<flush;
    cout<<" hData ("<<hData->Integral()<<")"<<endl;
  }
  
  Infomsgcolor("PerformFitWithCD", Form("Fit status: %i (0 means success!)", status), cyanTxt);
  
  Double_t value[ntemplates], error[ntemplates];//Fit results and errors
  Double_t yield[ntemplates], yielderror[ntemplates];//Fit results and errors
  
  //Drawing section for data and templates
  hData->SetLineColor(1);
  hData->SetTitle("Data");
  hData->SetName("Data");
  hData->DrawCopy("Ep");
  for(Int_t temp = 0; temp < ntemplates; temp++) htemplates[temp]->DrawCopy("same");
  
  if(status == 0){                       // check on fit status
    
    //Get the fit result
    TH1F* result = (TH1F*) fit->GeneratePrediction();
    //Get the fitted template
    TH1F* MCPred[ntemplates];
    for(Int_t temp = 0; temp < ntemplates; temp++){
      MCPred[temp] = (TH1F*)mc->At(temp)->Clone(Form("temp%i", temp));
      MCPred[temp]->SetLineColor(htemplates[temp]->GetLineColor());
    }
    
    //Drawing section after FIT
    result->SetTitle("Fit result");
    result->SetName("Fit result");
    result->DrawCopy("same");
    
    //Getting results section
    Double_t factor, dataIntegral, tempIntegral, sum = 0;
    
    dataIntegral = hData->GetSumOfWeights();
    TArrayD solutions(ntemplates);
    fit->GetSolutions(solutions);
    
    for(Int_t temp = 0; temp < ntemplates; temp++){
      value[temp] = solutions[temp];
      error[temp] = 1;
      tempIntegral = MCPred[temp]->GetSumOfWeights();
      sum += tempIntegral*value[temp];
      factor = dataIntegral/tempIntegral*value[temp];
      MCPred[temp]->SetLineStyle(2);
      MCPred[temp]->SetTitle(Form("Fitted Template %i unscaled", temp));
      MCPred[temp]->SetName(Form("Fitted Template %i unscaled", temp));
      //       MCPred[temp]->DrawCopy("same");
      MCPred[temp]->SetLineStyle(1);
      MCPred[temp]->Scale(factor);
      MCPred[temp]->SetTitle(Form("Fitted Template %i", temp));
      MCPred[temp]->SetName(Form("Fitted Template %i", temp));
      MCPred[temp]->DrawCopy("same");
      prediction->Add((TH1F *) MCPred[temp]->Clone(Form("prediction%i", temp)));
      Infomsgcolor("PerformFitWithCD", Form("Fraction %i: %f (%f) scaling factor: %f / %f * %f  = %f", temp, value[temp], error[temp], dataIntegral, tempIntegral, value[temp], factor), yellowTxt);
      
      if(range[0] >= range[1]) Fatalmsg("PerformFitWithCD", "Integration range badly defined");
      
      Int_t binlow = MCPred[temp]->FindBin((1-0.0001)*range[0]);
      Int_t binup = MCPred[temp]->FindBin((1-0.0001)*range[1]);
      if(temp == 0){
        Infomsgcolor("PerformFitWithCD", Form("The sum of the integrals of the %i functions from the fit is %f while the Integral of the data is %f, the difference is %f", ntemplates, sum, dataIntegral, sum-dataIntegral), yellowTxt);
        
        Infomsgcolor("PerformFitWithCD", Form("Now obtaining the total yield in the range [%f, %f]", MCPred[temp]->GetXaxis()->GetBinLowEdge(binlow), MCPred[temp]->GetXaxis()->GetBinUpEdge(binup)), cyanTxt);
      }
      yield[temp] = MCPred[temp]->IntegralAndError(binlow, binup, yielderror[temp]);
      Infomsgcolor("PerformFitWithCD", Form("Template %i from Integral: %f (%f)", temp, yield[temp], yielderror[temp]), yellowTxt);
    }
    
  }
  else{
    Warningmsg("PerformFitWithCD", Form("Fit not performed successfully so all fractions are set equal"));
    for(Int_t temp = 0; temp < ntemplates; temp++) yield[temp] = 1;
  }
  
  if(fit){
    //     delete fit;
    fit = 0;
  }
  
  Double_t Den = 0, eDen = 0, Num = 0, eNum = 0, Ratio = 0, RatioErr = 0;
  for(Int_t temp = 0; temp < ntemplates; temp++){
    Den += yield[temp];
    eDen += yielderror[temp]*yielderror[temp];
  }
  eDen = TMath::Sqrt(eDen);
  
  if(Den > 0){
    
    for(Int_t temp = 0; temp < ntemplates; temp++){
      Num = yield[temp];
      eNum = yielderror[temp];
      Ratio = Num/Den;
      
      if(Num>0) RatioErr = Ratio*TMath::Sqrt( eNum*eNum/(Num*Num) + eDen*eDen/(Den*Den));
      else RatioErr = 0;
      //       fraction[temp] = Ratio;
      //       fractionErr[temp] = RatioErr;
      fraction[temp] = yield[temp];
      fractionErr[temp] = yielderror[temp];
    }
  }
  
  TString f[2];
  for(Int_t temp = 0; temp < ntemplates; temp++){
    f[0] += Form(" Template %i = %f(%f)", temp, fraction[temp], fractionErr[temp]);
    f[1] += Form(" Template %i = %f(%f)", temp, value[temp], error[temp]);
  }
  
  Infomsgcolor("PerformFitWithCD", Form("Fractions from Integral: %s", f[0].Data()), cyanTxt);
  Infomsgcolor("PerformFitWithCD", Form("Fractions from fit     :  %s\n", f[1].Data()), cyanTxt);
  
  if(status == 0) return kTRUE;
  else return kFALSE;
  #else 
  return kFALSE;
  #endif
}

//_________________________________________________________________________________________________
Bool_t UseBinCounting(TH1F * hData, TObjArray* mc, Double_t *rangelol, TArrayD &fraction, TArrayD &fractionErr){
  const Int_t ntemplates = mc->GetEntries();
  
  const Double_t range[2] = {-240, 240};
  const Double_t sig = hData->Integral(hData->GetXaxis()->FindBin(range[0]), hData->GetXaxis()->FindBin(range[1]));
  const Double_t bkg_S = hData->Integral(hData->GetXaxis()->FindBin(2.*range[0]), hData->GetXaxis()->FindBin(range[0]));
  const Double_t bkg_D = hData->Integral(hData->GetXaxis()->FindBin(range[1]), hData->GetXaxis()->FindBin(2*range[1]));
  const Double_t bkg = (bkg_S/TMath::Abs(range[0]) + bkg_D/TMath::Abs(range[1]))*(range[1] - range[0]);
  
  for(Int_t i = 0; i < ntemplates; i++){
    fraction[i] = sig - bkg;
    fractionErr[i] = 0;
  }
  return kFALSE;
}
#endif
