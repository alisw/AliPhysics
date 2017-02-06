/**************************************************************************
 * Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
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

#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TVirtualFitter.h>
#include <TDatabasePDG.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>

#include "AliLog.h"
#include "AliVertexingHFUtils.h"
#include "AliHFInvMassFitter.h"

/// \cond CLASSIMP
ClassImp(AliHFInvMassFitter);
/// \endcond

/////////////////////////////////////////////////////////////
///
/// Implemenatation od AliHFInvMassFitter class for
/// the fit of invariant mass distribution of charm hadron candidates
/// reconstructed from their hadronic decays
///
/// Author: F. Prino
///   -> Simplified version of AliHFMassFitter class developed by C.Bianchin
///   -> With extra features:
///       - Polynomials with degree >2 as background fit function
///              (as implemented by A. Rossi in AliHFMassFitterVAR
///       - Possibility to incude a Gaussian function in the background
///              (for D+ -> KKpi decays in the D_s background)
///
/////////////////////////////////////////////////////////////

//__________________________________________________________________________
AliHFInvMassFitter::AliHFInvMassFitter() : 
  TNamed(),
  fHistoInvMass(0x0),
  fMinMass(0),
  fMaxMass(5),
  fTypeOfFit4Bkg(kExpo),
  fPolDegreeBkg(4),
  fCurPolDegreeBkg(-1),
  fMassParticle(1.864),
  fTypeOfFit4Sgn(kGaus),
  fMass(1.865),
  fMassErr(0.),
  fSigmaSgn(0.012),
  fSigmaSgnErr(0.),
  fFixedMean(kFALSE),
  fFixedSigma(kFALSE),
  fFixedRawYield(-1.),
  fNSigPars(3),
  fNBkgPars(2),
  fOnlySideBands(kFALSE),
  fFitOption("L,E"),
  fRawYield(0.),
  fRawYieldErr(0.),
  fSigFunc(0x0),
  fBkgFuncSb(0x0),
  fBkgFunc(0x0),
  fBkgFuncRef(0x0),
  fSecondPeak(kFALSE),
  fSecMass(-999.),
  fSecWidth(9999.),
  fFixSecMass(kFALSE),
  fFixSecWidth(kFALSE),
  fSecFunc(0x0),
  fFuncTot(0x0)
{
  /// default constructor
}

//__________________________________________________________________________
AliHFInvMassFitter::AliHFInvMassFitter(const TH1F *histoToFit, Double_t minvalue, Double_t maxvalue, Int_t fittypeb, Int_t fittypes):
  TNamed(),
  fHistoInvMass(0x0),
  fMinMass(minvalue),
  fMaxMass(maxvalue),
  fTypeOfFit4Bkg(fittypeb),
  fPolDegreeBkg(4),
  fCurPolDegreeBkg(-1),
  fMassParticle(1.864),
  fTypeOfFit4Sgn(fittypes),
  fMass(1.865),
  fMassErr(0.),
  fSigmaSgn(0.012),
  fSigmaSgnErr(0.),
  fFixedMean(kFALSE),
  fFixedSigma(kFALSE),
  fFixedRawYield(-1.),
  fNSigPars(3),
  fNBkgPars(2),
  fOnlySideBands(kFALSE),
  fFitOption("L,E"),
  fRawYield(0.),
  fRawYieldErr(0.),
  fSigFunc(0x0),
  fBkgFuncSb(0x0),
  fBkgFunc(0x0),
  fBkgFuncRef(0x0),
  fSecondPeak(kFALSE),
  fSecMass(-999.),
  fSecWidth(9999.),
  fFixSecMass(kFALSE),
  fFixSecWidth(kFALSE),
  fSecFunc(0x0),
  fFuncTot(0x0)
{
  /// standard constructor
  fHistoInvMass=(TH1F*)histoToFit->Clone("fHistoInvMass");
  fHistoInvMass->SetDirectory(0);
  SetNumberOfParams();
}
//_________________________________________________________________________
AliHFInvMassFitter::~AliHFInvMassFitter() {

  ///destructor

  delete fHistoInvMass;
  delete fSigFunc;
  delete fBkgFunc;
  delete fBkgFuncRef;
  delete fSecFunc;
  delete fFuncTot;
  
}
//__________________________________________________________________________
void AliHFInvMassFitter::SetNumberOfParams(){
  /// Configure number of parameters of fit functions
  ///

  switch (fTypeOfFit4Bkg) {
  case 0:
    fNBkgPars=2;
    break;
  case 1:
    fNBkgPars=2;
    break;
  case 2:
    fNBkgPars=3;
    break;
  case 3:
    fNBkgPars=1;
    break;
  case 4:
    fNBkgPars=2;	
    break;
  case 5:
    fNBkgPars=3;	
    break;
  case 6:
    fNBkgPars=fPolDegreeBkg+1;
    break;
  default:
    AliError("Error in computing fNBkgPars: check fTypeOfFit4Bkg");
    break;
  }
  switch (fTypeOfFit4Sgn) {
  case 0:
    fNSigPars=3;
    break;
  case 1:
    fNSigPars=5;
    break;
  default:
    AliError("Error in computing fNSigPars: check fTypeOfFit4Sgn");
    break;
  }

}
//__________________________________________________________________________
Bool_t AliHFInvMassFitter::MassFitter(Bool_t draw){
  /// Main function to fit the invariant mass distribution
  ///

  TVirtualFitter::SetDefaultFitter("Minuit");

  Double_t integralHisto=fHistoInvMass->Integral(fHistoInvMass->FindBin(fMinMass),fHistoInvMass->FindBin(fMaxMass),"width");

  fOnlySideBands = kTRUE;
  fBkgFuncSb = CreateBackgroundFitFunction("funcbkgsb",integralHisto);
  Int_t status=-1;
  if(fTypeOfFit4Bkg==6){
    if(PrepareHighPolFit(fBkgFuncSb)){
      fHistoInvMass->GetListOfFunctions()->Add(fBkgFuncSb);
      fHistoInvMass->GetFunction(fBkgFuncSb->GetName())->SetBit(1<<9,kTRUE);
      status=0;
    }
  }
  else status=fHistoInvMass->Fit("funcbkgsb",Form("R,%s,+,0",fFitOption.Data()));
  fBkgFuncSb->SetLineColor(kGray+1);
  if (status != 0){
    printf("Failed first fit with only background, minuit status = %d\n",status);
    return kFALSE;
  }

  fOnlySideBands = kFALSE;
  if(!fBkgFunc){
    fBkgFunc = CreateBackgroundFitFunction("funcbkg",integralHisto);
    for(Int_t ipar=0; ipar<fNBkgPars; ipar++) fBkgFunc->SetParameter(ipar,fBkgFuncSb->GetParameter(ipar));
  }
  fBkgFunc->SetLineColor(kGray+1);

  Double_t estimSignal=CheckForSignal(fMass,fSigmaSgn);
  if(estimSignal<0.){
    if(draw) DrawFit();
    return kTRUE;
  }
  if(!fBkgFuncRef){
    fBkgFuncRef = CreateBackgroundFitFunction("funcbkgrefit",integralHisto);
    for(Int_t ipar=0; ipar<fNBkgPars; ipar++) fBkgFuncRef->SetParameter(ipar,fBkgFunc->GetParameter(ipar));
  }
  fBkgFuncRef->SetLineColor(2);
  fSigFunc = CreateSignalFitFunction("fsigfit",estimSignal);
  if(fSecondPeak){
    Double_t estimSec=CheckForSignal(fSecMass,fSecWidth);
    fSecFunc = CreateSecondPeakFunction("fsecpeak",estimSec);
  }
  fFuncTot = CreateTotalFitFunction("funcmass");
  status=fHistoInvMass->Fit("funcmass",Form("R,%s,+,0",fFitOption.Data()));
  if (status != 0){
    printf("Failed fit with signal+background, minuit status = %d\n",status);
    return kFALSE;
  }
  for(Int_t ipar=0; ipar<fNBkgPars; ipar++) fBkgFuncRef->SetParameter(ipar,fFuncTot->GetParameter(ipar));
  for(Int_t ipar=0; ipar<fNSigPars; ipar++) fSigFunc->SetParameter(ipar,fFuncTot->GetParameter(ipar+fNBkgPars));
  fMass=fSigFunc->GetParameter(1);
  fMassErr=fSigFunc->GetParError(1);
  fSigmaSgn=fSigFunc->GetParameter(2);
  fSigmaSgnErr=fSigFunc->GetParError(2);
  fFuncTot->SetLineColor(4);
  fRawYield=fFuncTot->GetParameter(fNBkgPars)/fHistoInvMass->GetBinWidth(1);
  fRawYieldErr=fFuncTot->GetParError(fNBkgPars)/fHistoInvMass->GetBinWidth(1);
  if(draw) DrawFit();
  return kTRUE;
}
//______________________________________________________________________________
void AliHFInvMassFitter::DrawFit(){
  /// Steering method to draw the fit output
  ///

  TCanvas* c0=new TCanvas("c0");
  DrawHere(c0);
}
//______________________________________________________________________________
void AliHFInvMassFitter::DrawHere(TVirtualPad* c){
  /// Core method to draw the fit output
  ///

  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  c->cd();
  fHistoInvMass->GetXaxis()->SetRangeUser(fMinMass,fMaxMass);
  fHistoInvMass->SetMarkerStyle(20);
  fHistoInvMass->SetMinimum(0.);
  fHistoInvMass->Draw("PE");
  if(fBkgFunc) fBkgFunc->Draw("same");
  if(fBkgFuncRef) fBkgFuncRef->Draw("same");
  if(fFuncTot) fFuncTot->Draw("same");
  TPaveText *pinfos=new TPaveText(0.12,0.65,0.47,0.89,"NDC");
  TPaveText *pinfom=new TPaveText(0.6,0.7,1.,.87,"NDC");
  pinfos->SetBorderSize(0);
  pinfos->SetFillStyle(0);
  pinfom->SetBorderSize(0);
  pinfom->SetFillStyle(0);
  if(fFuncTot){
    pinfom->SetTextColor(kBlue);
    for(Int_t ipar=1; ipar<fNSigPars; ipar++){
      pinfom->AddText(Form("%s = %.3f #pm %.3f",fFuncTot->GetParName(ipar+fNBkgPars),fFuncTot->GetParameter(ipar+fNBkgPars),fFuncTot->GetParError(ipar+fNBkgPars)));
    }
    pinfom->Draw();
 
    Double_t nsigma=3;
    Double_t bkg,errbkg;
    Background(nsigma,bkg,errbkg);
    Double_t signif,errsignif;
    Significance(nsigma,signif,errsignif);

    pinfos->AddText(Form("S = %.0f #pm %.0f ",fRawYield,fRawYieldErr));
    pinfos->AddText(Form("B (%.0f#sigma) = %.0f #pm %.0f",nsigma,bkg,errbkg));
    pinfos->AddText(Form("S/B = (%.0f#sigma) %.4f ",nsigma,fRawYield/bkg)); 
    pinfos->AddText(Form("Signif (%.0f#sigma) = %.1f #pm %.1f ",nsigma,signif,errsignif));
    pinfos->Draw();
  }
  c->Update();
  return;
}
//______________________________________________________________________________
Double_t AliHFInvMassFitter::CheckForSignal(Double_t mean, Double_t sigma){
  /// Checks if there are signal counts above the background
  ///   in the invariant mass region of the peak

  Double_t minForSig=mean-4.*sigma;
  Double_t maxForSig=mean+4.*sigma;
  Int_t binForMinSig=fHistoInvMass->FindBin(minForSig);
  Int_t binForMaxSig=fHistoInvMass->FindBin(maxForSig);
  Double_t sum=0.;
  Double_t sumback=0.;
  fBkgFunc->Print();
  for(Int_t ibin=binForMinSig; ibin<=binForMaxSig; ibin++){
    sum+=fHistoInvMass->GetBinContent(ibin);
    sumback+=fBkgFunc->Eval(fHistoInvMass->GetBinCenter(ibin));
  }
  Double_t diffUnderPeak=(sum-sumback);
  printf("intUnderFunc=%f  intUnderHisto=%f   sestimSig=%f\n",sum,sumback,diffUnderPeak);
  if(diffUnderPeak/TMath::Sqrt(sum)<1.){
    printf("(Tot-Bkg)/sqrt(Tot)=%f ---> Likely no signal/\n",diffUnderPeak/TMath::Sqrt(sum));
    return -1;
  }
  return diffUnderPeak*fHistoInvMass->GetBinWidth(1);
}

//______________________________________________________________________________
TF1* AliHFInvMassFitter::CreateBackgroundFitFunction(TString fname, Double_t integral){
  /// Creates the background fit fucntion
  ///

  SetNumberOfParams();
  TF1* funcbkg =  new TF1(fname.Data(),this,&AliHFInvMassFitter::FitFunction4Bkg,fMinMass,fMaxMass,fNBkgPars,"AliHFInvMassFitter","FitFunction4Bkg");
  switch (fTypeOfFit4Bkg) {
  case 0: //gaus+expo
    funcbkg->SetParNames("BkgInt","Slope"); 
    funcbkg->SetParameters(integral,-2.); 
    break;
  case 1:
    funcbkg->SetParNames("BkgInt","Slope");
    funcbkg->SetParameters(integral,-100.); 
    break;
  case 2:
    funcbkg->SetParNames("BkgInt","Coef1","Coef2");
    funcbkg->SetParameters(integral,-10.,5);
    break;
  case 3:
    funcbkg->SetParNames("Const");
    funcbkg->SetParameter(0,0.);
    funcbkg->FixParameter(0,0.);
    break;
  case 4:     
    funcbkg->SetParNames("BkgInt","Coef1");
    funcbkg->SetParameters(integral,0.5);
    break;
  case 5:    
    funcbkg->SetParNames("BkgInt","Coef1","Coef2");
    funcbkg->SetParameters(integral,-10.,5.);
    break;
  case 6:
    for(Int_t j=0;j<fNBkgPars; j++){
      funcbkg->SetParName(j,Form("Coef%d",j));
      funcbkg->SetParameter(j,0);
    }
    funcbkg->SetParameter(0,integral);
    break;
  default:
    AliError(Form("Wrong choice of fTypeOfFit4Bkg (%d)",fTypeOfFit4Bkg));
    delete funcbkg;
    return 0x0;
    break;
  }
  funcbkg->SetLineColor(kBlue+3); 
  return funcbkg;
}
//______________________________________________________________________________
TF1* AliHFInvMassFitter::CreateSecondPeakFunction(TString fname, Double_t integsig){
  /// Creates a function for a gaussian peak in the background fit function
  /// Can be used e.g. to include the D+->KKpi peak in the D_s inv. mass fit

  TF1* funcsec =  new TF1(fname.Data(),this,&AliHFInvMassFitter::FitFunction4SecPeak,fMinMass,fMaxMass,3,"AliHFInvMassFitter","FitFunction4SecPeak");
  funcsec->SetParameter(0,integsig);
  funcsec->SetParameter(1,fSecMass);
  if(fFixSecMass) funcsec->FixParameter(1,fSecMass);
  funcsec->SetParameter(2,fSecWidth);
  if(fFixSecWidth) funcsec->FixParameter(2,fSecWidth);
  funcsec->SetParNames("SecPeakInt","SecMean","SecSigma");
  return funcsec;
}

//______________________________________________________________________________
TF1* AliHFInvMassFitter::CreateSignalFitFunction(TString fname, Double_t integsig){
  /// Creates the fit function fot he signal peak
  ///

  SetNumberOfParams();
  TF1* funcsig =  new TF1(fname.Data(),this,&AliHFInvMassFitter::FitFunction4Sgn,fMinMass,fMaxMass,fNSigPars,"AliHFInvMassFitter","FitFunction4Sgn");
  if(fTypeOfFit4Sgn==kGaus){
    funcsig->SetParameter(0,integsig);
    if(fFixedRawYield>-0.1) funcsig->FixParameter(0,fFixedRawYield);
    funcsig->SetParameter(1,fMass);
    if(fFixedMean) funcsig->FixParameter(1,fMass);
    funcsig->SetParameter(2,fSigmaSgn);
    if(fFixedSigma) funcsig->FixParameter(2,fSigmaSgn);
    funcsig->SetParNames("SgnInt","Mean","Sigma");
   }
  if(fTypeOfFit4Sgn==k2Gaus){
    funcsig->SetParameter(0,integsig);
    if(fFixedRawYield>-0.1) funcsig->FixParameter(0,fFixedRawYield);
    funcsig->SetParameter(1,fMass);
    if(fFixedMean) funcsig->FixParameter(1,fMass);
    funcsig->SetParameter(2,fSigmaSgn);
    funcsig->SetParLimits(2,0.004,0.05);
    if(fFixedSigma) funcsig->FixParameter(2,fSigmaSgn);
    funcsig->SetParameter(3,0.2);
    funcsig->SetParLimits(3,0.,1.);
    funcsig->SetParameter(4,fSigmaSgn);
    funcsig->SetParLimits(4,0.004,0.05);
    funcsig->SetParNames("SgnInt","Mean","Sigma1","Frac","Sigma2");
  }
  return funcsig;
}

//______________________________________________________________________________
TF1* AliHFInvMassFitter::CreateTotalFitFunction(TString fname){
  /// Creates the total fit fucntion (signal+background+possible second peak)
  ///

  SetNumberOfParams();
  Int_t nParSecPeak=0;
  if(fSecondPeak && fSecFunc) nParSecPeak=3;
  TF1* ftot=ftot=new TF1(fname.Data(),this,&AliHFInvMassFitter::FitFunction4Mass,fMinMass,fMaxMass,(fNBkgPars+nParSecPeak+fNSigPars),"AliHFInvMassFitter","FitFunction4Mass");
  for(Int_t ipar=0; ipar<fNBkgPars; ipar++){ 
    ftot->SetParameter(ipar,fBkgFunc->GetParameter(ipar));
    ftot->SetParName(ipar,fBkgFunc->GetParName(ipar));
  }
  printf("Copy parameters for signal = %d\n",fNSigPars);
  for(Int_t ipar=0; ipar<fNSigPars; ipar++){
    ftot->SetParameter(ipar+fNBkgPars,fSigFunc->GetParameter(ipar));
    ftot->SetParName(ipar+fNBkgPars,fSigFunc->GetParName(ipar));
    Double_t parmin,parmax;
    fSigFunc->GetParLimits(ipar,parmin,parmax);
    ftot->SetParLimits(ipar+fNBkgPars,parmin,parmax);
    printf("Par %d val %f vs %f\n",ipar,fSigFunc->GetParameter(ipar),ftot->GetParameter(ipar+fNBkgPars));
  }
  if(fSecondPeak && fSecFunc){
    for(Int_t ipar=0; ipar<nParSecPeak; ipar++){ 
      ftot->SetParameter(ipar+fNBkgPars+fNSigPars,fSecFunc->GetParameter(ipar));
      ftot->SetParName(ipar+fNBkgPars+fNSigPars,fSecFunc->GetParName(ipar));    
      Double_t parmin,parmax;
      fSecFunc->GetParLimits(ipar,parmin,parmax);
      ftot->SetParLimits(ipar+fNBkgPars+fNSigPars,parmin,parmax);
    }
  }
  return ftot;
}
//__________________________________________________________________________
Double_t AliHFInvMassFitter::FitFunction4Bkg (Double_t *x, Double_t *par){
  /// Fit function for the background
  ///

  Double_t maxDeltaM = 4.*fSigmaSgn;
  if(fOnlySideBands && TMath::Abs(x[0]-fMass) < maxDeltaM) {
    TF1::RejectPoint();
    return 0;
  }
  if(fOnlySideBands && fSecondPeak && TMath::Abs(x[0]-fSecMass) < (4.*fSecWidth)){
    TF1::RejectPoint();
    return 0;
  }
  Double_t total=0;

  switch (fTypeOfFit4Bkg){
  case 0:
    //exponential
    //exponential = A*exp(B*x) -> integral(exponential)=A/B*exp(B*x)](min,max)
    //-> A = B*integral/(exp(B*max)-exp(B*min)) where integral can be written
    //as integralTot- integralGaus (=par [2])
    //Par:
    // * [0] = integralBkg;
    // * [1] = B;
    //exponential = [1]*[0]/(exp([1]*max)-exp([1]*min))*exp([1]*x)
    total = par[0]*par[1]/(TMath::Exp(par[1]*fMaxMass)-TMath::Exp(par[1]*fMinMass))*TMath::Exp(par[1]*x[0]);
    //    AliInfo("Background function set to: exponential");
    break;
  case 1:
    //linear
    //y=a+b*x -> integral = a(max-min)+1/2*b*(max^2-min^2) -> a = (integral-1/2*b*(max^2-min^2))/(max-min)=integral/(max-min)-1/2*b*(max+min)
    // * [0] = integralBkg;
    // * [1] = b;
    total= par[0]/(fMaxMass-fMinMass)+par[1]*(x[0]-0.5*(fMaxMass+fMinMass));
    //    AliInfo("Background function set to: linear");
    break;
  case 2:
    //parabola
    //y=a+b*x+c*x**2 -> integral = a(max-min) + 1/2*b*(max^2-min^2) +
    //+ 1/3*c*(max^3-min^3) -> 
    //a = (integral-1/2*b*(max^2-min^2)-1/3*c*(max^3-min^3))/(max-min)
    // * [0] = integralBkg;
    // * [1] = b;
    // * [2] = c;
    total = par[0]/(fMaxMass-fMinMass)+par[1]*(x[0]-0.5*(fMaxMass+fMinMass))+par[2]*(x[0]*x[0]-1/3.*(fMaxMass*fMaxMass*fMaxMass-fMinMass*fMinMass*fMinMass)/(fMaxMass-fMinMass));
    //    AliInfo("Background function set to: polynomial");
    break;
  case 3:
    total=par[0];
    break;
  case 4:  
    //power function 
    //y=a(x-m_pi)^b -> integral = a/(b+1)*((max-m_pi)^(b+1)-(min-m_pi)^(b+1))
    //
    //a = integral*(b+1)/((max-m_pi)^(b+1)-(min-m_pi)^(b+1))
    // * [0] = integralBkg;
    // * [1] = b;
    // a(power function) = [0]*([1]+1)/((max-m_pi)^([1]+1)-(min-m_pi)^([1]+1))*(x-m_pi)^[1]
    {
    Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    total = par[0]*(par[1]+1.)/(TMath::Power(fMaxMass-mpi,par[1]+1.)-TMath::Power(fMinMass-mpi,par[1]+1.))*TMath::Power(x[0]-mpi,par[1]);
    //    AliInfo("Background function set to: powerlaw");
    }
    break;
  case 5:
   //power function wit exponential
    //y=a*Sqrt(x-m_pi)*exp(-b*(x-m_pi))  
    { 
    Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    total = par[1]*TMath::Sqrt(x[0] - mpi)*TMath::Exp(-1.*par[2]*(x[0]-mpi));
    //    AliInfo("Background function set to: wit exponential");
    } 
    break;
  case 6:
    // the following comment must be removed
    //     // pol 3, following convention for pol 2
    //     //y=a+b*x+c*x**2+d*x**3 -> integral = a(max-min) + 1/2*b*(max^2-min^2) +
    //     //+ 1/3*c*(max^3-min^3) + 1/4 d * (max^4-min^4) -> 
    //     //a = (integral-1/2*b*(max^2-min^2)-1/3*c*(max^3-min^3) - 1/4 d * (max^4-min^4) )/(max-min)
    //     // * [0] = integralBkg;
    //     // * [1] = b;
    //     // * [2] = c;
    //     // * [3] = d;
    {
      total=par[0];
      for(Int_t it=1;it<=fPolDegreeBkg;it++){
	total+=par[it]*TMath::Power(x[0]-fMassParticle,it)/TMath::Factorial(it);
      }
    }
    break;
  }
  return total;
}
//_________________________________________________________________________
Double_t AliHFInvMassFitter::FitFunction4Sgn (Double_t *x, Double_t *par){
  /// Fit function for the signal
  ///

  //  AliInfo("Signal function set to: Gaussian");
  Double_t sigval=0;
  switch (fTypeOfFit4Sgn){
  case 0:
    //gaussian = A/(sigma*sqrt(2*pi))*exp(-(x-mean)^2/2/sigma^2)
    //Par:
    // * [0] = integralSgn
    // * [1] = mean
    // * [2] = sigma
  //gaussian = [0]/TMath::Sqrt(2.*TMath::Pi())/[2]*exp[-(x-[1])*(x-[1])/(2*[2]*[2])]
    sigval=par[0]/TMath::Sqrt(2.*TMath::Pi())/par[2]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2]);
    break;
  case 1:
    //double gaussian = A/(sigma*sqrt(2*pi))*exp(-(x-mean)^2/2/sigma^2)
    //Par:
    // * [0] = integralSgn
    // * [1] = mean
    // * [2] = sigma1
    // * [3] = 2nd gaussian ratio
    // * [4] = deltaSigma
    //gaussian = [0]/TMath::Sqrt(2.*TMath::Pi())/[2]*exp[-(x-[1])*(x-[1])/(2*[2]*[2])]
    Double_t g1=(1.-par[3])/TMath::Sqrt(2.*TMath::Pi())/par[2]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2]);
    Double_t g2=par[3]/TMath::Sqrt(2.*TMath::Pi())/par[4]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[4]/par[4]);
    sigval=par[0]*(g1+g2);
    break;
  }
  return sigval;
}
//_________________________________________________________________________
Double_t AliHFInvMassFitter::FitFunction4SecPeak (Double_t *x, Double_t *par){
  /// Fit function for a second gaussian peak 
  /// To be used, e.g., for D+->KKpi in the Ds mass spectrum

  //gaussian = A/(sigma*sqrt(2*pi))*exp(-(x-mean)^2/2/sigma^2)
  //Par:
  // * [0] = integralSgn
  // * [1] = mean
  // * [2] = sigma
  Double_t secgaval=par[0]/TMath::Sqrt(2.*TMath::Pi())/par[2]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2]);
  return secgaval;
}
//_________________________________________________________________________
Double_t AliHFInvMassFitter::FitFunction4Mass(Double_t *x, Double_t *par){
  /// Total fit function (signal+background+possible second peak)
  ///
  
  Double_t bkg=FitFunction4Bkg(x,par);
  Double_t sig=FitFunction4Sgn(x,&par[fNBkgPars]);
  Double_t sec=0.;
  if(fSecondPeak) sec=FitFunction4SecPeak(x,&par[fNBkgPars+fNSigPars]);
  return bkg+sig+sec;
}

//_________________________________________________________________________
void AliHFInvMassFitter::Signal(Double_t nOfSigma,Double_t &signal,Double_t &errsignal) const {
  /// Return signal integral in mean +- n sigma
  ///

  Double_t min=fMass-nOfSigma*fSigmaSgn;
  Double_t max=fMass+nOfSigma*fSigmaSgn;
  Signal(min,max,signal,errsignal);
  return;
}

//_________________________________________________________________________
void AliHFInvMassFitter::Signal(Double_t min, Double_t max, Double_t &signal,Double_t &errsignal) const {
  /// Return signal integral in a range
  ///

  signal=fSigFunc->Integral(min, max)/(Double_t)fHistoInvMass->GetBinWidth(1);
  errsignal=(fRawYieldErr/fRawYield)*signal;/*assume relative error is the same as for total integral*/
  return;
}

//___________________________________________________________________________
void AliHFInvMassFitter::Background(Double_t nOfSigma,Double_t &background,Double_t &errbackground) const {
  /// Return background integral in mean +- n sigma
  ///

  Double_t min=fMass-nOfSigma*fSigmaSgn;
  Double_t max=fMass+nOfSigma*fSigmaSgn;
  Background(min,max,background,errbackground);
  return;
  
}
//___________________________________________________________________________
void AliHFInvMassFitter::Background(Double_t min, Double_t max, Double_t &background,Double_t &errbackground) const {
  /// Return background integral in a range
  ///

  TF1 *funcbkg=0x0;
  if(fBkgFuncRef) funcbkg=fBkgFuncRef;
  else if(fBkgFunc) funcbkg=fBkgFunc;
  if(!funcbkg){
    AliError("Bkg function not found!");
    return;
  }

  Double_t intB=funcbkg->GetParameter(0);
  Double_t intBerr=funcbkg->GetParError(0);

  //relative error evaluation: from histo

  Int_t leftBand=fHistoInvMass->FindBin(fMass-4*fSigmaSgn);
  Int_t rightBand=fHistoInvMass->FindBin(fMass+4*fSigmaSgn);
  intB=fHistoInvMass->Integral(1,leftBand)+fHistoInvMass->Integral(rightBand,fHistoInvMass->GetNbinsX());
  Double_t sum2=0;
  for(Int_t i=1;i<=leftBand;i++){
    sum2+=fHistoInvMass->GetBinError(i)*fHistoInvMass->GetBinError(i);
  }
  for(Int_t i=rightBand; i<=fHistoInvMass->GetNbinsX();i++){
    sum2+=fHistoInvMass->GetBinError(i)*fHistoInvMass->GetBinError(i);
  }

  intBerr=TMath::Sqrt(sum2);

  background=funcbkg->Integral(min,max)/(Double_t)fHistoInvMass->GetBinWidth(1);
  errbackground=intBerr/intB*background; 

  return;

}
//__________________________________________________________________________

void AliHFInvMassFitter::Significance(Double_t nOfSigma,Double_t &significance,Double_t &errsignificance) const  {
  /// Return significance in mean +- n sigma
  ///

  Double_t min=fMass-nOfSigma*fSigmaSgn;
  Double_t max=fMass+nOfSigma*fSigmaSgn;
  Significance(min, max, significance, errsignificance);

  return;
}

//__________________________________________________________________________

void AliHFInvMassFitter::Significance(Double_t min, Double_t max, Double_t &significance,Double_t &errsignificance) const {
  /// Return significance integral in a range
  ///

  Double_t background,errbackground;
  Background(min,max,background,errbackground);

  if (fRawYield+background <= 0.){
    significance=-1;
    errsignificance=0;
    return;
  }

  AliVertexingHFUtils::ComputeSignificance(fRawYield,fRawYieldErr,background,errbackground,significance,errsignificance);
  
  return;
}
//________________________________________________________________________
Bool_t AliHFInvMassFitter::PrepareHighPolFit(TF1 *fback){
  /// Perform intermediate fit steps up to fPolDegreeBkg-1
  /// in case of fit with a polynomial with degree > 2 (fTypeOfFit4Bkg=6)

  TH1F *hCp=(TH1F*)fHistoInvMass->Clone("htemp");
  Double_t estimatecent=0.5*(hCp->GetBinContent(hCp->FindBin(fMass-3.5*fSigmaSgn))+hCp->GetBinContent(hCp->FindBin(fMass+3.5*fSigmaSgn)));// just a first rough estimate
  Double_t estimateslope=(hCp->GetBinContent(hCp->FindBin(fMass+3.5*fSigmaSgn))-hCp->GetBinContent(hCp->FindBin(fMass-3.5*fSigmaSgn)))/(7*fSigmaSgn);// first rough estimate

  fCurPolDegreeBkg=2;
  TF1 *funcbkg,*funcPrev=0x0;
  while(fCurPolDegreeBkg<=fPolDegreeBkg){
    funcbkg = new TF1(Form("temp%d",fCurPolDegreeBkg),this,&AliHFInvMassFitter::BackFitFuncPolHelper,fMinMass,fMaxMass,fCurPolDegreeBkg+1,"AliHFInvMassFitter","BackFitFuncPolHelper");
    if(funcPrev){
      for(Int_t j=0;j<fCurPolDegreeBkg;j++){// now is +1 degree w.r.t. previous fit funct
	funcbkg->SetParameter(j,funcPrev->GetParameter(j));
      }
      delete funcPrev;
    }
    else{
      funcbkg->SetParameter(0,estimatecent);
      funcbkg->SetParameter(1,estimateslope);
    }
    printf(" ----- Pre-fit of background with pol degree %d ---\n",fCurPolDegreeBkg);
    hCp->Fit(funcbkg,"REMN","");
    funcPrev=(TF1*)funcbkg->Clone("ftemp");
    delete funcbkg;
    fCurPolDegreeBkg++;
  }

  for(Int_t j=0;j<=fPolDegreeBkg;j++){
    fback->SetParameter(j,funcPrev->GetParameter(j));
    fback->SetParError(j,funcPrev->GetParError(j));
  }
  printf(" ----- Final background fit with pol degree %d ---\n",fPolDegreeBkg);
  hCp->Fit(fback,"REMN","");// THIS IS JUST TO SET NOT ONLY THE PARAMETERS BUT ALSO chi2, etc...


  // The following lines might be useful for debugging
  //   TCanvas *cDebug=new TCanvas();
  //   cDebug->cd();
  //   hCp->Draw();
  //   TString strout=Form("Test%d.root",(Int_t)fhistoInvMass->GetBinContent(fhistoInvMass->FindBin(fMass)));
  //   cDebug->Print(strout.Data());
  // delete cDebug;

  delete funcPrev;
  delete hCp;
  return kTRUE;

}
// _______________________________________________________________________
Double_t AliHFInvMassFitter::BackFitFuncPolHelper(Double_t *x,Double_t *par){
  /// Helper function for polynomials with degree>2
  ///

  Double_t maxDeltaM = 4.*fSigmaSgn;
  if(fOnlySideBands && TMath::Abs(x[0]-fMass) < maxDeltaM) {
    TF1::RejectPoint();
    return 0;
  }
  Double_t back=par[0];
  for(Int_t it=1;it<=fCurPolDegreeBkg;it++){
    back+=par[it]*TMath::Power(x[0]-fMassParticle,it)/TMath::Factorial(it);
  }
  return back;
}
