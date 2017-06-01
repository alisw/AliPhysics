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
  fNParsSig(3),
  fNParsBkg(2),
  fOnlySideBands(kFALSE),
  fNSigma4SideBands(4.),
  fFitOption("L,E"),
  fRawYield(0.),
  fRawYieldErr(0.),
  fSigFunc(0x0),
  fBkgFuncSb(0x0),
  fBkgFunc(0x0),
  fBkgFuncRefit(0x0),
  fReflections(kFALSE),
  fNParsRfl(0),
  fRflOverSig(0),
  fFixRflOverSig(kFALSE),
  fHistoTemplRfl(0x0),
  fSmoothRfl(kFALSE),
  fRawYieldHelp(0),
  fRflFunc(0x0),
  fBkRFunc(0x0),
  fSecondPeak(kFALSE),
  fNParsSec(0),
  fSecMass(-999.),
  fSecWidth(9999.),
  fFixSecMass(kFALSE),
  fFixSecWidth(kFALSE),
  fSecFunc(0x0),
  fTotFunc(0x0)
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
  fNParsSig(3),
  fNParsBkg(2),
  fOnlySideBands(kFALSE),
  fNSigma4SideBands(4.),
  fFitOption("L,E"),
  fRawYield(0.),
  fRawYieldErr(0.),
  fSigFunc(0x0),
  fBkgFuncSb(0x0),
  fBkgFunc(0x0),
  fBkgFuncRefit(0x0),
  fReflections(kFALSE),
  fNParsRfl(0),
  fRflOverSig(0),
  fFixRflOverSig(kFALSE),
  fHistoTemplRfl(0x0),
  fSmoothRfl(kFALSE),
  fRawYieldHelp(0),
  fRflFunc(0x0),
  fBkRFunc(0x0),
  fSecondPeak(kFALSE),
  fNParsSec(0),
  fSecMass(-999.),
  fSecWidth(9999.),
  fFixSecMass(kFALSE),
  fFixSecWidth(kFALSE),
  fSecFunc(0x0),
  fTotFunc(0x0)
{
  /// standard constructor
  fHistoInvMass=(TH1F*)histoToFit->Clone("fHistoInvMass");
  fHistoInvMass->SetDirectory(0);
  SetNumberOfParams();
}
//_________________________________________________________________________
AliHFInvMassFitter::~AliHFInvMassFitter() {

  ///destructor

  delete fSigFunc;
  delete fBkgFunc;
  delete fBkgFuncSb;
  delete fBkgFuncRefit;
  delete fHistoInvMass;
  delete fRflFunc;
  delete fBkRFunc;
  delete fSecFunc;
  delete fTotFunc;
  delete fHistoTemplRfl;

}
//__________________________________________________________________________
void AliHFInvMassFitter::SetNumberOfParams(){
  /// Configure number of parameters of fit functions
  ///

  switch (fTypeOfFit4Bkg) {
  case 0:
    fNParsBkg=2;
    break;
  case 1:
    fNParsBkg=2;
    break;
  case 2:
    fNParsBkg=3;
    break;
  case 3:
    fNParsBkg=1;
    break;
  case 4:
    fNParsBkg=2;
    break;
  case 5:
    fNParsBkg=3;
    break;
  case 6:
    fNParsBkg=fPolDegreeBkg+1;
    break;
  default:
    AliError("Error in computing fNParsBkg: check fTypeOfFit4Bkg");
    break;
  }

  switch (fTypeOfFit4Sgn) {
  case 0:
    fNParsSig=3;
    break;
  case 1:
    fNParsSig=5;
    break;
  default:
    AliError("Error in computing fNParsSig: check fTypeOfFit4Sgn");
    break;
  }

  if(fReflections) fNParsRfl=1;
  else fNParsRfl=0;

  if(fSecondPeak) fNParsSec=3;
  else fNParsSec=0;

}
//__________________________________________________________________________
Int_t AliHFInvMassFitter::MassFitter(Bool_t draw){
  /// Main function to fit the invariant mass distribution
  /// returns 0 if the fit fails
  /// returns 1 if the fit succeeds
  /// returns 2 if there is no signal and the fit is performed with only background

  TVirtualFitter::SetDefaultFitter("Minuit");

  Double_t integralHisto=fHistoInvMass->Integral(fHistoInvMass->FindBin(fMinMass),fHistoInvMass->FindBin(fMaxMass),"width");

  fOnlySideBands = kTRUE;
  fBkgFuncSb = CreateBackgroundFitFunction("funcbkgsb",integralHisto);
  Int_t status=-1;
  printf("\n--- First fit with only background on the side bands - Exclusion region = %.2f sigma ---\n",fNSigma4SideBands);
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
    printf("   ---> Failed first fit with only background, minuit status = %d\n",status);
    return 0;
  }

  fOnlySideBands = kFALSE;
  if(!fBkgFunc){
    fBkgFunc = CreateBackgroundFitFunction("funcbkg",integralHisto);
    for(Int_t ipar=0; ipar<fNParsBkg; ipar++) fBkgFunc->SetParameter(ipar,fBkgFuncSb->GetParameter(ipar));
  }
  fBkgFunc->SetLineColor(kGray+1);

  printf("\n--- Estimate signal counts in the peak region ---\n");
  Double_t estimSignal=CheckForSignal(fMass,fSigmaSgn);
  Bool_t doFinalFit=kTRUE;
  if(estimSignal<0.){
    if(draw) DrawFit();
    estimSignal=0.;
    doFinalFit=kFALSE;
  } 

  fRawYieldHelp=estimSignal; // needed for reflection normalization
  if(!fBkgFuncRefit){
    fBkgFuncRefit = CreateBackgroundFitFunction("funcbkgrefit",integralHisto);
    for(Int_t ipar=0; ipar<fNParsBkg; ipar++) fBkgFuncRefit->SetParameter(ipar,fBkgFunc->GetParameter(ipar));
  }
  fBkgFuncRefit->SetLineColor(2);
  fSigFunc = CreateSignalFitFunction("fsigfit",estimSignal);
  if(fSecondPeak){
    printf("   ---> Final fit includes a second inv. mass peak\n");
    Double_t estimSec=CheckForSignal(fSecMass,fSecWidth);
    fSecFunc = CreateSecondPeakFunction("fsecpeak",estimSec);
  }
  if(fReflections){
    printf("   ---> Final fit includes reflections\n");
    fRflFunc = CreateReflectionFunction("freflect");
    fBkRFunc = CreateBackgroundPlusReflectionFunction("fbkgrfl");
  }
  fTotFunc = CreateTotalFitFunction("funcmass");

  if(doFinalFit){
    printf("\n--- Final fit with signal+background on the full range ---\n");
    status=fHistoInvMass->Fit("funcmass",Form("R,%s,+,0",fFitOption.Data()));
    if (status != 0){
      printf("   ---> Failed fit with signal+background, minuit status = %d\n",status);
      return 0;
    }
  }

  for(Int_t ipar=0; ipar<fNParsBkg; ipar++){
    fBkgFuncRefit->SetParameter(ipar,fTotFunc->GetParameter(ipar));
    fBkgFuncRefit->SetParError(ipar,fTotFunc->GetParError(ipar));
  }
  for(Int_t ipar=0; ipar<fNParsSig; ipar++){
    fSigFunc->SetParameter(ipar,fTotFunc->GetParameter(ipar+fNParsBkg));
    fSigFunc->SetParError(ipar,fTotFunc->GetParError(ipar+fNParsBkg));
  }
  if(fSecondPeak){
    for(Int_t ipar=0; ipar<fNParsSec; ipar++){
      fSecFunc->SetParameter(ipar,fTotFunc->GetParameter(ipar+fNParsBkg+fNParsSig));
      fSecFunc->SetParError(ipar,fTotFunc->GetParError(ipar+fNParsBkg+fNParsSig));
    }
    fSecFunc->SetLineColor(kMagenta+1);
    fSecFunc->SetLineStyle(3);
  }
  if(fReflections){
    for(Int_t ipar=0; ipar<fNParsRfl; ipar++){
      fRflFunc->SetParameter(ipar,fTotFunc->GetParameter(ipar+fNParsBkg+fNParsSig+fNParsSec));
      fRflFunc->SetParError(ipar,fTotFunc->GetParError(ipar+fNParsBkg+fNParsSig+fNParsSec));
      fBkRFunc->SetParameter(ipar+fNParsBkg,fTotFunc->GetParameter(ipar+fNParsBkg+fNParsSig+fNParsSec));
      fBkRFunc->SetParError(ipar+fNParsBkg,fTotFunc->GetParError(ipar+fNParsBkg+fNParsSig+fNParsSec));
    }
    for(Int_t ipar=0; ipar<fNParsBkg; ipar++){
      fBkRFunc->SetParameter(ipar,fTotFunc->GetParameter(ipar));
      fBkRFunc->SetParError(ipar,fTotFunc->GetParError(ipar));
    }
    fRflFunc->SetLineColor(kGreen+1);
    fBkRFunc->SetLineColor(kRed+1);
    fBkRFunc->SetLineStyle(7);
  }
  fMass=fSigFunc->GetParameter(1);
  fMassErr=fSigFunc->GetParError(1);
  fSigmaSgn=fSigFunc->GetParameter(2);
  fSigmaSgnErr=fSigFunc->GetParError(2);
  fTotFunc->SetLineColor(4);
  fRawYield=fTotFunc->GetParameter(fNParsBkg)/fHistoInvMass->GetBinWidth(1);
  fRawYieldErr=fTotFunc->GetParError(fNParsBkg)/fHistoInvMass->GetBinWidth(1);
  fRawYieldHelp=fRawYield;
  if(draw) DrawFit();
  if(doFinalFit) return 1;
  else return 2;
}
//______________________________________________________________________________
void AliHFInvMassFitter::DrawFit(){
  /// Steering method to draw the fit output
  ///

  TCanvas* c0=new TCanvas("c0");
  DrawHere(c0);
}
//______________________________________________________________________________
void AliHFInvMassFitter::DrawHere(TVirtualPad* c, Double_t nsigma,Int_t writeFitInfo){
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
  if(fBkgFuncRefit) fBkgFuncRefit->Draw("same");
  if(fBkRFunc) fBkRFunc->Draw("same");
  if(fRflFunc) fRflFunc->Draw("same");
  if(fSecFunc) fSecFunc->Draw("same");
  if(fTotFunc) fTotFunc->Draw("same");

  if(writeFitInfo > 0){
    TPaveText *pinfos=new TPaveText(0.12,0.65,0.47,0.89,"NDC");
    TPaveText *pinfom=new TPaveText(0.6,0.7,1.,.87,"NDC");
    pinfos->SetBorderSize(0);
    pinfos->SetFillStyle(0);
    pinfom->SetBorderSize(0);
    pinfom->SetFillStyle(0);
    if(fTotFunc){
      pinfom->SetTextColor(kBlue);
      for(Int_t ipar=1; ipar<fNParsSig; ipar++){
	pinfom->AddText(Form("%s = %.3f #pm %.3f",fTotFunc->GetParName(ipar+fNParsBkg),fTotFunc->GetParameter(ipar+fNParsBkg),fTotFunc->GetParError(ipar+fNParsBkg)));
      }
      pinfom->Draw();
      
      Double_t bkg,errbkg;
      Background(nsigma,bkg,errbkg);
      Double_t signif,errsignif;
      Significance(nsigma,signif,errsignif);
      
      pinfos->AddText(Form("S = %.0f #pm %.0f ",fRawYield,fRawYieldErr));
      pinfos->AddText(Form("B (%.0f#sigma) = %.0f #pm %.0f",nsigma,bkg,errbkg));
      pinfos->AddText(Form("S/B = (%.0f#sigma) %.4f ",nsigma,fRawYield/bkg));
      if(fRflFunc)  pinfos->AddText(Form("Refl/Sig =  %.3f #pm %.3f ",fRflFunc->GetParameter(0),fRflFunc->GetParError(0)));
      pinfos->AddText(Form("Signif (%.0f#sigma) = %.1f #pm %.1f ",nsigma,signif,errsignif));
      pinfos->Draw();
    }
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
  printf("   ---> IntegralUnderFitFunc=%f  IntegralUnderHisto=%f   EstimatedSignal=%f\n",sum,sumback,diffUnderPeak);
  if(diffUnderPeak/TMath::Sqrt(sum)<1.){
    printf("   ---> (Tot-Bkg)/sqrt(Tot)=%f ---> Likely no signal/\n",diffUnderPeak/TMath::Sqrt(sum));
    return -1;
  }
  return diffUnderPeak*fHistoInvMass->GetBinWidth(1);
}

//______________________________________________________________________________
TF1* AliHFInvMassFitter::CreateBackgroundFitFunction(TString fname, Double_t integral){
  /// Creates the background fit fucntion
  ///

  SetNumberOfParams();
  TF1* funcbkg =  new TF1(fname.Data(),this,&AliHFInvMassFitter::FitFunction4Bkg,fMinMass,fMaxMass,fNParsBkg,"AliHFInvMassFitter","FitFunction4Bkg");
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
    for(Int_t j=0;j<fNParsBkg; j++){
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
TF1* AliHFInvMassFitter::CreateReflectionFunction(TString fname){
  /// Creates a function for reflections contribution
  /// in the D0->Kpi inv. mass distribution
  TF1* funcrfl =  new TF1(fname.Data(),this,&AliHFInvMassFitter::FitFunction4Refl,fMinMass,fMaxMass,1,"AliHFInvMassFitter","FitFunction4Refl");
  funcrfl->SetParameter(0,fRflOverSig);
  funcrfl->SetParLimits(0,0.,1.);
  if(fFixRflOverSig) funcrfl->FixParameter(0,fRflOverSig);
  funcrfl->SetParNames("ReflOverS");
  return funcrfl;
}
//______________________________________________________________________________
TF1* AliHFInvMassFitter::CreateBackgroundPlusReflectionFunction(TString fname){
  /// Creates the function with sum of background and reflections
  ///
  SetNumberOfParams();
  Int_t totParams=fNParsBkg+fNParsRfl;
  TF1* fbr=new TF1(fname.Data(),this,&AliHFInvMassFitter::FitFunction4BkgAndRefl,fMinMass,fMaxMass,totParams,"AliHFInvMassFitter","FitFunction4BkgAndRefl");
  for(Int_t ipar=0; ipar<fNParsBkg; ipar++){
    fbr->SetParameter(ipar,fBkgFunc->GetParameter(ipar));
    fbr->SetParName(ipar,fBkgFunc->GetParName(ipar));
  }
  for(Int_t ipar=0; ipar<fNParsRfl; ipar++){
    fbr->SetParameter(ipar+fNParsBkg,fRflFunc->GetParameter(ipar));
    fbr->SetParName(ipar+fNParsBkg,fRflFunc->GetParName(ipar));
    // par limits not set because this function is not used for fitting but only for drawing
  }
  return fbr;
}
//______________________________________________________________________________
TF1* AliHFInvMassFitter::CreateSignalFitFunction(TString fname, Double_t integsig){
  /// Creates the fit function for the signal peak
  ///

  SetNumberOfParams();
  TF1* funcsig =  new TF1(fname.Data(),this,&AliHFInvMassFitter::FitFunction4Sgn,fMinMass,fMaxMass,fNParsSig,"AliHFInvMassFitter","FitFunction4Sgn");
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
  Int_t totParams=fNParsBkg+fNParsRfl+fNParsSec+fNParsSig;
  TF1* ftot=ftot=new TF1(fname.Data(),this,&AliHFInvMassFitter::FitFunction4Mass,fMinMass,fMaxMass,totParams,"AliHFInvMassFitter","FitFunction4Mass");
  for(Int_t ipar=0; ipar<fNParsBkg; ipar++){
    ftot->SetParameter(ipar,fBkgFunc->GetParameter(ipar));
    ftot->SetParName(ipar,fBkgFunc->GetParName(ipar));
  }
  for(Int_t ipar=0; ipar<fNParsSig; ipar++){
    ftot->SetParameter(ipar+fNParsBkg,fSigFunc->GetParameter(ipar));
    ftot->SetParName(ipar+fNParsBkg,fSigFunc->GetParName(ipar));
    Double_t parmin,parmax;
    fSigFunc->GetParLimits(ipar,parmin,parmax);
    ftot->SetParLimits(ipar+fNParsBkg,parmin,parmax);
  }
  if(fSecondPeak && fSecFunc){
    for(Int_t ipar=0; ipar<fNParsSec; ipar++){
      ftot->SetParameter(ipar+fNParsBkg+fNParsSig,fSecFunc->GetParameter(ipar));
      ftot->SetParName(ipar+fNParsBkg+fNParsSig,fSecFunc->GetParName(ipar));
      Double_t parmin,parmax;
      fSecFunc->GetParLimits(ipar,parmin,parmax);
      ftot->SetParLimits(ipar+fNParsBkg+fNParsSig,parmin,parmax);
    }
  }
  if(fReflections && fRflFunc){
    for(Int_t ipar=0; ipar<fNParsRfl; ipar++){
      ftot->SetParameter(ipar+fNParsBkg+fNParsSig+fNParsSec,fRflFunc->GetParameter(ipar));
      ftot->SetParName(ipar+fNParsBkg+fNParsSig+fNParsSec,fRflFunc->GetParName(ipar));
      Double_t parmin,parmax;
      fRflFunc->GetParLimits(ipar,parmin,parmax);
      ftot->SetParLimits(ipar+fNParsBkg+fNParsSig+fNParsSec,parmin,parmax);
    }
  }
  return ftot;
}
//__________________________________________________________________________
Double_t AliHFInvMassFitter::FitFunction4Bkg (Double_t *x, Double_t *par){
  /// Fit function for the background
  ///

  Double_t maxDeltaM = fNSigma4SideBands*fSigmaSgn;
  if(fOnlySideBands && TMath::Abs(x[0]-fMass) < maxDeltaM) {
    TF1::RejectPoint();
    return 0;
  }
  if(fOnlySideBands && fSecondPeak && TMath::Abs(x[0]-fSecMass) < (fNSigma4SideBands*fSecWidth)){
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
  fRawYieldHelp=par[0]/fHistoInvMass->GetBinWidth(1);
  return sigval;
}
//_________________________________________________________________________
Double_t AliHFInvMassFitter::FitFunction4Refl(Double_t *x,Double_t *par){
  /// Fit function for reflections:
  /// D0->Kpi decays with swapped mass assignment to pion and kaon decay tracks
  if(!fHistoTemplRfl) return 0;

  Int_t bin =fHistoTemplRfl->FindBin(x[0]);
  Double_t value=fHistoTemplRfl->GetBinContent(bin);
  Int_t binmin=fHistoTemplRfl->FindBin(fMinMass*1.00001);
  Int_t binmax=fHistoTemplRfl->FindBin(fMaxMass*0.99999);
  Double_t norm=fHistoTemplRfl->Integral(binmin,binmax)*fHistoTemplRfl->GetBinWidth(bin);
  if(TMath::Abs(value)<1.e-14&&fSmoothRfl){// very rough, assume a constant trend, much better would be a pol1 or pol2 over a broader range
    value+=fHistoTemplRfl->GetBinContent(bin-1)+fHistoTemplRfl->GetBinContent(bin+1);
    value/=3.;
  }
  return par[0]*value/norm*fRawYieldHelp*fHistoInvMass->GetBinWidth(1);
}
//_________________________________________________________________________
Double_t AliHFInvMassFitter::FitFunction4BkgAndRefl(Double_t *x, Double_t *par){
  /// Fit fucntion with the sum of background and reflections
  ///
  Double_t bkg=FitFunction4Bkg(x,par);
  Double_t refl=0;
  if(fReflections) refl=FitFunction4Refl(x,&par[fNParsBkg]);
  return bkg+refl;
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
  Double_t sig=FitFunction4Sgn(x,&par[fNParsBkg]);
  Double_t sec=0.;
  if(fSecondPeak) sec=FitFunction4SecPeak(x,&par[fNParsBkg+fNParsSig]);
  Double_t refl=0;
  if(fReflections) refl=FitFunction4Refl(x,&par[fNParsBkg+fNParsSig+fNParsSec]); 
  return bkg+sig+sec+refl;
}

//_________________________________________________________________________
void AliHFInvMassFitter::Signal(Double_t nOfSigma,Double_t &signal,Double_t &errsignal) const {
  /// Return signal integral in mean +- n sigma
  ///

  Double_t minMass=fMass-nOfSigma*fSigmaSgn;
  Double_t maxMass=fMass+nOfSigma*fSigmaSgn;
  Signal(minMass,maxMass,signal,errsignal);
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

  Double_t minMass=fMass-nOfSigma*fSigmaSgn;
  Double_t maxMass=fMass+nOfSigma*fSigmaSgn;
  Background(minMass,maxMass,background,errbackground);
  return;
  
}
//___________________________________________________________________________
void AliHFInvMassFitter::Background(Double_t min, Double_t max, Double_t &background,Double_t &errbackground) const {
  /// Return background integral in a range
  ///

  TF1 *funcbkg=0x0;
  if(fBkgFuncRefit) funcbkg=fBkgFuncRefit;
  else if(fBkgFunc) funcbkg=fBkgFunc;
  if(!funcbkg){
    AliError("Bkg function not found!");
    return;
  }

  Double_t intB=funcbkg->GetParameter(0);
  Double_t intBerr=funcbkg->GetParError(0);

  //relative error evaluation: from histo

  Int_t leftBand=fHistoInvMass->FindBin(fMass-fNSigma4SideBands*fSigmaSgn);
  Int_t rightBand=fHistoInvMass->FindBin(fMass+fNSigma4SideBands*fSigmaSgn);
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

  Double_t minMass=fMass-nOfSigma*fSigmaSgn;
  Double_t maxMass=fMass+nOfSigma*fSigmaSgn;
  Significance(minMass, maxMass, significance, errsignificance);

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
    printf("   ---> Pre-fit of background with pol degree %d ---\n",fCurPolDegreeBkg);
    hCp->Fit(funcbkg,"REMN","");
    funcPrev=(TF1*)funcbkg->Clone("ftemp");
    delete funcbkg;
    fCurPolDegreeBkg++;
  }

  for(Int_t j=0;j<=fPolDegreeBkg;j++){
    fback->SetParameter(j,funcPrev->GetParameter(j));
    fback->SetParError(j,funcPrev->GetParError(j));
  }
  printf("   ---> Final background fit with pol degree %d ---\n",fPolDegreeBkg);
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

  Double_t maxDeltaM = fNSigma4SideBands*fSigmaSgn;
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
// _______________________________________________________________________
TH1F* AliHFInvMassFitter::SetTemplateReflections(const TH1 *h, TString opt,Double_t minRange,Double_t maxRange){
  /// Method to create the reflection invariant mass distributions from MC templates
  /// option could be: 
  ///    "template"                use MC histograms
  ///    "1gaus" ot "singlegaus"   single gaussian function fit to MC templates
  ///    "2gaus" ot "doublegaus"   double gaussian function fit to MC templates
  ///    "pol3"                    3rd order polynomial fit to MC templates
  ///    "pol6"                    6th order polynomial fit to MC templates

  if(!h){
    fReflections=kFALSE;
    return 0x0;
  }
  fHistoTemplRfl=(TH1F*)h->Clone("hTemplRfl");
  opt.ToLower();
  fReflections=kTRUE;
  printf("\n--- Reflection templates from simulation ---\n");
  if(opt.Contains("templ")){
    printf("   ---> Reflection contribution using directly the histogram from simulation\n");
    return fHistoTemplRfl;
  }

  TF1 *f=0x0;
  Bool_t isPoissErr=kTRUE;
  Double_t xMinForFit=h->GetBinLowEdge(1);
  Double_t xMaxForFit=h->GetXaxis()->GetBinUpEdge(h->GetNbinsX());
  if(minRange>=0 && maxRange>=0){
    xMinForFit=TMath::Max(minRange,h->GetBinLowEdge(1));
    xMaxForFit=TMath::Min(maxRange,h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
  }
  if(opt.EqualTo("1gaus") || opt.EqualTo("singlegaus")){
    printf("   ---> Reflection contribution from single-Gaussian fit to histogram from simulation\n");
    f=new TF1("mygaus","gaus",xMinForFit,xMaxForFit);
    f->SetParameter(0,h->GetMaximum());
    //    f->SetParLimits(0,0,100.*h->Integral());
    f->SetParameter(1,1.865);
    f->SetParameter(2,0.050);
    fHistoTemplRfl->Fit(f,"REM","");//,h->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
  }
  else if(opt.EqualTo("2gaus") || opt.EqualTo("doublegaus")){
    printf("   ---> Reflection contribution from double-Gaussian fit to histogram from simulation\n");
    f=new TF1("my2gaus","[0]*([3]/( TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+(1.-[3])/( TMath::Sqrt(2.*TMath::Pi())*[5])*TMath::Exp(-(x-[4])*(x-[4])/(2.*[5]*[5])))",xMinForFit,xMaxForFit);
    f->SetParameter(0,h->GetMaximum());
    //    f->SetParLimits(0,0,100.*h->Integral());
    f->SetParLimits(3,0.,1.);
    f->SetParameter(3,0.5);
    f->SetParameter(1,1.84);
    f->SetParameter(2,0.050);
    f->SetParameter(4,1.88);
    f->SetParameter(5,0.050);
    fHistoTemplRfl->Fit(f,"REM","");//,h->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
  }
  else if(opt.EqualTo("pol3")){
    printf("   ---> Reflection contribution from pol3 fit to histogram from simulation\n");
    f=new TF1("mypol3","pol3",xMinForFit,xMaxForFit);
    f->SetParameter(0,h->GetMaximum());
    //    f->SetParLimits(0,0,100.*h->Integral());
    // Hard to initialize the other parameters...
    fHistoTemplRfl->Fit(f,"REM","");
    //    Printf("We USED %d POINTS in the Fit",f->GetNumberFitPoints());
  }
  else if(opt.EqualTo("pol6")){
    printf("   ---> Reflection contribution from pol6 fit to histogram from simulation\n");
    f=new TF1("mypol6","pol6",xMinForFit,xMaxForFit);
    f->SetParameter(0,h->GetMaximum());
    //    f->SetParLimits(0,0,100.*h->Integral());
    // Hard to initialize the other parameters...
    fHistoTemplRfl->Fit(f,"RLEMI","");//,h->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
  }
  else{
    // no good option passed
    printf("   ---> Bad option for reflection configuration -> reflections will not be included in the fit\n");
    fReflections=kFALSE;
    delete fHistoTemplRfl;
    fHistoTemplRfl=0x0;
    return 0x0;
  }

  // Fill fHistoTemplRfl with values of fit function
  if(f){
    for(Int_t j=1;j<=fHistoTemplRfl->GetNbinsX();j++){
      fHistoTemplRfl->SetBinContent(j,f->Integral(fHistoTemplRfl->GetBinLowEdge(j),fHistoTemplRfl->GetXaxis()->GetBinUpEdge(j))/fHistoTemplRfl->GetBinWidth(j));
      if(fHistoTemplRfl->GetBinContent(j)>=0.&&TMath::Abs(h->GetBinError(j)*h->GetBinError(j)-h->GetBinContent(j))>0.1*h->GetBinContent(j))isPoissErr=kFALSE;
    }
    for(Int_t j=1;j<=fHistoTemplRfl->GetNbinsX();j++){
      if(isPoissErr){
	fHistoTemplRfl->SetBinError(j,TMath::Sqrt(fHistoTemplRfl->GetBinContent(j)));
      }
      else fHistoTemplRfl->SetBinError(j,0.001*fHistoTemplRfl->GetBinContent(j));
    }
    fReflections=kTRUE;
    return fHistoTemplRfl;
  }else{
    printf("   ---> Fit to MC template for reflection failed -> reflections will not be included in the fit\n");
    fReflections=kFALSE;
    delete fHistoTemplRfl;
    fHistoTemplRfl=0x0;
    return 0x0;
  }
  return 0x0;
}
// _______________________________________________________________________
Double_t AliHFInvMassFitter::GetRawYieldBinCounting(Double_t& errRyBC, Double_t nOfSigma, Int_t option, Int_t pdgCode) const{
  /// Method to compute the signal using inv. mass histo bin counting 
  /// -> interface method to compute yield in nsigma range around peak
  /// pdgCode: if==411,421,413,413 or 4122: range defined based on PDG mass
  //           else (default) mean of gaussian fit

  Double_t massVal=fMass;
  switch (pdgCode) {
  case 411:
  case 421:
  case 431:
  case 4122:
    massVal=TDatabasePDG::Instance()->GetParticle(pdgCode)->Mass();
    break;
  case 413:
    massVal=TDatabasePDG::Instance()->GetParticle(pdgCode)->Mass();
    massVal-=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    break;
  default:
    massVal=fMass;
    break;
  }

  Double_t minMass=massVal-nOfSigma*fSigmaSgn;
  Double_t maxMass=massVal+nOfSigma*fSigmaSgn;
  return GetRawYieldBinCounting(errRyBC,minMass,maxMass,option);
}
// _______________________________________________________________________
Double_t AliHFInvMassFitter::GetRawYieldBinCounting(Double_t& errRyBC, Double_t minMass, Double_t maxMass, Int_t option) const{
  /// Method to compute the signal using inv. mass histo bin counting 
  /// after background subtraction from background fit function
  ///   option=0: background fit function from 1st fit step (only side bands)
  ///   option=1: background fit function from 2nd fit step (S+B)

  Int_t minBinSum=fHistoInvMass->FindBin(minMass);
  Int_t maxBinSum=fHistoInvMass->FindBin(maxMass);
  if(minBinSum<1){
    printf("Left range for bin counting smaller than allowed by histogram axis, setting it to the lower edge of the first histo bin\n");
    minBinSum=1;
  }
  if(maxBinSum>fHistoInvMass->GetNbinsX()){
    printf("Right range for bin counting larger than allowed by histogram axis, setting it to the upper edge of the last histo bin\n");
    maxBinSum=fHistoInvMass->GetNbinsX();
  }
  Double_t cntSig=0.;
  Double_t cntErr=0.;
  errRyBC=0;
  TF1* fbackground=fBkgFunc;
  if(option==1) fbackground=fBkgFuncRefit;
  if(!fbackground) return 0.;

  for(Int_t jb=minBinSum; jb<=maxBinSum; jb++){
    Double_t cntTot=fHistoInvMass->GetBinContent(jb);
    Double_t cntBkg=fbackground->Integral(fHistoInvMass->GetBinLowEdge(jb),fHistoInvMass->GetBinLowEdge(jb)+fHistoInvMass->GetBinWidth(jb))/fHistoInvMass->GetBinWidth(jb);
    //Double_t cntBkg=fbackground->Eval(fHistoInvMass->GetBinCenter(jb));
    cntSig+=(cntTot-cntBkg);
    cntErr+=(fHistoInvMass->GetBinError(jb)*fHistoInvMass->GetBinError(jb));
  }
  errRyBC=TMath::Sqrt(cntErr);
  return cntSig;
}


// _______________________________________________________________________
TH1F* AliHFInvMassFitter::GetResidualsAndPulls(TH1 *hPulls,TH1 *hResidualTrend,TH1 *hPullsTrend, Double_t minrange,Double_t maxrange){

  /// fill and return the residual and pull histos

  Int_t binmi=fHistoInvMass->FindBin(fMinMass*1.001);
  Int_t binma=fHistoInvMass->FindBin(fMaxMass*0.9999);
  if(maxrange>minrange){
    binmi=fHistoInvMass->FindBin(minrange*1.001);
    binma=fHistoInvMass->FindBin(maxrange*0.9999);
  }
  if(hResidualTrend){
    //fHistoInvMass->Copy(hResidualTrend);
    hResidualTrend->SetBins(fHistoInvMass->GetNbinsX(),fHistoInvMass->GetXaxis()->GetXmin(),fHistoInvMass->GetXaxis()->GetXmax());
    hResidualTrend->SetName(Form("%s_residualTrend",fHistoInvMass->GetName()));
    hResidualTrend->SetTitle(Form("%s  (Residuals)",fHistoInvMass->GetTitle()));
    hResidualTrend->SetMarkerStyle(20);
    hResidualTrend->SetMarkerSize(1.0);
    hResidualTrend->Reset();
  }
  if(hPullsTrend){
    hPullsTrend->SetBins(fHistoInvMass->GetNbinsX(),fHistoInvMass->GetXaxis()->GetXmin(),fHistoInvMass->GetXaxis()->GetXmax());
    hPullsTrend->Reset();
    hPullsTrend->SetName(Form("%s_pullTrend",fHistoInvMass->GetName()));
    hPullsTrend->SetTitle(Form("%s (Pulls)",fHistoInvMass->GetTitle()));
    hPullsTrend->SetMarkerStyle(20);
    hPullsTrend->SetMarkerSize(1.0);
  }
  if(hPulls){
    hPulls->SetName(Form("%s_pulls",fHistoInvMass->GetName()));
    hPulls->SetTitle(Form("%s ; Pulls",fHistoInvMass->GetTitle()));
    hPulls->SetBins(40,-10,10);
    hPulls->Reset();
  }

  Double_t res=-1.e-6,min=1.e+12,max=-1.e+12;
  TArrayD *arval=new TArrayD(binma-binmi+1);
  for(Int_t jst=1;jst<=fHistoInvMass->GetNbinsX();jst++){      
    
    res=fHistoInvMass->GetBinContent(jst)-fTotFunc->Integral(fHistoInvMass->GetBinLowEdge(jst),fHistoInvMass->GetBinLowEdge(jst)+fHistoInvMass->GetBinWidth(jst))/fHistoInvMass->GetBinWidth(jst);
    if(jst>=binmi&&jst<=binma){
      arval->AddAt(res,jst-binmi);
      if(res<min)min=res;
      if(res>max)max=res;
    }
    //      Printf("Res = %f from %f - %f",res,fHistoInvMass->GetBinContent(jst),fTotFunc->Integral(fHistoInvMass->GetBinLowEdge(jst),fHistoInvMass->GetBinLowEdge(jst)+fHistoInvMass->GetBinWidth(jst))/fHistoInvMass->GetBinWidth(jst));
    if(hResidualTrend){
      hResidualTrend->SetBinContent(jst,res);
      hResidualTrend->SetBinError(jst,fHistoInvMass->GetBinError(jst));
    }
    if(hPulls){
      if(jst>=binmi&&jst<=binma)hPulls->Fill(res/fHistoInvMass->GetBinError(jst));
    }    
    if(hPullsTrend){
      hPullsTrend->SetBinContent(jst,res/fHistoInvMass->GetBinError(jst));
      hPullsTrend->SetBinError(jst,0.0001);
    }
  }
  if(hResidualTrend)hResidualTrend->GetXaxis()->SetRange(binmi,binma);
  if(hPullsTrend){
    hPullsTrend->GetXaxis()->SetRange(binmi,binma);
    hPullsTrend->SetMinimum(-7);
    hPullsTrend->SetMaximum(+7);
  }
  if(TMath::Abs(min)>TMath::Abs(max))max=min;

  TH1F *hout=new TH1F(Form("%s_residuals",fHistoInvMass->GetName()),Form("%s ; residuals",fHistoInvMass->GetTitle()),25,-TMath::Abs(max)*1.5,TMath::Abs(max)*1.5);
  for(Int_t j=0;j<binma-binmi+1;j++){
    hout->Fill(arval->At(j));
  }
  hout->Sumw2();
  hout->Fit("gaus","LEM","",-TMath::Abs(max)*1.2,TMath::Abs(max)*1.2);

  if(hPulls){
    hPulls->Sumw2();
    hPulls->Fit("gaus","LEM","",-3,3);
  }
  delete arval;
  return hout;
}
// _______________________________________________________________________
void AliHFInvMassFitter::PrintFunctions(){
  /// dump the function parameters
  ///
  if(fBkgFunc){
    printf("--- Background function in 1st step fit to side bands ---\n");
    fBkgFunc->Print();
  }
  if(fTotFunc){
    printf("--- Total fit function from 2nd step fit ---\n");
    fTotFunc->Print();
  }
  if(fBkgFuncRefit){
    printf("--- Background function in 2nd step fit ---\n");
    fBkgFuncRefit->Print();
  }
  if(fRflFunc){
    printf("--- Reflections ---\n");
    fRflFunc->Print();
  }
  if(fBkRFunc){
    printf("--- Background + reflections ---\n");
    fBkRFunc->Print();
  }
  if(fSecFunc){
    printf("--- Additional Gaussian peak ---\n");
    fSecFunc->Print();
  }
}
