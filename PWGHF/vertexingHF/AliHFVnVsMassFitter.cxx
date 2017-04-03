#include "AliHFVnVsMassFitter.h"

#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TColor.h>
#include <TLegend.h>
#include <TList.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <TDatabasePDG.h>
#include <TPaveText.h>
#include "Fit/BinData.h"
#include "HFitInterface.h"
#include <vector>

#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliHFVnVsMassFitter);
/// \endcond

//________________________________________________________________
AliHFVnVsMassFitter::AliHFVnVsMassFitter()
  :TObject()
  ,fMassHisto(0x0)
  ,fVnVsMassHisto(0x0)
  ,fMassSgnFuncType(kGaus)
  ,fMassBkgFuncType(kExpo)
  ,fVnBkgFuncType(kLin)
  ,fMassFuncFromPrefit(0x0)
  ,fMassBkgFunc(0x0)
  ,fMassSgnFunc(0x0)
  ,fMassTotFunc(0x0)
  ,fVnBkgFuncSb(0x0)
  ,fVnBkgFunc(0x0)
  ,fVnTotFunc(0x0)
  ,fMassFitter(0x0)
  ,fMassMin(1.69)
  ,fMassMax(2.05)
  ,fVn(0.)
  ,fVnUncertainty(0.)
  ,fSigma(0.)
  ,fSigmaUncertainty(0.)
  ,fMean(0.)
  ,fMeanUncertainty(0.)
  ,fRawYield(0.)
  ,fRawYieldUncertainty(0.)
  ,fChiSquare(0.)
  ,fNDF(0)
  ,fProb(0.)
  ,fNSigmaForSB(3.)
  ,fSigmaInit(0.012)
  ,fMeanInit(1.870)
  ,fMeanFixedFromMassFit(kFALSE)
  ,fSigmaFixedFromMassFit(kFALSE)
  ,fMassParticle(1.870)
  ,fNParsMassSgn(3)
  ,fNParsMassBkg(2)
  ,fNParsVnBkg(2)
  ,fSigmaFixed(0)
  ,fMeanFixed(0)
  ,fPolDegreeBkg(3)
  ,fReflections(kFALSE)
  ,fNParsRfl(0)
  ,fRflOverSig(0.)
  ,fFixRflOverSig(kFALSE)
  ,fHistoTemplRfl(0x0)
  ,fHistoTemplRflInit(0x0)
  ,fMassRflFunc(0x0)
  ,fMassBkgRflFunc(0x0)
  ,fRflOpt("1gaus")
  ,fMinRefl(0.)
  ,fMaxRefl(0.)
  ,fSmoothRfl(kFALSE)
  ,fRawYieldHelp(0.)
  ,fSecondPeak(kFALSE)
  ,fMassSecPeakFunc(0x0)
  ,fNParsSec(0)
  ,fSecMass(-999.)
  ,fSecWidth(9999.)
  ,fFixSecMass(kFALSE)
  ,fFixSecWidth(kFALSE)
  ,fDoSecondPeakVn(kFALSE)
  ,fHarmonic(2) {

    //default constructor
}

//________________________________________________________________
AliHFVnVsMassFitter::AliHFVnVsMassFitter(TH1F* hMass, TH1F* hvn, Double_t min, Double_t max, Int_t funcMassBkg, Int_t funcMassSgn, Int_t funcvnBkg)
  :TObject()
  ,fMassSgnFuncType(funcMassSgn)
  ,fMassBkgFuncType(funcMassBkg)
  ,fVnBkgFuncType(funcvnBkg)
  ,fMassFuncFromPrefit(0x0)
  ,fMassBkgFunc(0x0)
  ,fMassSgnFunc(0x0)
  ,fMassTotFunc(0x0)
  ,fVnBkgFuncSb(0x0)
  ,fVnBkgFunc(0x0)
  ,fVnTotFunc(0x0)
  ,fMassFitter(0x0)
  ,fMassMin(min)
  ,fMassMax(max)
  ,fVn(0.)
  ,fVnUncertainty(0.)
  ,fSigma(0.)
  ,fSigmaUncertainty(0.)
  ,fMean(0.)
  ,fMeanUncertainty(0.)
  ,fRawYield(0.)
  ,fRawYieldUncertainty(0.)
  ,fChiSquare(0.)
  ,fNDF(0)
  ,fProb(0)
  ,fNSigmaForSB(3.)
  ,fSigmaInit(0.012)
  ,fMeanInit(1.870)
  ,fMeanFixedFromMassFit(kFALSE)
  ,fSigmaFixedFromMassFit(kFALSE)
  ,fMassParticle(1.870)
  ,fNParsMassSgn(3)
  ,fNParsMassBkg(2)
  ,fNParsVnBkg(2)
  ,fSigmaFixed(0)
  ,fMeanFixed(0)
  ,fPolDegreeBkg(3)
  ,fReflections(kFALSE)
  ,fNParsRfl(0)
  ,fRflOverSig(0.)
  ,fFixRflOverSig(kFALSE)
  ,fHistoTemplRfl(0x0)
  ,fHistoTemplRflInit(0x0)
  ,fMassRflFunc(0x0)
  ,fMassBkgRflFunc(0x0)
  ,fRflOpt("1gaus")
  ,fMinRefl(0.)
  ,fMaxRefl(0.)
  ,fSmoothRfl(kFALSE)
  ,fRawYieldHelp(0.)
  ,fSecondPeak(kFALSE)
  ,fMassSecPeakFunc(0x0)
  ,fNParsSec(0)
  ,fSecMass(-999.)
  ,fSecWidth(9999.)
  ,fFixSecMass(kFALSE)
  ,fFixSecWidth(kFALSE)
  ,fDoSecondPeakVn(kFALSE)
  ,fHarmonic(2) {

    //standard constructor
    fMassHisto = (TH1F*)hMass->Clone("fHistoInvMass");
    fMassHisto->SetDirectory(0);
    fVnVsMassHisto = (TH1F*)hvn->Clone(Form("fHistoV%dVsMass",fHarmonic));
    fVnVsMassHisto->SetDirectory(0);

    DefineNumberOfParameters();
}

//________________________________________________________________
AliHFVnVsMassFitter::~AliHFVnVsMassFitter() {

  //destructor
  if(fMassHisto)          delete fMassHisto;
  if(fVnVsMassHisto)      delete fVnVsMassHisto;
  if(fMassFuncFromPrefit) delete fMassFuncFromPrefit;
  if(fMassBkgFunc)        delete fMassBkgFunc;
  if(fMassBkgRflFunc)     delete fMassBkgRflFunc;
  if(fMassSgnFunc)        delete fMassSgnFunc;
  if(fMassTotFunc)        delete fMassTotFunc;
  if(fVnBkgFuncSb)        delete fVnBkgFuncSb;
  if(fVnBkgFunc)          delete fVnBkgFunc;
  if(fVnTotFunc)          delete fVnTotFunc;
  if(fMassFitter)         delete fMassFitter;
  if(fHistoTemplRfl)      delete fHistoTemplRfl;
  if(fHistoTemplRflInit)  delete fHistoTemplRflInit;
  if(fMassRflFunc)        delete fMassRflFunc;
  if(fMassSecPeakFunc)    delete fMassSecPeakFunc;
}

//________________________________________________________________
Bool_t AliHFVnVsMassFitter::SimultaneusFit(Bool_t drawFit) {

  if(!fMassHisto || !fVnVsMassHisto) {AliError("Histograms not set! Exit."); return kFALSE;}
  DefineNumberOfParameters();

  const Int_t nparsmass = fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl;
  Int_t NvnParsSgn = 1;
  if(fSecondPeak && fDoSecondPeakVn) {NvnParsSgn+=1;}
  const Int_t nparsvn = nparsmass+fNParsVnBkg+NvnParsSgn;

  Bool_t massprefit=MassPrefit();
  if(!massprefit) {AliError("Impossible to perform the mass prefit"); return kFALSE;}
  Bool_t vnprefit=VnSBPrefit();

  std::vector<Double_t> initpars;
  for(Int_t iBkgPar=0; iBkgPar<fNParsMassBkg; iBkgPar++) {
    initpars.push_back(fMassFuncFromPrefit->GetParameter(iBkgPar));
  }
  for(Int_t iSgnPar=0; iSgnPar<fNParsMassSgn; iSgnPar++) {
    initpars.push_back(fMassFuncFromPrefit->GetParameter(iSgnPar+fNParsMassBkg));
  }
  for(Int_t iSecPeakPar=0; iSecPeakPar<fNParsSec; iSecPeakPar++) {
    initpars.push_back(fMassFuncFromPrefit->GetParameter(iSecPeakPar+fNParsMassBkg+fNParsMassSgn));
  }
  for(Int_t iReflPar=0; iReflPar<fNParsRfl; iReflPar++) {
    initpars.push_back(fMassFuncFromPrefit->GetParameter(iReflPar+fNParsMassBkg+fNParsMassSgn+fNParsSec));
  }
  for(Int_t iVnBkgPar=0; iVnBkgPar<fNParsVnBkg; iVnBkgPar++) {
    if(vnprefit) {initpars.push_back(fVnBkgFuncSb->GetParameter(iVnBkgPar));}
    else {initpars.push_back(0.05);}
  }
  initpars.push_back(0.10); //initial parameter for signal vn
  if(fSecondPeak && fDoSecondPeakVn) {initpars.push_back(0.10);} //initial parameter for second peak vn

  fMassTotFunc = new TF1("fMassTotFunc",this,&AliHFVnVsMassFitter::MassFunc,fMassMin,fMassMax,nparsmass,"AliHFVnVsMassFitter","MassFunc");
  fVnTotFunc = new TF1("fVnTotFunc",this,&AliHFVnVsMassFitter::vnFunc,fMassMin,fMassMax,nparsvn,"AliHFVnVsMassFitter","vnFunc");
  SetParNames();

  ROOT::Math::WrappedMultiTF1 wfTotMass(*fMassTotFunc,1);
  ROOT::Math::WrappedMultiTF1 wfTotVn(*fVnTotFunc,1);

  // set data options and ranges
  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange rangeMass; //same range for two functions

  rangeMass.SetRange(fMassMin,fMassMax);
  ROOT::Fit::BinData dataMass(opt,rangeMass);
  ROOT::Fit::FillData(dataMass, fMassHisto);
  ROOT::Fit::BinData dataVn(opt,rangeMass);
  ROOT::Fit::FillData(dataVn, fVnVsMassHisto);

  //define the 2 chi squares
  ROOT::Fit::Chi2Function chi2Mass(dataMass, wfTotMass);
  ROOT::Fit::Chi2Function chi2Vn(dataVn, wfTotVn);

  //define the global chi square
  AliHFGlobalChi2 globalChi2(chi2Mass, chi2Vn);

  //define fitter
  ROOT::Fit::Fitter fitter;
  // create before the parameter settings in order to fix or set range on them
  fitter.Config().SetParamsSettings(nparsvn,initpars.data()); //set initial parameters from prefits
  if(fMeanFixed==2 || fMeanFixedFromMassFit) {fitter.Config().ParSettings(fNParsMassBkg+1).Fix();}
  if(fSigmaFixed==2 || fSigmaFixedFromMassFit) {fitter.Config().ParSettings(fNParsMassBkg+2).Fix();}
  if(fSecondPeak && fFixSecMass) {fitter.Config().ParSettings(fNParsMassBkg+fNParsMassSgn+1).Fix();}
  if(fSecondPeak && fFixSecWidth) {fitter.Config().ParSettings(fNParsMassBkg+fNParsMassSgn+2).Fix();}
  if(fReflections && fFixRflOverSig) {fitter.Config().ParSettings(fNParsMassBkg+fNParsMassSgn+fNParsSec).Fix();}

  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");
  for(Int_t iPar=0; iPar<nparsvn; iPar++) {fitter.Config().ParSettings(iPar).SetName(fVnTotFunc->GetParName(iPar));}
  // fit FCN function directly
  // (specify optionally data size and flag to indicate that is a chi2 fit
  fitter.FitFCN(nparsvn,globalChi2,0,dataMass.Size()+dataVn.Size(),kFALSE);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  //set parameters in every function
  fVnBkgFunc = new TF1("fVnBkgFunc",this,&AliHFVnVsMassFitter::vnBkgFunc,fMassMin,fMassMax,fNParsVnBkg,"AliHFVnVsMassFitter","vnBkgFunc");
  fMassBkgFunc = new TF1("fMassBkgFunc",this,&AliHFVnVsMassFitter::MassBkg,fMassMin,fMassMax,fNParsMassBkg,"AliHFVnVsMassFitter","MassBkg");
  fMassSgnFunc = new TF1("fMassSgnFunc",this,&AliHFVnVsMassFitter::MassSignal,fMassMin,fMassMax,fNParsMassSgn,"AliHFVnVsMassFitter","MassSignal");
  if(fReflections) {fMassRflFunc = new TF1("fMassRflFunc",this,&AliHFVnVsMassFitter::MassRfl,fMassMin,fMassMax,fNParsRfl,"AliHFVnVsMassFitter","MassRfl");}
  if(fReflections) {fMassBkgRflFunc = new TF1("fMassBkgRflFunc",this,&AliHFVnVsMassFitter::MassBkgRfl,fMassMin,fMassMax,fNParsMassBkg+fNParsRfl,"AliHFVnVsMassFitter","MassBkgRfl");}
  if(fSecondPeak) {fMassSecPeakFunc = new TF1("fMassSecPeakFunc",this,&AliHFVnVsMassFitter::MassSecondPeak,fMassMin,fMassMax,fNParsSec,"AliHFVnVsMassFitter","MassSecondPeak");}
  for(Int_t iPar=0; iPar<nparsvn; iPar++) {
    fVnTotFunc->SetParameter(iPar,result.Parameter(iPar));
    fVnTotFunc->SetParError(iPar,result.ParError(iPar));
    if(iPar<nparsmass) {
      fMassTotFunc->SetParameter(iPar,result.Parameter(iPar));
      fMassTotFunc->SetParError(iPar,result.ParError(iPar));
    }
    if(iPar>=nparsmass && iPar<nparsvn-NvnParsSgn) {
      fVnBkgFunc->SetParameter(iPar-nparsmass,result.Parameter(iPar));
      fVnBkgFunc->SetParError(iPar-nparsmass,result.ParError(iPar));
    }
    if(iPar>=fNParsMassBkg && iPar<fNParsMassBkg+fNParsMassSgn) {
      fMassSgnFunc->SetParameter(iPar-fNParsMassBkg,result.Parameter(iPar));
    }
    if(iPar<fNParsMassBkg) {
      fMassBkgFunc->SetParameter(iPar,result.Parameter(iPar));
      fMassBkgFunc->SetParError(iPar,result.ParError(iPar));
      if(fReflections) {
        fMassRflFunc->SetParameter(iPar,result.Parameter(iPar));
        fMassRflFunc->SetParError(iPar,result.ParError(iPar));
      }
    }
    if(fReflections && (iPar>=fNParsMassBkg+fNParsMassSgn+fNParsSec && iPar<fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl)) {
      fMassRflFunc->SetParameter(iPar-(fNParsMassBkg+fNParsMassSgn+fNParsSec),result.Parameter(iPar));
      fMassRflFunc->SetParError(iPar-(fNParsMassBkg+fNParsMassSgn+fNParsSec),result.ParError(iPar));
    }
    if(fSecondPeak && (iPar>=fNParsMassBkg+fNParsMassSgn && iPar<fNParsMassBkg+fNParsMassSgn+fNParsSec)) {
      fMassSecPeakFunc->SetParameter(iPar-(fNParsMassBkg+fNParsMassSgn),result.Parameter(iPar));
      fMassSecPeakFunc->SetParError(iPar-(fNParsMassBkg+fNParsMassSgn),result.ParError(iPar));
    }
  }
  if(drawFit) {DrawFit();}

  fVn = fVnTotFunc->GetParameter(fVnTotFunc->GetNpar()-NvnParsSgn);
  fVnUncertainty = fVnTotFunc->GetParError(fVnTotFunc->GetNpar()-NvnParsSgn);
  if(fDoSecondPeakVn) {
    fVnSecPeak = fVnTotFunc->GetParameter(fVnTotFunc->GetNpar()-1);
    fVnSecPeakUncertainty = fVnTotFunc->GetParError(fVnTotFunc->GetNpar()-1);
  }
  fRawYield = fVnTotFunc->GetParameter(fNParsMassBkg)/fMassHisto->GetBinWidth(10);
  fRawYieldUncertainty = fVnTotFunc->GetParError(fNParsMassBkg)/fMassHisto->GetBinWidth(10);
  fMean = fVnTotFunc->GetParameter(fNParsMassBkg+1);
  fMeanUncertainty = fVnTotFunc->GetParError(fNParsMassBkg+1);
  fSigma = fVnTotFunc->GetParameter(fNParsMassBkg+2);
  fSigmaUncertainty = fVnTotFunc->GetParError(fNParsMassBkg+2);
  fChiSquare = result.MinFcnValue();
  fNDF = result.Ndf();
  fProb = result.Prob();

  return kTRUE;
}

//______________________________________________________________________________
void AliHFVnVsMassFitter::DrawHere(TVirtualPad* c){
  /// Core method to draw the fit output

  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  c->Divide(1,2);

  c->cd(1);
  fMassHisto->SetTitle("");
  fMassHisto->SetMarkerStyle(20);
  fMassHisto->SetMarkerSize(1);
  fMassHisto->SetMarkerColor(kBlack);
  fMassHisto->SetLineColor(kBlack);
  fMassHisto->GetYaxis()->SetRangeUser(0.,fMassHisto->GetMaximum()*1.2);
  fMassHisto->GetXaxis()->SetRangeUser(fMassMin,fMassMax);
  fMassHisto->Draw("E");
  if(fMassFuncFromPrefit) {
    fMassFuncFromPrefit->SetLineColor(kGray+1);
    fMassFuncFromPrefit->SetLineStyle(7);
    fMassFuncFromPrefit->SetRange(fMassMin,fMassMax);
    fMassFuncFromPrefit->Draw("same");
  }
  if(fMassBkgFunc) {
    fMassBkgFunc->SetLineColor(kRed);
    fMassBkgFunc->SetRange(fMassMin,fMassMax);
    fMassBkgFunc->Draw("same");
  }
  if(fMassRflFunc) {
    fMassRflFunc->SetLineColor(kGreen+1);
    fMassRflFunc->SetRange(fMassMin,fMassMax);
    fMassRflFunc->Draw("same");
  }
  if(fMassBkgRflFunc) {
    fMassBkgRflFunc->SetLineColor(kRed+1);
    fMassBkgRflFunc->SetLineStyle(7);
    fMassBkgRflFunc->SetRange(fMassMin,fMassMax);
    fMassBkgRflFunc->Draw("same");
  }
  if(fMassSecPeakFunc) {
    fMassSecPeakFunc->SetLineColor(kMagenta+1);
    fMassSecPeakFunc->SetLineStyle(7);
    fMassSecPeakFunc->SetRange(fMassMin,fMassMax);
    fMassSecPeakFunc->Draw("same");
  }
  if(fMassTotFunc) {
    fMassTotFunc->SetLineColor(kBlue);
    fMassTotFunc->SetRange(fMassMin,fMassMax);
    fMassTotFunc->Draw("same");
  }
  TPaveText* massinfo = new TPaveText(0.45,0.7,1.,0.87,"NDC");
  massinfo->SetTextFont(42);
  massinfo->SetTextSize(0.05);
  massinfo->SetBorderSize(0);
  massinfo->SetFillStyle(0);
  massinfo->SetTextColor(kBlue);
  massinfo->AddText(Form("mean = %.3f #pm %.3f",fVnTotFunc->GetParameter(fNParsMassBkg+1),fVnTotFunc->GetParError(fNParsMassBkg+1)));
  massinfo->AddText(Form("sigma = %.3f #pm %.3f",fVnTotFunc->GetParameter(fNParsMassBkg+2),fVnTotFunc->GetParError(fNParsMassBkg+2)));
  if(fMassSgnFuncType==k2Gaus) {
    massinfo->AddText(Form("sigma2 = %.3f #pm %.3f",fVnTotFunc->GetParameter(fNParsMassBkg+3),fVnTotFunc->GetParError(fNParsMassBkg+3)));
  }
  massinfo->Draw("same");

  c->cd(2);
  fVnVsMassHisto->SetTitle("");
  fVnVsMassHisto->SetMarkerStyle(20);
  fVnVsMassHisto->SetMarkerSize(1);
  fVnVsMassHisto->SetMarkerColor(kBlack);
  fVnVsMassHisto->SetLineColor(kBlack);
  fVnVsMassHisto->GetYaxis()->SetRangeUser(fVnVsMassHisto->GetMinimum()-0.15,fVnVsMassHisto->GetMaximum()+0.20);
  fVnVsMassHisto->GetXaxis()->SetRangeUser(fMassMin,fMassMax);
  fVnBkgFuncSb->SetRange(fMassMin,fMassMax);
  fVnVsMassHisto->Draw("E");
  if(fVnBkgFuncSb) {
    fVnBkgFuncSb->SetLineColor(kGray+1);
    fVnBkgFuncSb->SetLineStyle(7);
    fVnBkgFuncSb->Draw("same");
  }
  if(fVnBkgFunc) {
    fVnBkgFunc->SetLineColor(kRed);
    fVnBkgFunc->SetRange(fMassMin,fMassMax);
    fVnBkgFunc->Draw("same");
  }
  if(fVnTotFunc) {
    fVnTotFunc->SetLineColor(kBlue);
    fVnTotFunc->SetRange(fMassMin,fMassMax);
    fVnTotFunc->Draw("same");
  }

  TPaveText* vninfo = new TPaveText(-0.45,0.7,1.,0.87,"NDC");
  vninfo->SetTextFont(42);
  vninfo->SetTextSize(0.05);
  vninfo->SetBorderSize(0);
  vninfo->SetFillStyle(0);
  Int_t NvnParsSgn = 1;
  if(fSecondPeak && fDoSecondPeakVn) {NvnParsSgn+=1;}
  vninfo->AddText(Form("#it{v}_{%d}^{sgn} = %.3f #pm %.3f",fHarmonic,fVnTotFunc->GetParameter(fVnTotFunc->GetNpar()-NvnParsSgn),fVnTotFunc->GetParError(fVnTotFunc->GetNpar()-NvnParsSgn)));
  if(fSecondPeak && fDoSecondPeakVn) {vninfo->AddText(Form("#it{v}_{%d}^{sec peak} = %.3f #pm %.3f",fHarmonic,fVnTotFunc->GetParameter(fVnTotFunc->GetNpar()-1),fVnTotFunc->GetParError(fVnTotFunc->GetNpar()-1)));}
  vninfo->AddText(Form("#chi^{2}/#it{ndf} = %.2f/%d",fChiSquare,fNDF));
  vninfo->Draw("same");

  c->Update();
}

//________________________________________________________________
Bool_t AliHFVnVsMassFitter::MassPrefit() {

  //define proper maxs and mins from histos
  Double_t tmpmin = TMath::Max(fMassHisto->GetBinLowEdge(1),fVnVsMassHisto->GetBinLowEdge(1));
  fMassMin=TMath::Max(fMassMin,tmpmin);
  Double_t tmpmax = TMath::Min(fMassHisto->GetBinLowEdge(fMassHisto->GetNbinsX())+fMassHisto->GetBinWidth(fMassHisto->GetNbinsX()),fVnVsMassHisto->GetBinLowEdge(fVnVsMassHisto->GetNbinsX())+fVnVsMassHisto->GetBinWidth(fVnVsMassHisto->GetNbinsX()));
  fMassMax=TMath::Min(fMassMax,tmpmax);

  fMassFitter = new AliHFInvMassFitter(fMassHisto,fMassMin,fMassMax,fMassBkgFuncType,fMassSgnFuncType);
  if(fSigmaFixed==1) fMassFitter->SetInitialGaussianSigma(fSigmaInit);
  else if(fSigmaFixed==2) fMassFitter->SetFixGaussianSigma(fSigmaInit);
  if(fMeanFixed==1) fMassFitter->SetInitialGaussianMean(fMeanInit);
  else if(fMeanFixed==2) fMassFitter->SetFixGaussianMean(fMeanInit);
  fMassFitter->SetUseLikelihoodFit();
  if(fMassBkgFuncType==kPoln) {fMassFitter->SetPolDegreeForBackgroundFit(fPolDegreeBkg);}
  if(fSecondPeak) {fMassFitter->IncludeSecondGausPeak(fSecMass,fFixSecMass,fSecWidth,fFixSecWidth);}
  if(fReflections) {
    fHistoTemplRfl = (TH1F*)fMassFitter->SetTemplateReflections(fHistoTemplRflInit,fRflOpt,fMinRefl,fMaxRefl);
    if(fRflOverSig>0) {fMassFitter->SetInitialReflOverS(fRflOverSig);}
    if(fFixRflOverSig) {fMassFitter->SetFixReflOverS(fRflOverSig);}
  }
  Bool_t status = fMassFitter->MassFitter(kFALSE);

  if(status) {
    fMassFuncFromPrefit = (TF1*)fMassFitter->GetMassFunc();
    fMassFuncFromPrefit->SetName("fMassFuncFromPrefit");
  }

  return status;
}

//________________________________________________________________
Bool_t AliHFVnVsMassFitter::VnSBPrefit() {

  Double_t mean = fMassFitter->GetMean();
  Double_t sigma = fMassFitter->GetSigma();
  const Int_t nMassBins = fVnVsMassHisto->GetNbinsX();
  Double_t SBbins[nMassBins];
  Int_t nSBbins=0;
  for(Int_t iBin=0; iBin<nMassBins; iBin++) {
    Double_t min = fVnVsMassHisto->GetBinLowEdge(iBin+1);
    Double_t max = fVnVsMassHisto->GetBinLowEdge(iBin+1)+fVnVsMassHisto->GetBinWidth(iBin+1);
    if(max<mean){
      if(max<(mean-fNSigmaForSB*sigma)) {SBbins[iBin]=1; nSBbins++;}
      else {SBbins[iBin]=0;}
    }
    if(min>=mean){
      if(min>(mean+fNSigmaForSB*sigma)) {SBbins[iBin]=1; nSBbins++;}
      else {SBbins[iBin]=0;}
    }
  }
  TGraphErrors* gVnVsMassSB = new TGraphErrors(nSBbins);
  for(Int_t iBin=0; iBin<nMassBins; iBin++) {
    if(SBbins[iBin]==1) {
      gVnVsMassSB->SetPoint(iBin,fVnVsMassHisto->GetBinCenter(iBin+1),fVnVsMassHisto->GetBinContent(iBin+1));
      gVnVsMassSB->SetPointError(iBin,fVnVsMassHisto->GetBinWidth(iBin+1)/2,fVnVsMassHisto->GetBinError(iBin+1));
    }
  }
  fVnBkgFuncSb = new TF1("fVnBkgFuncSb",this,&AliHFVnVsMassFitter::vnBkgFunc,fMassMin,fMassMax,fNParsVnBkg,"AliHFVnVsMassFitter","vnBkgFunc");
  switch(fVnBkgFuncType) {
    case 1:
      fVnBkgFuncSb->SetParName(0,"ConstVnBkg");
      fVnBkgFuncSb->SetParName(1,"SlopeVnBkg");
      break;
    case 2:
      fVnBkgFuncSb->SetParName(0,"ConstVnBkg");
      fVnBkgFuncSb->SetParName(1,"Coef1VnBkg");
      fVnBkgFuncSb->SetParName(2,"Coef2VnBkg");
      break;
    default:
      AliError("Error in setting signal par names: check fVnBkgFuncType");
      break;
  }
  gVnVsMassSB->Fit(fVnBkgFuncSb,"","",fMassMin,fMassMax);
  Bool_t status=kFALSE;
  if(fVnBkgFuncSb->GetChisquare()<1000) status=kTRUE;

  delete gVnVsMassSB;
  return status;
}

//________________________________________________________________
void AliHFVnVsMassFitter::DefineNumberOfParameters() {

  switch(fMassSgnFuncType) {
    case 0: //single gaus
      fNParsMassSgn=3;
      break;
    case 1: //double gaus
      fNParsMassSgn=5;
      break;
    default:
      AliError("Error in computing fMassSgnFuncType: check fMassSgnFuncType");
      break;
  }

  switch(fMassBkgFuncType) {
    case 0: //expo
      fNParsMassBkg=2;
      break;
    case 1: //lin
      fNParsMassBkg=2;
      break;
    case 2: //pol2
      fNParsMassBkg=3;
      break;
    case 3: //no bkg
      fNParsMassBkg=1;
      break;
    case 4: //power law
      fNParsMassBkg=2;
      break;
    case 5: //power expo
      fNParsMassBkg=3;
      break;
    case 6: //high degree pol
      fNParsMassBkg=fPolDegreeBkg+1;
      break;
    default:
      AliError("Error in computing fNParsMassBkg: check fMassBkgFuncType");
      break;
  }

  switch(fVnBkgFuncType) {
    case 1: //lin
      fNParsVnBkg=2;
      break;
    case 2: //pol2
      fNParsVnBkg=3;
      break;
    default:
      AliError("Error in computing fNParsVnBkg: check fVnBkgFuncType");
      break;
  }

  if(fReflections) fNParsRfl=1;
  else fNParsRfl=0;

  if(fSecondPeak) fNParsSec=3;
  else fNParsSec=0;
}

//________________________________________________________________
void AliHFVnVsMassFitter::SetParNames() {

  switch(fMassSgnFuncType) {
    case 0: //single gaus
      fVnTotFunc->SetParName(fNParsMassBkg,"SgnInt");
      fVnTotFunc->SetParName(fNParsMassBkg+1,"Mean");
      fVnTotFunc->SetParName(fNParsMassBkg+2,"Sigma");
      break;
    case 1: //double gaus
      fVnTotFunc->SetParName(fNParsMassBkg,"SgnInt");
      fVnTotFunc->SetParName(fNParsMassBkg+1,"Mean");
      fVnTotFunc->SetParName(fNParsMassBkg+2,"Sigma1");
      fVnTotFunc->SetParName(fNParsMassBkg+3,"Frac");
      fVnTotFunc->SetParName(fNParsMassBkg+4,"Sigma2");
      break;
    default:
      AliError("Error in setting signal par names: check fMassSgnFuncType");
      break;
  }

  switch(fMassBkgFuncType) {
    case 0: //expo
      fVnTotFunc->SetParName(0,"BkgInt");
      fVnTotFunc->SetParName(1,"Slope");
      break;
    case 1: //lin
      fVnTotFunc->SetParName(0,"BkgInt");
      fVnTotFunc->SetParName(1,"Slope");
      break;
    case 2: //pol2
      fVnTotFunc->SetParName(0,"BkgInt");
      fVnTotFunc->SetParName(1,"Coef1");
      fVnTotFunc->SetParName(2,"Coef1");
      break;
    case 3: //no bkg
      fVnTotFunc->SetParName(0,"Const");
      break;
    case 4: //power law
      fVnTotFunc->SetParName(0,"BkgInt");
      fVnTotFunc->SetParName(1,"Coef1");
      break;
    case 5: //power expo
      fVnTotFunc->SetParName(0,"BkgInt");
      fVnTotFunc->SetParName(1,"Coef1");
      fVnTotFunc->SetParName(2,"Coef2");
      break;
    case 6: //high degree pol
      fVnTotFunc->SetParName(0,"BkgInt");
      for(Int_t iPar=1; iPar<fNParsMassBkg; iPar++) {fVnTotFunc->SetParName(iPar,Form("Coef%d",iPar));}
    break;
    default:
      AliError("Error in setting signal par names: check fMassBkgFuncType");
      break;
  }

  for(Int_t iPar=0; iPar<fNParsVnBkg; iPar++) {fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl+iPar,fVnBkgFuncSb->GetParName(iPar));}

  if(fReflections) {fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+fNParsSec,"ReflOverS");}

  if(fSecondPeak) {
    fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn,"SecPeakInt");
    fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+1,"SecPeakMean");
    fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+2,"SecPeakSigma");
  }
  fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsVnBkg,Form("v%dSgn",fHarmonic));
  if(fSecondPeak && fDoSecondPeakVn) {fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsVnBkg+1,Form("v%dSecPeak",fHarmonic));}
}

//_________________________________________________________________________
void AliHFVnVsMassFitter::Signal(Double_t nOfSigma,Double_t &signal,Double_t &errsignal) const {
  /// Return signal integral in mean +- n sigma
  ///

  Double_t minMass=fMean-nOfSigma*fSigma;
  Double_t maxMass=fMean+nOfSigma*fSigma;
  Signal(minMass,maxMass,signal,errsignal);
  return;
}

//_________________________________________________________________________
void AliHFVnVsMassFitter::Signal(Double_t min, Double_t max, Double_t &signal,Double_t &errsignal) const {
  /// Return signal integral in a range
  ///
  if(!fMassSgnFunc) {signal=-1; errsignal=0; return;}

  signal=fMassSgnFunc->Integral(min, max)/(Double_t)fMassHisto->GetBinWidth(1);
  errsignal=(fRawYieldUncertainty/fRawYield)*signal;/*assume relative error is the same as for total integral*/

  return;
}

//___________________________________________________________________________
void AliHFVnVsMassFitter::Background(Double_t nOfSigma,Double_t &background,Double_t &errbackground) const {
  /// Return background integral in mean +- n sigma
  ///

  Double_t minMass=fMean-nOfSigma*fSigma;
  Double_t maxMass=fMean+nOfSigma*fSigma;
  Background(minMass,maxMass,background,errbackground);

  return;
}

//___________________________________________________________________________
void AliHFVnVsMassFitter::Background(Double_t min, Double_t max, Double_t &background,Double_t &errbackground) const {
  /// Return background integral in a range
  ///

  if(!fMassBkgFunc) {background=-1; errbackground=0; return;}

  Double_t intB=fMassBkgFunc->GetParameter(0);
  Double_t intBerr=fMassBkgFunc->GetParError(0);
  //relative error evaluation: from histo

  Int_t leftBand=fMassHisto->FindBin(fMean-4*fSigma);
  Int_t rightBand=fMassHisto->FindBin(fMean+4*fSigma);
  intB=fMassHisto->Integral(1,leftBand)+fMassHisto->Integral(rightBand,fMassHisto->GetNbinsX());
  Double_t sum2=0;
  for(Int_t iBin=1; iBin<=leftBand; iBin++){
    sum2+=fMassHisto->GetBinError(iBin)*fMassHisto->GetBinError(iBin);
  }
  for(Int_t iBin=rightBand; iBin<=fMassHisto->GetNbinsX(); iBin++){
    sum2+=fMassHisto->GetBinError(iBin)*fMassHisto->GetBinError(iBin);
  }

  intBerr=TMath::Sqrt(sum2);

  background=fMassBkgFunc->Integral(min,max)/(Double_t)fMassHisto->GetBinWidth(1);
  errbackground=intBerr/intB*background;

  return;
}

//__________________________________________________________________________

void AliHFVnVsMassFitter::Significance(Double_t nOfSigma,Double_t &significance,Double_t &errsignificance) const  {
  /// Return significance in mean +- n sigma
  ///

  Double_t minMass=fMean-nOfSigma*fSigma;
  Double_t maxMass=fMean+nOfSigma*fSigma;
  Significance(minMass, maxMass, significance, errsignificance);

  return;
}

//__________________________________________________________________________

void AliHFVnVsMassFitter::Significance(Double_t min, Double_t max, Double_t &significance,Double_t &errsignificance) const {
  /// Return significance integral in a range
  ///

  Double_t background,errbackground;
  Background(min,max,background,errbackground);

  if (fRawYield+background <= 0.){
    significance=-1;
    errsignificance=0;
    return;
  }

  AliVertexingHFUtils::ComputeSignificance(fRawYield,fRawYieldUncertainty,background,errbackground,significance,errsignificance);

  return;
}

//________________________________________________________________
Double_t AliHFVnVsMassFitter::GetGausPDF(Double_t x, Double_t mean, Double_t sigma) {

  return TMath::Gaus(x,mean,sigma,kTRUE);
}

//________________________________________________________________
Double_t AliHFVnVsMassFitter::GetExpoPDF(Double_t x, Double_t coeff) {

  return TMath::Exp(x/coeff)/(coeff*(TMath::Exp(fMassMax/coeff)-TMath::Exp(fMassMin/coeff)));
}

//________________________________________________________________
Double_t AliHFVnVsMassFitter::GetPolPDF(Double_t x, Double_t *pars, Int_t order, Bool_t isnorm) {

  switch(order) {
    case 0:
      if(isnorm) {return 1./(fMassMax-fMassMin);}
      else {return pars[0];}
      break;
    case 1:
      if(isnorm) {return (pars[0]-pars[1]/2*(fMassMax*fMassMax-fMassMin*fMassMin))/(fMassMax-fMassMin)+pars[1]*x;}
      else {return pars[0]+pars[1]*x;}
      break;
    case 2:
      if(isnorm) {return (pars[0]-pars[1]/2*(fMassMax*fMassMax-fMassMin*fMassMin)-pars[2]/3*(fMassMax*fMassMax*fMassMax-fMassMin*fMassMin*fMassMin))/(fMassMax-fMassMin)+pars[1]*x;}
      else {return pars[0]+pars[1]*x+pars[2]*x*x;}
  }
  return 0;
}

//________________________________________________________________
Double_t AliHFVnVsMassFitter::GetPowerFuncPDF(Double_t x, Double_t *pars) {

  Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  return pars[0]*(pars[1]+1.)/(TMath::Power(fMassMax-mpi,pars[1]+1.)-TMath::Power(fMassMin-mpi,pars[1]+1.))*TMath::Power(x-mpi,pars[1]);
}

//________________________________________________________________
Double_t AliHFVnVsMassFitter::GetPowerExpoPDF(Double_t x, Double_t *pars) {

  Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  return pars[0]*TMath::Sqrt(x - mpi)*TMath::Exp(-1.*pars[1]*(x-mpi));
}

//________________________________________________________________
Double_t AliHFVnVsMassFitter::GetHigherPolFuncPDF(Double_t x, Double_t *pars) {

  Double_t total=pars[0];
  for(Int_t iT=1; iT<=fPolDegreeBkg; iT++){
    total+=pars[iT]*TMath::Power(x-fMassParticle,iT)/TMath::Factorial(iT);
  }
  return total;
}

//________________________________________________________________
Double_t AliHFVnVsMassFitter::MassSignal(Double_t *m, Double_t *pars) {

  switch(fMassSgnFuncType) {
    case 0:
      return pars[0]*GetGausPDF(m[0],pars[1],pars[2]);
      break;
    case 1:
      return pars[0]*(pars[3]*GetGausPDF(m[0],pars[1],pars[2])+(1-pars[3])*GetGausPDF(m[0],pars[1],pars[4]));
      break;
  }
  fRawYieldHelp=pars[0]/fMassHisto->GetBinWidth(1);

  return 0;
}

//________________________________________________________________
Double_t AliHFVnVsMassFitter::MassBkg(Double_t *m, Double_t *pars) {

  switch(fMassBkgFuncType) {
    case 0: //exponential
      return pars[0]*GetExpoPDF(m[0],pars[1]);
      break;
    case 1: //linear
      return GetPolPDF(m[0],pars,1,kTRUE);
      break;
    case 2: //parabolic
      return GetPolPDF(m[0],pars,2,kTRUE);
      break;
    case 3: //constant
      return GetPolPDF(m[0],pars,0,kTRUE);
      break;
    case 4: //power law
      return GetPowerFuncPDF(m[0],pars);
      break;
    case 5: //power law expo
      return GetPowerExpoPDF(m[0],pars);
      break;
    case 6: //higher order (>=3) polinomial
      return GetHigherPolFuncPDF(m[0],pars);
      break;
  }
  return 0;
}

//_________________________________________________________________________
Double_t AliHFVnVsMassFitter::MassRfl(Double_t *m,Double_t *pars){
  /// Fit function for reflections:
  /// D0->Kpi decays with swapped mass assignment to pion and kaon decay tracks
  if(!fHistoTemplRfl) return 0;

  Int_t bin =fHistoTemplRfl->FindBin(m[0]);
  Double_t value=fHistoTemplRfl->GetBinContent(bin);
  Int_t binmin=fHistoTemplRfl->FindBin(fMassMin*1.00001);
  Int_t binmax=fHistoTemplRfl->FindBin(fMassMax*0.99999);
  Double_t norm=fHistoTemplRfl->Integral(binmin,binmax)*fHistoTemplRfl->GetBinWidth(bin);
  if(TMath::Abs(value)<1.e-14 && fSmoothRfl){// very rough, assume a constant trend, much better would be a pol1 or pol2 over a broader range
    value+=fHistoTemplRfl->GetBinContent(bin-1)+fHistoTemplRfl->GetBinContent(bin+1);
    value/=3.;
  }
  return pars[0]*value/norm*fRawYieldHelp*fMassHisto->GetBinWidth(1);
}

//_________________________________________________________________________
Double_t AliHFVnVsMassFitter::MassBkgRfl(Double_t *m,Double_t *pars){

  if(!fHistoTemplRfl) {return MassBkg(m,pars);}
  else {
    //bkg mass parameters
    const Int_t nBkgPars = fNParsMassBkg;
    Double_t bkgpars[nBkgPars];
    for(Int_t iPar=0; iPar<fNParsMassBkg; iPar++) {bkgpars[iPar] = pars[iPar];}
    //reflection parameters
    Double_t rflpars[1]; //maximum number of parameters for rfl = 1 for the implemented functions
    for(Int_t iPar=0; iPar<fNParsRfl; iPar++) {rflpars[iPar] = pars[iPar+fNParsMassBkg];}
    return MassBkg(m,bkgpars)+MassRfl(m,rflpars);
  }
}

//_________________________________________________________________________
Double_t AliHFVnVsMassFitter::MassSecondPeak(Double_t *m,Double_t *pars){
  /// Fit function for a second gaussian peak
  /// To be used, e.g., for D+->KKpi in the Ds mass spectrum

  return pars[0]*GetGausPDF(m[0],pars[1],pars[2]);
}

//________________________________________________________________
Double_t AliHFVnVsMassFitter::vnBkgFunc(Double_t *m, Double_t *pars) {

  switch(fVnBkgFuncType) {
    case 1: //linear
      return GetPolPDF(m[0],pars,1,kFALSE);
      break;
    case 2: //parabolic
      return GetPolPDF(m[0],pars,2,kFALSE);
      break;
  }
  return 0;
}

//________________________________________________________________
Double_t AliHFVnVsMassFitter::MassFunc(Double_t *m, Double_t *pars) {

  //bkg mass parameters
  const Int_t nBkgPars = fNParsMassBkg;
  Double_t bkgpars[nBkgPars];
  for(Int_t iPar=0; iPar<fNParsMassBkg; iPar++) {bkgpars[iPar] = pars[iPar];}
  //signal mass parameters
  Double_t sgnpars[5]; //maximum number of parameters for sgn = 5 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsMassSgn; iPar++) {sgnpars[iPar] = pars[iPar+fNParsMassBkg];}
  //second peak parameters
  Double_t secpeakpars[3]; //maximum number of parameters for second peak = 3 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsSec; iPar++) {secpeakpars[iPar] = pars[iPar+fNParsMassBkg+fNParsMassSgn];}
  //reflection parameters
  Double_t rflpars[1]; //maximum number of parameters for rfl = 1 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsRfl; iPar++) {rflpars[iPar] = pars[iPar+fNParsMassBkg+fNParsMassSgn+fNParsSec];}

  Double_t total = MassSignal(m,sgnpars)+MassBkg(m,bkgpars);
  if(fSecondPeak) {total += MassSecondPeak(m,secpeakpars);}
  if(fReflections) {total += MassRfl(m,rflpars);}

  return total;
}

//________________________________________________________________
Double_t AliHFVnVsMassFitter::vnFunc(Double_t *m, Double_t *pars) {

  //bkg mass parameters
  const Int_t nBkgPars = fNParsMassBkg;
  Double_t massbkgpars[nBkgPars];
  for(Int_t iPar=0; iPar<fNParsMassBkg; iPar++) {massbkgpars[iPar] = pars[iPar];}
  //signal mass parameters
  Double_t masssgnpars[5]; //maximum number of parameters for mass sgn = 5 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsMassSgn; iPar++) {masssgnpars[iPar] = pars[iPar+fNParsMassBkg];}
  //second peak parameters
  Double_t secpeakpars[3]; //maximum number of parameters for second peak = 3 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsSec; iPar++) {secpeakpars[iPar] = pars[iPar+fNParsMassBkg+fNParsMassSgn];}
  //reflection parameters
  Double_t rflpars[1]; //maximum number of parameters for rfl = 1 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsRfl; iPar++) {rflpars[iPar] = pars[iPar+fNParsMassBkg+fNParsMassSgn+fNParsSec];}
  //bkg vn parameters
  Double_t vnbkgpars[3]; //maximum number of parameters for vn bkg = 3 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsVnBkg; iPar++) {vnbkgpars[iPar] = pars[iPar+fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl];}
  //signal vn parameter
  Double_t vnSgn = pars[fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsVnBkg];
  //second peak vn parameter
  Double_t vnSecPeak = 0;
  if(fSecondPeak && fDoSecondPeakVn) {vnSecPeak = pars[fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsVnBkg+1];}


  Double_t vnBkg = vnBkgFunc(m,vnbkgpars);
  Double_t Sgn = MassSignal(m,masssgnpars);
  Double_t Bkg = MassBkg(m,massbkgpars);
  if(fReflections) {Bkg += MassRfl(m,rflpars);}
  Double_t SecPeak = 0;
  if(fSecondPeak) {SecPeak += MassSecondPeak(m,secpeakpars);}

  if(fSecondPeak && fDoSecondPeakVn) {return (vnSgn*Sgn+vnBkg*Bkg+vnSecPeak*SecPeak)/(Sgn+Bkg+SecPeak);}
  else {return (vnSgn*Sgn+vnBkg*(Bkg+SecPeak))/(Sgn+Bkg+SecPeak);}
}

//______________________________________________________________________________
void AliHFVnVsMassFitter::DrawFit(){
  /// Steering method to draw the fit output
  ///
  TCanvas* c0=new TCanvas("c0","",600,800);
  DrawHere(c0);
}
