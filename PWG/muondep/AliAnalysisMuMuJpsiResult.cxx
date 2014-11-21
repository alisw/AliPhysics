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

///
/// Class to hold results about J/psi
/// like number of of J/psi (before and after Acc x Eff correction),
/// Acc x Eff correction, Yield, RAB, etc...
///
/// A note on "naming conventions"
///
/// FitFunctionXY : denotes a function with a prototype double func(double*,double*)
/// which can be used in a fit. X = Background, Signal or Total for Background+Signal
/// Y is the functio name
///
///
/// author: Laurent Aphecetche (Subatech)
///

#include "AliAnalysisMuMuJpsiResult.h"

ClassImp(AliAnalysisMuMuJpsiResult)

#include "TF1.h"
#include "TProfile.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TH2.h"
#include "THashList.h"
#include "TLine.h"
#include "TList.h"
#include "TMap.h"
#include "TMath.h"
#include "TMethodCall.h"
#include "TObjArray.h"
#include "TParameter.h"
#include "AliAnalysisMuMuBinning.h"
#include "AliLog.h"
#include <map>

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"


namespace {
  
  //____________________________________________________________________________
  const std::map<std::string,Double_t>& MassMap()
  {
    /// a simple map of masses...
    static std::map<std::string,Double_t> massMap;
    // not super elegant, but that way we do not depend on TDatabasePDG and thus
    // can decide on our particle naming
    if (massMap.empty())
    {
      massMap["JPsi"] = 3.096916e+00;
      massMap["PsiP"] = 3.68609e+00;
      massMap["Upsilon"] = 9.46030e+00;
      massMap["UpsilonPrime"] = 1.00233e+01;
    }
    return massMap;
  }

  //____________________________________________________________________________
  Bool_t GetKeyValue(const TString& str, const char separator, TString& key, TString& value)
  {
    /// Get a key value pair, separated by separator character
    key=value="";
    if ( !str.CountChar(separator) ) return kFALSE;
    Int_t index = str.Index(separator);
    key = str(0,index);
    value = str(index+1,str.Length()-index-1);
    return kTRUE;
  }
    
  const TString kKeyFunc = "Func";
  const TString kKeyRange = "Range";
  const TString kKeyRebin = "Rebin";
  const TString kFitRangeLow = "FitRangeLow";
  const TString kFitRangeHigh = "FitRangeHigh";
  const TString kKeyCount = "Count";
  const TString kKeyHistoType = "HistoType";
  const TString kKeyTails = "Tails";
  const TString kKeySPsiP = "FSigmaPsiP"; //Factor to fix the psi' sigma to sigmaJPsi*SigmaPsiP (Usually factor SigmaPsiP = 1, 0.9 and 1.1)
  const TString kKeyMinvRS = "MinvRS"; // FIXME: not very correct since "MinvRS" is in AliAnalysisMuMu::GetParametersFromResult
  
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult::AliAnalysisMuMuJpsiResult(TRootIOCtor* /*io*/) :
AliAnalysisMuMuResult("",""),
fNofRuns(),
fNofTriggers(-1),
fHisto(0x0),
fBin(),
fTriggerClass(),
fEventSelection(),
fPairSelection(),
fCentralitySelection(),
fFitFunction(),
fFitRejectRangeLow(TMath::Limits<Double_t>::Max()),
fFitRejectRangeHigh(TMath::Limits<Double_t>::Max()),
fRejectFitPoints(kFALSE),
fParticle(""),
fMinvRS("")
{
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult::AliAnalysisMuMuJpsiResult(const char* particle,
                                                     const TH1& h,
                                                     const char* fitType)
:
AliAnalysisMuMuResult(fitType,""),
fNofRuns(1),
fNofTriggers(-1),
fHisto(0x0),
fBin(),
fTriggerClass(),
fEventSelection(),
fPairSelection(),
fCentralitySelection(),
fFitFunction(),
fFitRejectRangeLow(TMath::Limits<Double_t>::Max()),
fFitRejectRangeHigh(TMath::Limits<Double_t>::Max()),
fRejectFitPoints(kFALSE),
fParticle(particle),
fMinvRS("")
{
  SetHisto(h);
  
  DecodeFitType(fitType);
  
  TString name(fFitFunction);
  if ( !name.Contains("PSICB2") && !name.Contains("PSINA60NEW") && !name.Contains("PSICOUNT") ) //To avoid adding things to the name of simu results
  {
    Bool_t isMPt = kFALSE;
    
    if ( name.BeginsWith("PSI") ) name.ReplaceAll("PSIPSIPRIME","");
    else if ( name.BeginsWith("MPT") )
    {
      name.ReplaceAll("MPTPSIPSIPRIME","");
      isMPt = kTRUE;
    }
    name += "_";
    name += Form("%1.1f",GetValue(kFitRangeLow));
    name += "_";
    name += Form("%1.1f",GetValue(kFitRangeHigh));
    if ( !isMPt )
    {
      name += "_";
      name += Form("SP%1.1f",GetValue(kKeySPsiP));
    }
    else
    {
      name += Form("(Sig:%s)",fMinvRS.Data());
    }
  }
  SetName(name.Data());
  
  Int_t rebin = TMath::Nint(GetValue(kKeyRebin));
  
  if (rebin>0)
  {
    fHisto->Rebin(rebin);
  }
  
  if ( fHisto->GetEntries()<100 && !TString(GetName()).Contains(kKeyCount) )
  {
    // not enough statistics to perform a fit.
    Invalidate();
    std::cout << "Fit Excluded: Not enough statistics to perform a fit" << std::endl;
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult::AliAnalysisMuMuJpsiResult(const char* particle,
                                                     const TH1& h,
                                                     const char* triggerName,
                                                     const char* eventSelection,
                                                     const char* pairSelection,
                                                     const char* centSelection,
                                                     const AliAnalysisMuMuBinning::Range& bin)
:
AliAnalysisMuMuResult(Form("%s-%s-%s-%s",triggerName,eventSelection,pairSelection,centSelection),""),
fNofRuns(1),
fNofTriggers(-1),
fHisto(0x0),
fBin(bin),
fTriggerClass(triggerName),
fEventSelection(eventSelection),
fPairSelection(pairSelection),
fCentralitySelection(centSelection),
fFitFunction(),
fFitRejectRangeLow(TMath::Limits<Double_t>::Max()),
fFitRejectRangeHigh(TMath::Limits<Double_t>::Max()),
fRejectFitPoints(kFALSE),
fParticle(particle),
fMinvRS("")
{
  SetHisto(h);
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult::AliAnalysisMuMuJpsiResult(const AliAnalysisMuMuJpsiResult& rhs)
:
AliAnalysisMuMuResult(rhs),
fNofRuns(rhs.NofRuns()),
fNofTriggers(rhs.NofTriggers()),
fHisto(0x0),
fBin(rhs.Bin()),
fTriggerClass(rhs.fTriggerClass),
fEventSelection(rhs.fEventSelection),
fPairSelection(rhs.fPairSelection),
fCentralitySelection(rhs.fCentralitySelection),
fFitFunction(rhs.fFitFunction),
fFitRejectRangeLow(rhs.fFitRejectRangeLow),
fFitRejectRangeHigh(rhs.fFitRejectRangeHigh),
fRejectFitPoints(rhs.fRejectFitPoints),
fParticle(rhs.fParticle),
fMinvRS(rhs.fMinvRS)
{
  /// copy ctor
  /// Note that the mother is lost
  /// fKeys remains 0x0 so it will be recomputed if need be

  if ( rhs.fHisto )
  {
    fHisto = static_cast<TH1*>(rhs.fHisto->Clone());
  }  
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult& AliAnalysisMuMuJpsiResult::operator=(const AliAnalysisMuMuJpsiResult& rhs)
{
  /// Assignment operator
  
  if (this!=&rhs)
  {
    static_cast<AliAnalysisMuMuResult&>(*this) = static_cast<const AliAnalysisMuMuResult&>(rhs);
    delete fHisto;

    if ( rhs.fHisto )
    {
      fHisto = static_cast<TH1*>(rhs.fHisto->Clone());
    }
    
    fNofRuns = rhs.NofRuns();
    fNofTriggers = rhs.NofTriggers();
    fBin = rhs.Bin();
    fTriggerClass = rhs.fTriggerClass;
    fEventSelection = rhs.fEventSelection;
    fPairSelection = rhs.fPairSelection;
    fCentralitySelection = rhs.fCentralitySelection;
    fFitFunction = rhs.fFitFunction;
    fFitRejectRangeLow = rhs.fFitRejectRangeLow;
    fFitRejectRangeHigh = rhs.fFitRejectRangeHigh;
    fRejectFitPoints = rhs.fRejectFitPoints;
    fParticle = rhs.fParticle;
    fMinvRS = rhs.fMinvRS;

  }
  
  return *this;
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult::~AliAnalysisMuMuJpsiResult()
{
  // dtor
  delete fHisto;
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::AttachFunctionsToHisto(TF1* signal, TF1* bck, TF1* total,
                                                       Double_t xmin, Double_t xmax)
{
  /// Add some fit functions and some lines to the histogram
  if (signal)
  {
    signal->SetLineColor(1);
    signal->SetNpx(1000);
    fHisto->GetListOfFunctions()->Add(signal);
  }

  if ( bck )
  {
    bck->SetLineColor(2);
    bck->SetNpx(1000);
    fHisto->GetListOfFunctions()->Add(bck);
  }
  
  if ( total )
  {
    total->SetLineColor(4);
    total->SetNpx(1000);
    fHisto->GetListOfFunctions()->Add(total);
  }
  
  TLine* l1 = new TLine(xmin,0,xmin,fHisto->GetMaximum()*0.8);
  TLine* l2 = new TLine(xmax,0,xmax,fHisto->GetMaximum()*0.8);
  l1->SetLineColor(6);
  l2->SetLineColor(6);
  fHisto->GetListOfFunctions()->Add(l1);
  fHisto->GetListOfFunctions()->Add(l2);
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::AttachFunctionsToHisto(TF1* signal1, TF1* signal2, TF1* bck, TF1* total,
                                                       Double_t xmin, Double_t xmax)
{
  /// Add some fit functions and some lines to the histogram
  if (signal1)
  {
    signal1->SetLineColor(1);
    signal1->SetNpx(1000);
    fHisto->GetListOfFunctions()->Add(signal1);
  }
  
  if (signal2)
  {
    signal2->SetLineColor(3);
    signal2->SetNpx(1000);
    fHisto->GetListOfFunctions()->Add(signal2);
  }
  
  if ( bck )
  {
    bck->SetLineColor(2);
    bck->SetNpx(1000);
    fHisto->GetListOfFunctions()->Add(bck);
  }
  
  if ( total )
  {
    total->SetLineColor(4);
    total->SetNpx(1000);
    fHisto->GetListOfFunctions()->Add(total);
  }
  
  TLine* l1 = new TLine(xmin,0,xmin,fHisto->GetMaximum()*0.8);
  TLine* l2 = new TLine(xmax,0,xmax,fHisto->GetMaximum()*0.8);
  l1->SetLineColor(6);
  l2->SetLineColor(6);
  fHisto->GetListOfFunctions()->Add(l1);
  fHisto->GetListOfFunctions()->Add(l2);
}


//_____________________________________________________________________________
const AliAnalysisMuMuBinning::Range& AliAnalysisMuMuJpsiResult::Bin() const
{
  /// Get the bin of this result
  if ( !Mother() ) return fBin;
  else return Mother()->Bin();
}

//_____________________________________________________________________________
TObject* AliAnalysisMuMuJpsiResult::Clone(const char* /*newname*/) const
{
  /// Clone this result
  return new AliAnalysisMuMuJpsiResult(*this);
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuJpsiResult::Correct(const AliAnalysisMuMuJpsiResult& other, const char* particle, const char* subResultName)
{
  /// Assuming other has an AccxEff entry, correct this value by AccxEff of other
  
  if ( HasValue(Form("Nof%s",particle)) )
  {
    if (!other.HasValue(Form("AccEff%s",particle),subResultName))
    {
      AliError(Form("Cannot correct as I do not find the AccEff%s value (subResultName=%s)!",particle,subResultName));
      return kFALSE;
    }
    
    Double_t acc = other.GetValue(Form("AccEff%s",particle),subResultName);
    Double_t value = 0.0;
    
    if ( acc > 0 ) value = GetValue(Form("Nof%s",particle)) / acc;
    
    Double_t error = ErrorAB( GetValue(Form("Nof%s",particle)),
                              GetErrorStat(Form("Nof%s",particle)),
                              other.GetValue(Form("AccEff%s",particle),subResultName),
                              other.GetErrorStat(Form("AccEff%s",particle),subResultName) );
                                    
    Set(Form("CorrNof%s",particle),value,error*value);
    
    return kTRUE;
  }
  else
  {
    AliError(Form("Result does not have Nof%s : cannot correct it !",particle));
  }
  return kFALSE;
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::Draw(Option_t* opt)
{
  /// short cut method to draw our internal histogram
  if (fHisto) fHisto->Draw(opt);
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::CountParticle(const TH1& hminv, const char* particle, Double_t sigma)
{
  /// Count the number of entries in an invariant mass range
  
  const std::map<std::string,Double_t>& m = MassMap();
  
  std::map<std::string,Double_t>::const_iterator it = m.find(particle);
  
  if ( it == m.end() )
  {
    AliErrorClass(Form("Don't know the mass of particle %s",particle));
    return 0.0;
  }
  
  Double_t mass = it->second;
  
  if ( sigma < 0 )
  {
    AliDebugClass(1,Form("Oups. Got a sigma of %e for particle %s !",sigma,particle));
    return hminv.Integral();
  }
  
  const TAxis* x = hminv.GetXaxis();

  Int_t b1 = x->FindFixBin(mass-sigma);
  Int_t b2 = x->FindFixBin(mass+sigma);
  
  AliDebugClass(1,Form("hminv getentries %e integral %e",hminv.GetEntries(),hminv.Integral(b1,b2)));
  
  return hminv.Integral(b1,b2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundLin(Double_t *x, Double_t *par)
{
  // Linear function for Bkg 2 params
  
  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
  return par[0]+par[1]*x[0];
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2(Double_t *x, Double_t *par)
{
  // pol2 3 params
  
  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol3(Double_t *x, Double_t *par)
{
  // pol2 4 params
  
  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0];
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4(Double_t *x, Double_t *par)
{
  // pol2 5 params
  
  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0];
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp(Double_t *x, Double_t *par)
{
  // pol2 x exp : 4 params
  
  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
//  return par[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0])*TMath::Exp(par[4]/x[0]);
  return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0])*TMath::Exp(par[3]*x[0]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Exp(Double_t *x, Double_t *par)
{
  // pol4 x exp : 6 params
  
  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
//  return par[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]+par[5]*x[0]*x[0]*x[0]*x[0])*TMath::Exp(par[6]/x[0]);
//  return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0])*TMath::Exp(par[5]/x[0]);
  return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0])*TMath::Exp(par[5]*x[0]);

}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG(Double_t *x, Double_t *par)
{
  // gaussian variable width : 4 params
  
  if (fRejectFitPoints && x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
  Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
  return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended(Double_t *x,Double_t *par)
{
  // Extended Crystal Ball : 7 parameters
  //
  // par[0] = Normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = alpha
  // par[4] = n
  // par[5] = alpha'
  // par[6] = n'
  
  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  
  if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
    return par[0]*(exp(-0.5*t*t));
  }
  
  if (t < -absAlpha) //left tail
  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;
    return par[0]*(a/TMath::Power(b - t, par[4]));
  }
  
  if (t >= absAlpha2) //right tail
  {
    
    Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
    Double_t d = par[6]/absAlpha2 - absAlpha2;
    return par[0]*(c/TMath::Power(d + t, par[6]));
  }
  
  return 0. ;
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionNA60New(Double_t *x,Double_t *par)
{
  // New Formulation of NA60 : 11 parameters
  //
  // par[0] = Normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = p1Left
  // par[4] = p2Left
  // par[5] = p3Left
  // par[6] = p1Right
  // par[7] = p2Right
  // par[8] = p3Right
  // par[9] = alphaLeft
  // par[10] = alphaRight

  
  const Double_t t = (x[0]-par[1])/par[2];
  
  Double_t sigmaRatio;
  if( t < par[9] ) sigmaRatio = ( 1.0 + TMath::Power( par[3]*(par[9]-t), par[4]-par[5]*TMath::Sqrt(par[9] - t) ) );
  else if( t >= par[9] && t < par[10] ) sigmaRatio = 1;
  else if( t >= par[10] ) sigmaRatio = ( 1.0 + TMath::Power( par[6]*(t-par[10]), par[7]-par[8]*TMath::Sqrt(t - par[10]) ) );
  
  return par[0]*TMath::Exp( -(1/2.)*TMath::Power(t/sigmaRatio,2.));
  
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewVWG(Double_t *x, Double_t *par)
{
  /// 2 NA60 (new) + pol2 x exp
  /// width of the second NA60 related to the first (free) one.
  
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[11] = {
    par[15],
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609
    par[7],
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
  };
  return FitFunctionBackgroundVWG(x, par) + FitFunctionNA60New(x, &par[4]) + FitFunctionNA60New(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewPol2Exp(Double_t *x, Double_t *par)
{
  /// 2 NA60 (new) + pol2 x exp
  /// width of the second NA60 related to the first (free) one.
  
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[11] = {
    par[15],
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor,  // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
  };
  return FitFunctionBackgroundPol2Exp(x, par) + FitFunctionNA60New(x, &par[4]) + FitFunctionNA60New(x, par2);
}


//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewPol4Exp(Double_t *x, Double_t *par)
{
  /// 2 NA60 (new) + pol4 x exp
  /// width of the second NA60 related to the first (free) one.
  
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[11] = {
    par[17],
    3.68609+(par[7]-3.096916)/3.096916*3.68609,
    par[8]*SPsiPFactor, // /3.096916*3.68609,
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
    par[15],
    par[16],
  };
  return FitFunctionBackgroundPol4Exp(x, par) + FitFunctionNA60New(x, &par[6]) + FitFunctionNA60New(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Lin(Double_t *x, Double_t *par)
{
  /// 2 extended crystal balls + Pol1
  /// width of the second CB related to the first (free) one.
  
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[7] = {
    par[9],
    3.68609+(par[3]-3.096916)/3.096916*3.68609,
    par[4]*SPsiPFactor, // /3.096916*3.68609,
    par[5],
    par[6],
    par[7],
    par[8]
  };
  return FitFunctionBackgroundLin(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[2]) + FitFunctionSignalCrystalBallExtended(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol2Exp(Double_t *x, Double_t *par)
{
  /// 2 extended crystal balls + pol2 x exp
  /// width of the second CB related to the first (free) one.
  
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[7] = {
    par[11],
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };
  return FitFunctionBackgroundPol2Exp(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionSignalCrystalBallExtended(x, par2);
}

////____________________________________________________________________________
//Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol2Exp(Double_t *x, Double_t *par)
//{
//  /// 2 extended crystal balls + pol2 x exp
//  /// width of the second CB related to the first (free) one.
//
//  Double_t par2[7] = {
//    par[12],
//    3.68609+(par[6]-3.096916)/3.096916*3.68609,
//    par[7]/3.096916*3.68609,
//    par[8],
//    par[9],
//    par[10],
//    par[11]
//  };
//  return FitFunctionBackgroundPol2Exp(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[5]) + FitFunctionSignalCrystalBallExtended(x, par2);
//}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol4Exp(Double_t *x, Double_t *par)
{
  /// 2 extended crystal balls + pol4 x exp
  /// width of the second CB related to the first (free) one.
  
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[7] = {
    par[13],
    3.68609+(par[7]-3.096916)/3.096916*3.68609,
    par[8]*SPsiPFactor, // /3.096916*3.68609,
    par[9],
    par[10],
    par[11],
    par[12]
  };
  return FitFunctionBackgroundPol4Exp(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[6]) + FitFunctionSignalCrystalBallExtended(x, par2);
}

////____________________________________________________________________________
//Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol4Exp(Double_t *x, Double_t *par)
//{
//  /// 2 extended crystal balls + pol4 x exp
//  /// width of the second CB related to the first (free) one.
//  
//  Double_t par2[7] = {
//    par[14],
//    3.68609+(par[8]-3.096916)/3.096916*3.68609,
//    par[9]/3.096916*3.68609,
//    par[10],
//    par[11],
//    par[12],
//    par[13]
//  };
//  return FitFunctionBackgroundPol4Exp(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[7]) + FitFunctionSignalCrystalBallExtended(x, par2);
//}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG(Double_t *x, Double_t *par)
{
  /// 2 extended crystal balls + VWG
  /// width of the second CB related to the first (free) one.
  
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[7] = {
    par[11],
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };
  return FitFunctionBackgroundVWG(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionSignalCrystalBallExtended(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWGINDEPTAILS(Double_t *x, Double_t *par)
{
  /// 2 extended crystal balls + pol2 x exp
  /// The tail parameters are independent but the sPsiP and mPsiP are fixed to the one of the JPsi
  
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[7] = {
    par[11],
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[12],
    par[13],
    par[14],
    par[15]
  };
  return FitFunctionBackgroundVWG(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionSignalCrystalBallExtended(x, par2);

//  return FitFunctionBackgroundVWG(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionSignalCrystalBallExtended(x, &par[11]);
}

//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::alphaCB2VWG(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[4])/(FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionBackgroundVWG(x,par));
}

//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::alphaCB2POL2EXP(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[4])/(FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionBackgroundPol2Exp(x,par));
}

//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::alphaNA60NEWVWG(Double_t*x, Double_t* par)
{
  return FitFunctionNA60New(x, &par[4])/(FitFunctionNA60New(x, &par[4]) + FitFunctionBackgroundVWG(x,par));
}

//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::alphaNA60NEWPOL2EXP(Double_t*x, Double_t* par)
{
  return FitFunctionNA60New(x, &par[4])/(FitFunctionNA60New(x, &par[4]) + FitFunctionBackgroundPol2Exp(x,par));
}


//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2Lin(Double_t* x, Double_t* par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi
  return alphaCB2VWG(x,par)*par[12] + (1. - alphaCB2VWG(x,par))*FitFunctionBackgroundPol2(x,&par[13]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2Lin(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[11] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[11], //kPsi'
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };
  
  
  return alphaCB2VWG(x,par)*par[12] + alphaCB2VWG(x,par2)*par[15] + (1. - alphaCB2VWG(x,par) - alphaCB2VWG(x,par2))*FitFunctionBackgroundPol2(x,&par[13]);
}


//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2VWGPOL2(Double_t* x, Double_t* par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi
  return alphaCB2VWG(x,par)*par[12] + (1. - alphaCB2VWG(x,par))*FitFunctionBackgroundPol2(x,&par[13]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL2(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[11] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[11], //kPsi'
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };
  
  
  return alphaCB2VWG(x,par)*par[12] + alphaCB2VWG(x,par2)*par[16] + (1. - alphaCB2VWG(x,par) - alphaCB2VWG(x,par2))*FitFunctionBackgroundPol2(x,&par[13]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL2EXP(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[11] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[11], //kPsi'
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };
  
  
  return alphaCB2VWG(x,par)*par[12] + alphaCB2VWG(x,par2)*par[17] + (1. - alphaCB2VWG(x,par) - alphaCB2VWG(x,par2))*FitFunctionBackgroundPol2Exp(x,&par[13]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2POL2EXPPOL2(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[11] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[11], //kPsi'
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };
  
  
  return alphaCB2POL2EXP(x,par)*par[12] + alphaCB2POL2EXP(x,par2)*par[16] + (1. - alphaCB2POL2EXP(x,par) - alphaCB2POL2EXP(x,par2))*FitFunctionBackgroundPol2(x,&par[13]);

}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2POL2EXPPOL2EXP(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[11] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[11], //kPsi'
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };
  
  
  return alphaCB2POL2EXP(x,par)*par[12] + alphaCB2POL2EXP(x,par2)*par[17] + (1. - alphaCB2POL2EXP(x,par) - alphaCB2POL2EXP(x,par2))*FitFunctionBackgroundPol2Exp(x,&par[13]);
  
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2NA60NEWVWGPOL2(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[15] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[15],
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
  };

  return alphaNA60NEWVWG(x,par)*par[16] + alphaNA60NEWVWG(x,par2)*par[20] + (1. - alphaNA60NEWVWG(x,par) - alphaNA60NEWVWG(x,par2))*FitFunctionBackgroundPol2(x,&par[17]);
  
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2NA60NEWVWGPOL2EXP(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[15] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[15],
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
  };
  
  
  return alphaNA60NEWVWG(x,par)*par[16] + alphaNA60NEWVWG(x,par2)*par[21] + (1. - alphaNA60NEWVWG(x,par) - alphaNA60NEWVWG(x,par2))*FitFunctionBackgroundPol2Exp(x,&par[17]);
  
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2NA60NEWPOL2EXPPOL2(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[15] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[15],
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
  };  
  
  return alphaNA60NEWPOL2EXP(x,par)*par[16] + alphaNA60NEWPOL2EXP(x,par2)*par[20] + (1. - alphaNA60NEWPOL2EXP(x,par) - alphaNA60NEWPOL2EXP(x,par2))*FitFunctionBackgroundPol2(x,&par[17]);
  
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2NA60NEWPOL2EXPPOL2EXP(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[15] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[15],
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
  };  
  
  return alphaNA60NEWPOL2EXP(x,par)*par[16] + alphaNA60NEWPOL2EXP(x,par2)*par[21] + (1. - alphaNA60NEWPOL2EXP(x,par) - alphaNA60NEWPOL2EXP(x,par2))*FitFunctionBackgroundPol2Exp(x,&par[17]);
  
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL3(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[11] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[11], //kPsi'
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };
  
  
  return alphaCB2VWG(x,par)*par[12] + alphaCB2VWG(x,par2)*par[17] + (1. - alphaCB2VWG(x,par) - alphaCB2VWG(x,par2))*FitFunctionBackgroundPol2(x,&par[13]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL4(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[11] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[11], //kPsi'
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };
  
  
  return alphaCB2VWG(x,par)*par[12] + alphaCB2VWG(x,par2)*par[18] + (1. - alphaCB2VWG(x,par) - alphaCB2VWG(x,par2))*FitFunctionBackgroundPol4(x,&par[13]);
}

//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL2INDEPTAILS(Double_t* x, Double_t* par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP with independent tails
  Double_t SPsiPFactor = GetValue(kKeySPsiP);
  
  Double_t par2[11] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[11], //kPsi'
   3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[12],
    par[13],
    par[14],
    par[15]
  };

  return alphaCB2VWG(x,par)*par[16] + alphaCB2VWG(x,par2)*par[20] + (1. - alphaCB2VWG(x,par) - alphaCB2VWG(x,par2))*FitFunctionBackgroundPol2(x,&par[17]);
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSICOUNT()
{
  /// Simple counting of the number of j/psi ...
  Double_t n = CountParticle(*fHisto,Form("%s",GetParticle()));
  Set(Form("Nof%s",GetParticle()),n,TMath::Sqrt(n));
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSICB2()
{
  /// Fit using 1 extended crystal balls (pure signal)
  TString particleName(GetParticle());
  
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t alphaLow = GetValue(Form("al%s",particleName.Data()));
  Double_t nLow = GetValue(Form("nl%s",particleName.Data()));
  Double_t alphaUp = GetValue(Form("au%s",particleName.Data()));
  Double_t nUp = GetValue(Form("nu%s",particleName.Data()));
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  TString msg;
  
  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  
  AliDebug(1,Form("Fit with jpsi CB2 %s",msg.Data()));
  
  // Extended Crystal Ball : 7 parameters
  //
  // par[0] = Normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = alpha
  // par[4] = n
  // par[5] = alpha'
  // par[6] = n'
  
  TF1* fitTotal = new TF1("signal",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fitRangeLow,fitRangeHigh,7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  
  fitTotal->SetParNames(Form("k%s",particleName.Data()),Form("m%s",particleName.Data()),Form("s%s",particleName.Data()),Form("al%s",particleName.Data()),
  //                                 0                                1                              2                                 3
                        Form("nl%s",particleName.Data()),Form("au%s",particleName.Data()),Form("nu%s",particleName.Data()));
  //                                   4                                 5                              6
  
  fitTotal->SetParameter(0, fHisto->GetMaximum()); // norm
  
  if (particleName.Contains("JPsi"))
  {
    fitTotal->SetParameter(1, 3.1);
    fitTotal->SetParLimits(1, 3.08, 3.2);
    fitTotal->SetParameter(2, 0.08);
    fitTotal->SetParLimits(2, 0.05, 0.15);
  }
  else if (particleName.Contains("PsiP"))
  {
    fitTotal->SetParameter(1, 3.7);
    fitTotal->SetParLimits(1, 3.63, 3.72);
    fitTotal->SetParameter(2, 0.08);
    fitTotal->SetParLimits(2, 0.05, 0.15);
  }
  else AliError(Form("Could not set initial fit parameters for particle %s: The fit might not converge",particleName.Data()));
  
  SetParameter(fitTotal,3,alphaLow,0.9,0.1,10.0);
  SetParameter(fitTotal,4,nLow,5.0,0.0,10.0);
  SetParameter(fitTotal,5,alphaUp,2.0,0.1,10.0);
  SetParameter(fitTotal,6,nUp,3.0,0.0,10.0);
  
//  fitTotal->FixParameter(3, alphaLow);
//  fitTotal->FixParameter(4, nLow);
//  fitTotal->FixParameter(5, alphaUp);
//  fitTotal->FixParameter(6, nUp);
  
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,"SER","");
  
  // Check parameters...
  if (
      StrongCorrelation(fitResult,fitTotal,3,4,2) ||
      StrongCorrelation(fitResult,fitTotal,5,6,2) ||
      WrongParameter(fitTotal,3,1) ||
      WrongParameter(fitTotal,4,0)  ||
      WrongParameter(fitTotal,5,1)  ||
      WrongParameter(fitTotal,6,0)
      )
  {
    // try again...
    fitResult = fHisto->Fit(fitTotal,"SER","");
  }
  
  TF1* signal = new TF1("signal",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","SignalCB2");
  
  signal->SetParameters(fitTotal->GetParameter(0),fitTotal->GetParameter(1),fitTotal->GetParameter(2),fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6));

  Set("FitResult",static_cast<int>(fitResult)*1.0,0.0);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);
  
  Set(Form("m%s",particleName.Data()),fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set(Form("s%s",particleName.Data()),fitTotal->GetParameter(2),fitTotal->GetParError(2));
  
  Set(Form("al%s",particleName.Data()),fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set(Form("nl%s",particleName.Data()),fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set(Form("au%s",particleName.Data()),fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set(Form("nu%s",particleName.Data()),fitTotal->GetParameter(6),fitTotal->GetParError(6));
  
  AttachFunctionsToHisto(signal,0x0,fitTotal,fitRangeLow,fitRangeHigh);
  
  Double_t a = fHisto->GetXaxis()->GetXmin();
  Double_t b = fHisto->GetXaxis()->GetXmax();
  double njpsi = fitTotal->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr = fitTotal->IntegralError(a,b)/fHisto->GetBinWidth(1);
  
  Set(Form("Nof%s",particleName.Data()),njpsi,nerr);
  
  double m = GetValue(Form("m%s",particleName.Data()));
  double s = GetValue(Form("s%s",particleName.Data()));
  double al = GetValue(Form("al%s",particleName.Data()));
  double au = GetValue(Form("au%s",particleName.Data()));
  
  double njpsiCore = fitTotal->Integral(m-al*s,m+au*s)/fHisto->GetBinWidth(1);
  double nerrCore = fitTotal->IntegralError(m-al*s,m+au*s)/fHisto->GetBinWidth(1);
  
  double njpsiTailL = fitTotal->Integral(a,m-al*s)/fHisto->GetBinWidth(1);
  double nerrTailL = fitTotal->IntegralError(a,m-al*s)/fHisto->GetBinWidth(1);
  
  double njpsiTailR = fitTotal->Integral(m+au*s,b)/fHisto->GetBinWidth(1);
  double nerrTailR = fitTotal->IntegralError(m+au*s,b)/fHisto->GetBinWidth(1);
  
  Set(Form("Nof%sCore",particleName.Data()),njpsiCore,nerrCore);
  Set(Form("Nof%sTailL",particleName.Data()),njpsiTailL,nerrTailL);
  Set(Form("Nof%sTailR",particleName.Data()),njpsiTailR,nerrTailR);
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSINA60NEW()
{
  /// Fit using 1 new NA60 (pure signal)
  
  TString particleName(GetParticle());
  
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t p1Left = GetValue(Form("p1L%s",particleName.Data()));
  Double_t p2Left = GetValue(Form("p2L%s",particleName.Data()));
  Double_t p3Left = GetValue(Form("p3L%s",particleName.Data()));
  Double_t p1Right = GetValue(Form("p1R%s",particleName.Data()));
  Double_t p2Right = GetValue(Form("p2R%s",particleName.Data()));
  Double_t p3Right = GetValue(Form("p3R%s",particleName.Data()));
  
  Double_t alphaLeft = GetValue(Form("aL%s",particleName.Data()));
  Double_t alphaRight = GetValue(Form("aR%s",particleName.Data()));
  
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  TString msg;
  
  if (IsValidValue(p1Left)) msg += TString::Format("p1L=%e ",p1Left);
  if (IsValidValue(p2Left)) msg += TString::Format("p2L=%e ",p2Left);
  if (IsValidValue(p3Left)) msg += TString::Format("p3L=%e ",p3Left);
  if (IsValidValue(p1Right)) msg += TString::Format("p1R=%e ",p1Right);
  if (IsValidValue(p2Right)) msg += TString::Format("p2R=%e ",p2Right);
  if (IsValidValue(p3Right)) msg += TString::Format("p3R=%e ",p3Right);
  
  if (IsValidValue(alphaLeft)) msg += TString::Format("aL=%e ",alphaLeft);
  if (IsValidValue(alphaRight)) msg += TString::Format("aR=%e ",alphaRight);
  
  AliDebug(1,Form("Fit with jpsi NA60 new %s",msg.Data()));
  
  // New NA60 : 11 parameters
  //
  // par[0] = Normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = p1Left
  // par[4] = p2Left
  // par[5] = p3Left
  // par[6] = p1Right
  // par[7] = p2Right
  // par[8] = p3Right
  // par[9] = alphaLeft
  // par[10] = alphaRight
  
  TF1* fitTotal = new TF1("fitSignal",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fitRangeLow,fitRangeHigh,11,"AliAnalysisMuMuJpsiResult","FitFunctionNA60New");
  
  fitTotal->SetParNames(Form("k%s",particleName.Data()),Form("m%s",particleName.Data()),Form("s%s",particleName.Data()),Form("p1L%s",particleName.Data()),
        //                                 0                                1                              2                                 3
                        Form("p2L%s",particleName.Data()),Form("p3L%s",particleName.Data()),Form("p1R%s",particleName.Data()),Form("p2R%s",particleName.Data()),
        //                                   4                                 5                              6                            7
                        Form("p3R%s",particleName.Data()),Form("aL%s",particleName.Data()),Form("aR%s",particleName.Data()));
        //                                   8                                 9                              10
  
  fitTotal->SetParameter(0, fHisto->GetMaximum()); // norm
  
  if (particleName.Contains("JPsi"))
  {
    fitTotal->SetParameter(1, 3.1);
    fitTotal->SetParLimits(1, 3.08, 3.2);
    fitTotal->SetParameter(2, 0.08);
    fitTotal->SetParLimits(2, 0.05, 0.15);
  }
  else if (particleName.Contains("PsiP"))
  {
    fitTotal->SetParameter(1, 3.7);
    fitTotal->SetParLimits(1, 3.63, 3.72);
    fitTotal->SetParameter(2, 0.08);
    fitTotal->SetParLimits(2, 0.05, 0.15);
  }
  else AliError(Form("Could not set initial fit parameters for particle %s: The fit might not converge",particleName.Data()));
  
  SetParameter(fitTotal,3,p1Left,0.02,0.001,2.0);
  SetParameter(fitTotal,4,p2Left,0.4,0.2,0.6);
  SetParameter(fitTotal,5,p3Left,0.2,0.05,0.4);
  SetParameter(fitTotal,6,p1Right,0.2,0.001,0.4);
  SetParameter(fitTotal,7,p2Right,1.0,0.0,1.4);
  SetParameter(fitTotal,8,p3Right,0.02,0.005,0.4);
  
  SetParameter(fitTotal,9,alphaLeft,0.0,-1.5,1.5);
  SetParameter(fitTotal,10,alphaRight,2.3,2.0,2.5);
  
//  fitTotal->SetParameter(3, 0.2);
//  fitTotal->SetParameter(4, 2.0);
//  fitTotal->SetParameter(5, 1.);
//  fitTotal->SetParameter(6, 0.1);
//  fitTotal->SetParameter(7, 2.4);
//  fitTotal->SetParameter(8, 1.1);
//  
//  fitTotal->SetParameter(9, 1.0);
//  fitTotal->SetParameter(10, 1.0);
  
//  SetParameter(fitTotal,3,p1Left,0.1,-1.E5,1.E5);
//  SetParameter(fitTotal,4,p2Left,1.6,-1.E5,1.E5);
//  SetParameter(fitTotal,5,p3Left,0.06,-1.E5,1.E5);
//  SetParameter(fitTotal,6,p1Right,0.1,-1.E5,1.E5);
//  SetParameter(fitTotal,7,p2Right,1.5,-1.E5,1.E5);
//  SetParameter(fitTotal,8,p3Right,0.1,-1.E5,1.E5);
//  
//  SetParameter(fitTotal,9,alphaLeft,1.,-1.E8,1.E5);
//  SetParameter(fitTotal,10,alphaRight,1.0,-1.E5,1.E5);
  
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,"SER","");
  
  // Check parameters...
//  if (
//      StrongCorrelation(fitResult,fitTotal,3,4,2) ||
//      StrongCorrelation(fitResult,fitTotal,5,6,2) ||
//      WrongParameter(fitTotal,3,1) ||
//      WrongParameter(fitTotal,4,0)  ||
//      WrongParameter(fitTotal,5,1)  ||
//      WrongParameter(fitTotal,6,0)
//      )
//  {
//    // try again...
//    fitResult = fHisto->Fit(fitTotal,"SER","");
//  }
  TF1* signal = new TF1("signal",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","SignalNA60New");
  
  signal->SetParameters(fitTotal->GetParameter(0),fitTotal->GetParameter(1),fitTotal->GetParameter(2),fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitTotal->GetParameter(7),fitTotal->GetParameter(8),fitTotal->GetParameter(9));
  signal->SetParameter(10,fitTotal->GetParameter(10));
  
  Set("FitResult",static_cast<int>(fitResult)*1.0,0.0);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);
  
  Set(Form("m%s",particleName.Data()),fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set(Form("s%s",particleName.Data()),fitTotal->GetParameter(2),fitTotal->GetParError(2));
  
  Set(Form("p1L%s",particleName.Data()),fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set(Form("p2L%s",particleName.Data()),fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set(Form("p3L%s",particleName.Data()),fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set(Form("p1R%s",particleName.Data()),fitTotal->GetParameter(6),fitTotal->GetParError(6));
  Set(Form("p2R%s",particleName.Data()),fitTotal->GetParameter(7),fitTotal->GetParError(7));
  Set(Form("p3R%s",particleName.Data()),fitTotal->GetParameter(8),fitTotal->GetParError(8));
  
  Set(Form("aL%s",particleName.Data()),fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set(Form("aR%s",particleName.Data()),fitTotal->GetParameter(10),fitTotal->GetParError(10));
  
  AttachFunctionsToHisto(signal,0x0,fitTotal,fitRangeLow,fitRangeHigh);
  
  Double_t a = fHisto->GetXaxis()->GetXmin();
  Double_t b = fHisto->GetXaxis()->GetXmax();
  double njpsi = fitTotal->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr = fitTotal->IntegralError(a,b)/fHisto->GetBinWidth(1);
  
  Set(Form("Nof%s",particleName.Data()),njpsi,nerr);
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2VWG()
{
  /// Fit using 2 extended crystal balls (signal) + variable width gaussian (background)
  
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  Double_t paramSPsiP = GetValue("FSigmaPsiP");
  
  TString msg;
  
  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  
  AliDebug(1,Form("Fit with jpsi + psiprime VWG %s",msg.Data()));
  
//  std::cout << "Tails parameters: " << msg.Data() << std::endl;

  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,12,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");

  fitTotal->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi",
//                        0      1       2       3       4      5      6
                        "alJPsi","nlJPsi","auJPsi","nuJPsi");
//                         7        8        9        10
  fitTotal->SetParName(11, "kPsiP");
//                            11

  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG");

  const char* fitOption = "SER"; //We can add NO to avoid plotting
  
#if 0
  bck->SetParameter(0,fHisto->GetMaximum());
  bck->SetParameter(1,3);
  bck->SetParameter(2,10);
  bck->SetParameter(3,10);

  bck->SetParLimits(1, 0., 100.);
  bck->SetParLimits(2, 0., 100.);
  bck->SetParLimits(3, 0., 100.);
  
  SetFitRejectRange(2.7,3.5);
  
  fHisto->Fit(bck,fitOption,"");
  
  for ( Int_t i = 0; i < 4; ++i )
  {
    Double_t a,b;
    bck->GetParLimits(i,a,b);
    fitTotal->SetParameter(i,bck->GetParameter(i));
    fitTotal->SetParLimits(i,a,b);
  }
#endif
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG,1.7,6.,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG");
  
  Int_t bin = fHisto->FindBin(0.26);
  
  bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
  
  SetFitRejectRange(2.7,4.0);
  
  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,"SR");
  
  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  
  if ( static_cast<int>(fitResultInit) )
  {
    bin = fHisto->FindBin(0.82);
    bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
    fitResultInit = fHisto->Fit(bckInit,"SR");
  }
  else if ( bckInit->GetParameter(0) < 0 )
  {
    bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
  }
  
  SetFitRejectRange();
  
  for ( Int_t i = 0; i < 4; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }
  
  

//  fitTotal->SetParameter(0, fHisto->GetMaximum()); // kVWG
//  fitTotal->SetParameter(1, 1.9); // mVWG
//  
//  fitTotal->SetParameter(2, 0.5); // sVWG1
//  fitTotal->SetParLimits(2, 0., 100.);
//  
//  fitTotal->SetParameter(3, 0.3); // sVWG2
//  fitTotal->SetParLimits(3, 0., 100.);
  
  bin = fHisto->FindBin(3.09);
  fitTotal->SetParameter(4, fHisto->GetBinContent(bin)); // norm
  
  fitTotal->SetParameter(5, 3.1); // mean
  fitTotal->SetParLimits(5, 3.0, 3.2);
  
  fitTotal->SetParameter(6, 0.08); // sigma
  fitTotal->SetParLimits(6, 0.05, 0.09);
  
  if ( IsValidValue(alphaLow) )
  {
    fitTotal->FixParameter(7, alphaLow);
  }
  else
  {
    fitTotal->SetParameter(7,0.9);
    fitTotal->SetParLimits(7,0.1,10.0);
  }
  
  if ( IsValidValue(nLow) ) 
  {
    fitTotal->FixParameter(8, nLow);
  }
  else
  {
    fitTotal->SetParameter(8,5.0);
    fitTotal->SetParLimits(8,0.0,10.0);
  }
  
  if ( IsValidValue(alphaUp) )
  {
    fitTotal->FixParameter(9, alphaUp);
  }
  else
  {
    fitTotal->SetParameter(9, 2.0);
    fitTotal->SetParLimits(9,0.1,10.0);
  }
  
  if ( IsValidValue(nUp) )
  {
    fitTotal->FixParameter(10, nUp);    
  }
  else
  {
    fitTotal->SetParameter(10,3.0);
    fitTotal->SetParLimits(10,0.0,10.0);
  }
  
  fitTotal->SetParameter(11, 10.); //kPsi'
  fitTotal->SetParLimits(11, 0.,100000.);
  
  SetFitRejectRange();
  
//  std::cout << fitTotal->GetParameter(0) << std::endl; //Just a xcheck

  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
  
//  std::cout << fitTotal->GetParameter(0) << " ?= " << fitResult->Parameter(0) << std::endl; //Just a xcheck
  
  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
   std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
  
//  std::cout << fitTotal->GetParameter(0) << std::endl; //Just a xcheck
  
  if ( static_cast<int>(fitResult) )
  {
    if ( 0.5*fitTotal->GetParameter(11) <= fitTotal->GetParError(11) ) //kPsi'
    {
      std::cout << "//-------Refitting again (setting Psi'norm=0)" << std::endl;
       bin = fHisto->FindBin(3.68);
      fitTotal->SetParameter(11, fHisto->GetBinContent(bin)/2.);
    }
    
    if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) )
    {
      std::cout << "//-------Refitting again (setting kVWG=kVWG/2)" << std::endl;
      
      fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
    }
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  
  if ( static_cast<int>(fitResult) )
  {
    std::cout << "//-------Refitting again (setting kVWG=kVWG*2)" << std::endl;
    
    fitTotal->SetParameter(0, fHisto->GetMaximum()*2.); // kVWG
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  
  if ( static_cast<int>(fitResult) )
  {
    std::cout << "//-------Refitting again (setting kVWG=kVWG/2)" << std::endl;
    
    fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    std::cout << "//-------Refitting again (setting Bkg norm= Bkg norm/2)" << std::endl;
//    fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//  }
  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
//    fitTotal->SetParameter(11, 5.); // //kPsi'
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//    
//    std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//  }
//  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    std::cout << "//-------Refitting again (refitting Bkg)" << std::endl;
//    
//    for ( Int_t i = 0; i < 4; ++i )
//    {
//      bck->SetParameter(i, fitTotal->GetParameter(i));
//    }
//    
//    SetFitRejectRange(2.7,3.5);
//    
//    fHisto->Fit(bck,"R");
//    
//    SetFitRejectRange();
//    
//    for ( Int_t i = 0; i < 4; ++i )
//    {
//      fitTotal->SetParameter(i, bck->GetParameter(i));
//    }
//    
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//    std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
//  }
//  if ( static_cast<int>(fitResult) )
//  {
//    fitTotal->SetParameter(0, fHisto->GetMaximum()); // kVWG
//    fitTotal->SetParameter(1, 1.9); // mVWG
//    
//    fitTotal->SetParameter(2, 0.5); // sVWG1
//    fitTotal->SetParLimits(2, 0., 100.);
//    
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//    std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
//
//  }
  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitTotal->SetParameter(i, bckInit->GetParameter(i));
    }
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    
    if ( static_cast<int>(fitResult) ) std::cout << "//-------Cannot fit properly, try something else..." << std::endl;
  }
  
  delete bckInit;
  
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);

  Set("mJPsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set("sJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  
  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalJPsi->SetParameters(fitTotal->GetParameter(4),
                     fitTotal->GetParameter(5),
                     fitTotal->GetParameter(6),
                     fitTotal->GetParameter(7),
                     fitTotal->GetParameter(8),
                     fitTotal->GetParameter(9),
                     fitTotal->GetParameter(10));
  
  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalPsiP->SetParameters(fitTotal->GetParameter(11),
                        3.68609+(fitTotal->GetParameter(5)-3.096916)/3.096916*3.68609,
                        fitTotal->GetParameter(6)*paramSPsiP, // /3.096916*3.68609,
                        fitTotal->GetParameter(7),
                        fitTotal->GetParameter(8),
                        fitTotal->GetParameter(9),
                        fitTotal->GetParameter(10));
  
  bck->SetParameters(fitTotal->GetParameter(0),
                     fitTotal->GetParameter(1),
                     fitTotal->GetParameter(2),
                     fitTotal->GetParameter(3));
  
  
  Set("kVWG",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("mVWG",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("sVWG1",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("sVWG2",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("kJPsi",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  
  Set("alJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
  Set("nlJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("auJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set("nuJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  
//  Set("alPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
//  Set("nlPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
//  Set("auPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
//  Set("nuPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("kPsiP",fitTotal->GetParameter(11),fitTotal->GetParError(11));


  SetFitRejectRange();
  
  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);
  
  
  Double_t cbParameters[7];
  Double_t covarianceMatrix[7][7];
  
  for ( int ix = 0; ix < 7; ++ix )
  {
    cbParameters[ix] = fitTotal->GetParameter(ix+4);
  }
  
  for ( int iy = 0; iy < 7; ++iy )
  {
    for ( int ix = 0; ix < 7; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+4,iy+4);
    }
  }
  
  Double_t a = fHisto->GetXaxis()->GetXmin();
  Double_t b = fHisto->GetXaxis()->GetXmax();
  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr = signalJPsi->IntegralError(a,b,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);

  Set("NofJPsi",njpsi,nerr);
  
  double m = GetValue("mJPsi");
  double s = GetValue("sJPsi");
  double njpsi3s = signalJPsi->Integral(m-3*s,m+3*s)/fHisto->GetBinWidth(1);
  double nerr3s = signalJPsi->IntegralError(m-3*s,m+3*s,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi3s",njpsi3s,nerr3s);
  
  //Computation of bin significance and signal over background
  
  Double_t bkgParameters[4];
  Double_t bkgcovarianceMatrix[4][4];
  
  for ( int ix = 0; ix < 4; ++ix )
  {
    bkgParameters[ix] = fitTotal->GetParameter(ix);
  }
  
  for ( int iy = 0; iy < 4; ++iy )
  {
    for ( int ix = 0; ix < 4; ++ix )
    {
      bkgcovarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix,iy);
    }
  }
  
  double nbck3s = bck->Integral(m-3.*s,m+3.*s)/fHisto->GetBinWidth(1);
  double nbck3sErr = bck->IntegralError(m-3.*s,m+3.*s,&bkgParameters[0],&bkgcovarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  double sOverB3s = njpsi3s / nbck3s;
  double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));
  
  Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);
  
  double sig = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
  double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
                               TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );
  
  Set("Significance3s",sig,sigErr);

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2VWGINDEPTAILS()
{
  /// Fit using 2 extended crystal balls with independent tails (signal) + variable width gaussian (background)
  
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  
  Double_t alphaLowP = GetValue("alPsiP");
  Double_t nLowP = GetValue("nlPsiP");
  Double_t alphaUpP = GetValue("auPsiP");
  Double_t nUpP = GetValue("nuPsiP");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  Double_t paramSPsiP = GetValue("FSigmaPsiP");
  
  TString msg;
  
  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  
  if (IsValidValue(alphaLowP)) msg += TString::Format("alphaLowP=%e ",alphaLowP);
  if (IsValidValue(nLowP)) msg += TString::Format("nLowP=%e ",nLowP);
  if (IsValidValue(alphaUpP)) msg += TString::Format("alphaUpP=%e ",alphaUpP);
  if (IsValidValue(nUpP)) msg += TString::Format("nUpP=%e ",nUpP);
  
  AliDebug(1,Form("Fit with jpsi + psiprime VWG %s",msg.Data()));
  
  // Add fit with indep tails
  
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWGINDEPTAILS,fitRangeLow,fitRangeHigh,16,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWGINDEPTAILS");
  
  fitTotal->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi",
  //                        0      1       2       3       4      5      6
                        "alJPsi","nlJPsi","auJPsi","nuJPsi");
  //                       7        8        9        10
//  fitTotal->SetParName(11, "kPsiP");
//  //                            11
//  fitTotal->SetParName(12, "mPsiP");
//  //                            12
//  fitTotal->SetParName(13, "sPsiP");
//  //                            13
//  fitTotal->SetParName(14, "alPsiP");
//  //                            14
//  fitTotal->SetParName(15, "nlPsiP");
//  //                            15
//  fitTotal->SetParName(16, "auPsiP");
//  //                            16
//  fitTotal->SetParName(17, "nuPsiP");
//  //                            17
  
  fitTotal->SetParName(11, "kPsiP");
  //                            11
   fitTotal->SetParName(12, "alPsiP");
  //                            12
  fitTotal->SetParName(13, "nlPsiP");
  //                            13
  fitTotal->SetParName(14, "auPsiP");
  //                            14
  fitTotal->SetParName(15, "nuPsiP");
  //                            15
  
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG");
  
  const char* fitOption = "SER";
  
#if 0
  bck->SetParameter(0,fHisto->GetMaximum());
  bck->SetParameter(1,3);
  bck->SetParameter(2,10);
  bck->SetParameter(3,10);
  
  bck->SetParLimits(1, 0., 100.);
  bck->SetParLimits(2, 0., 100.);
  bck->SetParLimits(3, 0., 100.);
  
  SetFitRejectRange(2.7,3.5);
  
  fHisto->Fit(bck,fitOption,"");
  
  for ( Int_t i = 0; i < 4; ++i )
  {
    Double_t a,b;
    bck->GetParLimits(i,a,b);
    fitTotal->SetParameter(i,bck->GetParameter(i));
    fitTotal->SetParLimits(i,a,b);
  }
#endif
  
  Int_t bin = fHisto->FindBin(fitRangeLow);
  fitTotal->SetParameter(0, fHisto->GetBinContent(bin)); // kVWG
  fitTotal->SetParameter(1, 1.9); // mVWG
  
  fitTotal->SetParameter(2, 0.5); // sVWG1
  fitTotal->SetParLimits(2, 0., 100.);
  
  fitTotal->SetParameter(3, 0.3); // sVWG2
  fitTotal->SetParLimits(3, 0., 100.);
  
  bin = fHisto->FindBin(3.09);
  fitTotal->SetParameter(4, fHisto->GetBinContent(bin)); // norm(kJPsi)
  
  fitTotal->SetParameter(5, 3.1); // mean
  fitTotal->SetParLimits(5, 3.0, 3.2);
  
  fitTotal->SetParameter(6, 0.08); // sigma
  fitTotal->SetParLimits(6, 0.03, 0.15);
  
  if ( IsValidValue(alphaLow) )
  {
    fitTotal->FixParameter(7, alphaLow);
  }
  else
  {
    fitTotal->SetParameter(7,0.9);
    fitTotal->SetParLimits(7,0.1,10.0);
  }
  
  if ( IsValidValue(nLow) )
  {
    fitTotal->FixParameter(8, nLow);
  }
  else
  {
    fitTotal->SetParameter(8,5.0);
    fitTotal->SetParLimits(8,0.0,10.0);
  }
  
  if ( IsValidValue(alphaUp) )
  {
    fitTotal->FixParameter(9, alphaUp);
  }
  else
  {
    fitTotal->SetParameter(9, 2.0);
    fitTotal->SetParLimits(9,0.1,10.0);
  }
  
  if ( IsValidValue(nUp) )
  {
    fitTotal->FixParameter(10, nUp);
  }
  else
  {
    fitTotal->SetParameter(10,3.0);
    fitTotal->SetParLimits(10,0.0,10.0);
  }
  
  bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(11, fHisto->GetBinContent(bin)); //kPsiP
  fitTotal->SetParLimits(11, 0., fHisto->GetBinContent(bin)*2); //kPsiP
  
//  fitTotal->SetParameter(12, 3.7); // mean PsiP
//  fitTotal->SetParLimits(12, 3.6, 3.71);
//  
//  fitTotal->SetParameter(13, 0.08/3.096916*3.68609); // sigma PsiP
//  fitTotal->SetParLimits(13, 0.03/3.096916*3.68609, 0.15/3.096916*3.68609);
  
  if ( IsValidValue(alphaLowP) )
  {
    fitTotal->FixParameter(12, alphaLowP);
  }
  else
  {
    fitTotal->SetParameter(12,0.9);
    fitTotal->SetParLimits(12,0.1,10.0);
  }
  
  if ( IsValidValue(nLowP) )
  {
    fitTotal->FixParameter(13, nLowP);
  }
  else
  {
    fitTotal->SetParameter(13,5.0);
    fitTotal->SetParLimits(13,0.0,10.0);
  }
  
  if ( IsValidValue(alphaUpP) )
  {
    fitTotal->FixParameter(14, alphaUpP);
  }
  else
  {
    fitTotal->SetParameter(14, 2.0);
    fitTotal->SetParLimits(14,0.1,10.0);
  }
  
  if ( IsValidValue(nUpP) )
  {
    fitTotal->FixParameter(15, nUpP);
  }
  else
  {
    fitTotal->SetParameter(15,3.0);
    fitTotal->SetParLimits(15,0.0,10.0);
  }
  
  //  SetFitRejectRange();
  
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
  
  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    if ( fitTotal->GetParameter(11) <= fitTotal->GetParError(11) ) //kPsi'
    {
      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
      
      fitTotal->SetParameter(11, fHisto->GetBinContent(bin)/2.);
    }
    
    if ( fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) //kVWG
    {
      std::cout << "//-------Refitting again (setting VWG norm= VWG norm /2)" << std::endl;
      bin = fHisto->FindBin(fitRangeLow);
      
      fitTotal->SetParameter(0, fHisto->GetBinContent(bin)/2.);
    }
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  
  if ( static_cast<int>(fitResult) )
  {
    if ( fitTotal->GetParameter(11) <= fitTotal->GetParError(11) ) //kPsi'
    {
      std::cout << "//-------Refitting again (setting Psi'norm=0)" << std::endl;
      
      fitTotal->FixParameter(11, 0.);
    }
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  
  if ( static_cast<int>(fitResult) )
  {
    std::cout << "//-------Refitting again (setting kVWG=kVWG/2)" << std::endl;
    
    fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  
  if ( static_cast<int>(fitResult) )
  {
    std::cout << "//-------Cannot fit properly, try something else..." << std::endl;
  }
  
  
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);
  
  Set("mJPsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set("sJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  
  Set("mPsiP",3.68609+(fitTotal->GetParameter(5)-3.096916)/3.096916*3.68609,fitTotal->GetParError(5)/3.096916*3.68609);
  Set("sPsiP",fitTotal->GetParameter(6)*paramSPsiP,fitTotal->GetParError(6)*paramSPsiP);
  
  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalJPsi->SetParameters(fitTotal->GetParameter(4),
                        fitTotal->GetParameter(5),
                        fitTotal->GetParameter(6),
                        fitTotal->GetParameter(7),
                        fitTotal->GetParameter(8),
                        fitTotal->GetParameter(9),
                        fitTotal->GetParameter(10));
  
  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalPsiP->SetParameters(fitTotal->GetParameter(11),
                            3.68609+(fitTotal->GetParameter(5)-3.096916)/3.096916*3.68609,
                            fitTotal->GetParameter(6)*paramSPsiP, // /3.096916*3.68609,
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13),
                            fitTotal->GetParameter(14),
                            fitTotal->GetParameter(15));
  
  bck->SetParameters(fitTotal->GetParameter(0),
                     fitTotal->GetParameter(1),
                     fitTotal->GetParameter(2),
                     fitTotal->GetParameter(3));
  
  
  Set("kVWG",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("mVWG",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("sVWG1",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("sVWG2",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  
  Set("kJPsi",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set("alJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
  Set("nlJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("auJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set("nuJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  
  Set("kPsiP",fitTotal->GetParameter(11),fitTotal->GetParError(11));
  Set("alPsiP",fitTotal->GetParameter(12),fitTotal->GetParError(12));
  Set("nlPsiP",fitTotal->GetParameter(13),fitTotal->GetParError(13));
  Set("auPsiP",fitTotal->GetParameter(14),fitTotal->GetParError(14));
  Set("nuPsiP",fitTotal->GetParameter(15),fitTotal->GetParError(15));
 
  
  
//  SetFitRejectRange();
  
  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);
  
  
  Double_t cbParameters[7];
  Double_t covarianceMatrix[7][7];
  
  for ( int ix = 0; ix < 7; ++ix )
  {
    cbParameters[ix] = fitTotal->GetParameter(ix+4);
  }
  
  for ( int iy = 0; iy < 7; ++iy )
  {
    for ( int ix = 0; ix < 7; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+4,iy+4);
    }
  }
  
  Double_t a = fHisto->GetXaxis()->GetXmin();
  Double_t b = fHisto->GetXaxis()->GetXmax();
  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr = signalJPsi->IntegralError(a,b,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi",njpsi,nerr);
  
  double m = GetValue("mJPsi");
  double s = GetValue("sJPsi");
  double njpsi3s = signalJPsi->Integral(m-3*s,m+3*s)/fHisto->GetBinWidth(1);
  double nerr3s = signalJPsi->IntegralError(m-3*s,m+3*s,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi3s",njpsi3s,nerr3s);
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2POL2EXP()
{
  /// Fit using 2 extended crystal balls (signal) + pol2 x exp (background)
  /// 13 parameters
  
  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  Double_t paramSPsiP = GetValue("FSigmaPsiP");
  
  TString msg;
  
  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  
  AliDebug(1,Form("Fit with jpsi + psiprime Pol2 x Exp %s",msg.Data()));
  
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol2Exp,fitRangeLow,fitRangeHigh,12,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2Pol2Exp");
  
  fitTotal->SetParNames("pol0","pol1","pol2","exp","kJPsi","mJPsi","sJPsi","alJPsi",
  //                      0      1       2     3      4      5       6       7
                        "nlJPsi","auJPsi","nuJPsi");
  //                          8       9        10
  fitTotal->SetParName(11,"kPsiP");
  //                            11
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  
  
  //___
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,1.7,6.,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  
  Int_t bin = fHisto->FindBin(0.26);
  
  bckInit->SetParameters(fHisto->GetBinContent(bin),-fHisto->GetBinContent(bin)/3.,100.,0.05);//fHisto->GetBinContent(bin)
  
  SetFitRejectRange(2.7,4.0);
  
  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,"SR");
  
  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  
  if ( static_cast<int>(fitResultInit) )
  {
    bin = fHisto->FindBin(0.82);
    bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
    fitResultInit = fHisto->Fit(bckInit,"SR");
  }
  
  SetFitRejectRange();
  
  for ( Int_t i = 0; i < 4; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }
  
  delete bckInit;
  //__
  //
  //  bck->SetParameters(fHisto->GetMaximum(),0.05,0.05,0.05,1.);
  //
  //  SetFitRejectRange(2.7,3.5);
  //
  //  fHisto->Fit(bck,"R");
  //
  //  SetFitRejectRange();
  //
  //  for ( Int_t i = 0; i < 5; ++i )
  //  {
  //    fitTotal->SetParameter(i, bck->GetParameter(i));
  //  }
  
  bin = fHisto->FindBin(3.09);
  fitTotal->SetParameter(4, fHisto->GetBinContent(bin)); // norm(kJPsi)
  
  fitTotal->SetParameter(5, 3.1);
  fitTotal->SetParLimits(5, 3.07, 3.2);
  fitTotal->SetParameter(6, 0.08);
  fitTotal->SetParLimits(6, 0.05, 0.15);
  
  if ( IsValidValue(alphaLow) )
  {
    fitTotal->FixParameter(7, alphaLow);
  }
  else
  {
    fitTotal->SetParameter(7,0.9);
    fitTotal->SetParLimits(7,0.1,10.0);
  }
  
  if ( IsValidValue(nLow) )
  {
    fitTotal->FixParameter(8, nLow);
  }
  else
  {
    fitTotal->SetParameter(8,5.0);
    fitTotal->SetParLimits(8,0.0,10.0);
  }
  
  if ( IsValidValue(alphaUp) )
  {
    fitTotal->FixParameter(9, alphaUp);
  }
  else
  {
    fitTotal->SetParameter(9, 2.0);
    fitTotal->SetParLimits(9,0.1,10.0);
  }
  
  if ( IsValidValue(nUp) )
  {
    fitTotal->FixParameter(10, nUp);
  }
  else
  {
    fitTotal->SetParameter(10,3.0);
    fitTotal->SetParLimits(10,0.0,10.0);
  }
  
  fitTotal->SetParameter(11, 10.);
  
  const char* fitOption = "SER";
  
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
  
  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    if ( 0.5*fitTotal->GetParameter(11) <= fitTotal->GetParError(11) ) //kPsi'
    {
      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
      
      bin = fHisto->FindBin(3.68);
      fitTotal->SetParameter(11, fHisto->GetBinContent(bin)/2.);
    }
    
    if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) //kPol2Exp
    {
      std::cout << "//-------Refitting again (setting kPol2Exp norm= kPol2Exp norm /2)" << std::endl;
      bin = fHisto->FindBin(fitRangeLow);
      
      fitTotal->SetParameter(0, fHisto->GetBinContent(bin)/2.);
    }
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  
  //  if ( fitResult->CovMatrixStatus() != 3 )
  //  {
  //    std::cout << "//-------Refitting again (setting Bkg norm= Bkg norm/2)" << std::endl;
  //    fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
  //    fitResult = fHisto->Fit(fitTotal,fitOption,"");
  //  }
  
  if ( fitResult->CovMatrixStatus() != 3 )
  {
    if ( 0.5*fitTotal->GetParameter(11) <= fitTotal->GetParError(11) ) //kPsi'
    {
      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
      
      bin = fHisto->FindBin(3.68);
      fitTotal->FixParameter(11, 0.);
    }
    
    else
    {
      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
      
      fHisto->FindBin(3.68);
      fitTotal->SetParameter(11, fHisto->GetBinContent(bin)/2.); // kVWG
    }
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    
    std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
  }
  
  //  if ( fitResult->CovMatrixStatus() != 3 )
  //  {
  //    std::cout << "//-------Refitting again (refitting Bkg)" << std::endl;
  //
  //    if( fitTotal->GetParameter(9) >= 0.09 || fitTotal->GetParameter(9) <= 0.06 || fitTotal->GetParameter(6) >= 3.11 || fitTotal->GetParameter(6) <= 3.08)
  //    {
  //      bck->SetParameters(0.01,fHisto->GetMaximum(),0.01,0.05,0.05,0.05,1.);
  //      //      for ( Int_t i = 0; i < 7; ++i )
  //      //      {
  //      //        bck->SetParameter(i, 0.);
  //      //      }
  //    }
  //    else
  //    {
  //      for ( Int_t i = 0; i < 7; ++i )
  //      {
  //        bck->SetParameter(i, fitTotal->GetParameter(i));
  //      }
  //    }
  
  //
  //    for ( Int_t i = 0; i < 5; ++i )
  //    {
  //      bck->SetParameter(i, fitTotal->GetParameter(i));
  //    }
  //
  //    SetFitRejectRange(2.7,3.5);
  //
  //    fHisto->Fit(bck,"R");
  //
  //    SetFitRejectRange();
  //
  //    for ( Int_t i = 0; i < 5; ++i )
  //    {
  //      fitTotal->SetParameter(i, bck->GetParameter(i));
  //    }
  //
  //    fitResult = fHisto->Fit(fitTotal,fitOption,"");
  //  }
  //  std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
  
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);
  
  Set("mJPsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set("sJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  
  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalJPsi->SetParameters(fitTotal->GetParameter(4),
                            fitTotal->GetParameter(5),
                            fitTotal->GetParameter(6),
                            fitTotal->GetParameter(7),
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10));
  
  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalPsiP->SetParameters(fitTotal->GetParameter(11),
                            3.68609+(fitTotal->GetParameter(5)-3.096916)/3.096916*3.68609,
                            fitTotal->GetParameter(6)*paramSPsiP, // /3.096916*3.68609,
                            fitTotal->GetParameter(7),
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10));
  
  for ( Int_t i = 0; i < 4; ++i )
  {
    bck->SetParameter(i, fitTotal->GetParameter(i));
  }
  
//  Set("kPol2Exp",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("pol0",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("pol1",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("pol2",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("exp",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("kJPsi",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  
  Set("alJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
  Set("nlJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("auJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set("nuJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  
  Set("kPsiP",fitTotal->GetParameter(11),fitTotal->GetParError(11));
  
  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);
  
  Double_t cbParameters[7];
  Double_t covarianceMatrix[7][7];
  
  for ( int ix = 0; ix < 7; ++ix )
  {
    cbParameters[ix] = fitTotal->GetParameter(ix+4);
  }
  
  for ( int iy = 0; iy < 7; ++iy )
  {
    for ( int ix = 0; ix < 7; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+4,iy+4);
    }
  }
  
  Double_t a = fHisto->GetXaxis()->GetXmin();
  Double_t b = fHisto->GetXaxis()->GetXmax();
  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr = signalJPsi->IntegralError(a,b,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi",njpsi,nerr);
  
  double m = GetValue("mJPsi");
  double s = GetValue("sJPsi");
  double njpsi3s = signalJPsi->Integral(m-3*s,m+3*s)/fHisto->GetBinWidth(1);
  double nerr3s = signalJPsi->IntegralError(m-3*s,m+3*s,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi3s",njpsi3s,nerr3s);
  
  //Computation of bin significance and signal over background
  
  Double_t bkgParameters[4];
  Double_t bkgcovarianceMatrix[4][4];
  
  for ( int ix = 0; ix < 4; ++ix )
  {
    bkgParameters[ix] = fitTotal->GetParameter(ix);
  }
  
  for ( int iy = 0; iy < 4; ++iy )
  {
    for ( int ix = 0; ix < 4; ++ix )
    {
      bkgcovarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix,iy);
    }
  }
  
  double nbck3s = bck->Integral(m-3.*s,m+3.*s)/fHisto->GetBinWidth(1);
  double nbck3sErr = bck->IntegralError(m-3.*s,m+3.*s,&bkgParameters[0],&bkgcovarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  double sOverB3s = njpsi3s / nbck3s;
  double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));
  
  Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);
  
  double sig = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
  double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
                              TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );
  
  Set("Significance3s",sig,sigErr);
}

////_____________________________________________________________________________
//void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2POL2EXP()
//{
//  /// Fit using 2 extended crystal balls (signal) + pol2 x exp (background)
//  /// 13 parameters
//  
//  Double_t alphaLow = GetValue("alJPsi");
//  Double_t nLow = GetValue("nlJPsi");
//  Double_t alphaUp = GetValue("auJPsi");
//  Double_t nUp = GetValue("nuJPsi");
//  Double_t fitRangeLow = GetValue(kFitRangeLow);
//  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
//  
//  TString msg;
//  
//  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
//  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
//  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
//  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
//  
//  AliDebug(1,Form("Fit with jpsi + psiprime Pol2 x Exp %s",msg.Data()));
//  
//  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol2Exp,fitRangeLow,fitRangeHigh,13,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2Pol2Exp");
//  
//  fitTotal->SetParNames("kPol2Exp","pol0","pol1","pol2","exp","kJPsi","mJPsi","sJPsi",
////                           0      1       2      3     4      5       6       7 
//                        "alJPsi","nlJPsi","auJPsi");
////                          8       9        10
//  fitTotal->SetParName(11,"nuJPsi");
////                            11
//  fitTotal->SetParName(12, "kPsiP");
////                            12
//  
//  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
//  
//  
//  //___
//  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,1.7,6.,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
//  
//  Int_t bin = fHisto->FindBin(0.26);
//  
//  bckInit->SetParameters(-0.1,1.,0.5,0.3,1.);//fHisto->GetBinContent(bin)
//  
//  SetFitRejectRange(2.7,4.0);
//  
//  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,"SR");
//  
//  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
//  
//  if ( static_cast<int>(fitResultInit) )
//  {
//    bin = fHisto->FindBin(0.82);
//    bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
//    fitResultInit = fHisto->Fit(bckInit,"SR");
//  }
//  
//  SetFitRejectRange();
//  
//  for ( Int_t i = 0; i < 5; ++i )
//  {
//    fitTotal->SetParameter(i, bckInit->GetParameter(i));
//  }
//  
//  delete bckInit;
////__
////  
////  bck->SetParameters(fHisto->GetMaximum(),0.05,0.05,0.05,1.);
////  
////  SetFitRejectRange(2.7,3.5);
////  
////  fHisto->Fit(bck,"R");
////
////  SetFitRejectRange();
////
////  for ( Int_t i = 0; i < 5; ++i )
////  {
////    fitTotal->SetParameter(i, bck->GetParameter(i));
////  }
//  
//  bin = fHisto->FindBin(3.09);
//  fitTotal->SetParameter(5, fHisto->GetBinContent(bin)); // norm(kJPsi)
//  
//  fitTotal->SetParameter(6, 3.1);
//  fitTotal->SetParLimits(6, 3.07, 3.2);
//  fitTotal->SetParameter(7, 0.08);
//  fitTotal->SetParLimits(7, 0.05, 0.15);
//  
//  if ( IsValidValue(alphaLow) )
//  {
//    fitTotal->FixParameter(8, alphaLow);
//  }
//  else
//  {
//    fitTotal->SetParameter(8,0.9);
//    fitTotal->SetParLimits(8,0.1,10.0);
//  }
//  
//  if ( IsValidValue(nLow) )
//  {
//    fitTotal->FixParameter(9, nLow);
//  }
//  else
//  {
//    fitTotal->SetParameter(9,5.0);
//    fitTotal->SetParLimits(9,0.0,10.0);
//  }
//  
//  if ( IsValidValue(alphaUp) )
//  {
//    fitTotal->FixParameter(10, alphaUp);
//  }
//  else
//  {
//    fitTotal->SetParameter(10, 2.0);
//    fitTotal->SetParLimits(10,0.1,10.0);
//  }
//  
//  if ( IsValidValue(nUp) )
//  {
//    fitTotal->FixParameter(11, nUp);
//  }
//  else
//  {
//    fitTotal->SetParameter(11,3.0);
//    fitTotal->SetParLimits(11,0.0,10.0);
//  }
//  
//  fitTotal->SetParameter(12, 10.);
//  
//  const char* fitOption = "SER";
//  
//  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
//  
//  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
//  std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//  
//  if ( static_cast<int>(fitResult) )
//  {
//    if ( fitTotal->GetParameter(12) <= fitTotal->GetParError(12) ) //kPsi'
//    {
//      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
//      
//      bin = fHisto->FindBin(3.68);
//      fitTotal->SetParameter(12, fHisto->GetBinContent(bin)/2.);
//    }
//    
//    if ( fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) //kPol2Exp
//    {
//      std::cout << "//-------Refitting again (setting kPol2Exp norm= kPol2Exp norm /2)" << std::endl;
//      bin = fHisto->FindBin(fitRangeLow);
//      
//      fitTotal->SetParameter(0, fHisto->GetBinContent(bin)/2.);
//    }
//    
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
//  }
//  
////  if ( fitResult->CovMatrixStatus() != 3 )
////  {
////    std::cout << "//-------Refitting again (setting Bkg norm= Bkg norm/2)" << std::endl;
////    fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
////    fitResult = fHisto->Fit(fitTotal,fitOption,"");
////  }
//  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    if ( fitTotal->GetParameter(12) <= fitTotal->GetParError(12) ) //kPsi'
//    {
//      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
//      
//      bin = fHisto->FindBin(3.68);
//      fitTotal->FixParameter(12, 0.);
//    }
//
//    else
//    {
//      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
//      
//      fHisto->FindBin(3.68);
//      fitTotal->SetParameter(12, fHisto->GetBinContent(bin)/2.); // kVWG
//    }
//    
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//    
//     std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//  }
//  
////  if ( fitResult->CovMatrixStatus() != 3 )
////  {
////    std::cout << "//-------Refitting again (refitting Bkg)" << std::endl;
////    
////    if( fitTotal->GetParameter(9) >= 0.09 || fitTotal->GetParameter(9) <= 0.06 || fitTotal->GetParameter(6) >= 3.11 || fitTotal->GetParameter(6) <= 3.08)
////    {
////      bck->SetParameters(0.01,fHisto->GetMaximum(),0.01,0.05,0.05,0.05,1.);
////      //      for ( Int_t i = 0; i < 7; ++i )
////      //      {
////      //        bck->SetParameter(i, 0.);
////      //      }
////    }
////    else
////    {
////      for ( Int_t i = 0; i < 7; ++i )
////      {
////        bck->SetParameter(i, fitTotal->GetParameter(i));
////      }
////    }
//
////    
////    for ( Int_t i = 0; i < 5; ++i )
////    {
////      bck->SetParameter(i, fitTotal->GetParameter(i));
////    }
////  
////    SetFitRejectRange(2.7,3.5);
////    
////    fHisto->Fit(bck,"R");
////    
////    SetFitRejectRange();
////    
////    for ( Int_t i = 0; i < 5; ++i )
////    {
////      fitTotal->SetParameter(i, bck->GetParameter(i));
////    }
////    
////    fitResult = fHisto->Fit(fitTotal,fitOption,"");
////  }
////  std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//  
//  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
//  Set("FitNDF",fitTotal->GetNDF(),0.0);
//  
//  Set("mJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
//  Set("sJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
//  
//  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
//  signalJPsi->SetParameters(fitTotal->GetParameter(5),
//                            fitTotal->GetParameter(6),
//                            fitTotal->GetParameter(7),
//                            fitTotal->GetParameter(8),
//                            fitTotal->GetParameter(9),
//                            fitTotal->GetParameter(10),
//                            fitTotal->GetParameter(11));
//  
//  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
//  signalPsiP->SetParameters(fitTotal->GetParameter(12),
//                            3.68609+(fitTotal->GetParameter(6)-3.096916)/3.096916*3.68609,
//                            fitTotal->GetParameter(7)/3.096916*3.68609,
//                            fitTotal->GetParameter(8),
//                            fitTotal->GetParameter(9),
//                            fitTotal->GetParameter(10),
//                            fitTotal->GetParameter(11));
//  
//  for ( Int_t i = 0; i < 5; ++i )
//  {
//    bck->SetParameter(i, fitTotal->GetParameter(i));
//  }
//
//  Set("kPol2Exp",fitTotal->GetParameter(0),fitTotal->GetParError(0));
//  Set("pol0",fitTotal->GetParameter(1),fitTotal->GetParError(1));
//  Set("pol1",fitTotal->GetParameter(2),fitTotal->GetParError(2));
//  Set("pol2",fitTotal->GetParameter(3),fitTotal->GetParError(3));
//  Set("exp",fitTotal->GetParameter(4),fitTotal->GetParError(4));
//
//  Set("alJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
//  Set("nlJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
//  Set("auJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
//  Set("nuJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));
//  
//  Set("kPsiP",fitTotal->GetParameter(12),fitTotal->GetParError(12));
//  
//  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);
//  
//  Double_t cbParameters[7];
//  Double_t covarianceMatrix[7][7];
//  
//  for ( int ix = 0; ix < 7; ++ix )
//  {
//    cbParameters[ix] = fitTotal->GetParameter(ix+5);
//  }
//  
//  for ( int iy = 0; iy < 7; ++iy )
//  {
//    for ( int ix = 0; ix < 7; ++ix )
//    {
//      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+5,iy+5);
//    }
//  }
//  
//  Double_t a = fHisto->GetXaxis()->GetXmin();
//  Double_t b = fHisto->GetXaxis()->GetXmax();
//  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
//  double nerr = signalJPsi->IntegralError(a,b,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
//  
//  Set("NofJPsi",njpsi,nerr);
//}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2POL4EXP()
{
  /// Fit using 2 extended crystal balls (signal) + pol4 x exp (background)
  /// 15 parameters
  
  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  Double_t paramSPsiP = GetValue("FSigmaPsiP");
  
  TString msg;
  
  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  
  AliDebug(1,Form("Fit with jpsi + psiprime Pol4 x Exp %s",msg.Data()));
  
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol4Exp,fitRangeLow,fitRangeHigh,14,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2Pol4Exp");
  
  fitTotal->SetParNames("pol0","pol1","pol2","pol3","pol4","exp","kJPsi",
  //                      0       1      2      3     4       5     6
                        "mJPsi","sJPsi","alJPsi","nlJPsi");
  //                       7      8       9         10
  fitTotal->SetParName(11,"auJPsi");
  //                        11
  fitTotal->SetParName(12,"nuJPsi");
  //                         12
  fitTotal->SetParName(13,"kPsiP");
  //                         13

  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Exp,fitRangeLow,fitRangeHigh,6,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Exp");
  
  //    bck->SetParameters(fHisto->GetMaximum(),0.01,0.01,0.01,0.01,0.01,1.);
  //
  //  SetFitRejectRange(2.7,3.5);
  //
  //  fHisto->Fit(bck,"R");
  //
  //  SetFitRejectRange();
  //
  //  for ( Int_t i = 0; i < 7; ++i )
  //  {
  //    fitTotal->SetParameter(i, bck->GetParameter(i));
  //  }
  
  //___
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Exp,1.6,7.,6,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Exp");
  
  Int_t bin = fHisto->FindBin(1.6);
  
  bckInit->SetParameters(fHisto->GetBinContent(bin),-fHisto->GetBinContent(bin),fHisto->GetBinContent(bin)/2.,-fHisto->GetBinContent(bin)/10.,fHisto->GetBinContent(bin)/100.,-2.);
  
  SetFitRejectRange(2.6,4.0);
  
  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,"SR");
  
  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  
  if ( static_cast<int>(fitResultInit) )
  {
    bin = fHisto->FindBin(0.82);
    bckInit->SetParameters(0.,1.,1.,1.,2.,0.5);
    fitResultInit = fHisto->Fit(bckInit,"SR");
  }
  
  SetFitRejectRange();
  
  for ( Int_t i = 0; i < 6; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }
  
  delete bckInit;
  //___
  
  //  bck->SetRange(fitRangeLow,fitRangeHigh);
  
  bin = fHisto->FindBin(3.09);
  fitTotal->SetParameter(6, fHisto->GetBinContent(bin)); // norm(kJPsi)
  
  fitTotal->SetParameter(7, 3.1); // mean
  fitTotal->SetParLimits(7, 3.0, 3.2);
  
  fitTotal->SetParameter(8, 0.08); // sigma
  fitTotal->SetParLimits(8, 0.03, 0.15);
  
  if ( IsValidValue(alphaLow) )
  {
    fitTotal->FixParameter(9, alphaLow);
  }
  else
  {
    fitTotal->SetParameter(9,0.9);
    fitTotal->SetParLimits(9,0.1,10.0);
  }
  
  if ( IsValidValue(nLow) )
  {
    fitTotal->FixParameter(10, nLow);
  }
  else
  {
    fitTotal->SetParameter(10,5.0);
    fitTotal->SetParLimits(10,0.0,10.0);
  }
  
  if ( IsValidValue(alphaUp) )
  {
    fitTotal->FixParameter(11, alphaUp);
  }
  else
  {
    fitTotal->SetParameter(11, 2.0);
    fitTotal->SetParLimits(11,0.1,10.0);
  }
  
  if ( IsValidValue(nUp) )
  {
    fitTotal->FixParameter(12, nUp);
  }
  else
  {
    fitTotal->SetParameter(12,3.0);
    fitTotal->SetParLimits(12,0.0,10.0);
  }
  
  bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(13, fHisto->GetBinContent(bin)); //kPsiP
  fitTotal->SetParLimits(13, 0., fHisto->GetBinContent(bin)*2); //kPsiP
  
  //  fitTotal->SetParameter(12, 3.7); // mean PsiP
  //  fitTotal->SetParLimits(12, 3.6, 3.71);
  //
  //  fitTotal->SetParameter(13, 0.08/3.096916*3.68609); // sigma PsiP
  //  fitTotal->SetParLimits(13, 0.03/3.096916*3.68609, 0.15/3.096916*3.68609);
  
  const char* fitOption = "SER";
  
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
  
  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    if ( fitTotal->GetParameter(13) <= fitTotal->GetParError(13) ) //kPsi'
    {
      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
      
      fitTotal->SetParameter(13, fHisto->GetBinContent(bin)/2.);
    }
    
    if ( fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) //kPol4Exp
    {
      std::cout << "//-------Refitting again (setting kPol4Exp norm= kPol4Exp norm /2)" << std::endl;
      bin = fHisto->FindBin(fitRangeLow);
      
      fitTotal->SetParameter(0, fHisto->GetBinContent(bin)/2.);
    }
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  
  if ( fitResult->CovMatrixStatus() != 3 )
  {
    std::cout << "//-------Refitting again (setting Bkg norm= Bkg norm/2)" << std::endl;
    fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
  }
  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    bin = fHisto->FindBin(3.68);
//    std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
//    fitTotal->SetParameter(13, fHisto->GetBinContent(bin)/2.); // kVWG
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//    
//    std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//  }
  
  //  if ( fitResult->CovMatrixStatus() != 3 )
  //  {
  //    std::cout << "//-------Refitting again (refitting Bkg)" << std::endl;
  //    if( fitTotal->GetParameter(9) >= 0.09 || fitTotal->GetParameter(9) <= 0.06 )
  //    {
  //      bck->SetParameters(0.01,fHisto->GetMaximum(),0.01,0.01,0.01,0.01,1.);
  ////      for ( Int_t i = 0; i < 7; ++i )
  ////      {
  ////        bck->SetParameter(i, 0.);
  ////      }
  //    }
  //    else
  //    {
  //      for ( Int_t i = 0; i < 7; ++i )
  //      {
  //        bck->SetParameter(i, fitTotal->GetParameter(i));
  //      }
  //    }
  //
  //    SetFitRejectRange(2.7,3.5);
  //
  //    fHisto->Fit(bck,"R");
  //
  //    SetFitRejectRange();
  //
  //    for ( Int_t i = 0; i < 7; ++i )
  //    {
  //      fitTotal->SetParameter(i, bck->GetParameter(i));
  //    }
  //
  //    fitResult = fHisto->Fit(fitTotal,fitOption,"");
  //  }
  //  std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
  
  
  
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);
  
  Set("mJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
  Set("sJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  
  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalJPsi->SetParameters(fitTotal->GetParameter(6),
                            fitTotal->GetParameter(7),
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12));
  
  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalPsiP->SetParameters(fitTotal->GetParameter(13),
                            3.68609+(fitTotal->GetParameter(7)-3.096916)/3.096916*3.68609,
                            fitTotal->GetParameter(8)*paramSPsiP, // /3.096916*3.68609,
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12));
  
  for ( Int_t i = 0; i < 7; ++i )
  {
    bck->SetParameter(i, fitTotal->GetParameter(i));
  }
  
  Set("alJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set("nlJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("auJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));
  Set("nuJPsi",fitTotal->GetParameter(12),fitTotal->GetParError(12));
  
//  Set("kPol4Exp",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("pol0",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("pol1",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("pol2",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("pol3",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("pol4",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set("exp",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  
  Set("kJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  Set("kPsiP",fitTotal->GetParameter(13),fitTotal->GetParError(13));
  
  
  
  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);
  
  Double_t cbParameters[7];
  Double_t covarianceMatrix[7][7];
  
  for ( int ix = 0; ix < 7; ++ix )
  {
    cbParameters[ix] = fitTotal->GetParameter(ix+6);
  }
  
  for ( int iy = 0; iy < 7; ++iy )
  {
    for ( int ix = 0; ix < 7; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+6,iy+6);
    }
  }
  
  Double_t a = fHisto->GetXaxis()->GetXmin();
  Double_t b = fHisto->GetXaxis()->GetXmax();
  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr = signalJPsi->IntegralError(a,b,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi",njpsi,nerr);
  
  double m = GetValue("mJPsi");
  double s = GetValue("sJPsi");
  double njpsi3s = signalJPsi->Integral(m-3*s,m+3*s)/fHisto->GetBinWidth(1);
  double nerr3s = signalJPsi->IntegralError(m-3*s,m+3*s,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi3s",njpsi3s,nerr3s);
  
  //Computation of bin significance and signal over background
  
  Double_t bkgParameters[6];
  Double_t bkgcovarianceMatrix[6][6];
  
  for ( int ix = 0; ix < 6; ++ix )
  {
    bkgParameters[ix] = fitTotal->GetParameter(ix);
  }
  
  for ( int iy = 0; iy < 6; ++iy )
  {
    for ( int ix = 0; ix < 6; ++ix )
    {
      bkgcovarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix,iy);
    }
  }
  
  double nbck3s = bck->Integral(m-3.*s,m+3.*s)/fHisto->GetBinWidth(1);
  double nbck3sErr = bck->IntegralError(m-3.*s,m+3.*s,&bkgParameters[0],&bkgcovarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  double sOverB3s = njpsi3s / nbck3s;
  double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));
  
  Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);
  
  double sig = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
  double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
                              TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );
  
  Set("Significance3s",sig,sigErr);

  
}

//_____________________________________________________________________________
//void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2POL4EXP()
//{
//  /// Fit using 2 extended crystal balls (signal) + pol4 x exp (background)
//  /// 15 parameters
//  
//  Double_t alphaLow = GetValue("alJPsi");
//  Double_t nLow = GetValue("nlJPsi");
//  Double_t alphaUp = GetValue("auJPsi");
//  Double_t nUp = GetValue("nuJPsi");
//  Double_t fitRangeLow = GetValue(kFitRangeLow);
//  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
//  
//  TString msg;
//  
//  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
//  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
//  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
//  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
//  
//  AliDebug(1,Form("Fit with jpsi + psiprime Pol4 x Exp %s",msg.Data()));
//  
//  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol4Exp,fitRangeLow,fitRangeHigh,15,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2Pol4Exp");
//  
//  fitTotal->SetParNames("kPol4Exp","pol0","pol1","pol2","pol3","pol4","exp","kJPsi",
////                          0         1      2      3     4       5     6      7    
//                        "mJPsi","sJPsi","alJPsi");
//  //                      8       9         10    
//  fitTotal->SetParName(11,"nlJPsi");
//  //                        11
//  fitTotal->SetParName(12,"auJPsi");
//  //                         12
//  fitTotal->SetParName(13,"nuJPsi");
//  //                         13
//  fitTotal->SetParName(14,"kPsiP");
//  //                        14
//  
//  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Exp,fitRangeLow,fitRangeHigh,7,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Exp");
//  
////    bck->SetParameters(fHisto->GetMaximum(),0.01,0.01,0.01,0.01,0.01,1.);
////  
////  SetFitRejectRange(2.7,3.5);
////
////  fHisto->Fit(bck,"R");
////  
////  SetFitRejectRange();
////  
////  for ( Int_t i = 0; i < 7; ++i )
////  {
////    fitTotal->SetParameter(i, bck->GetParameter(i));
////  }
//  
//  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Exp,1.6,7.,6,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Exp");
//  
//  Int_t bin = fHisto->FindBin(1.6);
//  
//  bckInit->SetParameters(fHisto->GetBinContent(bin),-fHisto->GetBinContent(bin)/2.,fHisto->GetBinContent(bin)/2.,-fHisto->GetBinContent(bin)/20,fHisto->GetBinContent(bin/100.),-2.);
//  
//  SetFitRejectRange(2.6,4.0);
//  
//  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,"SR");
//  
//  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
//  
//  if ( static_cast<int>(fitResultInit) )
//  {
//    bin = fHisto->FindBin(0.82);
//    bckInit->SetParameters(0.,1.,1.,1.,2.,0.5);
//    fitResultInit = fHisto->Fit(bckInit,"SR");
//  }
//  
//  SetFitRejectRange();
//  
//  for ( Int_t i = 0; i < 6; ++i )
//  {
//    fitTotal->SetParameter(i, bckInit->GetParameter(i));
//  }
//
//  
//  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Exp,1.2,6.,7,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Exp");
//  
//  Int_t bin = fHisto->FindBin(1.2);
//  
//  bckInit->SetParameters(-1.,fHisto->GetBinContent(bin),-fHisto->GetBinContent(bin),fHisto->GetBinContent(bin)/10,-fHisto->GetBinContent(bin/10),10.,1.);
//  
//  SetFitRejectRange(2.7,4.0);
//  
//  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,"SR");
//  
//  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
//  
//  if ( static_cast<int>(fitResultInit) )
//  {
//    bin = fHisto->FindBin(0.82);
//    bckInit->SetParameters(0.,1.,1.,1.,2.,0.5,0.3);
//    fitResultInit = fHisto->Fit(bckInit,"SR");
//  }
//
//  SetFitRejectRange();
//  
//  for ( Int_t i = 0; i < 7; ++i )
//  {
//    fitTotal->SetParameter(i, bckInit->GetParameter(i));
//  }
//  
//  delete bckInit;
//
//
////  bck->SetRange(fitRangeLow,fitRangeHigh);
//  
//  bin = fHisto->FindBin(3.09);
//  fitTotal->SetParameter(7, fHisto->GetBinContent(bin)); // norm(kJPsi)
//  
//  fitTotal->SetParameter(8, 3.1); // mean
//  fitTotal->SetParLimits(8, 3.0, 3.2);
//  
//  fitTotal->SetParameter(9, 0.08); // sigma
//  fitTotal->SetParLimits(9, 0.03, 0.15);
//  
//  if ( IsValidValue(alphaLow) )
//  {
//    fitTotal->FixParameter(10, alphaLow);
//  }
//  else
//  {
//    fitTotal->SetParameter(10,0.9);
//    fitTotal->SetParLimits(10,0.1,10.0);
//  }
//  
//  if ( IsValidValue(nLow) )
//  {
//    fitTotal->FixParameter(11, nLow);
//  }
//  else
//  {
//    fitTotal->SetParameter(11,5.0);
//    fitTotal->SetParLimits(11,0.0,10.0);
//  }
//  
//  if ( IsValidValue(alphaUp) )
//  {
//    fitTotal->FixParameter(12, alphaUp);
//  }
//  else
//  {
//    fitTotal->SetParameter(12, 2.0);
//    fitTotal->SetParLimits(12,0.1,10.0);
//  }
//  
//  if ( IsValidValue(nUp) )
//  {
//    fitTotal->FixParameter(13, nUp);
//  }
//  else
//  {
//    fitTotal->SetParameter(13,3.0);
//    fitTotal->SetParLimits(13,0.0,10.0);
//  }
//  
//  bin = fHisto->FindBin(3.68);
//  fitTotal->SetParameter(14, fHisto->GetBinContent(bin)); //kPsiP
//  fitTotal->SetParLimits(14, 0., fHisto->GetBinContent(bin)*2); //kPsiP
//  
//  //  fitTotal->SetParameter(12, 3.7); // mean PsiP
//  //  fitTotal->SetParLimits(12, 3.6, 3.71);
//  //
//  //  fitTotal->SetParameter(13, 0.08/3.096916*3.68609); // sigma PsiP
//  //  fitTotal->SetParLimits(13, 0.03/3.096916*3.68609, 0.15/3.096916*3.68609);
//  
//  const char* fitOption = "SER";
//  
//  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
//  
//  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
//   std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//  
//  if ( static_cast<int>(fitResult) )
//  {
//    if ( fitTotal->GetParameter(14) <= fitTotal->GetParError(14) ) //kPsi'
//    {
//      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
//      
//      fitTotal->SetParameter(14, fHisto->GetBinContent(bin)/2.);
//    }
//    
//    if ( fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) //kPol4Exp
//    {
//      std::cout << "//-------Refitting again (setting kPol4Exp norm= kPol4Exp norm /2)" << std::endl;
//      bin = fHisto->FindBin(fitRangeLow);
//      
//      fitTotal->SetParameter(0, fHisto->GetBinContent(bin)/2.);
//    }
//    
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
//  }
//  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    std::cout << "//-------Refitting again (setting Bkg norm= Bkg norm/2)" << std::endl;
//    fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//  }
//  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    bin = fHisto->FindBin(3.68);
//    std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
//    fitTotal->SetParameter(14, fHisto->GetBinContent(bin)/2.); // kVWG
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//    
//    std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//  }
//  
////  if ( fitResult->CovMatrixStatus() != 3 )
////  {
////    std::cout << "//-------Refitting again (refitting Bkg)" << std::endl;
////    if( fitTotal->GetParameter(9) >= 0.09 || fitTotal->GetParameter(9) <= 0.06 )
////    {
////      bck->SetParameters(0.01,fHisto->GetMaximum(),0.01,0.01,0.01,0.01,1.);
//////      for ( Int_t i = 0; i < 7; ++i )
//////      {
//////        bck->SetParameter(i, 0.);
//////      }
////    }
////    else
////    {
////      for ( Int_t i = 0; i < 7; ++i )
////      {
////        bck->SetParameter(i, fitTotal->GetParameter(i));
////      }
////    }
////    
////    SetFitRejectRange(2.7,3.5);
////    
////    fHisto->Fit(bck,"R");
////    
////    SetFitRejectRange();
////    
////    for ( Int_t i = 0; i < 7; ++i )
////    {
////      fitTotal->SetParameter(i, bck->GetParameter(i));
////    }
////    
////    fitResult = fHisto->Fit(fitTotal,fitOption,"");
////  }
////  std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//
//
//  
//  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
//  Set("FitNDF",fitTotal->GetNDF(),0.0);
//  
//  Set("mJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
//  Set("sJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
//  
//  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
//  signalJPsi->SetParameters(fitTotal->GetParameter(7),
//                            fitTotal->GetParameter(8),
//                            fitTotal->GetParameter(9),
//                            fitTotal->GetParameter(10),
//                            fitTotal->GetParameter(11),
//                            fitTotal->GetParameter(12),
//                            fitTotal->GetParameter(13));
//  
//  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
//  signalPsiP->SetParameters(fitTotal->GetParameter(14),
//                            3.68609+(fitTotal->GetParameter(8)-3.096916)/3.096916*3.68609,
//                            fitTotal->GetParameter(9)/3.096916*3.68609,
//                            fitTotal->GetParameter(10),
//                            fitTotal->GetParameter(11),
//                            fitTotal->GetParameter(12),
//                            fitTotal->GetParameter(13));
//
//  for ( Int_t i = 0; i < 7; ++i )
//  {
//    bck->SetParameter(i, fitTotal->GetParameter(i));
//  }
//  
//  Set("alJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
//  Set("nlJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));
//  Set("auJPsi",fitTotal->GetParameter(12),fitTotal->GetParError(12));
//  Set("nuJPsi",fitTotal->GetParameter(13),fitTotal->GetParError(13));
//  
//  Set("kPol4Exp",fitTotal->GetParameter(0),fitTotal->GetParError(0));
//  Set("pol0",fitTotal->GetParameter(1),fitTotal->GetParError(1));
//  Set("pol1",fitTotal->GetParameter(2),fitTotal->GetParError(2));
//  Set("pol2",fitTotal->GetParameter(3),fitTotal->GetParError(3));
//  Set("pol3",fitTotal->GetParameter(4),fitTotal->GetParError(4));
//  Set("pol4",fitTotal->GetParameter(5),fitTotal->GetParError(5));
//  Set("exp",fitTotal->GetParameter(6),fitTotal->GetParError(6));
//  
//  Set("kJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
//  Set("kPsiP",fitTotal->GetParameter(14),fitTotal->GetParError(14));
//  
// 
//  
//  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);
//  
//  Double_t cbParameters[7];
//  Double_t covarianceMatrix[7][7];
//  
//  for ( int ix = 0; ix < 7; ++ix )
//  {
//    cbParameters[ix] = fitTotal->GetParameter(ix+7);
//  }
//  
//  for ( int iy = 0; iy < 7; ++iy )
//  {
//    for ( int ix = 0; ix < 7; ++ix )
//    {
//      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+7,iy+7);
//    }
//  }
//  
//  Double_t a = fHisto->GetXaxis()->GetXmin();
//  Double_t b = fHisto->GetXaxis()->GetXmax();
//  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
//  double nerr = signalJPsi->IntegralError(a,b,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
//  
//  Set("NofJPsi",njpsi,nerr);
//
//}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMENA60NEWVWG()
{
  /// Fit using 2 NA60(new) (signal) + variable width gaussian (background)
  
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");
  
  Double_t paramSPsiP = GetValue("FSigmaPsiP");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  TString msg;
  
  if (IsValidValue(p1Left)) msg += TString::Format("p1L=%e ",p1Left);
  if (IsValidValue(p2Left)) msg += TString::Format("p2L=%e ",p2Left);
  if (IsValidValue(p3Left)) msg += TString::Format("p3L=%e ",p3Left);
  if (IsValidValue(p1Right)) msg += TString::Format("p1R=%e ",p1Right);
  if (IsValidValue(p2Right)) msg += TString::Format("p2R=%e ",p2Right);
  if (IsValidValue(p3Right)) msg += TString::Format("p3R=%e ",p3Right);
  
  if (IsValidValue(alphaLeft)) msg += TString::Format("aL=%e ",alphaLeft);
  if (IsValidValue(alphaRight)) msg += TString::Format("aR=%e ",alphaRight);
  
  AliDebug(1,Form("Fit with jpsi + psiprime NA60 new and VWG %s",msg.Data()));
  
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewVWG,fitRangeLow,fitRangeHigh,16,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoNA60NewVWG");
  
  fitTotal->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi",
  //                       0      1      2       3       4       5       6
                        "p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi");
  //                        7         8         9        10
  fitTotal->SetParName(11, "p2RJPsi");
  //                           11
  fitTotal->SetParName(12, "p3RJPsi");
  //                           12
  fitTotal->SetParName(13, "aLJPsi");
  //                           13
  fitTotal->SetParName(14, "aRJPsi");
  //                           14
  fitTotal->SetParName(15, "kPsiP");
  //                           15
  
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG");
  
  
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG,1.7,6.,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG");
  
  Int_t bin = fHisto->FindBin(0.26);
  
  bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
  
  SetFitRejectRange(2.7,4.0);
  
  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,"SR");
  
  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  
  if ( static_cast<int>(fitResultInit) )
  {
    bin = fHisto->FindBin(0.82);
    bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
    fitResultInit = fHisto->Fit(bckInit,"SR");
  }
  else if ( bckInit->GetParameter(0) < 0 )
  {
    bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
  }
  
  SetFitRejectRange();
  
  for ( Int_t i = 0; i < 4; ++i )
  {
    TString name(GetName());
    if (name.Contains("NA60NEWVWG_2.2_4.7_SP1.0"))
    {
      fitTotal->SetParameter(i, bckInit->GetParameter(i)*0.8);
    }
    else fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }
  
  delete bckInit;

  
  
//  fitTotal->SetParameter(0, fHisto->GetMaximum()); // kVWG
//  fitTotal->SetParameter(1, 1.9); // mVWG
//  
//  fitTotal->SetParameter(2, 0.5); // sVWG1
//  fitTotal->SetParLimits(2, 0., 100.);
//  
//  fitTotal->SetParameter(3, 0.3); // sVWG2
//  fitTotal->SetParLimits(3, 0., 100.);
//  
//  Int_t bin = fHisto->FindBin(3.09);
//  fitTotal->SetParameter(4, fHisto->GetBinContent(bin)); // norm
  
  fitTotal->SetParameter(5, 3.1); // mean
  fitTotal->SetParLimits(5, 3.0, 3.2);
  
  fitTotal->SetParameter(6, 0.08); // sigma
  fitTotal->SetParLimits(6, 0.05, 0.15);
  
  fitTotal->FixParameter(7, p1Left);
  fitTotal->FixParameter(8, p2Left);
  fitTotal->FixParameter(9, p3Left);
  fitTotal->FixParameter(10, p1Right);
  fitTotal->FixParameter(11, p2Right);
  fitTotal->FixParameter(12, p3Right);
  
  fitTotal->FixParameter(13, alphaLeft);
  fitTotal->FixParameter(14, alphaRight);
  
  bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(15, fHisto->GetBinContent(bin)); //kPsiP
  fitTotal->SetParLimits(15, 0., fHisto->GetBinContent(bin)*2); //kPsiP
  
  const char* fitOption = "SER";
  
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
  
  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
   std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
  
  // Check parameters...
  //  if (
  //      StrongCorrelation(fitResult,fitTotal,3,4,2) ||
  //      StrongCorrelation(fitResult,fitTotal,5,6,2) ||
  //      WrongParameter(fitTotal,3,1) ||
  //      WrongParameter(fitTotal,4,0)  ||
  //      WrongParameter(fitTotal,5,1)  ||
  //      WrongParameter(fitTotal,6,0)
  //      )
  //  {
  //    // try again...
  //    fitResult = fHisto->Fit(fitTotal,"SER","");
  //  }
  
  if ( static_cast<int>(fitResult) )
  {
    if ( 0.5*fitTotal->GetParameter(15) <= fitTotal->GetParError(15) ) //kPsi'
    {
      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
      
      fitTotal->SetParameter(15, fHisto->GetBinContent(bin)/2.);
    }
    
    if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) //kVWG
    {
      std::cout << "//-------Refitting again (setting kVWG norm= kVWG norm /2)" << std::endl;
      bin = fHisto->FindBin(fitRangeLow);
      
      fitTotal->SetParameter(0, fHisto->GetBinContent(bin)/2.);
    }
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  
  if ( static_cast<int>(fitResult) )
  {
    std::cout << "//-------Refitting again (setting kVWG=kVWG*2)" << std::endl;
    
    fitTotal->SetParameter(0, fHisto->GetMaximum()*2.); // kVWG
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  
  if ( static_cast<int>(fitResult) )
  {
    std::cout << "//-------Refitting again (setting kVWG=kVWG/2)" << std::endl;
    
    fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    std::cout << "//-------Refitting again (setting Bkg norm= Bkg norm/2)" << std::endl;
//    fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//  }
//  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    bin = fHisto->FindBin(3.68);
//    std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
//    fitTotal->SetParameter(15, fHisto->GetBinContent(bin)/2.); // kVWG
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//    
//    std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//  }
//  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    std::cout << "//-------Refitting again (refitting Bkg)" << std::endl;
//    
//    for ( Int_t i = 0; i < 4; ++i )
//    {
//      bck->SetParameter(i, fitTotal->GetParameter(i));
//    }
//    
//    SetFitRejectRange(2.7,3.5);
//    
//    fHisto->Fit(bck,"R");
//    
//    SetFitRejectRange();
//    
//    for ( Int_t i = 0; i < 4; ++i )
//    {
//      fitTotal->SetParameter(i, bck->GetParameter(i));
//    }
//    
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//  }
//  std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
//  if ( static_cast<int>(fitResult) )
//  {
//    fitTotal->SetParameter(0, fHisto->GetMaximum()); // kVWG
//    fitTotal->SetParameter(1, 1.9); // mVWG
//    
//    fitTotal->SetParameter(2, 0.5); // sVWG1
//    fitTotal->SetParLimits(2, 0., 100.);
//    
//    fitTotal->SetParameter(3, 0.3); // sVWG2
//    fitTotal->SetParLimits(3, 0., 100.);
//
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//    
//    std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
//    
//    if ( static_cast<int>(fitResult) )
//    {
//      if ( fitTotal->GetParameter(15) <= fitTotal->GetParError(15) ) //kPsi'
//      {
//        std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
//        
//        fitTotal->SetParameter(15, fHisto->GetBinContent(bin)/2.);
//      }
//      
//      if ( fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) //kVWG
//      {
//        std::cout << "//-------Refitting again (setting kVWG norm= kVWG norm /2)" << std::endl;
//        bin = fHisto->FindBin(fitRangeLow);
//        
//        fitTotal->SetParameter(0, fHisto->GetBinContent(bin)/2.);
//      }
//      
//      fitResult = fHisto->Fit(fitTotal,fitOption,"");
//      std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
//    }
//
//  }
  
  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","FitFunctionNA60New");
  
  signalJPsi->SetParameters(fitTotal->GetParameter(4),
                            fitTotal->GetParameter(5),
                            fitTotal->GetParameter(6),
                            fitTotal->GetParameter(7),
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13));
  
  signalJPsi->SetParameter(10,fitTotal->GetParameter(14));
  
  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","FitFunctionNA60New");
  
  signalPsiP->SetParameters(fitTotal->GetParameter(15),
                            3.68609+(fitTotal->GetParameter(5)-3.096916)/3.096916*3.68609,
                            fitTotal->GetParameter(6)*paramSPsiP, // /3.096916*3.68609,
                            fitTotal->GetParameter(7),
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13));
  
  signalPsiP->SetParameter(10,fitTotal->GetParameter(14));
  
  bck->SetParameter(0,fitTotal->GetParameter(0));
  bck->SetParameter(1,fitTotal->GetParameter(1));
  bck->SetParameter(2,fitTotal->GetParameter(2));
  bck->SetParameter(3,fitTotal->GetParameter(3));
  
  Set("FitResult",static_cast<int>(fitResult)*1.0,0.0);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);
  
  Set("kVWG",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("mVWG",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("sVWG1",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("sVWG2",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  
  Set("kJPsi",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set("kPsiP",fitTotal->GetParameter(15),fitTotal->GetParError(15));
  
  Set("mJPsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set("sJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  
  Set("p1LJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
  Set("p2LJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("p3LJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set("p1RJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("p2RJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));
  Set("p3RJPsi",fitTotal->GetParameter(12),fitTotal->GetParError(12));
  
  Set("aLJPsi",fitTotal->GetParameter(13),fitTotal->GetParError(13));
  Set("aRJPsi",fitTotal->GetParameter(14),fitTotal->GetParError(14));
  
  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);
  
  Double_t na60Parameters[11];
  Double_t covarianceMatrix[11][11];
  
  for ( int ix = 0; ix < 11; ++ix )
  {
    na60Parameters[ix] = fitTotal->GetParameter(ix+4);
  }
  
  for ( int iy = 0; iy < 11; ++iy )
  {
    for ( int ix = 0; ix < 11; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+4,iy+4);
    }
  }

  
  Double_t a = fHisto->GetXaxis()->GetXmin();
  Double_t b = fHisto->GetXaxis()->GetXmax();
  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr = signalJPsi->IntegralError(a,b,&na60Parameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi",njpsi,nerr);
  
  double m = GetValue("mJPsi");
  double s = GetValue("sJPsi");
  double njpsi3s = signalJPsi->Integral(m-3*s,m+3*s)/fHisto->GetBinWidth(1);
  double nerr3s = signalJPsi->IntegralError(m-3*s,m+3*s,&na60Parameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi3s",njpsi3s,nerr3s);
  
  //Computation of bin significance and signal over background
  
  Double_t bkgParameters[4];
  Double_t bkgcovarianceMatrix[4][4];
  
  for ( int ix = 0; ix < 4; ++ix )
  {
    bkgParameters[ix] = fitTotal->GetParameter(ix);
  }
  
  for ( int iy = 0; iy < 4; ++iy )
  {
    for ( int ix = 0; ix < 4; ++ix )
    {
      bkgcovarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix,iy);
    }
  }
  
  double nbck3s = bck->Integral(m-3.*s,m+3.*s)/fHisto->GetBinWidth(1);
  double nbck3sErr = bck->IntegralError(m-3.*s,m+3.*s,&bkgParameters[0],&bkgcovarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  double sOverB3s = njpsi3s / nbck3s;
  double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));
  
  Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);
  
  double sig = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
  double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
                              TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );
  
  Set("Significance3s",sig,sigErr);

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMENA60NEWPOL2EXP()
{
  /// Fit using 2 NA60(new) (signal) + pol2 x exp (background)
  
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");
  
  Double_t paramSPsiP = GetValue("FSigmaPsiP");  
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  TString msg;
  
  if (IsValidValue(p1Left)) msg += TString::Format("p1L=%e ",p1Left);
  if (IsValidValue(p2Left)) msg += TString::Format("p2L=%e ",p2Left);
  if (IsValidValue(p3Left)) msg += TString::Format("p3L=%e ",p3Left);
  if (IsValidValue(p1Right)) msg += TString::Format("p1R=%e ",p1Right);
  if (IsValidValue(p2Right)) msg += TString::Format("p2R=%e ",p2Right);
  if (IsValidValue(p3Right)) msg += TString::Format("p3R=%e ",p3Right);
  
  if (IsValidValue(alphaLeft)) msg += TString::Format("aL=%e ",alphaLeft);
  if (IsValidValue(alphaRight)) msg += TString::Format("aR=%e ",alphaRight);
  
  AliDebug(1,Form("Fit with jpsi + psiprime NA60 new and pol4 x exp %s",msg.Data()));
  
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewPol2Exp,fitRangeLow,fitRangeHigh,16,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoNA60NewPol2Exp");
  
  fitTotal->SetParNames("pol0","pol1","pol2","exp","kJPsi","mJPsi","sJPsi",
  //                      0      1       2      3     4      5       6
                        "p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi");
  //                        7         8         9        10
                      
  fitTotal->SetParName(11, "p2RJPsi");
  //                           11
  fitTotal->SetParName(12, "p3RJPsi");
  //                           12
  fitTotal->SetParName(13, "aLJPsi");
  //                           13
  fitTotal->SetParName(14, "aRJPsi");
  //                           14
  fitTotal->SetParName(15, "kPsiP");
  //                           15
  
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  
  //___
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,1.7,6.,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  
  Int_t bin = fHisto->FindBin(0.26);
  
  bckInit->SetParameters(fHisto->GetBinContent(bin),-fHisto->GetBinContent(bin)/3.,100.,0.05);//fHisto->GetBinContent(bin)
  
  SetFitRejectRange(2.7,4.0);
  
  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,"SR");
  
  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  
  if ( static_cast<int>(fitResultInit) )
  {
    bin = fHisto->FindBin(0.82);
    bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
    fitResultInit = fHisto->Fit(bckInit,"SR");
  }
  
  SetFitRejectRange();
  
  for ( Int_t i = 0; i < 4; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }
  
  delete bckInit;
  //__

  
//  bck->SetParameters(fHisto->GetMaximum(),0.,0.,0.,1.);
//  
//  SetFitRejectRange(2.7,3.5);
//  
//  fHisto->Fit(bck,"R");
//  
//  SetFitRejectRange();
//  
//  for ( Int_t i = 0; i < 5; ++i )
//  {
//    fitTotal->SetParameter(i, bck->GetParameter(i));
//  }

  
  bin = fHisto->FindBin(3.09);
  fitTotal->SetParameter(4, fHisto->GetBinContent(bin)); // kJPsi
  
  fitTotal->SetParameter(5, 3.1); // mean
  fitTotal->SetParLimits(5, 3.0, 3.2);
  
  fitTotal->SetParameter(6, 0.08); // sigma
  fitTotal->SetParLimits(6, 0.05, 0.15);
  
  fitTotal->FixParameter(7, p1Left);
  fitTotal->FixParameter(8, p2Left);
  fitTotal->FixParameter(9, p3Left);
  fitTotal->FixParameter(10, p1Right);
  fitTotal->FixParameter(11, p2Right);
  fitTotal->FixParameter(12, p3Right);
  
  fitTotal->FixParameter(13, alphaLeft);
  fitTotal->FixParameter(14, alphaRight);
  
  bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(15, fHisto->GetBinContent(bin)); //kPsiP
  fitTotal->SetParLimits(15, 0., fHisto->GetBinContent(bin)*2); //kPsiP
  
  const char* fitOption = "SER";
  
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
  
  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
   std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
  
  // Check parameters...
  //  if (
  //      StrongCorrelation(fitResult,fitTotal,3,4,2) ||
  //      StrongCorrelation(fitResult,fitTotal,5,6,2) ||
  //      WrongParameter(fitTotal,3,1) ||
  //      WrongParameter(fitTotal,4,0)  ||
  //      WrongParameter(fitTotal,5,1)  ||
  //      WrongParameter(fitTotal,6,0)
  //      )
  //  {
  //    // try again...
  //    fitResult = fHisto->Fit(fitTotal,"SER","");
  //  }
  
  if ( static_cast<int>(fitResult) )
  {
    if ( 0.5*fitTotal->GetParameter(15) <= fitTotal->GetParError(15) ) //kPsi'
    {
      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
      bin = fHisto->FindBin(3.68);
      fitTotal->SetParameter(15, fHisto->GetBinContent(bin)/2.);
    }
    
    if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) //kPol2Exp
    {
      std::cout << "//-------Refitting again (setting kPol2Exp norm= kPol2Exp norm /2)" << std::endl;
      bin = fHisto->FindBin(fitRangeLow);
      
      fitTotal->SetParameter(0, fHisto->GetBinContent(bin)/2.);
    }
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  
  if ( fitResult->CovMatrixStatus() != 3 )
  {
    if ( 0.5*fitTotal->GetParameter(15) <= fitTotal->GetParError(15) ) //kPsi'
    {
      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
      
      bin = fHisto->FindBin(3.68);
      fitTotal->FixParameter(15, 0.);
    }
    
    else
    {
      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
      
      fHisto->FindBin(3.68);
      fitTotal->SetParameter(15, fHisto->GetBinContent(bin)/2.);
    }
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    
    std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
  }

  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    std::cout << "//-------Refitting again (setting Bkg norm= Bkg norm/2)" << std::endl;
//    fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//  }
//  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    bin = fHisto->FindBin(3.68);
//    std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
//    fitTotal->SetParameter(16, fHisto->GetBinContent(bin)/2.); // kVWG
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//    
//    std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//  }
//  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    std::cout << "//-------Refitting again (refitting Bkg)" << std::endl;
//    
//    for ( Int_t i = 0; i < 5; ++i )
//    {
//      bck->SetParameter(i, fitTotal->GetParameter(i));
//    }
//    
//    SetFitRejectRange(2.7,3.5);
//    
//    fHisto->Fit(bck,"R");
//    
//    SetFitRejectRange();
//    
//    for ( Int_t i = 0; i < 5; ++i )
//    {
//      fitTotal->SetParameter(i, bck->GetParameter(i));
//    }
//    
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//  }
//  std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;


  
  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","SignalNA60New");
  
  signalJPsi->SetParameters(fitTotal->GetParameter(4),
                            fitTotal->GetParameter(5),
                            fitTotal->GetParameter(6),
                            fitTotal->GetParameter(7),
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13));
  
  signalJPsi->SetParameter(10,fitTotal->GetParameter(14));
  
  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","SignalNA60New");
  
  signalPsiP->SetParameters(fitTotal->GetParameter(15),
                            3.68609+(fitTotal->GetParameter(5)-3.096916)/3.096916*3.68609,
                            fitTotal->GetParameter(6)*paramSPsiP, // /3.096916*3.68609,
                            fitTotal->GetParameter(7),
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13));
  
  signalPsiP->SetParameter(10,fitTotal->GetParameter(14));
  
  for ( Int_t i = 0; i < 4; ++i )
  {
    bck->SetParameter(i, fitTotal->GetParameter(i));
  }

  
  Set("FitResult",static_cast<int>(fitResult)*1.0,0.0);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);
  
//  Set("kPol2Exp",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("pol0",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("pol1",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("pol2",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("exp",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  
  Set("kJPsi",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set("kPsiP",fitTotal->GetParameter(15),fitTotal->GetParError(15));
  
  Set("mJPsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set("sJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  
  Set("p1LJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
  Set("p2LJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("p3LJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set("p1RJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("p2RJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));
  Set("p3RJPsi",fitTotal->GetParameter(12),fitTotal->GetParError(12));
  
  Set("aLJPsi",fitTotal->GetParameter(13),fitTotal->GetParError(13));
  Set("aRJPsi",fitTotal->GetParameter(14),fitTotal->GetParError(14));
  
  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);
  
  Double_t na60Parameters[11];
  Double_t covarianceMatrix[11][11];
  
  for ( int ix = 0; ix < 11; ++ix )
  {
    na60Parameters[ix] = fitTotal->GetParameter(ix+4);
  }
  
  for ( int iy = 0; iy < 11; ++iy )
  {
    for ( int ix = 0; ix < 11; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+4,iy+4);
    }
  }
  
  Double_t a = fHisto->GetXaxis()->GetXmin();
  Double_t b = fHisto->GetXaxis()->GetXmax();
  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr = signalJPsi->IntegralError(a,b,&na60Parameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi",njpsi,nerr);
  
  double m = GetValue("mJPsi");
  double s = GetValue("sJPsi");
  double njpsi3s = signalJPsi->Integral(m-3*s,m+3*s)/fHisto->GetBinWidth(1);
  double nerr3s = signalJPsi->IntegralError(m-3*s,m+3*s,&na60Parameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi3s",njpsi3s,nerr3s);
  
  //Computation of bin significance and signal over background
  
  Double_t bkgParameters[4];
  Double_t bkgcovarianceMatrix[4][4];
  
  for ( int ix = 0; ix < 4; ++ix )
  {
    bkgParameters[ix] = fitTotal->GetParameter(ix);
  }
  
  for ( int iy = 0; iy < 4; ++iy )
  {
    for ( int ix = 0; ix < 4; ++ix )
    {
      bkgcovarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix,iy);
    }
  }
  
  double nbck3s = bck->Integral(m-3.*s,m+3.*s)/fHisto->GetBinWidth(1);
  double nbck3sErr = bck->IntegralError(m-3.*s,m+3.*s,&bkgParameters[0],&bkgcovarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  double sOverB3s = njpsi3s / nbck3s;
  double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));
  
  Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);
  
  double sig = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
  double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
                              TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );
  
  Set("Significance3s",sig,sigErr);

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMENA60NEWPOL4EXP()
{
  /// Fit using 2 NA60(new) (signal) + pol4 x exp (background)
  
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");
  
  Double_t paramSPsiP = GetValue("FSigmaPsiP");
    
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  TString msg;
  
  if (IsValidValue(p1Left)) msg += TString::Format("p1L=%e ",p1Left);
  if (IsValidValue(p2Left)) msg += TString::Format("p2L=%e ",p2Left);
  if (IsValidValue(p3Left)) msg += TString::Format("p3L=%e ",p3Left);
  if (IsValidValue(p1Right)) msg += TString::Format("p1R=%e ",p1Right);
  if (IsValidValue(p2Right)) msg += TString::Format("p2R=%e ",p2Right);
  if (IsValidValue(p3Right)) msg += TString::Format("p3R=%e ",p3Right);
  
  if (IsValidValue(alphaLeft)) msg += TString::Format("aL=%e ",alphaLeft);
  if (IsValidValue(alphaRight)) msg += TString::Format("aR=%e ",alphaRight);
  
  AliDebug(1,Form("Fit with jpsi + psiprime NA60 new and pol4 x exp %s",msg.Data()));
  
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewPol4Exp,fitRangeLow,fitRangeHigh,18,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoNA60NewPol4Exp");
  
  fitTotal->SetParNames("pol0","pol1","pol2","pol3","pol4","exp","kJPsi",
  //                        0    1       2      3      4     5      6
                        "mJPsi","sJPsi","p1LJPsi","p2LJPsi");
  //                       7       8       9       10                      
  
  fitTotal->SetParName(11, "p3LJPsi");
  //                           11
  fitTotal->SetParName(12, "p1RJPsi");
  //                           12
  fitTotal->SetParName(13, "p2RJPsi");
  //                           13
  fitTotal->SetParName(14, "p3RJPsi");
  //                           14
  fitTotal->SetParName(15, "aLJPsi");
  //                           15
  fitTotal->SetParName(16, "aRJPsi");
  //                           16
  fitTotal->SetParName(17, "kPsiP");
  //                           17
  
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Exp,fitRangeLow,fitRangeHigh,6,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Exp");
  
  //___
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Exp,1.6,7.,6,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Exp");
  
  Int_t bin = fHisto->FindBin(1.6);
  
  bckInit->SetParameters(fHisto->GetBinContent(bin),-fHisto->GetBinContent(bin),fHisto->GetBinContent(bin)/2.,-fHisto->GetBinContent(bin)/10.,fHisto->GetBinContent(bin)/100.,-2.);
  
  SetFitRejectRange(2.6,4.0);
  
  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,"SR");
  
  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  
  if ( static_cast<int>(fitResultInit) )
  {
    bin = fHisto->FindBin(0.82);
    bckInit->SetParameters(0.,1.,1.,1.,2.,0.5);
    fitResultInit = fHisto->Fit(bckInit,"SR");
  }
  
  SetFitRejectRange();
  
  for ( Int_t i = 0; i < 6; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }
  
  delete bckInit;
  //___

//  
//  bck->SetParameters(fHisto->GetMaximum(),0.,0.,0.,0.,0.,1.);
//  
//  SetFitRejectRange(2.7,3.5);
//  
//  fHisto->Fit(bck,"R");
//  
//  SetFitRejectRange();
//  
//  for ( Int_t i = 0; i < 7; ++i )
//  {
//    fitTotal->SetParameter(i, bck->GetParameter(i));
//  }
  
  
  bin = fHisto->FindBin(3.09);
  fitTotal->SetParameter(6, fHisto->GetBinContent(bin)); // kJPsi
  
  fitTotal->SetParameter(7, 3.1); // mean
  fitTotal->SetParLimits(7, 3.0, 3.2);
  
  fitTotal->SetParameter(8, 0.08); // sigma
  fitTotal->SetParLimits(8, 0.05, 0.15);
  
  fitTotal->FixParameter(9, p1Left);
  fitTotal->FixParameter(10, p2Left);
  fitTotal->FixParameter(11, p3Left);
  fitTotal->FixParameter(12, p1Right);
  fitTotal->FixParameter(13, p2Right);
  fitTotal->FixParameter(14, p3Right);
  
  fitTotal->FixParameter(15, alphaLeft);
  fitTotal->FixParameter(16, alphaRight);
  
  bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(17, fHisto->GetBinContent(bin)); //kPsiP
  fitTotal->SetParLimits(17, 0., fHisto->GetBinContent(bin)*2); //kPsiP
  
  const char* fitOption = "SER";
  
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
  
  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
   std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
  
  // Check parameters...
  //  if (
  //      StrongCorrelation(fitResult,fitTotal,3,4,2) ||
  //      StrongCorrelation(fitResult,fitTotal,5,6,2) ||
  //      WrongParameter(fitTotal,3,1) ||
  //      WrongParameter(fitTotal,4,0)  ||
  //      WrongParameter(fitTotal,5,1)  ||
  //      WrongParameter(fitTotal,6,0)
  //      )
  //  {
  //    // try again...
  //    fitResult = fHisto->Fit(fitTotal,"SER","");
  //  }
  
  
  if ( static_cast<int>(fitResult) )
  {
    if ( 0.5*fitTotal->GetParameter(17) <= fitTotal->GetParError(17) ) //kPsi'
    {
      std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
      
      fitTotal->SetParameter(17, fHisto->GetBinContent(bin)/2.);
    }
    
    if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) //kPol4Exp
    {
      std::cout << "//-------Refitting again (setting kPol4Exp norm= kPol4Exp norm /2)" << std::endl;
      bin = fHisto->FindBin(fitRangeLow);
      
      fitTotal->SetParameter(0, fHisto->GetBinContent(bin)/2.);
    }
    
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  if ( fitResult->CovMatrixStatus() != 3 )
  {
    std::cout << "//-------Refitting again (setting Bkg norm= Bkg norm/2)" << std::endl;
    fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
  }
  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    bin = fHisto->FindBin(3.68);
//    std::cout << "//-------Refitting again (setting Psi'norm= Psi'norm/2)" << std::endl;
//    fitTotal->SetParameter(18, fHisto->GetBinContent(bin)/2.); // kVWG
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//    
//    std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;
//  }
//  
//  if ( fitResult->CovMatrixStatus() != 3 )
//  {
//    std::cout << "//-------Refitting again (refitting Bkg)" << std::endl;
//    
//    for ( Int_t i = 0; i < 7; ++i )
//    {
//      bck->SetParameter(i, fitTotal->GetParameter(i));
//    }
//    
//    SetFitRejectRange(2.7,3.5);
//    
//    fHisto->Fit(bck,"R");
//    
//    SetFitRejectRange();
//    
//    for ( Int_t i = 0; i < 7; ++i )
//    {
//      fitTotal->SetParameter(i, bck->GetParameter(i));
//    }
//    
//    fitResult = fHisto->Fit(fitTotal,fitOption,"");
//  }
//  std::cout << "CovMatrixStatus=" << fitResult->CovMatrixStatus() << std::endl;

  
  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","SignalNA60New");
  
  signalJPsi->SetParameters(fitTotal->GetParameter(6),
                            fitTotal->GetParameter(7),
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13),
                            fitTotal->GetParameter(14),
                            fitTotal->GetParameter(15));
  
  signalJPsi->SetParameter(10,fitTotal->GetParameter(16));
  
  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","SignalNA60New");
  
  signalPsiP->SetParameters(fitTotal->GetParameter(17),
                            3.68609+(fitTotal->GetParameter(7)-3.096916)/3.096916*3.68609,
                            fitTotal->GetParameter(8)*paramSPsiP, // /3.096916*3.68609,
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13),
                            fitTotal->GetParameter(14),
                            fitTotal->GetParameter(15));
  
  signalPsiP->SetParameter(10,fitTotal->GetParameter(16));
  
  for ( Int_t i = 0; i < 6; ++i )
  {
    bck->SetParameter(i, fitTotal->GetParameter(i));
  }

  
  Set("FitResult",static_cast<int>(fitResult)*1.0,0.0);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);
  
//  Set("kPol4Exp",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("pol0",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("pol1",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("pol2",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("pol3",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("pol4",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set("exp",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  
  Set("kJPsi",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set("kPsiP",fitTotal->GetParameter(17),fitTotal->GetParError(17));
  
  Set("mJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
  Set("sJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  
  Set("p1LJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set("p2LJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("p3LJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));
  Set("p1RJPsi",fitTotal->GetParameter(12),fitTotal->GetParError(12));
  Set("p2RJPsi",fitTotal->GetParameter(13),fitTotal->GetParError(13));
  Set("p3RJPsi",fitTotal->GetParameter(14),fitTotal->GetParError(14));
  
  Set("aLJPsi",fitTotal->GetParameter(15),fitTotal->GetParError(15));
  Set("aRJPsi",fitTotal->GetParameter(16),fitTotal->GetParError(16));
  
  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);
  
  Double_t na60Parameters[11];
  Double_t covarianceMatrix[11][11];
  
  for ( int ix = 0; ix < 11; ++ix )
  {
    na60Parameters[ix] = fitTotal->GetParameter(ix+6);
  }
  
  for ( int iy = 0; iy < 11; ++iy )
  {
    for ( int ix = 0; ix < 11; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+6,iy+6);
    }
  }
  
  Double_t a = fHisto->GetXaxis()->GetXmin();
  Double_t b = fHisto->GetXaxis()->GetXmax();
  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr = signalJPsi->IntegralError(a,b,&na60Parameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi",njpsi,nerr);
  
  double m = GetValue("mJPsi");
  double s = GetValue("sJPsi");
  double njpsi3s = signalJPsi->Integral(m-3*s,m+3*s)/fHisto->GetBinWidth(1);
  double nerr3s = signalJPsi->IntegralError(m-3*s,m+3*s,&na60Parameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  Set("NofJPsi3s",njpsi3s,nerr3s);
  
  //Computation of bin significance and signal over background
  
  Double_t bkgParameters[6];
  Double_t bkgcovarianceMatrix[6][6];
  
  for ( int ix = 0; ix < 6; ++ix )
  {
    bkgParameters[ix] = fitTotal->GetParameter(ix);
  }
  
  for ( int iy = 0; iy < 6; ++iy )
  {
    for ( int ix = 0; ix < 6; ++ix )
    {
      bkgcovarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix,iy);
    }
  }
  
  double nbck3s = bck->Integral(m-3.*s,m+3.*s)/fHisto->GetBinWidth(1);
  double nbck3sErr = bck->IntegralError(m-3.*s,m+3.*s,&bkgParameters[0],&bkgcovarianceMatrix[0][0])/fHisto->GetBinWidth(1);
  
  double sOverB3s = njpsi3s / nbck3s;
  double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));
  
  Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);
  
  double sig = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
  double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
                              TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );
  
  Set("Significance3s",sig,sigErr);

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2VWG_BKGMPTPOL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2 
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  
  Double_t kVWG = GetValue("kVWG");
  Double_t mVWG = GetValue("mVWG");
  Double_t sVWG1 = GetValue("sVWG1");
  Double_t sVWG2 = GetValue("sVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
//  TString msg;
//  
//  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
//  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
//  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
//  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
//  
//  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
//  
//  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
//    
//TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  
//  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }
  
  //_____________
  Int_t minBin = p->FindBin(fitRangeLow);
  Int_t maxBin = p->FindBin(fitRangeHigh);
  
  for ( Int_t i = minBin ; i <= maxBin ; i++ )
  {
    if ( p->GetBinEffectiveEntries(i) < 10 )
    {
      Double_t effEntries(0.),sumErr(0.);
      for ( Int_t j = i - 5 ; j < i ; j++ )
      {
        if ( j <= 0 ) continue;
//        if ( j > p->GetNbinsX() ) break;
        
        effEntries += p->GetBinEffectiveEntries(j);
        sumErr += p->GetBinEffectiveEntries(j)*p->GetBinError(j);
      }
      
      Double_t meanErr = sumErr/effEntries;
      
      if ( p->GetBinError(i) < meanErr/2. )
      {
        std::cout << "Resetting bin " << i << " error" <<std::endl;
        p->SetBinError(i,meanErr);
      }
    }
  }
  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(3.,1.,0.);
  
  bck->SetParLimits(0, 1.,8.0);
//  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
//  
//  bck->SetParameters(1,0,0,0);
  
  SetFitRejectRange(2.3,4.0);
  
  p->Fit(bck,"SER","",fitRangeLow,fitRangeHigh);
  
  SetFitRejectRange();
  
  
  
  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL2,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2VWGPOL2");
  
  fitMeanpt->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanpt->SetParName(11,"kPsiP");
  fitMeanpt->SetParName(12,"<pt>JPsi");
  fitMeanpt->SetParName(13,"<pt>BG0");
  fitMeanpt->SetParName(14,"<pt>BG1");
  fitMeanpt->SetParName(15,"<pt>BG2");
  fitMeanpt->SetParName(16,"<pt>PsiP");
  
  fitMeanpt->FixParameter(0,kVWG);
  fitMeanpt->FixParameter(1,mVWG);
  fitMeanpt->FixParameter(2,sVWG1);
  fitMeanpt->FixParameter(3,sVWG2);
  fitMeanpt->FixParameter(4,kJPsi);
  fitMeanpt->FixParameter(5,mJPsi);
  fitMeanpt->FixParameter(6,sJPsi);
  fitMeanpt->FixParameter(7,alphaLow);
  fitMeanpt->FixParameter(8,nLow);
  fitMeanpt->FixParameter(9,alphaUp);
  fitMeanpt->FixParameter(10,nUp);
  fitMeanpt->FixParameter(11,kPsiP);
  
  fitMeanpt->SetParameter(12, 3.);
  fitMeanpt->SetParLimits(12, 1.0,5.);
  
  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }

//  fitMeanpt->SetParameter(13, 3.);
//  fitMeanpt->SetParLimits(13, 0.5,10.);
//  
//  fitMeanpt->SetParameter(14, 0.2);
//  //  fitMeanpt->SetParLimits(14, 0.1,0.2);
//  
//  fitMeanpt->SetParameter(15, 0.1);
//  //  fitMeanpt->SetParLimits(15, 0.,1.);
  
  fitMeanpt->SetParameter(16, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(16, 0.,psipPtLim);
  
  
//  TProfile::Approximate();
  
  const char* fitOption = "SER"; //+";//SER
  
  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");
  
  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(16) <= fitMeanpt->GetParError(16) || fitMeanpt->GetParError(16) >= 0.75*fitMeanpt->GetParameter(16) || (fitMeanpt->GetParameter(16)/psipPtLim > 0.9) )
  {
    fitMeanpt->FixParameter(16, bck->Eval(3.68));
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }

  
//  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(12),fitMeanpt->GetParError(12));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2VWG_BKGMPTPOL2EXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  
  Double_t kVWG = GetValue("kVWG");
  Double_t mVWG = GetValue("mVWG");
  Double_t sVWG1 = GetValue("sVWG1");
  Double_t sVWG2 = GetValue("sVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  //  TString msg;
  //
  //  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  //  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  //  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  //  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  //
  //  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
  //
  //  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
  //
  //TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  
  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }
  
  //_____________
  Int_t minBin = p->FindBin(fitRangeLow);
  Int_t maxBin = p->FindBin(fitRangeHigh);
  
  for ( Int_t i = minBin ; i <= maxBin ; i++ )
  {
    if ( p->GetBinEffectiveEntries(i) < 10 )
    {
      Double_t effEntries(0.),sumErr(0.);
      for ( Int_t j = i - 5 ; j < i ; j++ )
      {
        if ( j <= 0 ) continue;
        //        if ( j > p->GetNbinsX() ) break;
        
        effEntries += p->GetBinEffectiveEntries(j);
        sumErr += p->GetBinEffectiveEntries(j)*p->GetBinError(j);
      }
      
      Double_t meanErr = sumErr/effEntries;
      
      if ( p->GetBinError(i) < meanErr/2. )
      {
        std::cout << "Resetting bin " << i << " error" <<std::endl;
        p->SetBinError(i,meanErr);
      }
    }
  }
  //_____________
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  
  bck->SetParameters(3.,-2.,0.4,-0.0);
  
  bck->SetParLimits(0, 1.,10.0);
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  //
  //  bck->SetParameters(1,0,0,0);
  
  SetFitRejectRange(2.3,4.0);
  
  p->Fit(bck,"SER","",fitRangeLow,fitRangeHigh);
  
  SetFitRejectRange();
  
  
  
  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL2EXP,fitRangeLow,fitRangeHigh,18,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2VWGPOL2EXP");
  
  fitMeanpt->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanpt->SetParName(11,"kPsiP");
  fitMeanpt->SetParName(12,"<pt>JPsi");
  fitMeanpt->SetParName(13,"<pt>BG0");
  fitMeanpt->SetParName(14,"<pt>BG1");
  fitMeanpt->SetParName(15,"<pt>BG2");
  fitMeanpt->SetParName(16,"<pt>BGEXP");
  fitMeanpt->SetParName(17,"<pt>PsiP");
  
  fitMeanpt->FixParameter(0,kVWG);
  fitMeanpt->FixParameter(1,mVWG);
  fitMeanpt->FixParameter(2,sVWG1);
  fitMeanpt->FixParameter(3,sVWG2);
  fitMeanpt->FixParameter(4,kJPsi);
  fitMeanpt->FixParameter(5,mJPsi);
  fitMeanpt->FixParameter(6,sJPsi);
  fitMeanpt->FixParameter(7,alphaLow);
  fitMeanpt->FixParameter(8,nLow);
  fitMeanpt->FixParameter(9,alphaUp);
  fitMeanpt->FixParameter(10,nUp);
  fitMeanpt->FixParameter(11,kPsiP);
  
  fitMeanpt->SetParameter(12, 3.);
  fitMeanpt->SetParLimits(12, 1.0,5.);
  
  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }
  
  //  fitMeanpt->SetParameter(13, 3.);
  //  fitMeanpt->SetParLimits(13, 0.5,10.);
  //
  //  fitMeanpt->SetParameter(14, 0.2);
  //  //  fitMeanpt->SetParLimits(14, 0.1,0.2);
  //
  //  fitMeanpt->SetParameter(15, 0.1);
  //  //  fitMeanpt->SetParLimits(15, 0.,1.);
  
  fitMeanpt->SetParameter(17, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(17, 0.,psipPtLim);
  
  
  //  TProfile::Approximate();
  
  const char* fitOption = "SER"; //+";//SER
  
  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");
  
  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(17) <= fitMeanpt->GetParError(17) || fitMeanpt->GetParError(17) >= 0.75*fitMeanpt->GetParameter(17) || (fitMeanpt->GetParameter(17)/psipPtLim > 0.9) )
  {
    fitMeanpt->FixParameter(17, bck->Eval(3.68));
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  
  if ( fitMeanpt->GetParameter(13) <= fitMeanpt->GetParError(13) || fitMeanpt->GetParError(13) >= 0.75*fitMeanpt->GetParameter(13) )
  {
    fitMeanpt->SetParameter(13, 2.);
    
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15),fitMeanpt->GetParameter(16));
  
  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(12),fitMeanpt->GetParError(12));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(17),fitMeanpt->GetParError(17));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  
  Double_t pol0 = GetValue("pol0");
  Double_t pol1 = GetValue("pol1");
  Double_t pol2 = GetValue("pol2");
  Double_t exp = GetValue("exp");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  //  TString msg;
  //
  //  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  //  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  //  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  //  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  //
  //  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
  //
  //  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
  //
  //TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  
  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }
  
  //_____________
  Int_t minBin = p->FindBin(fitRangeLow);
  Int_t maxBin = p->FindBin(fitRangeHigh);
  
  for ( Int_t i = minBin ; i <= maxBin ; i++ )
  {
    if ( p->GetBinEffectiveEntries(i) < 10 )
    {
      Double_t effEntries(0.),sumErr(0.);
      for ( Int_t j = i - 5 ; j < i ; j++ )
      {
        if ( j <= 0 ) continue;
        //        if ( j > p->GetNbinsX() ) break;
        
        effEntries += p->GetBinEffectiveEntries(j);
        sumErr += p->GetBinEffectiveEntries(j)*p->GetBinError(j);
      }
      
      Double_t meanErr = sumErr/effEntries;
      
      if ( p->GetBinError(i) < meanErr/2. )
      {
        std::cout << "Resetting bin " << i << " error" <<std::endl;
        p->SetBinError(i,meanErr);
      }
    }
  }
  //_____________
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(3.,1.,0.);
  
  bck->SetParLimits(0, 1.,8.0);
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  //
  //  bck->SetParameters(1,0,0,0);
  
  SetFitRejectRange(2.3,4.0);
  
  p->Fit(bck,"SER","",fitRangeLow,fitRangeHigh);
  
  SetFitRejectRange();
  
  
  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2POL2EXPPOL2,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2POL2EXPPOL2");
  
  fitMeanpt->SetParNames("pol0","pol1","pol2","exp","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanpt->SetParName(11,"kPsiP");
  fitMeanpt->SetParName(12,"<pt>JPsi");
  fitMeanpt->SetParName(13,"<pt>BG0");
  fitMeanpt->SetParName(14,"<pt>BG1");
  fitMeanpt->SetParName(15,"<pt>BG2");
  fitMeanpt->SetParName(16,"<pt>PsiP");
  
  fitMeanpt->FixParameter(0,pol0);
  fitMeanpt->FixParameter(1,pol1);
  fitMeanpt->FixParameter(2,pol2);
  fitMeanpt->FixParameter(3,exp);
  fitMeanpt->FixParameter(4,kJPsi);
  fitMeanpt->FixParameter(5,mJPsi);
  fitMeanpt->FixParameter(6,sJPsi);
  fitMeanpt->FixParameter(7,alphaLow);
  fitMeanpt->FixParameter(8,nLow);
  fitMeanpt->FixParameter(9,alphaUp);
  fitMeanpt->FixParameter(10,nUp);
  fitMeanpt->FixParameter(11,kPsiP);
  
  fitMeanpt->SetParameter(12, 3.);
  fitMeanpt->SetParLimits(12, 1.0,5.);
  
  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }
  
  //  fitMeanpt->SetParameter(13, 3.);
  //  fitMeanpt->SetParLimits(13, 0.5,10.);
  //
  //  fitMeanpt->SetParameter(14, 0.2);
  //  //  fitMeanpt->SetParLimits(14, 0.1,0.2);
  //
  //  fitMeanpt->SetParameter(15, 0.1);
  //  //  fitMeanpt->SetParLimits(15, 0.,1.);
  
  fitMeanpt->SetParameter(16, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(16, 0.,psipPtLim);
  
  
  //  TProfile::Approximate();
  
  const char* fitOption = "SER"; //+";//SER
  
  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");
  
  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(16) <= fitMeanpt->GetParError(16) || fitMeanpt->GetParError(16) >= 0.75*fitMeanpt->GetParameter(16) || (fitMeanpt->GetParameter(16)/psipPtLim > 0.9) )
  {
    fitMeanpt->FixParameter(16, bck->Eval(3.68));
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15));
  
  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(12),fitMeanpt->GetParError(12));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2EXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  
  Double_t pol0 = GetValue("pol0");
  Double_t pol1 = GetValue("pol1");
  Double_t pol2 = GetValue("pol2");
  Double_t exp = GetValue("exp");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  //  TString msg;
  //
  //  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  //  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  //  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  //  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  //
  //  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
  //
  //  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
  //
  //TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  
  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }
  
  //_____________
  Int_t minBin = p->FindBin(fitRangeLow);
  Int_t maxBin = p->FindBin(fitRangeHigh);
  
  for ( Int_t i = minBin ; i <= maxBin ; i++ )
  {
    if ( p->GetBinEffectiveEntries(i) < 10 )
    {
      Double_t effEntries(0.),sumErr(0.);
      for ( Int_t j = i - 5 ; j < i ; j++ )
      {
        if ( j <= 0 ) continue;
        //        if ( j > p->GetNbinsX() ) break;
        
        effEntries += p->GetBinEffectiveEntries(j);
        sumErr += p->GetBinEffectiveEntries(j)*p->GetBinError(j);
      }
      
      Double_t meanErr = sumErr/effEntries;
      
      if ( p->GetBinError(i) < meanErr/2. )
      {
        std::cout << "Resetting bin " << i << " error" <<std::endl;
        p->SetBinError(i,meanErr);
      }
    }
  }
  //_____________
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  
  bck->SetParameters(3.,-1.,0.4,-0.1);
  
  bck->SetParLimits(0, 1.,10.0);
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  //
  //  bck->SetParameters(1,0,0,0);
  
  SetFitRejectRange(2.3,4.0);
  
  p->Fit(bck,"SER","",fitRangeLow,fitRangeHigh);
  
  SetFitRejectRange();
  
  
  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2POL2EXPPOL2EXP,fitRangeLow,fitRangeHigh,18,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2POL2EXPPOL2EXP");
  
  fitMeanpt->SetParNames("pol0","pol1","pol2","exp","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanpt->SetParName(11,"kPsiP");
  fitMeanpt->SetParName(12,"<pt>JPsi");
  fitMeanpt->SetParName(13,"<pt>BG0");
  fitMeanpt->SetParName(14,"<pt>BG1");
  fitMeanpt->SetParName(15,"<pt>BG2");
  fitMeanpt->SetParName(16,"<pt>BGEXP");
  fitMeanpt->SetParName(17,"<pt>PsiP");
  
  fitMeanpt->FixParameter(0,pol0);
  fitMeanpt->FixParameter(1,pol1);
  fitMeanpt->FixParameter(2,pol2);
  fitMeanpt->FixParameter(3,exp);
  fitMeanpt->FixParameter(4,kJPsi);
  fitMeanpt->FixParameter(5,mJPsi);
  fitMeanpt->FixParameter(6,sJPsi);
  fitMeanpt->FixParameter(7,alphaLow);
  fitMeanpt->FixParameter(8,nLow);
  fitMeanpt->FixParameter(9,alphaUp);
  fitMeanpt->FixParameter(10,nUp);
  fitMeanpt->FixParameter(11,kPsiP);
  
  fitMeanpt->SetParameter(12, 3.);
  fitMeanpt->SetParLimits(12, 1.0,5.);
  
  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }
  
  //  fitMeanpt->SetParameter(13, 3.);
  //  fitMeanpt->SetParLimits(13, 0.5,10.);
  //
  //  fitMeanpt->SetParameter(14, 0.2);
  //  //  fitMeanpt->SetParLimits(14, 0.1,0.2);
  //
  //  fitMeanpt->SetParameter(15, 0.1);
  //  //  fitMeanpt->SetParLimits(15, 0.,1.);
  
  fitMeanpt->SetParameter(17, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(17, 0.,psipPtLim);
  
  
  //  TProfile::Approximate();
  
  const char* fitOption = "SER"; //+";//SER
  
  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");
  
  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(17) <= fitMeanpt->GetParError(17) || fitMeanpt->GetParError(17) >= 0.75*fitMeanpt->GetParameter(17) || (fitMeanpt->GetParameter(17)/psipPtLim > 0.9) )
  {
    fitMeanpt->FixParameter(17, bck->Eval(3.68));
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( fitMeanpt->GetParameter(13) <= fitMeanpt->GetParError(13) || fitMeanpt->GetParError(13) >= 0.75*fitMeanpt->GetParameter(13) )
  {
    fitMeanpt->SetParameter(13, 2.);
    
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15),fitMeanpt->GetParameter(16));
  
  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(12),fitMeanpt->GetParError(12));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(17),fitMeanpt->GetParError(17));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");
  
  Double_t kVWG = GetValue("kVWG");
  Double_t mVWG = GetValue("mVWG");
  Double_t sVWG1 = GetValue("sVWG1");
  Double_t sVWG2 = GetValue("sVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  //  TString msg;
  //
  //  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  //  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  //  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  //  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  //
  //  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
  //
  //  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
  //
  //TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  
  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }
  
  //_____________
  Int_t minBin = p->FindBin(fitRangeLow);
  Int_t maxBin = p->FindBin(fitRangeHigh);
  
  for ( Int_t i = minBin ; i <= maxBin ; i++ )
  {
    if ( p->GetBinEffectiveEntries(i) < 10 )
    {
      Double_t effEntries(0.),sumErr(0.);
      for ( Int_t j = i - 5 ; j < i ; j++ )
      {
        if ( j <= 0 ) continue;
        //        if ( j > p->GetNbinsX() ) break;
        
        effEntries += p->GetBinEffectiveEntries(j);
        sumErr += p->GetBinEffectiveEntries(j)*p->GetBinError(j);
      }
      
      Double_t meanErr = sumErr/effEntries;
      
      if ( p->GetBinError(i) < meanErr/2. )
      {
        std::cout << "Resetting bin " << i << " error" <<std::endl;
        p->SetBinError(i,meanErr);
      }
    }
  }
  //_____________
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(3.,1.,0.);
  
  bck->SetParLimits(0, 1.,8.0);
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  //
  //  bck->SetParameters(1,0,0,0);
  
  SetFitRejectRange(2.3,4.0);
  
  p->Fit(bck,"SER","",fitRangeLow,fitRangeHigh);
  
  SetFitRejectRange();
  
  
  
  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2NA60NEWVWGPOL2,fitRangeLow,fitRangeHigh,21,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2NA60NEWVWGPOL2");

  fitMeanpt->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi");
  fitMeanpt->SetParName(11,"p2RJPsi");
  fitMeanpt->SetParName(12,"p3RJPsi");
  fitMeanpt->SetParName(13,"aLJPsi");
  fitMeanpt->SetParName(14,"aRJPsi");
  fitMeanpt->SetParName(15,"kPsiP");
  fitMeanpt->SetParName(16,"<pt>JPsi");
  fitMeanpt->SetParName(17,"<pt>BG0");
  fitMeanpt->SetParName(18,"<pt>BG1");
  fitMeanpt->SetParName(19,"<pt>BG2");
  fitMeanpt->SetParName(20,"<pt>PsiP");
  
  fitMeanpt->FixParameter(0,kVWG);
  fitMeanpt->FixParameter(1,mVWG);
  fitMeanpt->FixParameter(2,sVWG1);
  fitMeanpt->FixParameter(3,sVWG2);
  fitMeanpt->FixParameter(4,kJPsi);
  fitMeanpt->FixParameter(5,mJPsi);
  fitMeanpt->FixParameter(6,sJPsi);
  fitMeanpt->FixParameter(7,p1Left);
  fitMeanpt->FixParameter(8,p2Left);
  fitMeanpt->FixParameter(9,p3Left);
  fitMeanpt->FixParameter(10,p1Right);
  fitMeanpt->FixParameter(11,p2Right);
  fitMeanpt->FixParameter(12,p3Right);
  fitMeanpt->FixParameter(13,alphaLeft);
  fitMeanpt->FixParameter(14,alphaRight);
  fitMeanpt->FixParameter(15,kPsiP);
  
  fitMeanpt->SetParameter(16, 3.);
  fitMeanpt->SetParLimits(16, 1.0,5.);

  
  for ( Int_t i = 0; i < 3; ++i )
  {
//    TString name(GetName());
//    if (name.Contains("NA60NEWVWG_BKGMPTPOL2_2.0_5.0(Sig:2.2_4.7_SP0.9)") || name.Contains("NA60NEWVWG_BKGMPTPOL2_2.2_4.7(Sig:2.2_4.7_SP0.9)") )
//    {
//      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)/2.);
//      fitMeanpt->SetParameter(16, 4.);
////      if (fitMeanpt->GetParameter(3) < 2.e-06) fitMeanpt->FixParameter(3,3e-6);
//
//    }
//    else
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i));
  }
  
  //  fitMeanpt->SetParameter(13, 3.);
  //  fitMeanpt->SetParLimits(13, 0.5,10.);
  //
  //  fitMeanpt->SetParameter(14, 0.2);
  //  //  fitMeanpt->SetParLimits(14, 0.1,0.2);
  //
  //  fitMeanpt->SetParameter(15, 0.1);
  //  //  fitMeanpt->SetParLimits(15, 0.,1.);
  
  fitMeanpt->SetParameter(20, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(20, 0.,psipPtLim);
  
  
  //  TProfile::Approximate();
  
  const char* fitOption = "SER"; //+";//SER
  
  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");
  
  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(20) <= fitMeanpt->GetParError(20) || fitMeanpt->GetParError(20) >= 0.75*fitMeanpt->GetParameter(20) || (fitMeanpt->GetParameter(20)/psipPtLim > 0.9) )
  {
    fitMeanpt->FixParameter(20, bck->Eval(3.68));
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(fitMeanpt->GetParameter(17),fitMeanpt->GetParameter(18),fitMeanpt->GetParameter(19));
  
  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(20),fitMeanpt->GetParError(20));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2EXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");
  
  Double_t kVWG = GetValue("kVWG");
  Double_t mVWG = GetValue("mVWG");
  Double_t sVWG1 = GetValue("sVWG1");
  Double_t sVWG2 = GetValue("sVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  //  TString msg;
  //
  //  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  //  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  //  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  //  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  //
  //  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
  //
  //  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
  //
  //TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  
  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }
  
  //_____________
  Int_t minBin = p->FindBin(fitRangeLow);
  Int_t maxBin = p->FindBin(fitRangeHigh);
  
  for ( Int_t i = minBin ; i <= maxBin ; i++ )
  {
    if ( p->GetBinEffectiveEntries(i) < 10 )
    {
      Double_t effEntries(0.),sumErr(0.);
      for ( Int_t j = i - 5 ; j < i ; j++ )
      {
        if ( j <= 0 ) continue;
        //        if ( j > p->GetNbinsX() ) break;
        
        effEntries += p->GetBinEffectiveEntries(j);
        sumErr += p->GetBinEffectiveEntries(j)*p->GetBinError(j);
      }
      
      Double_t meanErr = sumErr/effEntries;
      
      if ( p->GetBinError(i) < meanErr/2. )
      {
        std::cout << "Resetting bin " << i << " error" <<std::endl;
        p->SetBinError(i,meanErr);
      }
    }
  }
  //_____________
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  
  bck->SetParameters(2.5,-1.,0.2,-0.05);
  
  bck->SetParLimits(0, 1.,10.0);
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  //
  //  bck->SetParameters(1,0,0,0);
  
  SetFitRejectRange(2.3,4.0);
  
  p->Fit(bck,"SER","",fitRangeLow,fitRangeHigh);
  
  SetFitRejectRange();
  
  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2NA60NEWVWGPOL2EXP,fitRangeLow,fitRangeHigh,22,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2NA60NEWVWGPOL2EXP");
  
  fitMeanpt->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi");
  fitMeanpt->SetParName(11,"p2RJPsi");
  fitMeanpt->SetParName(12,"p3RJPsi");
  fitMeanpt->SetParName(13,"aLJPsi");
  fitMeanpt->SetParName(14,"aRJPsi");
  fitMeanpt->SetParName(15,"kPsiP");
  fitMeanpt->SetParName(16,"<pt>JPsi");
  fitMeanpt->SetParName(17,"<pt>BG0");
  fitMeanpt->SetParName(18,"<pt>BG1");
  fitMeanpt->SetParName(19,"<pt>BG2");
  fitMeanpt->SetParName(20,"<pt>BGEXP");
  fitMeanpt->SetParName(21,"<pt>PsiP");
  
  fitMeanpt->FixParameter(0,kVWG);
  fitMeanpt->FixParameter(1,mVWG);
  fitMeanpt->FixParameter(2,sVWG1);
  fitMeanpt->FixParameter(3,sVWG2);
  fitMeanpt->FixParameter(4,kJPsi);
  fitMeanpt->FixParameter(5,mJPsi);
  fitMeanpt->FixParameter(6,sJPsi);
  fitMeanpt->FixParameter(7,p1Left);
  fitMeanpt->FixParameter(8,p2Left);
  fitMeanpt->FixParameter(9,p3Left);
  fitMeanpt->FixParameter(10,p1Right);
  fitMeanpt->FixParameter(11,p2Right);
  fitMeanpt->FixParameter(12,p3Right);
  fitMeanpt->FixParameter(13,alphaLeft);
  fitMeanpt->FixParameter(14,alphaRight);
  fitMeanpt->FixParameter(15,kPsiP);
  
  fitMeanpt->SetParameter(16, 3.);
  fitMeanpt->SetParLimits(16, 1.0,5.);
  
  for ( Int_t i = 0; i < 4; ++i )
  {
//    TString name(GetName());
//    if ( name.Contains("NA60NEWVWG_BKGMPTPOL2EXP_2.0_5.0(Sig:2.2_4.7_SP0.9)") ||name.Contains("NA60NEWVWG_BKGMPTPOL2EXP_2.2_4.7(Sig:2.2_4.7_SP0.9)") )
//    {
//      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)/2.);
//      fitMeanpt->SetParameter(16, 4.);
////      if (fitMeanpt->GetParameter(3) < 2.e-06) fitMeanpt->FixParameter(3,3e-4);
//      
//    }
//    else
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i));
//    fitMeanpt->SetParameter(i + 17, bck->GetParameter(i));
  }
  
  fitMeanpt->SetParameter(21, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(21, 0.,psipPtLim);
  
  
  //  TProfile::Approximate();
  
  const char* fitOption = "SER"; //+";//SER
  
  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");
  
  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(21) <= fitMeanpt->GetParError(21) || fitMeanpt->GetParError(21) >= 0.75*fitMeanpt->GetParameter(21) || (fitMeanpt->GetParameter(21)/psipPtLim > 0.9) )
  {
    fitMeanpt->FixParameter(21, bck->Eval(3.68));
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( fitMeanpt->GetParameter(17) <= fitMeanpt->GetParError(17) || fitMeanpt->GetParError(17) >= 0.75*fitMeanpt->GetParameter(13) )
  {
    fitMeanpt->SetParameter(17, 2.);
    
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }

  
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(fitMeanpt->GetParameter(17),fitMeanpt->GetParameter(18),fitMeanpt->GetParameter(19),fitMeanpt->GetParameter(20));
  
  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(21),fitMeanpt->GetParError(21));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");
  
  Double_t pol0 = GetValue("pol0");
  Double_t pol1 = GetValue("pol1");
  Double_t pol2 = GetValue("pol2");
  Double_t exp = GetValue("exp");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  //  TString msg;
  //
  //  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  //  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  //  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  //  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  //
  //  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
  //
  //  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
  //
  //TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  
  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }
  
  //_____________
  Int_t minBin = p->FindBin(fitRangeLow);
  Int_t maxBin = p->FindBin(fitRangeHigh);
  
  for ( Int_t i = minBin ; i <= maxBin ; i++ )
  {
    if ( p->GetBinEffectiveEntries(i) < 10 )
    {
      Double_t effEntries(0.),sumErr(0.);
      for ( Int_t j = i - 5 ; j < i ; j++ )
      {
        if ( j <= 0 ) continue;
        //        if ( j > p->GetNbinsX() ) break;
        
        effEntries += p->GetBinEffectiveEntries(j);
        sumErr += p->GetBinEffectiveEntries(j)*p->GetBinError(j);
      }
      
      Double_t meanErr = sumErr/effEntries;
      
      if ( p->GetBinError(i) < meanErr/2. )
      {
        std::cout << "Resetting bin " << i << " error" <<std::endl;
        p->SetBinError(i,meanErr);
      }
    }
  }
  //_____________
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(3.,1.,0.);
  
  bck->SetParLimits(0, 1.,8.0);
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  //
  //  bck->SetParameters(1,0,0,0);
  
  SetFitRejectRange(2.3,4.0);
  
  p->Fit(bck,"SER","",fitRangeLow,fitRangeHigh);
  
  SetFitRejectRange();

  
  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2NA60NEWPOL2EXPPOL2,fitRangeLow,fitRangeHigh,21,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2NA60NEWPOL2EXPPOL2");
  
  fitMeanpt->SetParNames("pol0","pol1","pol2","exp","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi");
  fitMeanpt->SetParName(11,"p2RJPsi");
  fitMeanpt->SetParName(12,"p3RJPsi");
  fitMeanpt->SetParName(13,"aLJPsi");
  fitMeanpt->SetParName(14,"aRJPsi");
  fitMeanpt->SetParName(15,"kPsiP");
  fitMeanpt->SetParName(16,"<pt>JPsi");
  fitMeanpt->SetParName(17,"<pt>BG0");
  fitMeanpt->SetParName(18,"<pt>BG1");
  fitMeanpt->SetParName(19,"<pt>BG2");
  fitMeanpt->SetParName(20,"<pt>PsiP");
  
  fitMeanpt->FixParameter(0,pol0);
  fitMeanpt->FixParameter(1,pol1);
  fitMeanpt->FixParameter(2,pol2);
  fitMeanpt->FixParameter(3,exp);
  fitMeanpt->FixParameter(4,kJPsi);
  fitMeanpt->FixParameter(5,mJPsi);
  fitMeanpt->FixParameter(6,sJPsi);
  fitMeanpt->FixParameter(7,p1Left);
  fitMeanpt->FixParameter(8,p2Left);
  fitMeanpt->FixParameter(9,p3Left);
  fitMeanpt->FixParameter(10,p1Right);
  fitMeanpt->FixParameter(11,p2Right);
  fitMeanpt->FixParameter(12,p3Right);
  fitMeanpt->FixParameter(13,alphaLeft);
  fitMeanpt->FixParameter(14,alphaRight);
  
  fitMeanpt->FixParameter(15,kPsiP);
  
  fitMeanpt->SetParameter(16, 3.);
  fitMeanpt->SetParLimits(16, 1.0,5.);
  
  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanpt->SetParameter(i + 17, bck->GetParameter(i));
  }
  
  //  fitMeanpt->SetParameter(13, 3.);
  //  fitMeanpt->SetParLimits(13, 0.5,10.);
  //
  //  fitMeanpt->SetParameter(14, 0.2);
  //  //  fitMeanpt->SetParLimits(14, 0.1,0.2);
  //
  //  fitMeanpt->SetParameter(15, 0.1);
  //  //  fitMeanpt->SetParLimits(15, 0.,1.);
  
  fitMeanpt->SetParameter(20, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(20, 0.,psipPtLim);
  
  
  //  TProfile::Approximate();
  
  const char* fitOption = "SER"; //+";//SER
  
  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");
  
  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(20) <= fitMeanpt->GetParError(20) || fitMeanpt->GetParError(20) >= 0.75*fitMeanpt->GetParameter(20) || (fitMeanpt->GetParameter(20)/psipPtLim > 0.9) )
  {
    fitMeanpt->FixParameter(20, bck->Eval(3.68));
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(fitMeanpt->GetParameter(17),fitMeanpt->GetParameter(18),fitMeanpt->GetParameter(19));
  
  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(20),fitMeanpt->GetParError(20));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2EXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");
  
  Double_t pol0 = GetValue("pol0");
  Double_t pol1 = GetValue("pol1");
  Double_t pol2 = GetValue("pol2");
  Double_t exp = GetValue("exp");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  //  TString msg;
  //
  //  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  //  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  //  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  //  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  //
  //  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
  //
  //  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
  //
  //TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  
  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }
  
  //_____________
  Int_t minBin = p->FindBin(fitRangeLow);
  Int_t maxBin = p->FindBin(fitRangeHigh);
  
  for ( Int_t i = minBin ; i <= maxBin ; i++ )
  {
    if ( p->GetBinEffectiveEntries(i) < 10 )
    {
      Double_t effEntries(0.),sumErr(0.);
      for ( Int_t j = i - 5 ; j < i ; j++ )
      {
        if ( j <= 0 ) continue;
        //        if ( j > p->GetNbinsX() ) break;
        
        effEntries += p->GetBinEffectiveEntries(j);
        sumErr += p->GetBinEffectiveEntries(j)*p->GetBinError(j);
      }
      
      Double_t meanErr = sumErr/effEntries;
      
      if ( p->GetBinError(i) < meanErr/2. )
      {
        std::cout << "Resetting bin " << i << " error" <<std::endl;
        p->SetBinError(i,meanErr);
      }
    }
  }
  //_____________
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  
  bck->SetParameters(2.8,-0.5,0.2,0.05);
  
  bck->SetParLimits(0, 1.,10.0);
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  //
  //  bck->SetParameters(1,0,0,0);
  
  SetFitRejectRange(2.3,4.0);
  
  p->Fit(bck,"SER","",fitRangeLow,fitRangeHigh);
  
  SetFitRejectRange();  
  
  
  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2NA60NEWPOL2EXPPOL2EXP,fitRangeLow,fitRangeHigh,22,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2NA60NEWPOL2EXPPOL2EXP");
  
  fitMeanpt->SetParNames("pol0","pol1","pol2","exp","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi");
  fitMeanpt->SetParName(11,"p2RJPsi");
  fitMeanpt->SetParName(12,"p3RJPsi");
  fitMeanpt->SetParName(13,"aLJPsi");
  fitMeanpt->SetParName(14,"aRJPsi");
  fitMeanpt->SetParName(15,"kPsiP");
  fitMeanpt->SetParName(16,"<pt>JPsi");
  fitMeanpt->SetParName(17,"<pt>BG0");
  fitMeanpt->SetParName(18,"<pt>BG1");
  fitMeanpt->SetParName(19,"<pt>BG2");
  fitMeanpt->SetParName(20,"<pt>BGEXP");
  fitMeanpt->SetParName(21,"<pt>PsiP");
  
  fitMeanpt->FixParameter(0,pol0);
  fitMeanpt->FixParameter(1,pol1);
  fitMeanpt->FixParameter(2,pol2);
  fitMeanpt->FixParameter(3,exp);
  fitMeanpt->FixParameter(4,kJPsi);
  fitMeanpt->FixParameter(5,mJPsi);
  fitMeanpt->FixParameter(6,sJPsi);
  fitMeanpt->FixParameter(7,p1Left);
  fitMeanpt->FixParameter(8,p2Left);
  fitMeanpt->FixParameter(9,p3Left);
  fitMeanpt->FixParameter(10,p1Right);
  fitMeanpt->FixParameter(11,p2Right);
  fitMeanpt->FixParameter(12,p3Right);
  fitMeanpt->FixParameter(13,alphaLeft);
  fitMeanpt->FixParameter(14,alphaRight);
  
  fitMeanpt->FixParameter(15,kPsiP);
  
  fitMeanpt->SetParameter(16, 3.);
  fitMeanpt->SetParLimits(16, 1.0,5.);
  
  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanpt->SetParameter(i + 17, bck->GetParameter(i));
  }
  
  //  fitMeanpt->SetParameter(13, 3.);
  //  fitMeanpt->SetParLimits(13, 0.5,10.);
  //
  //  fitMeanpt->SetParameter(14, 0.2);
  //  //  fitMeanpt->SetParLimits(14, 0.1,0.2);
  //
  //  fitMeanpt->SetParameter(15, 0.1);
  //  //  fitMeanpt->SetParLimits(15, 0.,1.);
  
  fitMeanpt->SetParameter(21, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(21, 0.,psipPtLim);

  
  
  //  TProfile::Approximate();
  
  const char* fitOption = "SER"; //+";//SER
  
  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");
  
  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(21) <= fitMeanpt->GetParError(21) || fitMeanpt->GetParError(21) >= 0.75*fitMeanpt->GetParameter(21) || (fitMeanpt->GetParameter(21)/psipPtLim > 0.9) )
  {
    fitMeanpt->FixParameter(21, bck->Eval(3.68));
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( fitMeanpt->GetParameter(17) <= fitMeanpt->GetParError(17) || fitMeanpt->GetParError(17) >= 0.75*fitMeanpt->GetParameter(13) )
  {
    fitMeanpt->SetParameter(17, 2.);
    
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(fitMeanpt->GetParameter(17),fitMeanpt->GetParameter(18),fitMeanpt->GetParameter(19),fitMeanpt->GetParameter(20));
  
  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(21),fitMeanpt->GetParError(21));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2VWG_BKGMPTPOL3()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  
  Double_t kVWG = GetValue("kVWG");
  Double_t mVWG = GetValue("mVWG");
  Double_t sVWG1 = GetValue("sVWG1");
  Double_t sVWG2 = GetValue("sVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  //  TString msg;
  //
  //  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  //  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  //  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  //  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  //
  //  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
  //
  //  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
  //
  //TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  
  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }
  
  //_____________
  Int_t minBin = p->FindBin(fitRangeLow);
  Int_t maxBin = p->FindBin(fitRangeHigh);
  
  for ( Int_t i = minBin ; i <= maxBin ; i++ )
  {
    if ( p->GetBinEffectiveEntries(i) < 10 )
    {
      Double_t effEntries(0.),sumErr(0.);
      for ( Int_t j = i - 5 ; j < i ; j++ )
      {
        if ( j <= 0 ) continue;
        //        if ( j > p->GetNbinsX() ) break;
        
        effEntries += p->GetBinEffectiveEntries(j);
        sumErr += p->GetBinEffectiveEntries(j)*p->GetBinError(j);
      }
      
      Double_t meanErr = sumErr/effEntries;
      
      if ( p->GetBinError(i) < meanErr/2. )
      {
        std::cout << "Resetting bin " << i << " error" <<std::endl;
        p->SetBinError(i,meanErr);
      }
    }
  }
  //_____________
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol3,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol3");
  
  bck->SetParameters(3.,1.,-0.4,0.05);
  
  bck->SetParLimits(0, 1.,8.0);
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  //
  //  bck->SetParameters(1,0,0,0);
  
  SetFitRejectRange(2.3,4.0);
  
  p->Fit(bck,"SER","",fitRangeLow,fitRangeHigh);
  
  SetFitRejectRange();
  
  
  
  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL3,fitRangeLow,fitRangeHigh,18,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2VWGPOL3");
  
  fitMeanpt->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanpt->SetParName(11,"kPsiP");
  fitMeanpt->SetParName(12,"<pt>JPsi");
  fitMeanpt->SetParName(13,"<pt>BG0");
  fitMeanpt->SetParName(14,"<pt>BG1");
  fitMeanpt->SetParName(15,"<pt>BG2");
  fitMeanpt->SetParName(16,"<pt>BG3");
  fitMeanpt->SetParName(17,"<pt>PsiP");
  
  fitMeanpt->FixParameter(0,kVWG);
  fitMeanpt->FixParameter(1,mVWG);
  fitMeanpt->FixParameter(2,sVWG1);
  fitMeanpt->FixParameter(3,sVWG2);
  fitMeanpt->FixParameter(4,kJPsi);
  fitMeanpt->FixParameter(5,mJPsi);
  fitMeanpt->FixParameter(6,sJPsi);
  fitMeanpt->FixParameter(7,alphaLow);
  fitMeanpt->FixParameter(8,nLow);
  fitMeanpt->FixParameter(9,alphaUp);
  fitMeanpt->FixParameter(10,nUp);
  fitMeanpt->FixParameter(11,kPsiP);
  
  fitMeanpt->SetParameter(12, 3.);
  fitMeanpt->SetParLimits(12, 1.0,5.);
  
  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }
  
  //  fitMeanpt->SetParameter(13, 3.);
  //  fitMeanpt->SetParLimits(13, 0.5,10.);
  //
  //  fitMeanpt->SetParameter(14, 0.2);
  //  //  fitMeanpt->SetParLimits(14, 0.1,0.2);
  //
  //  fitMeanpt->SetParameter(15, 0.1);
  //  //  fitMeanpt->SetParLimits(15, 0.,1.);
  
  fitMeanpt->SetParameter(17, 3.);
  fitMeanpt->SetParLimits(17, 0.,10.);
  
  
  //  TProfile::Approximate();
  
  const char* fitOption = "SER"; //+";//SER
  
  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");
  
  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(17) <= fitMeanpt->GetParError(17) || fitMeanpt->GetParError(17) >= 0.75*fitMeanpt->GetParameter(17) )
  {
    fitMeanpt->FixParameter(17, 0.);
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15),fitMeanpt->GetParameter(16));
  
  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(12),fitMeanpt->GetParError(12));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(17),fitMeanpt->GetParError(17));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2VWG_BKGMPTLIN()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  
  Double_t kVWG = GetValue("kVWG");
  Double_t mVWG = GetValue("mVWG");
  Double_t sVWG1 = GetValue("sVWG1");
  Double_t sVWG2 = GetValue("sVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  //  TString msg;
  //
  //  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  //  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  //  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  //  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  //
  //  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
  //
  //  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
  //
  //TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  
  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundLin,fitRangeLow,fitRangeHigh,2,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundLin");
  
  bck->SetParameters(3.,0.);
  
  bck->SetParLimits(0, 2.0,4.0);
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  //
  //  bck->SetParameters(1,0,0,0);
  
  SetFitRejectRange(2.0,4.0);
  
  fHisto->Fit(bck,"SER","",fitRangeLow,fitRangeHigh);
  
  SetFitRejectRange();
  
  
  
  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2Lin,fitRangeLow,fitRangeHigh,16,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2Lin");
  
  fitMeanpt->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanpt->SetParName(11,"kPsiP");
  fitMeanpt->SetParName(12,"<pt>JPsi");
  fitMeanpt->SetParName(13,"<pt>BG0");
  fitMeanpt->SetParName(14,"<pt>BG1");
  fitMeanpt->SetParName(15,"<pt>PsiP");
  
  fitMeanpt->FixParameter(0,kVWG);
  fitMeanpt->FixParameter(1,mVWG);
  fitMeanpt->FixParameter(2,sVWG1);
  fitMeanpt->FixParameter(3,sVWG2);
  fitMeanpt->FixParameter(4,kJPsi);
  fitMeanpt->FixParameter(5,mJPsi);
  fitMeanpt->FixParameter(6,sJPsi);
  fitMeanpt->FixParameter(7,alphaLow);
  fitMeanpt->FixParameter(8,nLow);
  fitMeanpt->FixParameter(9,alphaUp);
  fitMeanpt->FixParameter(10,nUp);
  fitMeanpt->FixParameter(11,kPsiP);
  
  fitMeanpt->SetParameter(12, 3.);
  fitMeanpt->SetParLimits(12, 2.0,4.);
  
  for ( Int_t i = 0; i < 2; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }
  
  //  fitMeanpt->SetParameter(13, 3.);
  //  fitMeanpt->SetParLimits(13, 0.5,10.);
  //
  //  fitMeanpt->SetParameter(14, 0.2);
  //  //  fitMeanpt->SetParLimits(14, 0.1,0.2);
  //
  //  fitMeanpt->SetParameter(15, 0.1);
  //  //  fitMeanpt->SetParLimits(15, 0.,1.);
  
  fitMeanpt->SetParameter(15, 3.);
  fitMeanpt->SetParLimits(15, 0.5,6.);
  
  
  //  TProfile::Approximate();
  
  const char* fitOption = "SER"; //+";
  
  TFitResultPtr fitResult = fHisto->Fit(fitMeanpt,fitOption,"");
  
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14));
  
  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(12),fitMeanpt->GetParError(12));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(15),fitMeanpt->GetParError(15));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2VWGINDEPTAILS_BKGMPTPOL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  
  Double_t alphaLowP = GetValue("alPsiP");
  Double_t nLowP = GetValue("nlPsiP");
  Double_t alphaUpP = GetValue("auPsiP");
  Double_t nUpP = GetValue("nuPsiP");
  
  Double_t kVWG = GetValue("kVWG");
  Double_t mVWG = GetValue("mVWG");
  Double_t sVWG1 = GetValue("sVWG1");
  Double_t sVWG2 = GetValue("sVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
//  Double_t mPsiP = GetValue("mPsiP");
//  Double_t sPsiP = GetValue("sPsiP");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  //  TString msg;
  //
  //  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  //  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  //  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  //  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  //
  //  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
  //
  //  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
  //
  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  
  
  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  
  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL2INDEPTAILS,fitRangeLow,fitRangeHigh,21,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2VWGPOL2INDEPTAILS");
  
  fitMeanpt->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanpt->SetParName(11,"kPsiP");
//  fitMeanpt->SetParName(12,"mPsiP");
//  fitMeanpt->SetParName(13,"sPsiP");
  fitMeanpt->SetParName(12,"alPsiP");
  fitMeanpt->SetParName(13,"nlPsiP");
  fitMeanpt->SetParName(14,"auPsiP");
  fitMeanpt->SetParName(15,"nuPsiP");
  fitMeanpt->SetParName(16,"<pt>JPsi");//12
  fitMeanpt->SetParName(17,"<pt>BG0");//13
  fitMeanpt->SetParName(18,"<pt>BG1");//14
  fitMeanpt->SetParName(19,"<pt>BG2");//15
  fitMeanpt->SetParName(20,"<pt>PsiP");//16
  
  fitMeanpt->FixParameter(0,kVWG);
  fitMeanpt->FixParameter(1,mVWG);
  fitMeanpt->FixParameter(2,sVWG1);
  fitMeanpt->FixParameter(3,sVWG2);
  fitMeanpt->FixParameter(4,kJPsi);
  fitMeanpt->FixParameter(5,mJPsi);
  fitMeanpt->FixParameter(6,sJPsi);
  fitMeanpt->FixParameter(7,alphaLow);
  fitMeanpt->FixParameter(8,nLow);
  fitMeanpt->FixParameter(9,alphaUp);
  fitMeanpt->FixParameter(10,nUp);
  fitMeanpt->FixParameter(11,kPsiP);
//  fitMeanpt->FixParameter(12,mPsiP);
//  fitMeanpt->FixParameter(13,sPsiP);
  fitMeanpt->FixParameter(12,alphaLowP);
  fitMeanpt->FixParameter(13,nLowP);
  fitMeanpt->FixParameter(14,alphaUpP);
  fitMeanpt->FixParameter(15,nUpP);
  
  fitMeanpt->SetParameter(16, 3.);
  fitMeanpt->SetParLimits(16, 2.,3.);
  
  fitMeanpt->SetParameter(17, 1.);
  //  fitMeanpt->SetParLimits(13, 0.01,1.5);
  
  fitMeanpt->SetParameter(18, 0.2);
  //  fitMeanpt->SetParLimits(14, 0.1,0.2);
  
  fitMeanpt->SetParameter(19, 0.1);
  //  fitMeanpt->SetParLimits(15, 0.,1.);
  
  fitMeanpt->SetParameter(20, 3.);
  fitMeanpt->SetParLimits(20, 2.5,8.);
  
  
  
  const char* fitOption = "SER"; //+";
  
  TFitResultPtr fitResult = fHisto->Fit(fitMeanpt,fitOption,"");
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  bck->SetParameters(fitMeanpt->GetParameter(17),fitMeanpt->GetParameter(18),fitMeanpt->GetParameter(19));
  
  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);
  
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(20),fitMeanpt->GetParError(20));
  
}

////_____________________________________________________________________________
//void AliAnalysisMuMuJpsiResult::FitMPT2CB2VWG_BKGMPTPOL2EXP()
//{
//  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
//  fHisto->GetListOfFunctions()->Delete();
//  
//  Double_t alphaLow = GetValue("alJPsi");
//  Double_t nLow = GetValue("nlJPsi");
//  Double_t alphaUp = GetValue("auJPsi");
//  Double_t nUp = GetValue("nuJPsi");
//  
//  Double_t kVWG = GetValue("kVWG");
//  Double_t mVWG = GetValue("mVWG");
//  Double_t sVWG1 = GetValue("sVWG1");
//  Double_t sVWG2 = GetValue("sVWG2");
//  Double_t kJPsi = GetValue("kJPsi");
//  Double_t kPsiP = GetValue("kPsiP");
//  Double_t mJPsi = GetValue("mJPsi");
//  Double_t sJPsi = GetValue("sJPsi");
//  Double_t NofJPsi = GetValue("NofJPsi");
//  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
//  
//  Double_t fitRangeLow = GetValue(kFitRangeLow);
//  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
//  
//  //  TString msg;
//  //
//  //  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
//  //  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
//  //  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
//  //  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
//  //
//  //  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
//  //
//  //  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
//  //
//  //TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
//  
//  
//  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
//  
//  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL2EXP,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2VWGPOL2EXP");
//  
//  fitMeanpt->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
//  fitMeanpt->SetParName(11,"kPsiP");
//  fitMeanpt->SetParName(12,"<pt>JPsi");
//  fitMeanpt->SetParName(13,"<pt>BG0");
//  fitMeanpt->SetParName(14,"<pt>BG1");
//  fitMeanpt->SetParName(15,"<pt>BG2");
//  fitMeanpt->SetParName(16,"<pt>PsiP");
//  
//  fitMeanpt->FixParameter(0,kVWG);
//  fitMeanpt->FixParameter(1,mVWG);
//  fitMeanpt->FixParameter(2,sVWG1);
//  fitMeanpt->FixParameter(3,sVWG2);
//  fitMeanpt->FixParameter(4,kJPsi);
//  fitMeanpt->FixParameter(5,mJPsi);
//  fitMeanpt->FixParameter(6,sJPsi);
//  fitMeanpt->FixParameter(7,alphaLow);
//  fitMeanpt->FixParameter(8,nLow);
//  fitMeanpt->FixParameter(9,alphaUp);
//  fitMeanpt->FixParameter(10,nUp);
//  fitMeanpt->FixParameter(11,kPsiP);
//  
//  fitMeanpt->SetParameter(12, 3.);
//  fitMeanpt->SetParLimits(12, 2.0,4.);
//  
//  fitMeanpt->SetParameter(13, 3.);
//  fitMeanpt->SetParLimits(13, 3.,10.);
//  
//  fitMeanpt->SetParameter(14, 0.2);
//  //  fitMeanpt->SetParLimits(14, 0.1,0.2);
//  
//  fitMeanpt->SetParameter(15, 0.1);
//  //  fitMeanpt->SetParLimits(15, 0.,1.);
//  
//  fitMeanpt->SetParameter(16, 3.);
//  fitMeanpt->SetParLimits(16, 1.5,6.);
//  
//  
//  TProfile::Approximate();
//  
//  const char* fitOption = "SER"; //+";
//  
//  TFitResultPtr fitResult = fHisto->Fit(fitMeanpt,fitOption,"");
//  
//  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
//  
//  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15));
//  
//  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);
//  
//  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
//  Set("MeanPtJPsi",fitMeanpt->GetParameter(12),fitMeanpt->GetParError(12));
//  Set("MeanPtPsiP",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
//  
//}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2VWG_BKGMPTPOL4()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();
  
  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  
  Double_t kVWG = GetValue("kVWG");
  Double_t mVWG = GetValue("mVWG");
  Double_t sVWG1 = GetValue("sVWG1");
  Double_t sVWG2 = GetValue("sVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");
  
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  
  //  TString msg;
  //
  //  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  //  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  //  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  //  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  //
  //  AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
  //
  //  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
  //
  //TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
  
  
  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }
  
  //_____________
  Int_t minBin = p->FindBin(fitRangeLow);
  Int_t maxBin = p->FindBin(fitRangeHigh);
  
  for ( Int_t i = minBin ; i <= maxBin ; i++ )
  {
    if ( p->GetBinEffectiveEntries(i) < 10 )
    {
      Double_t effEntries(0.),sumErr(0.);
      for ( Int_t j = i - 5 ; j < i ; j++ )
      {
        if ( j <= 0 ) continue;
        //        if ( j > p->GetNbinsX() ) break;
        
        effEntries += p->GetBinEffectiveEntries(j);
        sumErr += p->GetBinEffectiveEntries(j)*p->GetBinError(j);
      }
      
      Double_t meanErr = sumErr/effEntries;
      
      if ( p->GetBinError(i) < meanErr/2. )
      {
        std::cout << "Resetting bin " << i << " error" <<std::endl;
        p->SetBinError(i,meanErr);
      }
    }
  }
  //_____________
  
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4");
  
  bck->SetParameters(3.,-0.5,-0.5,1.5,-0.01);
  
  bck->SetParLimits(0, 1.,8.0);
  
  SetFitRejectRange(2.7,4.0);
  
  p->Fit(bck,"SER","");
  
  SetFitRejectRange();

  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL4,fitRangeLow,fitRangeHigh,19,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2VWGPOL4");
  
  fitMeanpt->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanpt->SetParName(11,"kPsiP");
  fitMeanpt->SetParName(12,"<pt>JPsi");
  fitMeanpt->SetParName(13,"<pt>BG0");
  fitMeanpt->SetParName(14,"<pt>BG1");
  fitMeanpt->SetParName(15,"<pt>BG2");
  fitMeanpt->SetParName(16,"<pt>BG3");
  fitMeanpt->SetParName(17,"<pt>BG4");
  fitMeanpt->SetParName(18,"<pt>PsiP");
  
  fitMeanpt->FixParameter(0,kVWG);
  fitMeanpt->FixParameter(1,mVWG);
  fitMeanpt->FixParameter(2,sVWG1);
  fitMeanpt->FixParameter(3,sVWG2);
  fitMeanpt->FixParameter(4,kJPsi);
  fitMeanpt->FixParameter(5,mJPsi);
  fitMeanpt->FixParameter(6,sJPsi);
  fitMeanpt->FixParameter(7,alphaLow);
  fitMeanpt->FixParameter(8,nLow);
  fitMeanpt->FixParameter(9,alphaUp);
  fitMeanpt->FixParameter(10,nUp);
  fitMeanpt->FixParameter(11,kPsiP);
  
  fitMeanpt->SetParameter(12, 3.);
  fitMeanpt->SetParLimits(12, 2.0,4.);
  
  for ( Int_t i = 0; i < 5; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }

  fitMeanpt->SetParameter(18, 3.);
  fitMeanpt->SetParLimits(18, 0.,10.);
     
  const char* fitOption = "SER"; //+";//SER
  
  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");
  
  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  
  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(18) <= fitMeanpt->GetParError(18) || fitMeanpt->GetParError(18) >= 0.75*fitMeanpt->GetParameter(18) )
  {
    fitMeanpt->FixParameter(18, 0.);
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");
    
    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }

  for ( Int_t i = 0; i < 5; ++i )
  {
    bck->SetParameter(i, fitMeanpt->GetParameter(i+13));
  }
  
  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(12),fitMeanpt->GetParError(12));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
 
}

//Int_t iparMinv[12] ={0,1,2,3,4,5,6,7,8,9,10,11};
//Int_t iparMpt[17] ={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
//
//struct GlobalChi2 {
//  GlobalChi2(  vector< ROOT::Fit::Chi2Function * > & fVec)
//  {
//    for(unsigned int i = 0; i < fVec.size(); i++){
//      fChi2_Vec.push_back( fVec[i] );
//    }
//  }
//  
//  // parameter vector is first background (in common 1 and 2)
//  // and then is signal (only in 2)
//  double operator() (const double *par) const {
//    vector< vector<double> > pVec;
//    vector< double > dummyVec;
//    for (int i = 0; i < 2; ++i) dummyVec.push_back(par[iparMinv[i] ]);
//    pVec.push_back(dummyVec);
//    
//    dummyVec.clear();
//    for (int i = 0; i < 5; ++i) dummyVec.push_back(par[iparMpt[i] ]);
//    pVec.push_back(dummyVec);
//    
//    double val = 0;
//    for( size_t i = 0; i < fChi2_Vec.size(); i++ ){
//      val += (*fChi2_Vec[i])(&(pVec[i][0]));
//    }
//    
//    return val;
//  }
//  
//  vector< const ROOT::Math::IMultiGenFunction * > fChi2_Vec;
//};

////_____________________________________________________________________________
//void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECOMB_CB2VWG_MPTCB2VWG_BKGMPTPOL2()
//{
//  
//  //============== Take the initial parameters if any
//  fHisto->GetListOfFunctions()->Delete(); //Minv Histo
//  
//  Double_t alphaLow = GetValue("alJPsi"); //Initial parameters for Minv Fit
//  Double_t nLow = GetValue("nlJPsi");
//  Double_t alphaUp = GetValue("auJPsi");
//  Double_t nUp = GetValue("nuJPsi");
//  Double_t fitRangeLow = GetValue(kFitRangeLow);
//  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
//  Double_t paramSPsiP = GetValue("FSigmaPsiP");
//  
//  TString msg;
//  
//  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
//  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
//  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
//  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
//  
//  AliDebug(1,Form("Fit with jpsi + psiprime VWG %s",msg.Data()));
//  
//  //==============
//  
//  //============== Definition of the Minv and Mpt fitting funtions
//     //_____Minv 
//  TF1* fitTotalMinv = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,12,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");
//  
//  fitTotalMinv->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi",
//  //                        0      1       2       3       4      5      6
//                        "alJPsi","nlJPsi","auJPsi","nuJPsi");
//  //                         7        8        9        10
//  fitTotalMinv->SetParName(11, "kPsiP");
//  //                             11
//  
//  
//  TF1* bckMinv = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG");
//  
//  const char* fitOption = "SER"; //We can add NO to avoid plotting  
//  
//  Int_t bin = fHisto->FindBin(0.26);
//  
//  bckMinv->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
//  
//  SetFitRejectRange(2.7,4.0);
//  
//  TFitResultPtr fitResultInit = fHisto->Fit(bckMinv,"SR");
//  
//  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
//  
//  if ( static_cast<int>(fitResultInit) )
//  {
//    bin = fHisto->FindBin(0.82);
//    bckMinv->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
//    fitResultInit = fHisto->Fit(bckMinv,"SR");
//  }
//  else if ( bckMinv->GetParameter(0) < 0 )
//  {
//    bckMinv->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
//  }
//  
//  SetFitRejectRange();
//  
//  TH1* fHistoMpt(0x0);
//      //______Mpt
//  TProfile* p(0x0);
//  if ( fHistoMpt->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHistoMpt);
//  else
//  {
//    AliError("Mean pt histo has to be a TProfile");
//    return;
//  }
//  
//                 //_____________ To assing a correct error to bins with low nEntries
//  Int_t minBin = p->FindBin(fitRangeLow);
//  Int_t maxBin = p->FindBin(fitRangeHigh);
//  
//  for ( Int_t i = minBin ; i <= maxBin ; i++ )
//  {
//    if ( p->GetBinEffectiveEntries(i) < 10 )
//    {
//      Double_t effEntries(0.),sumErr(0.);
//      for ( Int_t j = i - 5 ; j < i ; j++ )
//      {
//        if ( j <= 0 ) continue;
//        
//        effEntries += p->GetBinEffectiveEntries(j);
//        sumErr += p->GetBinEffectiveEntries(j)*p->GetBinError(j);
//      }
//      
//      Double_t meanErr = sumErr/effEntries;
//      
//      if ( p->GetBinError(i) < meanErr/2. )
//      {
//        std::cout << "Resetting bin " << i << " error" <<std::endl;
//        p->SetBinError(i,meanErr);
//      }
//    }
//  }
//               //_____________
//  
//  TF1* bckMpt = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");
//  
//  bckMpt->SetParameters(3.,1.,0.);
//  
//  bckMpt->SetParLimits(0, 1.,8.0);
//  
//  SetFitRejectRange(2.3,4.0);
//  
//  p->Fit(bckMpt,"SER","",fitRangeLow,fitRangeHigh);
//  
//  SetFitRejectRange();
//  
//  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL2,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2VWGPOL2");
//  
//  fitMeanpt->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
//  fitMeanpt->SetParName(11,"kPsiP");
//  fitMeanpt->SetParName(12,"<pt>JPsi");
//  fitMeanpt->SetParName(13,"<pt>BG0");
//  fitMeanpt->SetParName(14,"<pt>BG1");
//  fitMeanpt->SetParName(15,"<pt>BG2");
//  fitMeanpt->SetParName(16,"<pt>PsiP");
//  //___________
//  
//  //============== End of definition of the Minv and Mpt fitting funtions
//  
//  
//  
//  //============== Set the initial parameters in the Minv and Mpt fitting funtions
//  Double_t parMinv[12] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//  Double_t parMpt[17] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//  
//  
//  for ( Int_t i = 0; i < 4; ++i )
//  {
//    fitTotalMinv->SetParameter(i, bckMinv->GetParameter(i));
//    fitMeanpt->SetParameter(i, bckMinv->GetParameter(i));
//    parMinv[i] = bckMinv->GetParameter(i);
//    parMpt[i] = bckMinv->GetParameter(i);
//  }
//  
//  bin = fHisto->FindBin(3.09);
//  fitTotalMinv->SetParameter(4, fHisto->GetBinContent(bin)); // norm
//  fitMeanpt->SetParameter(4, fHisto->GetBinContent(bin)); // norm
//  parMinv[4] = fHisto->GetBinContent(bin);
//  parMpt[4] = fHisto->GetBinContent(bin);
//  
//  fitTotalMinv->SetParameter(5, 3.1); // mean
//  fitTotalMinv->SetParLimits(5, 3.0, 3.2);
//  fitMeanpt->SetParameter(5, 3.1); // mean
//  fitMeanpt->SetParLimits(5, 3.0, 3.2);
//  parMinv[5] = 3.1;
//  parMpt[5] = 3.1;
//  
//  fitTotalMinv->SetParameter(6, 0.08); // sigma
//  fitTotalMinv->SetParLimits(6, 0.05, 0.09);
//  fitMeanpt->SetParameter(6, 0.08); // sigma
//  fitMeanpt->SetParLimits(6, 0.05, 0.09);
//  parMinv[6] = 0.08;
//  parMpt[6] = 0.08;
//  
//  
//  fitTotalMinv->SetParameter(7,0.9);
//  fitTotalMinv->SetParLimits(7,0.1,10.0);
//  fitMeanpt->SetParameter(7,0.9);
//  fitMeanpt->SetParLimits(7,0.1,10.0);
//  parMinv[7] = 0.9;
//  parMpt[7] = 0.9;
//  
//  fitTotalMinv->SetParameter(8,5.0);
//  fitTotalMinv->SetParLimits(8,0.0,10.0);
//  fitMeanpt->SetParameter(8,5.0);
//  fitMeanpt->SetParLimits(8,0.0,10.0);
//  parMinv[8] = 5.0;
//  parMpt[8] = 5.0;
//  
//  
//  fitTotalMinv->SetParameter(9, 2.0);
//  fitTotalMinv->SetParLimits(9,0.1,10.0);
//  fitMeanpt->SetParameter(9, 2.0);
//  fitMeanpt->SetParLimits(9,0.1,10.0);
//  parMinv[9] = 2.0;
//  parMpt[9] = 2.0;
//  
//  fitTotalMinv->SetParameter(10,3.0);
//  fitTotalMinv->SetParLimits(10,0.0,10.0);
//  fitMeanpt->SetParameter(10,3.0);
//  fitMeanpt->SetParLimits(10,0.0,10.0);
//  parMinv[10] = 3.0;
//  parMpt[10] = 3.0;
//  
//  
//  fitTotalMinv->SetParameter(11, 10.); //kPsi'
//  fitTotalMinv->SetParLimits(11, 0.,100000.);
//  fitMeanpt->SetParameter(11, 10.); //kPsi'
//  fitMeanpt->SetParLimits(11, 0.,100000.);
//  parMinv[11] = 10.;
//  parMpt[11] = 10.;
//  
//  
//  fitMeanpt->SetParameter(12, 3.);
//  fitMeanpt->SetParLimits(12, 1.0,5.);
//  parMpt[12] = 3.;
//  
//  for ( Int_t i = 0; i < 3; ++i )
//  {
//    fitMeanpt->SetParameter(i + 13, bckMpt->GetParameter(i));
//    parMpt[i + 13] = bckMpt->GetParameter(i);
//  }
//  
//  fitMeanpt->SetParameter(16, 3.);
//  parMpt[16] = 3.;
//  Double_t psipPtLim = 10.;
//  fitMeanpt->SetParLimits(16, 0.,psipPtLim);
//
//  //================== End of setting Minv and Mpt initial parameters
//
//  //================== Chi2 definition for the global fit
//  
//  
//
////  struct GlobalChi2 {
////    GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
////               ROOT::Math::IMultiGenFunction & f2) :
////    fChi2_1(&f1), fChi2_2(&f2) {}
////    
////    // parameter vector is first background (in common 1 and 2)
////    // and then is signal (only in 2)
////    double operator() (const double *par) const {
////      double p1[12];
////      for (int i = 0; i < 12; ++i) p1[i] = par[iparMinv[i] ];
////      
////      double p2[17];
////      for (int i = 0; i < 17; ++i) p2[i] = par[iparMpt[i] ];
////      
////      return (*fChi2_1)(p1) + (*fChi2_2)(p2);
////    }
////    
////    const  ROOT::Math::IMultiGenFunction * fChi2_1;
////    const  ROOT::Math::IMultiGenFunction * fChi2_2;
////  };
//  //================== End of chi2 definition for the global fit
//  
//  
//  
//  //================== Combined fit
//  ROOT::Math::WrappedMultiTF1 wfMinv(*fitTotalMinv,1);
//  ROOT::Math::WrappedMultiTF1 wfMpt(*fitMeanpt,1);
//  
//  ROOT::Fit::DataOptions opt;
//  ROOT::Fit::DataRange rangeMinv;
//  // set the data range
//  rangeMinv.SetRange(fitRangeLow,fitRangeHigh);
//  ROOT::Fit::BinData dataMinv(opt,rangeMinv);
//  ROOT::Fit::FillData(dataMinv, fHisto);
//  
//  ROOT::Fit::DataRange rangeMpt;
//  rangeMpt.SetRange(fitRangeLow,fitRangeHigh);
//  ROOT::Fit::BinData dataMpt(opt,rangeMpt);
//  ROOT::Fit::FillData(dataMpt, fHistoMpt);
//  
////  ROOT::Fit::Chi2Function chi2_Minv(dataMinv, wfMinv);
////  ROOT::Fit::Chi2Function chi2_Mpt(dataMpt, wfMpt);
//  
//  vector< ROOT::Fit::Chi2Function * > chi2_vec;
//  chi2_vec.push_back( new ROOT::Fit::Chi2Function(dataMinv, wfMinv) );
//  chi2_vec.push_back( new ROOT::Fit::Chi2Function(dataMpt, wfMpt) );
//  
//  GlobalChi2 globalChi2(chi2_vec);
//
////  GlobalChi2 globalChi2(chi2_Minv, chi2_Mpt);
//  
//  ROOT::Fit::Fitter fitter;
//
////  const int Npar = 6;
////  double par0[Npar] = { 5,5,-0.1,100, 30,10};
////  
////  // create before the parameter settings in order to fix or set range on them
////  fitter.Config().SetParamsSettings(6,par0);
//  // fix 5-th parameter
////  fitter.Config().ParSettings(4).Fix();
////  // set limits on the third and 4-th parameter
////  fitter.Config().ParSettings(2).SetLimits(-10,-1.E-4);
////  fitter.Config().ParSettings(3).SetLimits(0,10000);
////  fitter.Config().ParSettings(3).SetStepSize(5);
//  
//  fitter.Config().MinimizerOptions().SetPrintLevel(0);
//  fitter.Config().SetMinimizer("Minuit2","Migrad");
//  
//  // fit FCN function directly
//  // (specify optionally data size and flag to indicate that is a chi2 fit)
//  fitter.FitFCN(6,globalChi2,0,dataMinv.Size()+dataMpt.Size(),true);
//  ROOT::Fit::FitResult result = fitter.Result();
//  result.Print(std::cout);
//  
//  TCanvas * c1 = new TCanvas("Simfit","Simultaneous fit of two histograms",
//                             10,10,700,700);
//  c1->Divide(1,2);
//  c1->cd(1);
//  gStyle->SetOptFit(1111);
//  
//  fitTotalMinv->SetFitResult( result, iparMinv);
//  fitTotalMinv->SetRange(rangeMinv().first, rangeMinv().second);
//  fitTotalMinv->SetLineColor(kBlue);
////  hB->GetListOfFunctions()->Add(fB);
//  fitTotalMinv->Draw();
//  
//  c1->cd(2);
//  fitMeanpt->SetFitResult( result, iparMpt);
//  fitMeanpt->SetRange(rangeMpt().first, rangeMpt().second);
//  fitMeanpt->SetLineColor(kRed);
////  fitMeanpt->GetListOfFunctions()->Add(fSB);
//  fitMeanpt->Draw();
//
////==============
//}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuJpsiResult::AddFit(const char* fitType)
{
  // Add a fit to this result
  
  if ( !fHisto ) return kFALSE;

  TH1* histo = static_cast<TH1*>(fHisto->Clone(fitType));
  
  AliAnalysisMuMuJpsiResult* r = new AliAnalysisMuMuJpsiResult(GetParticle(),*histo,fitType);

  if ( !r->IsValid() )
  {
    delete r;
    return kFALSE;
  }
  
  TMethodCall callEnv;

  TString fittingMethod(r->GetFitFunctionMethodName().Data());
  
  std::cout << "+Using fitting method " << fittingMethod.Data() << "..." << std::endl;
  std::cout << "" << std::endl;
  
  callEnv.InitWithPrototype(IsA(),fittingMethod.Data(),"");
  
  if (callEnv.IsValid())
  {
    callEnv.Execute(r);
  }
  else
  {
    AliError(Form("Could not get the method %s",fittingMethod.Data()));
    delete r;
    return kFALSE;
  }
  
//  Float_t lpar[] = { -1.0, -1.0, -1.0, -1.0 }; // free tails by default
//
//  Bool_t ok(kTRUE);
//
//  if ( fitFunction.Contains("CB2") )
//  {
//    if ( extra.Length() )
//    {
//      AliDebug(1,Form("sFitType=%s",fitType));
//      
//      sscanf(extra,"ALPHALOW%fNLOW%fALPHAUP%fNUP%f",
//             &lpar[0],&lpar[1],&lpar[2],&lpar[3]);
//      
//      AliDebug(1,Form("PSILOW ALPHALOW=%f NLOW=%f ALPHAUP=%f NUP=%f",lpar[0],lpar[1],lpar[2],lpar[3]));
//    }
//    
//    if ( lpar[0] == 0.0 && lpar[1] == 0.0 && lpar[0] == 0.0 && lpar[1] == 0.0 )
//    {
//      AliError("Cannot work with zero tails !");
//      ok = kFALSE;
//    }
//  }
  
//  if ( ok == kTRUE )
//  {
//    
//    if ( fitFunction=="PSICB2" )
//    {
//      r = FitPSICB2(*hminv,fitMinvMin,fitMinvMax);
//    }
//    else if ( fitFunction == "PSIPSIPRIMECB2VWG")
//    {
//      r = FitPSIPSIPRIMECB2VWG(*hminv,fitMinvMin,fitMinvMax,lpar[0],lpar[1],lpar[2],lpar[3]);
//    }
////    else if ( fitFunction == "PSIPSIPRIMECB2POL2EXP")
////    {
////      r = FitPSIPSIPRIMECB2POL2EXP(*hminv,fitMinvMin,fitMinvMax,lpar[0],lpar[1],lpar[2],lpar[3]);
////    }
//    else if ( fitFunction == "PSILOWMCTAILS" )
//    {
//      if ( npar!= 4 )
//      {
//        AliError("Cannot use PSILOWMCTAILS without being given the MC tails !");
//        delete hminv;
//        return kFALSE;
//      }
//      r = FitPSIPSIPRIMECB2VWG(*hminv,fitMinvMin,fitMinvMax,par[0],par[1],par[2],par[3]);
//      if (r)
//      {
//        r->SetAlias("MCTAILS");
//      }
//    }
//    else if ( fitFunction == "COUNTJPSI" )
//    {
//      r = new AliAnalysisMuMuJpsiResult(*hminv,"COUNTJPSI");
//      Double_t n = CountParticle(*hminv,"Jpsi");
//      r->Set("NofJPsi",n,TMath::Sqrt(n));
//    }
//  
//  }
  
  if ( r->IsValid() )
  {
    StdoutToAliDebug(1,r->Print(););
    r->SetBin(Bin());
    r->SetNofTriggers(NofTriggers());
    r->SetNofRuns(NofRuns());

    Bool_t adoptOK = AdoptSubResult(r);
    if ( adoptOK ) std::cout << "Subresult " << r->GetName() << " adopted in " << GetName() <<  std::endl;
    else AliError(Form("Could not adopt subresult %s",r->GetName()));
  }
  else
  {
    delete r;
    r=0x0;
  }
  
  return (r!=0x0);
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::DecodeFitType(const char* fitType)
{
  /// decode the string containing all the information for the fit to be done
  ///
  /// the fittype is a combination of key=value pairs separated by ":"
  ///     +func=Name fo the function used in the fit(See ? to read the naming convention)
  ///     +rebin=Rebin of the histo(1,2...)
  ///     +histoType=The type of histo we wan to fit(minv,mpt ???or minv&mpt for a combined fit???)
  ///     +tails=from where we take the particle tails(mctails,mctailsJPsi&PsiP to use an indep sets of tails for each particle, predefined tails)
  ///     +alJpsi , nlPsiP ... = Fit parameters (they carry the parameter name and the particle name)
  ///
  /// e.g.
  /// func=gaus:rebin=2:range=2.3;4.7;nsigmas=3
  /// func=vwm+cb2:range=2;5;alphalow=5.2;alphahigh=6.2
  ///
  /// except from the func key, all the other ones have default values
  /// if the key is different from func,rebin,range,histoType it is assumed
  /// to be the value for a parameter of the fit function.
  
  TString resultName("");

  // default values
  TString fitFunction;
  TString histoType("minv");
  TString tails("");
  Int_t rebin=1;
  Double_t fitMinvMin=2.0;
  Double_t fitMinvMax=5.0;
  Double_t paramSPsiP= 3.68609/3.096916;
  
  TString sFitType(fitType);
  
  if (!sFitType.Contains(kKeyFunc,TString::kIgnoreCase)) return;
  
  TObjArray* parts = sFitType.Tokenize(":");
  TObjString* str;
  TIter next(parts);
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TString key,value;
    Bool_t ok = GetKeyValue(str->String(),'=',key,value);
    if (!ok)
    {
      AliErrorClass(Form("Invalid key=value pair %s",str->String().Data()));
      continue;
    }
    
    if ( key.CompareTo(kKeyFunc,TString::kIgnoreCase) == 0 )
    {
      fitFunction = value;
      resultName += fitFunction;
    }
    else if ( key.CompareTo(kKeyRange,TString::kIgnoreCase) == 0 )
    {
      TString xmin,xmax;
      if (GetKeyValue(value,';',xmin,xmax))
      {
        fitMinvMin = xmin.Atof();
        fitMinvMax = xmax.Atof();
        AliInfoClass(Form("xmin=%e xmax=%e",fitMinvMin,fitMinvMax));
        
        resultName += "_";
        resultName += xmin.Data();
        resultName += "-";
        resultName += xmax.Data();
      }
      else
      {
        AliErrorClass(Form("Improper range specification %s",value.Data()));
        continue;
      }
    }
    else if ( key.CompareTo(kKeyRebin,TString::kIgnoreCase) == 0 )
    {
      rebin = value.Atoi();
    }
    else if ( key.CompareTo(kKeyHistoType,TString::kIgnoreCase) == 0 ) //FIXME::Is really necesary to save the histoType? I think I dont use it
    {
      histoType = value;
   
      if ( histoType.CompareTo("minv",TString::kIgnoreCase) == 0 ) Set(kKeyHistoType,0.,0.0); //histoType=0 means minv histo
      else if ( histoType.CompareTo("mpt",TString::kIgnoreCase) == 0 ) Set(kKeyHistoType,1.,0.0); //histoType=1 means mpt histo
      else if ( histoType.CompareTo("minv&mpt",TString::kIgnoreCase) == 0 ) Set(kKeyHistoType,2.,0.0); //histoType=1 means combined fit minv and mpt
      else
      {
        AliErrorClass(Form("Improper histoType specification %s",value.Data()));
        continue;
      }
     
    }
    else if ( key.CompareTo(kKeyTails,TString::kIgnoreCase) == 0 )
    {
      tails = value;
      if ( tails.CompareTo("mctails",TString::kIgnoreCase) == 0 ) Set(kKeyTails,0.,0.0);
      else if ( tails.CompareTo("mctailsJPsi&PsiP",TString::kIgnoreCase) == 0 ) Set(kKeyTails,1.,0.0);
      else if ( tails.CompareTo("",TString::kIgnoreCase) == 0 ) Set(kKeyTails,2.,0.0); // Predefined tails
      else
      {
        AliErrorClass(Form("Improper tails specification %s",value.Data()));
        continue;
      }
    }
    else if ( key.CompareTo(kKeySPsiP,TString::kIgnoreCase) == 0 )
    {
      paramSPsiP = value.Atof();
    }
    else if ( key.CompareTo(kKeyMinvRS,TString::kIgnoreCase) == 0 )
    {
      fMinvRS = value.Data();
    }
    else
    {
      Set(key.Data(),value.Atof(),0.0);
    }
  }

  if ( fitFunction.CountChar('-') )
  {
    AliError(Form("Invalid fit function name : %s",fitFunction.Data()));
    Invalidate();
  }
  else
  {
  
    fFitFunction = fitFunction;

    Set(kKeySPsiP,paramSPsiP,0.0);
    Set(kKeyRebin,rebin,0.0);
    Set(kFitRangeLow,fitMinvMin,0.0);
    Set(kFitRangeHigh,fitMinvMax,0.0);
  }
  
  delete parts;
}

//_____________________________________________________________________________
TString AliAnalysisMuMuJpsiResult::GetFitFunctionMethodName() const
{
    /// Get the name of the function used to fit this result (if any)
  TString name(FitFunctionName());
  if ( name.CountChar('-') )
  {
    return "";
  }
  if ( GetValue("tails") == 1. ) name += "INDEPTAILS";
  return TString::Format("Fit%s",name.Data());
}

//_____________________________________________________________________________
Long64_t AliAnalysisMuMuJpsiResult::Merge(TCollection* list)
{
  /// Merge method
  ///
  /// Merge a list of AliAnalysisMuMuJpsiResult objects with this
  /// Returns the number of merged objects (including this).
  ///
  /// Note that the merging is to be understood here as an weighted mean operation
  ///
  /// FIXME ! (compared to base class Merge, should only Minv merging ?)
  
  AliError("Implement me !");
  if (!list) return 0;
  
  if (list->IsEmpty()) return 1;
  
  return 0;
  
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuJpsiResult::NofRuns() const
{
  /// Get the number of runs
  if ( !Mother() ) return fNofRuns;
  else return Mother()->NofRuns();
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuJpsiResult::NofTriggers() const
{
  /// Get the number of triggers
  
  if ( !Mother() ) return fNofTriggers;
  else return Mother()->NofTriggers();
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::Print(Option_t* opt) const
{
  /// printout

  std::cout << Form("NRUNS %d - NTRIGGER %10d - %s",
                                      NofRuns(),
                                      NofTriggers(),
                                      fBin.AsString().Data());
  
  AliAnalysisMuMuResult::Print(opt);
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::PrintValue(const char* key, const char* opt, Double_t value, Double_t errorStat,
                                           Double_t rms) const
{
  /// exclude the particles with zero stat
  
  const std::map<std::string,Double_t>& m = MassMap();

  for( std::map<std::string,Double_t>::const_iterator it = m.begin(); it != m.end(); ++it )
  {
    TString particle(it->first.c_str());
    
    if (TString(key).Contains(particle.Data()))
    {
      if ( GetValue("Nof%s",particle.Data()) <= 0.0 ) return;
    }
  }
  
  AliAnalysisMuMuResult::PrintValue(key,opt,value,errorStat,rms);

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::PrintParticle(const char* particle, const char* opt) const
{
  /// Print all information about one particule type
  
  Double_t npart = GetValue(Form("Nof%s",particle));
  if (npart<=0) return;
  
  
  std::cout << opt << Form("\t%s",particle) << std::endl;
  
  //  Double_t npartError = GetErrorStat(Form("Nof%s",particle));
//  std::cout << opt << Form("\t\t%20s %9.2f +- %5.2f","Count",npart,npartError) << std::endl;
  
  TIter next(Keys());
  TObjString* key;
  
  while ( ( key = static_cast<TObjString*>(next()) ) )
  {
    if ( !key->String().Contains(particle) ) continue;
    
    PrintValue(key->String(),opt,GetValue(key->String()),GetErrorStat(key->String()),GetRMS(key->String()));
  } 
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::SetBin(const AliAnalysisMuMuBinning::Range& bin)
{
  /// Set the bin
  
  if (!Mother()) fBin = bin;
  else Mother()->SetBin(bin);
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::SetFitRejectRange(Double_t a, Double_t b)
{
  /// Set a range the fit function(s) can ignore

  fRejectFitPoints = kFALSE;

  fFitRejectRangeLow = a;
  fFitRejectRangeHigh = b;
  if ( a <= TMath::Limits<Double_t>::Max() && b <= TMath::Limits<Double_t>::Max() )
  {
    fRejectFitPoints = kTRUE;
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::SetHisto(const TH1& h)
{
  /// Set the spectrum to be fitted.
  static UInt_t n(0);
  
  delete fHisto;
  fHisto = static_cast<TH1*>(h.Clone(Form("AliAnalysisMuMuJpsiResultHisto%u",n++)));
  fHisto->SetDirectory(0);
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::SetNofInputParticles(const TH1& hminv)
{
  /// Set the number of input particle from the invariant mass spectra
  
  static const char* particleNames[] = { "JPsi" , "PsiP", "Upsilon","UpsilonPrime" };

  const std::map<std::string,Double_t>& m = MassMap();

  for ( Int_t i = 0; i < 4; ++i )
  {
    std::map<std::string,Double_t>::const_iterator it = m.find(particleNames[i]);

    Double_t sigma(-1.0);

    if (it != m.end() )
    {
      sigma = it->second*0.1;
    }

    Double_t n = CountParticle(hminv,particleNames[i],sigma);

    AliDebug(1,Form("i=%d particle %s n %e",i,particleNames[i],n));
    
    if ( n > 0 )
    {
      SetNofInputParticles(particleNames[i],TMath::Nint(n));
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::SetNofInputParticles(const char* particle, int n)
{
  /// Set the number of input particles (so it is a MC result)
  /// and (re)compute the AccxEff values
  
  Set(Form("NofInput%s",particle),n,TMath::Sqrt(n));
  
  if (n<=0)
  {
    Set(Form("AccEff%s",particle),0,0);
    return;
  }
  
  Double_t npart = GetValue(Form("Nof%s",particle));
  Double_t npartErr  = GetErrorStat(Form("Nof%s",particle));
  Double_t ninput = GetValue(Form("NofInput%s",particle));
  Double_t ninputErr = GetErrorStat(Form("NofInput%s",particle));
  
  Set(Form("AccEff%s",particle),
      npart/ninput,
      (npart/ninput)*ErrorAB(npart,npartErr,ninput,ninputErr));
  
  TIter next(SubResults());
  AliAnalysisMuMuJpsiResult* r;
  
  while ( ( r = static_cast<AliAnalysisMuMuJpsiResult*>(next())) )
  {
    r->Set(Form("NofInput%s",particle),n,TMath::Sqrt(n));

    npart = r->GetValue(Form("Nof%s",particle));
    npartErr = r->GetErrorStat(Form("Nof%s",particle));
    
    r->Set(Form("AccEff%s",particle),
           npart/ninput,
           (npart/ninput)*ErrorAB(npart,npartErr,ninput,ninputErr));

  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::SetNofRuns(Int_t n)
{
  if ( !Mother() ) fNofRuns=n;
  else Mother()->SetNofRuns(n);
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::SetNofTriggers(Int_t n)
{
  if ( !Mother() ) fNofTriggers=n;
  else Mother()->SetNofTriggers(n);
}


//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::SetParameter(TF1* func, Int_t npar, Double_t fix, Double_t initialValue,
                                             Double_t min, Double_t max) const
{
  /// Fix one parameter or set its initial value and limits
  
  if ( IsValidValue(fix) )
  {
    func->FixParameter(npar, fix);
  }
  else
  {
    func->SetParameter(npar,initialValue);
    func->SetParLimits(npar,min,max);
  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuJpsiResult::StrongCorrelation(TFitResultPtr& r,
                                                    TF1* fitFunction,
                                                    Int_t npar1,
                                                    Int_t npar2,
                                                    Double_t fixValueIfWrong)
{
  // return kTRUE if the npar1-th and npar2-th parameters of the fit function
  // are too strongly correlated,
  // and in that case fix the npar2-th parameter's value to fixValueIfWrong
  
  Bool_t strongCorrelation = TMath::Abs(r->GetCorrelationMatrix()(npar1,npar2)) > 0.90;
  
  if ( strongCorrelation )
  {
    fitFunction->FixParameter(npar2,fixValueIfWrong);
    return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuJpsiResult::WrongParameter(TF1* fitFunction, Int_t npar,
                                                 Double_t fixValueIfWrong)
{
  // return kTRUE if npar-th parameter of fit function has a big error,
  // and in that case fix the parameter's value to fixValueIfWrong
  
  Bool_t wrong = (fitFunction->GetParError(npar) > 0.8*TMath::Abs(fitFunction->GetParameter(npar)));
  
  AliWarning(Form("npar %d error %e val %e wrong %d",
                  npar,fitFunction->GetParError(npar),
                  fitFunction->GetParameter(npar),wrong));
  
  if ( wrong )
  {
    AliWarning(Form("Fixing parameter %d of %s to %e",
                    npar,fitFunction->GetName(),fixValueIfWrong));
    
    fitFunction->FixParameter(npar,fixValueIfWrong);
    return kTRUE;
  }
  return kFALSE;
}
