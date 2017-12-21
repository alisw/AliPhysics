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
#include <iostream>

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TMinuit.h"
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
      massMap["JPsi"]         = 3.096916e+00;
      massMap["PsiP"]         = 3.68609e+00;
      massMap["Upsilon"]      = 9.46030e+00;
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

  const TString kKeyFunc      = "Func";
  const TString kKeyRange     = "Range";
  const TString kKeyRebin     = "Rebin";
  const TString kFitRangeLow  = "FitRangeLow";
  const TString kFitRangeHigh = "FitRangeHigh";
  const TString kKeyCount     = "Count";
  const TString kKeyHistoType = "HistoType";
  const TString kKeyTails     = "Tails";
  const TString kKeyWeight    = "Weight";
  const TString kKeySPsiP     = "FSigmaPsiP"; //Factor to fix the psi' sigma to sigmaJPsi*SigmaPsiP (Usually factor SigmaPsiP = 1, 0.9 and 1.1)
  const TString kKeyMinvRS    = "MinvRS"; // FIXME: not very correct since "MinvRS" is in AliAnalysisMuMu::GetParametersFromResult

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

  SetWeight(GetValue(kKeyWeight));

  AliDebug(1,Form("weight = %f\n",Weight()));

  TString name(fFitFunction);
  if ( !name.Contains("PSICB2") && !name.Contains("PSINA60NEW") && !name.Contains("PSICOUNT") ) //To avoid adding things to the name of simu results
  {
    Bool_t isMPt = kFALSE;

    if ( name.BeginsWith("PSI") ) name.ReplaceAll("PSIPSIPRIME","");
    else if ( name.BeginsWith("MPT") ) {
      name.ReplaceAll("MPTPSIPSIPRIME","");
      isMPt = kTRUE;
    }
    else if ( name.BeginsWith("MV2") ) {
      name.ReplaceAll("MV2PSIPSIPRIME","");
      isMPt = kTRUE;
    }
    if(( TString(fitType).Contains("AccEff") )){
      name +="AccEff";
    }

    name += "_";
    name += Form("%1.1f",GetValue(kFitRangeLow));
    name += "_";
    name += Form("%1.1f",GetValue(kFitRangeHigh));
    if(GetValue(kKeyWeight)!=1.){
      name += "_";
      name += Form("Weight=%1.1f",GetValue(kKeyWeight));
    }

    if ( !isMPt ) {
      name += "_";
      name += Form("SP%1.1f",GetValue(kKeySPsiP));
    }
    else name += Form("(Sig:%s)",fMinvRS.Data());
  }
  SetName(name.Data());

  Int_t rebin = TMath::Nint(GetValue(kKeyRebin));

  if (rebin>0) fHisto->Rebin(rebin);

  if ( fHisto->GetEntries()<100 && !TString(GetName()).Contains(kKeyCount) ){
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

    fNofRuns             = rhs.NofRuns();
    fNofTriggers         = rhs.NofTriggers();
    fBin                 = rhs.Bin();
    fTriggerClass        = rhs.fTriggerClass;
    fEventSelection      = rhs.fEventSelection;
    fPairSelection       = rhs.fPairSelection;
    fCentralitySelection = rhs.fCentralitySelection;
    fFitFunction         = rhs.fFitFunction;
    fFitRejectRangeLow   = rhs.fFitRejectRangeLow;
    fFitRejectRangeHigh  = rhs.fFitRejectRangeHigh;
    fRejectFitPoints     = rhs.fRejectFitPoints;
    fParticle            = rhs.fParticle;
    fMinvRS              = rhs.fMinvRS;

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
    signal->SetNpx(8000);
    fHisto->GetListOfFunctions()->Add(signal);
  }

  if ( bck )
  {
    bck->SetLineColor(2);
    bck->SetNpx(8000);
    fHisto->GetListOfFunctions()->Add(bck);
  }

  if ( total )
  {
    total->SetLineColor(4);
    total->SetNpx(8000);
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
    signal2->SetNpx(8000);
    fHisto->GetListOfFunctions()->Add(signal2);
  }

  if ( bck )
  {
    bck->SetLineColor(2);
    bck->SetNpx(8000);
    fHisto->GetListOfFunctions()->Add(bck);
  }

  if ( total )
  {
    total->SetLineColor(4);
    total->SetNpx(8000);
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

  if ( HasValue(Form("Nof%s",particle)) ){
    //Protection
    if (!other.HasValue(Form("AccEff%s",particle),subResultName)){
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


     AliDebug(1,Form("Nof%s = %f +/- %f",particle,GetValue(Form("Nof%s",particle)),GetErrorStat(Form("Nof%s",particle))));
     AliDebug(1,Form("AccEff%s = %f +/- %f",particle,other.GetValue(Form("AccEff%s",particle),subResultName),other.GetErrorStat(Form("AccEff%s",particle),subResultName)));
     AliDebug(1,Form("CorrNof%s = %f +/- %f",particle,value,error*value));

    return kTRUE;
  }

  else AliError(Form("Result does not have Nof%s : cannot correct it !",particle));
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
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol1Pol2(Double_t *x, Double_t *par)
{
  // Linear function for Bkg 2 params

  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
  else TF1::RejectPoint(kFALSE);
  return (par[0]*x[0] + par[1] )/(par[2]*x[0]*x[0] + par[3]*x[0] + par[4]) ;
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Pol3(Double_t *x, Double_t *par)
{
  // Linear function for Bkg 7 params

  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
  else TF1::RejectPoint(kFALSE);
  return (  par[0]*x[0]*x[0] + x[0]*par[1] + par[2]  )/( par[3]*x[0]*x[0]*x[0]+ par[4]*x[0]*x[0] + par[5]*x[0] + par[6]) ;
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Pol3V2(Double_t *x, Double_t *par)
{
  // Linear function for Bkg 6 params

  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
  else TF1::RejectPoint(kFALSE);
  return (  par[0]*x[0]*x[0] + x[0]*par[1] + par[2]  )/( (par[3] + x[0])*(par[4] + x[0])*(par[5] + x[0])) ;
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
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Cheb(Double_t *x, Double_t *par)
{
  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
  Double_t xmin = 2.2;
  Double_t xmax = 4.7;
  double xx = (2.0 * x[0] - xmin -xmax)/(xmax-xmin);
  const int order = 4;
  Double_t fT[order+1] = {0};
  if (order == 1) return par[0];
  if (order == 2) return par[0] + xx*par[1];
  // build the polynomials
  fT[0] = 1;
  fT[1] = xx;
  for (int i = 1; i< order; ++i) {
    fT[i+1] =  2 *xx * fT[i] - fT[i-1];
  }
  double sum = par[0]*fT[0];
  for (int i = 1; i<= order; ++i) {
    sum += par[i] * fT[i];
  }
  return sum;
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPolExp(Double_t *x, Double_t *par)
{
  // pol2 x exp : 4 params

  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
//  return par[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0])*TMath::Exp(par[4]/x[0]);
  return (par[0]+par[1]*x[0])*TMath::Exp(par[2]*x[0]);
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
  else TF1::RejectPoint(kFALSE);
  Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
  return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG2(Double_t *x, Double_t *par)
{
  // gaussian variable width : 5 params

  if (fRejectFitPoints && x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
  else TF1::RejectPoint(kFALSE);
  Double_t sigma = par[2] + par[3]*((x[0]-par[1])/par[1]) + par[4]*((x[0]-par[1])/par[1])*((x[0]-par[1])/par[1]);
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

  Double_t sigmaRatio(0.);
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
    par[5]+(3.68609-3.096916),
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
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewVWG2(Double_t *x, Double_t *par)
{
  /// 2 NA60 (new) + pol2 x exp
  /// width of the second NA60 related to the first (free) one.

  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[11] = {
    par[16],
    par[6]+(3.68609-3.096916),
    par[7]*SPsiPFactor, // /3.096916*3.68609
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
    par[15],
  };
  return FitFunctionBackgroundVWG2(x, par) + FitFunctionNA60New(x, &par[5]) + FitFunctionNA60New(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewPol1Pol2(Double_t *x, Double_t *par)
{
  /// 2 NA60 (new) + pol2 x exp
  /// width of the second NA60 related to the first (free) one.

  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[11] = {
    par[16],
    par[6]+(3.68609-3.096916),
    par[7]*SPsiPFactor, // /3.096916*3.68609
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
    par[15],
  };
  return FitFunctionBackgroundPol1Pol2(x, par) + FitFunctionNA60New(x, &par[5]) + FitFunctionNA60New(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewPol2Pol3(Double_t *x, Double_t *par)
{
  /// 2 NA60 (new) + pol2 x exp
  /// width of the second NA60 related to the first (free) one.

  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[11] = {
    par[18],
    par[8]+(3.68609-3.096916),
    par[9]*SPsiPFactor, // /3.096916*3.68609
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
    par[15],
    par[16],
    par[17],
  };
  return FitFunctionBackgroundPol2Pol3(x, par) + FitFunctionNA60New(x, &par[7]) + FitFunctionNA60New(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewPol2Exp(Double_t *x, Double_t *par)
{
  /// 2 NA60 (new) + pol2 x exp
  /// width of the second NA60 related to the first (free) one.

  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[11] = {
    par[15],
    par[5]+(3.68609-3.096916),
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
    par[7]+(3.68609-3.096916),
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
    par[3]+(3.68609-3.096916),
    par[4]*SPsiPFactor, // /3.096916*3.68609,
    par[5],
    par[6],
    par[7],
    par[8]
  };
  return FitFunctionBackgroundLin(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[2]) + FitFunctionSignalCrystalBallExtended(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol1Pol2(Double_t *x, Double_t *par)
{
  /// 2 extended crystal balls + Pol1
  /// width of the second CB related to the first (free) one.

  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[7] = {
    par[12],
    par[6]+(3.68609-3.096916),
    par[7]*SPsiPFactor, // /3.096916*3.68609,
    par[8],
    par[9],
    par[10],
    par[11]
  };
  return FitFunctionBackgroundPol1Pol2(x, par)   + FitFunctionSignalCrystalBallExtended(x, &par[5]) + FitFunctionSignalCrystalBallExtended(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol2Pol3(Double_t *x, Double_t *par)
{
  /// 2 extended crystal balls + Pol1
  /// width of the second CB related to the first (free) one.

  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[7] = {
    par[14],
    par[8]+(3.68609-3.096916),
    par[9]*SPsiPFactor, // /3.096916*3.68609,
    par[10],
    par[11],
    par[12],
    par[13]
  };
  return FitFunctionBackgroundPol2Pol3(x, par)   + FitFunctionSignalCrystalBallExtended(x, &par[7]) + FitFunctionSignalCrystalBallExtended(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol2Pol3V2(Double_t *x, Double_t *par)
{
  /// 2 extended crystal balls + Pol2/pol3
  /// width of the second CB related to the first (free) one.

  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[7] = {
    par[13],
    par[7]+(3.68609-3.096916),
    par[9]*SPsiPFactor, // /3.096916*3.68609,
    par[9],
    par[10],
    par[11],
    par[12]
  };
  return FitFunctionBackgroundPol2Pol3V2(x, par)   + FitFunctionSignalCrystalBallExtended(x, &par[6]) + FitFunctionSignalCrystalBallExtended(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol2Exp(Double_t *x, Double_t *par)
{
  /// 2 extended crystal balls + pol2 x exp
  /// width of the second CB related to the first (free) one.

  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[7] = {
    par[11],
    par[5]+(3.68609-3.096916),
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
//    par[6]+(3.68609-3.096916),
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
    par[7]+(3.68609-3.096916),
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
//    par[8]+(3.68609-3.096916),
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
    par[5]+(3.68609-3.096916),
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };
  return FitFunctionBackgroundVWG(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionSignalCrystalBallExtended(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG2(Double_t *x, Double_t *par)
{
  /// 2 extended crystal balls + VWG2
  /// width of the second CB related to the first (free) one.

  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[7] = {
    par[12],
    par[6]+(3.68609-3.096916),
    par[7]*SPsiPFactor, // /3.096916*3.68609,
    par[8],
    par[9],
    par[10],
    par[11]
  };
  return FitFunctionBackgroundVWG2(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[5]) + FitFunctionSignalCrystalBallExtended(x, par2);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWGINDEPTAILS(Double_t *x, Double_t *par)
{
  /// 2 extended crystal balls + pol2 x exp
  /// The tail parameters are independent but the sPsiP and mPsiP are fixed to the one of the JPsi

  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[7] = {
    par[11],
    par[5]+(3.68609-3.096916),
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
Double_t AliAnalysisMuMuJpsiResult::hFunction(Double_t*x, Double_t* par)
{
  return par[0] + par[1]*(x[0]-3.096916) + TMath::Power(par[2]*(x[0]-3.096916),2) + TMath::Power(par[3]*(x[0]-3.096916),3) + TMath::Power(par[4]*(x[0]-3.096916),4) + TMath::Power(par[5]*(x[0]-3.096916),5);
}

//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::alphaCB2VWG(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[4])/(FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionBackgroundVWG(x,par));
}
//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::alphaCB2VWG2(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[5])/(FitFunctionSignalCrystalBallExtended(x, &par[5]) + FitFunctionBackgroundVWG2(x,par));
}
//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::alphaCB2POL1POL2(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[5])/(FitFunctionSignalCrystalBallExtended(x, &par[5]) + FitFunctionBackgroundPol1Pol2(x,par));
}
//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::alphaCB2POL2POL3(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[7])/(FitFunctionSignalCrystalBallExtended(x, &par[7]) + FitFunctionBackgroundPol2Pol3(x,par));
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
Double_t AliAnalysisMuMuJpsiResult::alphaNA60NEWPOL1POL2(Double_t*x, Double_t* par)
{
  return FitFunctionNA60New(x, &par[5])/(FitFunctionNA60New(x, &par[5]) + FitFunctionBackgroundPol1Pol2(x,par));
}

//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::alphaNA60NEWPOL2EXP(Double_t*x, Double_t* par)
{
  return FitFunctionNA60New(x, &par[4])/(FitFunctionNA60New(x, &par[4]) + FitFunctionBackgroundPol2Exp(x,par));
}
//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::alphaNA60NEWVWG2(Double_t*x, Double_t* par)
{
  return FitFunctionNA60New(x, &par[5])/(FitFunctionNA60New(x, &par[5]) + FitFunctionBackgroundVWG2(x,par));
}

//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::alphaNA60NEWPOL2POL3(Double_t*x, Double_t* par)
{
  return FitFunctionNA60New(x, &par[7])/(FitFunctionNA60New(x, &par[7]) + FitFunctionBackgroundPol2Pol3(x,par));
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
    par[5]+(3.68609-3.096916),
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
  AliDebug(2,Form("pol2 parameters in FitFunctionMeanPtSCB2VWGPOL2 : %f,%f,%f",par[13],par[14],par[15]));
  return alphaCB2VWG(x,par)*par[12] + (1. - alphaCB2VWG(x,par))*FitFunctionBackgroundPol2(x,&par[13]);
}
//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2VWG2POL2(Double_t* x, Double_t* par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi
  AliDebug(2,Form("pol2 parameters in FitFunctionMeanPtSCB2VWG2POL2 : %f,%f,%f",par[14],par[15],par[16]));
  return alphaCB2VWG2(x,par)*par[13] + (1. - alphaCB2VWG2(x,par))*FitFunctionBackgroundPol2(x,&par[14]);
}
//------------------------------------------------------------------------------
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2VWG2POLEXP(Double_t* x, Double_t* par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi
  AliDebug(2,Form("pol2 parameters in FitFunctionMeanPtSCB2VWG2POL2 : %f,%f,%f",par[14],par[15],par[16]));
  return alphaCB2VWG2(x,par)*par[13] + (1. - alphaCB2VWG2(x,par))*FitFunctionBackgroundPolExp(x,&par[14]);
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2VWG2POL2EXP(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi
  AliDebug(2,Form("pol2 parameters in FitFunctionMeanPtSCB2VWG2POL2EXP : %f,%f,%f,%f",par[14],par[15],par[16],par[17]));

  return alphaCB2VWG2(x,par)*par[13] + (1. - alphaCB2VWG2(x,par))*FitFunctionBackgroundPol2Exp(x,&par[14]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2VWG2POL4(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi
  AliDebug(2,Form("pol2 parameters in FitFunctionMeanPtSCB2VWG2POL4 : %f,%f,%f,%f,%f",par[14],par[15],par[16],par[17],par[18]));

  return alphaCB2VWG2(x,par)*par[13] + (1. - alphaCB2VWG2(x,par))*FitFunctionBackgroundPol4(x,&par[14]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2VWG2POL4Cheb(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi
  AliDebug(2,Form("pol2 parameters in FitFunctionMeanPtSCB2VWG2POL4Cheb : %f,%f,%f,%f,%f",par[14],par[15],par[16],par[17],par[18]));

  return alphaCB2VWG2(x,par)*par[13] + (1. - alphaCB2VWG2(x,par))*FitFunctionBackgroundPol4Cheb(x,&par[14]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2POL2POL3POL2(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi
  AliDebug(2,Form("pol2 parameters in FitFunctionMeanPtSCB2POL2POL3POL2 : %f,%f,%f",par[16],par[17],par[18]));

  return alphaCB2POL2POL3(x,par)*par[15] + (1. - alphaCB2POL2POL3(x,par))*FitFunctionBackgroundPol2(x,&par[16]);
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2POL2POL3_POLEXP(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi
  AliDebug(2,Form("pol2 parameters in FitFunctionMeanPtSCB2POL2POL3POL2 : %f,%f,%f",par[16],par[17],par[18]));

  return alphaCB2POL2POL3(x,par)*par[15] + (1. - alphaCB2POL2POL3(x,par))*FitFunctionBackgroundPolExp(x,&par[16]);
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2POL2POL3POL2EXP(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi
  AliDebug(2,Form("pol2exp parameters in FitFunctionMeanPtSCB2POL2POL3POL2EXP : %f,%f,%f,%f",par[16],par[17],par[18],par[19]));

  return alphaCB2POL2POL3(x,par)*par[15] + (1. - alphaCB2POL2POL3(x,par))*FitFunctionBackgroundPol2Exp(x,&par[16]);
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2POL2POL3POL4(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi
  AliDebug(2,Form("pol4 parameters in FitFunctionMeanPtSCB2POL2POL3POL4 : %f,%f,%f,%f,%f",par[16],par[17],par[18],par[19],par[20]));

  return alphaCB2POL2POL3(x,par)*par[15] + (1. - alphaCB2POL2POL3(x,par))*FitFunctionBackgroundPol4(x,&par[16]);
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2POL2POL3POL4Cheb(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi
  AliDebug(2,Form("pol4 parameters in FitFunctionMeanPtSCB2POL2POL3POL4 : %f,%f,%f,%f,%f",par[16],par[17],par[18],par[19],par[20]));

  return alphaCB2POL2POL3(x,par)*par[15] + (1. - alphaCB2POL2POL3(x,par))*FitFunctionBackgroundPol4Cheb(x,&par[16]);
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
    par[5]+(3.68609-3.096916),
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };


  return alphaCB2VWG(x,par)*par[12] + alphaCB2VWG(x,par2)*par[16] + (1. - alphaCB2VWG(x,par) - alphaCB2VWG(x,par2))*FitFunctionBackgroundPol2(x,&par[13]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2POL1POL2POL2(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[12] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[4],
    par[12], //kPsi'
    par[6]+(3.68609-3.096916),
    par[7]*SPsiPFactor, // /3.096916*3.68609,
    par[8],
    par[9],
    par[10],
    par[11]
  };


  return alphaCB2POL1POL2(x,par)*par[13] + alphaCB2POL1POL2(x,par2)*par[17] + (1. - alphaCB2POL1POL2(x,par) - alphaCB2POL1POL2(x,par2))*FitFunctionBackgroundPol2(x,&par[14]);
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
    par[5]+(3.68609-3.096916),
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };


  return alphaCB2VWG(x,par)*par[12] + alphaCB2VWG(x,par2)*par[17] + (1. - alphaCB2VWG(x,par) - alphaCB2VWG(x,par2))*FitFunctionBackgroundPol2Exp(x,&par[13]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2POL1POL2POL2EXP(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[12] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[4],
    par[12], //kPsi'
    par[6]+(3.68609-3.096916),
    par[7]*SPsiPFactor, // /3.096916*3.68609,
    par[8],
    par[9],
    par[10],
    par[11]
  };


  return alphaCB2POL1POL2(x,par)*par[13] + alphaCB2POL1POL2(x,par2)*par[18] + (1. - alphaCB2POL1POL2(x,par) - alphaCB2POL1POL2(x,par2))*FitFunctionBackgroundPol2Exp(x,&par[14]);
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
    par[5]+(3.68609-3.096916),
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
    par[5]+(3.68609-3.096916),
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
    par[5]+(3.68609-3.096916),
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
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWVWG2POLEXP(Double_t *x, Double_t *par)
{
 AliDebug(2,Form("pol2 parameters in FitFunctionMeanPtSNA60NEWVWG2POLEXP : %f,%f,%f",par[18],par[19],par[20]));
  return alphaNA60NEWVWG2(x,par)*par[17] + (1. - alphaNA60NEWVWG2(x,par))*FitFunctionBackgroundPolExp(x,&par[18]);
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWVWG2POL2(Double_t *x, Double_t *par)
{
 AliDebug(2,Form("pol2 parameters in FitFunctionMeanPtSNA60NEWVWG2POL2 : %f,%f,%f",par[18],par[19],par[20]));
  return alphaNA60NEWVWG2(x,par)*par[17] + (1. - alphaNA60NEWVWG2(x,par))*FitFunctionBackgroundPol2(x,&par[18]);
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWVWG2POL2EXP(Double_t *x, Double_t *par)
{
 AliDebug(2,Form("pol2exp parameters in FitFunctionMeanPtSNA60NEWVWG2POL2EXP : %f,%f,%f,%f",par[18],par[19],par[20],par[21]));
  return alphaNA60NEWVWG2(x,par)*par[17] + (1. - alphaNA60NEWVWG2(x,par))*FitFunctionBackgroundPol2Exp(x,&par[18]);
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWVWG2POL4(Double_t *x, Double_t *par)
{
 AliDebug(2,Form("pol4 parameters in FitFunctionMeanPtSNA60NEWVWG2POL4 : %f,%f,%f,%f,%f",par[18],par[19],par[20],par[21],par[22]));
  return alphaNA60NEWVWG2(x,par)*par[17] + (1. - alphaNA60NEWVWG2(x,par))*FitFunctionBackgroundPol4(x,&par[18]);
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWVWG2POL4Cheb(Double_t *x, Double_t *par)
{
 AliDebug(2,Form("pol4 parameters in FitFunctionMeanPtSNA60NEWVWG2POL4 : %f,%f,%f,%f,%f",par[18],par[19],par[20],par[21],par[22]));
  return alphaNA60NEWVWG2(x,par)*par[17] + (1. - alphaNA60NEWVWG2(x,par))*FitFunctionBackgroundPol4Cheb(x,&par[18]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWPOL2POL3_POLEXP(Double_t *x, Double_t *par)
{
 AliDebug(2,Form("pol2 parameters in FitFunctionMeanPtSNA60NEWPOL2POL3_POLEXP : %f,%f,%f",par[20],par[21],par[22]));
  return alphaNA60NEWPOL2POL3(x,par)*par[19] + (1. - alphaNA60NEWPOL2POL3(x,par))*FitFunctionBackgroundPolExp(x,&par[20]);
}

//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWPOL2POL3POL2(Double_t *x, Double_t *par)
{
 AliDebug(2,Form("pol2 parameters in FitFunctionMeanPtSNA60NEWPOL2POL3POL2 : %f,%f,%f",par[20],par[21],par[22]));
  return alphaNA60NEWPOL2POL3(x,par)*par[19] + (1. - alphaNA60NEWPOL2POL3(x,par))*FitFunctionBackgroundPol2(x,&par[20]);
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWPOL2POL3POL2EXP(Double_t *x, Double_t *par)
{
 AliDebug(2,Form("pol2exp parameters in FitFunctionMeanPtSNA60NEWPOL2POL3POL2EXP : %f,%f,%f,%f",par[20],par[21],par[22],par[23]));
  return alphaNA60NEWPOL2POL3(x,par)*par[19] + (1. - alphaNA60NEWPOL2POL3(x,par))*FitFunctionBackgroundPol2Exp(x,&par[20]);
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWPOL2POL3POL4(Double_t *x, Double_t *par)
{
 AliDebug(2,Form("pol4 parameters in FitFunctionMeanPtSNA60NEWPOL2POL3POL4 : %f,%f,%f,%f,%f",par[20],par[21],par[22],par[23],par[24]));
  return alphaNA60NEWPOL2POL3(x,par)*par[19] + (1. - alphaNA60NEWPOL2POL3(x,par))*FitFunctionBackgroundPol4(x,&par[20]);
}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWPOL2POL3POL4Cheb(Double_t *x, Double_t *par)
{
 AliDebug(2,Form("pol4 parameters in FitFunctionMeanPtSNA60NEWPOL2POL3POL4 : %f,%f,%f,%f,%f",par[20],par[21],par[22],par[23],par[24]));
  return alphaNA60NEWPOL2POL3(x,par)*par[19] + (1. - alphaNA60NEWPOL2POL3(x,par))*FitFunctionBackgroundPol4Cheb(x,&par[20]);
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
    par[5]+(3.68609-3.096916),
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
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2NA60NEWPOL1POL2POL2(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[16] = {
    par[0],//a
    par[1],//b
    par[2],//a'
    par[3],//b'
    par[4],//c'
    par[16],//kPsiP
    par[6]+(3.68609-3.096916),//mPsiP
    par[7]*SPsiPFactor, // /3.096916*3.68609,
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
    par[15],
  };

  return alphaNA60NEWPOL1POL2(x,par)*par[17] + alphaNA60NEWPOL1POL2(x,par2)*par[21] + (1. - alphaNA60NEWPOL1POL2(x,par) - alphaNA60NEWPOL1POL2(x,par2))*FitFunctionBackgroundPol2(x,&par[18]);

}
//____________________________________________________________________________
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2NA60NEWPOL1POL2POL2EXP(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = GetValue(kKeySPsiP);

  Double_t par2[16] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[4],
    par[16],
    par[6]+(3.68609-3.096916),
    par[7]*SPsiPFactor, // /3.096916*3.68609,
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
    par[15],
  };


  return alphaNA60NEWPOL1POL2(x,par)*par[17] + alphaNA60NEWPOL1POL2(x,par2)*par[22] + (1. - alphaNA60NEWPOL1POL2(x,par) - alphaNA60NEWPOL1POL2(x,par2))*FitFunctionBackgroundPol2Exp(x,&par[18]);

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
    par[5]+(3.68609-3.096916),
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
    par[5]+(3.68609-3.096916),
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
Double_t AliAnalysisMuMuJpsiResult::FitFunctionMeanPtHFunction(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi mean pt resolution test

  Double_t xlim[2] = {2.5,3.29};

  if(x[0] <= xlim[0])       return hFunction(&xlim[0],par) + par[6]*(3.096916-xlim[0]);
  else if (x[0] >= xlim[1]) return hFunction(&xlim[1],par) + par[7]*(3.096916-xlim[1]);
  else                      return hFunction(x,par);
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
    par[5]+(3.68609-3.096916),
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
    par[5]+(3.68609-3.096916),
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
   par[5]+(3.68609-3.096916),
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

  Double_t alphaLow     = GetValue(Form("al%s",particleName.Data()));
  Double_t nLow         = GetValue(Form("nl%s",particleName.Data()));
  Double_t alphaUp      = GetValue(Form("au%s",particleName.Data()));
  Double_t nUp          = GetValue(Form("nu%s",particleName.Data()));
  Double_t fitRangeLow  = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  TString msg;

  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  if (IsValidValue(nLow)) msg     += TString::Format("nLow=%e ",nLow);
  if (IsValidValue(alphaUp)) msg  += TString::Format("alphaUp=%e ",alphaUp);
  if (IsValidValue(nUp)) msg      += TString::Format("nUp=%e ",nUp);

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
    fitTotal->SetParLimits(2, 0.03, 0.2);
  }
  else if (particleName.Contains("PsiP"))
  {
    fitTotal->SetParameter(1, 3.7);
    fitTotal->SetParLimits(1, 3.63, 3.72);
    fitTotal->SetParameter(2, 0.08);
    fitTotal->SetParLimits(2, 0.03, 0.2);
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

  TFitResultPtr fitResult = fHisto->Fit(fitTotal,"SERL","");

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
    fitResult = fHisto->Fit(fitTotal,"SERL","");
  }

  TF1* signal = new TF1("signal",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","SignalCB2");

  signal->SetParameters(fitTotal->GetParameter(0),fitTotal->GetParameter(1),fitTotal->GetParameter(2),fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6));

  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
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

  // new TCanvas;
  // fHisto->DrawCopy();

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
    fitTotal->SetParLimits(2, 0.03, 0.2);
  }
  else if (particleName.Contains("PsiP"))
  {
    fitTotal->SetParameter(1, 3.7);
    fitTotal->SetParLimits(1, 3.63, 3.72);
    fitTotal->SetParameter(2, 0.08);
    fitTotal->SetParLimits(2, 0.03, 0.2);
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

  TFitResultPtr fitResult = fHisto->Fit(fitTotal,"SERL","");

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

  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
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

// //_____________________________________________________________________________
// void AliAnalysisMuMuJpsiResult::FitPSICB2VWG()
// {
//     /// Fit using 1 extended crystal ball (signal) + variable width gaussian (background)
//     //Refresh
//     fHisto->GetListOfFunctions()->Delete();
//     //Check if values given in MuMuConfig
//     Double_t alphaLow = GetValue("alJPsi");
//     Double_t nLow = GetValue("nlJPsi");
//     Double_t alphaUp = GetValue("auJPsi");
//     Double_t nUp = GetValue("nuJPsi");
//     Double_t fitRangeLow = GetValue(kFitRangeLow);
//     Double_t fitRangeHigh = GetValue(kFitRangeHigh);

//     TString msg;
//     //Output Message
//     if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
//     if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
//     if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
//     if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);

//     AliDebug(1,Form("Fit with jpsi VWG %s",msg.Data()));

//     //  std::cout << "Tails parameters: " << msg.Data() << std::endl;
//     //Create pointer and configure signal+bck fit functions
//     TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalOneCB2VWG,fitRangeLow,fitRangeHigh,11,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");

//     fitTotal->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi",
//     //                        0      1       2       3       4      5      6
//                           "alJPsi","nlJPsi","auJPsi","nuJPsi");
//     //                         7        8        9        10


//     // Create pointer for bck function
//     TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG");

//     const char* fitOption = "SERLI"; //We can add NO to avoid plotting

// //#if 0
// //    bck->SetParameter(0,fHisto->GetMaximum());
// //    bck->SetParameter(1,3);
// //    bck->SetParameter(2,10);
// //    bck->SetParameter(3,10);
// //
// //    bck->SetParLimits(1, 0., 100.);
// //    bck->SetParLimits(2, 0., 100.);
// //    bck->SetParLimits(3, 0., 100.);
// //
// //    SetFitRejectRange(2.7,3.5);
// //
// //    fHisto->Fit(bck,fitOption,"");
// //
// //    for ( Int_t i = 0; i < 4; ++i )
// //        {
// //        Double_t a,b;
// //        bck->GetParLimits(i,a,b);
// //        fitTotal->SetParameter(i,bck->GetParameter(i));
// //        fitTotal->SetParLimits(i,a,b);
// //        }
// //#endif
//     // Create pointer for bck fitting
//     TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG,1.,6.,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG");
//     // To set first parameter
//     Int_t bin = fHisto->FindBin(0.26);
//     // Setting parameters
//     bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
//     // Reject region ot J/Psi
//     SetFitRejectRange(2.7,4.0);
//     // Store fit results in a TFitResult
//     TFitResultPtr fitResultInit = fHisto->Fit(bckInit,"SRL");

//     std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
//     // if success a new fit is proceed
//     if ( static_cast<int>(fitResultInit) )
//         {
//         bin = fHisto->FindBin(0.82);
//         bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
//         fitResultInit = fHisto->Fit(bckInit,"SRL");
//         }
//     else if ( bckInit->GetParameter(0) < 0 )
//         {
//         bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);
//         }
//     // ??
//     SetFitRejectRange();
//     // Set bck paramaters from precedent fit to fitTotal
//     for ( Int_t i = 0; i < 4; ++i )
//         {
//         fitTotal->SetParameter(i, bckInit->GetParameter(i));
//         }

//     // Set J/Psi parameters to FitTotal
//     bin = fHisto->FindBin(3.10);
//     fitTotal->SetParameter(4, fHisto->GetBinContent(bin)); // norm

//     fitTotal->SetParameter(5, 3.12); // mean
//     fitTotal->SetParLimits(5, 3.0, 3.2);

//     fitTotal->SetParameter(6, 0.10); // sigma
//     fitTotal->SetParLimits(6, 0.05, 0.09);
//     // Set parameters from MuMuConfig if present
//     if ( IsValidValue(alphaLow) )
//         {
//         fitTotal->FixParameter(7, alphaLow);
//         }
//     else
//         {
//         fitTotal->SetParameter(7,0.9);
//         fitTotal->SetParLimits(7,0.1,10.0);
//         }

//     if ( IsValidValue(nLow) )
//         {
//         fitTotal->FixParameter(8, nLow);
//         }
//     else
//         {
//         fitTotal->SetParameter(8,5.0);
//         fitTotal->SetParLimits(8,0.0,10.0);
//         }

//     if ( IsValidValue(alphaUp) )
//         {
//         fitTotal->FixParameter(9, alphaUp);
//         }
//     else
//         {
//         fitTotal->SetParameter(9, 2.0);
//         fitTotal->SetParLimits(9,0.1,10.0);
//         }

//     if ( IsValidValue(nUp) )
//         {
//         fitTotal->FixParameter(10, nUp);
//         }
//     else
//         {
//         fitTotal->SetParameter(10,3.0);
//         fitTotal->SetParLimits(10,0.0,10.0);
//         }


//     SetFitRejectRange();

//     //  std::cout << fitTotal->GetParameter(0) << std::endl; //Just a xcheck

//     TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");

//     //  std::cout << fitTotal->GetParameter(0) << " ?= " << fitResult->Parameter(0) << std::endl; //Just a xcheck
//     // Output message
//     std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
//     std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

//     //  std::cout << fitTotal->GetParameter(0) << std::endl; //Just a xcheck
//     // Refit in case of error on kVWG to high
//     if ( static_cast<int>(fitResult) )
//         {
//         if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) )
//             {
//             std::cout << "//-------Refitting again (setting kVWG=kVWG/2)" << std::endl;

//             fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG
//             }

//         fitResult = fHisto->Fit(fitTotal,fitOption,"");
//         std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
//         }
//     // Refit in case of KkVWG to low
//     if ( static_cast<int>(fitResult) )
//         {
//         std::cout << "//-------Refitting again (setting kVWG=kVWG*2)" << std::endl;

//         fitTotal->SetParameter(0, fHisto->GetMaximum()*2.); // kVWG

//         fitResult = fHisto->Fit(fitTotal,fitOption,"");
//         std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
//         }
//     // Idem
//     if ( static_cast<int>(fitResult) )
//         {
//         std::cout << "//-------Refitting again (setting kVWG=kVWG/2)" << std::endl;

//         fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG

//         fitResult = fHisto->Fit(fitTotal,fitOption,"");
//         std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
//         }
//     // Set parameters to fit Total and fit with it
//     if ( static_cast<int>(fitResult) )
//         {
//         for ( Int_t i = 0; i < 4; ++i )
//             {
//             fitTotal->SetParameter(i, bckInit->GetParameter(i));
//             }
//         fitResult = fHisto->Fit(fitTotal,fitOption,"");

//         if ( static_cast<int>(fitResult) ) std::cout << "//-------Cannot fit properly, try something else..." << std::endl;
//         }

//     delete bckInit;
//     // Set things ...
//     Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
//     Set("FitNDF",fitTotal->GetNDF(),0.0);

//     Set("mJPsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));
//     Set("sJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
//     // Create JPsi function and set parameters
//     TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
//     signalJPsi->SetParameters(fitTotal->GetParameter(4),
//                               fitTotal->GetParameter(5),
//                               fitTotal->GetParameter(6),
//                               fitTotal->GetParameter(7),
//                               fitTotal->GetParameter(8),
//                               fitTotal->GetParameter(9),
//                               fitTotal->GetParameter(10));
//     // Set bck parameter
//     bck->SetParameters(fitTotal->GetParameter(0),
//                        fitTotal->GetParameter(1),
//                        fitTotal->GetParameter(2),
//                        fitTotal->GetParameter(3));

//     // Set value + error for each parameters
//     Set("kVWG",fitTotal->GetParameter(0),fitTotal->GetParError(0));
//     Set("mVWG",fitTotal->GetParameter(1),fitTotal->GetParError(1));
//     Set("sVWG1",fitTotal->GetParameter(2),fitTotal->GetParError(2));
//     Set("sVWG2",fitTotal->GetParameter(3),fitTotal->GetParError(3));
//     Set("kJPsi",fitTotal->GetParameter(4),fitTotal->GetParError(4));

//     Set("alJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
//     Set("nlJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
//     Set("auJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
//     Set("nuJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));

//     SetFitRejectRange();

//     AttachFunctionsToHisto(signalJPsi,bck,fitTotal,fitRangeLow,fitRangeHigh);


//     Double_t cbParameters[7];
//     Double_t covarianceMatrix[7][7];

//     for ( int ix = 0; ix < 7; ++ix )
//         {
//         cbParameters[ix] = fitTotal->GetParameter(ix+4);
//         }

//     for ( int iy = 0; iy < 7; ++iy )
//         {
//         for ( int ix = 0; ix < 7; ++ix )
//             {
//             covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+4,iy+4);
//             }
//         }

//     Double_t a = fHisto->GetXaxis()->GetXmin();
//     Double_t b = fHisto->GetXaxis()->GetXmax();
//     double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
//     double nerr = signalJPsi->IntegralError(a,b,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);

//     Set("NofJPsi",njpsi,nerr);

//     double m = GetValue("mJPsi");
//     double s = GetValue("sJPsi");
//     double njpsi3s = signalJPsi->Integral(m-3*s,m+3*s)/fHisto->GetBinWidth(1);
//     double nerr3s = signalJPsi->IntegralError(m-3*s,m+3*s,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);

//     Set("NofJPsi3s",njpsi3s,nerr3s);

//     //Computation of bin significance and signal over background

//     Double_t bkgParameters[4];
//     Double_t bkgcovarianceMatrix[4][4];

//     for ( int ix = 0; ix < 4; ++ix )
//         {
//         bkgParameters[ix] = fitTotal->GetParameter(ix);
//         }

//     for ( int iy = 0; iy < 4; ++iy )
//         {
//         for ( int ix = 0; ix < 4; ++ix )
//             {
//             bkgcovarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix,iy);
//             }
//         }

//     double nbck3s = bck->Integral(m-3.*s,m+3.*s)/fHisto->GetBinWidth(1);
//     double nbck3sErr = bck->IntegralError(m-3.*s,m+3.*s,&bkgParameters[0],&bkgcovarianceMatrix[0][0])/fHisto->GetBinWidth(1);

//     double sOverB3s = njpsi3s / nbck3s;
//     double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));

//     Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);

//     double sig = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
//     double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
//                                 TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );

//     Set("Significance3s",sig,sigErr);

// }



//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2VWG()
{
  /// Fit using 2 extended crystal balls (signal) + variable width gaussian (background)

  fHisto->GetListOfFunctions()->Delete();

  TString histoName       = fHisto->GetTitle();
  TString sfitOption      = histoName.Contains("Corrected") ? "0SERL" : "NO0SRLM"; //SERL
  const char* fitOption   = sfitOption.Data(); //We can add NO to avoid plotting
  const char* fitOptionBg = "NOSR"; //We can add NO to avoid plotting //SER


  //__________ Get tails parameters, fitting range and SigmaPsiP
  Double_t alphaLow     = GetValue("alJPsi");
  Double_t nLow         = GetValue("nlJPsi");
  Double_t alphaUp      = GetValue("auJPsi");
  Double_t nUp          = GetValue("nuJPsi");
  Double_t fitRangeLow  = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  Double_t paramSPsiP   = GetValue("FSigmaPsiP");
  Double_t meanJPsi     = GetValue("meanJPsi");
  Double_t sigmaJPsi    = GetValue("sigmaJPsi");
  Double_t binNormJPsi  = GetValue("binNormJPsi");
  Double_t binNormPsiP  = GetValue("binNormPsiP");

  TString msg;

  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  if (IsValidValue(nLow)) msg     += TString::Format("nLow=%e ",nLow);
  if (IsValidValue(alphaUp)) msg  += TString::Format("alphaUp=%e ",alphaUp);
  if (IsValidValue(nUp)) msg      += TString::Format("nUp=%e ",nUp);
  //__________


  AliDebug(1,Form("Fit with jpsi + psiprime VWG %s",msg.Data()));


  //__________ Define the function to fit the spectrum, and the background just for plotting
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG,fitRangeLow,fitRangeHigh,12,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG");

  fitTotal->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi",
//                        0      1       2       3       4      5      6
                        "alJPsi","nlJPsi","auJPsi","nuJPsi");
//                         7        8        9        10
  fitTotal->SetParName(11, "kPsiP");
//                            11

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG");

  //___________


  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG,1.8,10.,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG");

  Int_t bin = fHisto->FindBin(0.27);

  // Get param from FitType or FitSingle if any
  Double_t par0 = GetValue("kVWG");
  Double_t par1 = GetValue("mVWG");
  Double_t par2 = GetValue("sVWG1");
  Double_t par3 = GetValue("sVWG2");

  if (! IsValidValue(par0) ) par0 = fHisto->GetBinContent(bin);
  if (! IsValidValue(par1) ) par1 = 2.;
  if (! IsValidValue(par2) ) par2 = 0.5;
  if (! IsValidValue(par3) ) par3 = 0.3;

  bckInit->SetParameters(par0,par1,par2,par3);

  SetFitRejectRange(2.3,3.4);

  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOptionBg);

  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

  //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundVWG",fitOptionBg);
  //___________

  SetFitRejectRange();
  //____________


  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 4; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }

  if(IsValidValue(binNormJPsi)) bin = fHisto->FindBin(binNormJPsi); //
  else                          bin = fHisto->FindBin(3.09);
  fitTotal->SetParameter(4, fHisto->GetBinContent(bin)); // norm

  if(IsValidValue(meanJPsi))fitTotal->SetParameter(5,meanJPsi); // mean
  else                      fitTotal->SetParameter(5, 3.16); // mean

  fitTotal->SetParLimits(5, 2.95, 3.2);

  if(IsValidValue(sigmaJPsi)) fitTotal->SetParameter(5,sigmaJPsi);
  else                        fitTotal->SetParameter(6, 0.08); // sigma
  fitTotal->SetParLimits(6, 0.03, 0.2);

  if ( IsValidValue(alphaLow) )
  {
    fitTotal->FixParameter(7, alphaLow);
  }
  else
  {
    fitTotal->SetParameter(7,0.9);
    fitTotal->SetParLimits(7,0.8,1.2);
  }

  if ( IsValidValue(nLow) )
  {
    fitTotal->FixParameter(8, nLow);
  }
  else
  {
    fitTotal->SetParameter(8,3.0);
    fitTotal->SetParLimits(8,1.0,5.0);
  }

  if ( IsValidValue(alphaUp) )
  {
    fitTotal->FixParameter(9, alphaUp);
  }
  else
  {
    fitTotal->SetParameter(9, 2.0);
    fitTotal->SetParLimits(9,1.0,3.0);
  }

  if ( IsValidValue(nUp) )
  {
    fitTotal->FixParameter(10, nUp);
  }
  else
  {
    fitTotal->SetParameter(10,2.0);
    fitTotal->SetParLimits(10,1.0,2.3);
  }

  if(IsValidValue(binNormJPsi)) bin = fHisto->FindBin(binNormPsiP); //
  else bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(11, fHisto->GetBinContent(bin)*0.5); //kPsi'
  fitTotal->SetParLimits(11, 0.,fHisto->GetBinContent(bin));
  //______________


  //_____________First fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");

  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  //___________


  //___________ Further attempts to fit if the first one fails
  if ( ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 ) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,11,3);
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  //___________


  delete bckInit; //Delete the initial background funtion

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
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
                        fitTotal->GetParameter(5) + (3.68609-3.096916),
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

  double npsip = signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerrpsip = (fitTotal->GetParError(11)/fitTotal->GetParameter(11))*signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);

  Set("NofPsiP",npsip,nerrpsip);

  double mpsip = GetValue("mJPsi")+ (3.68609 - 3.096916);
  double spsip = GetValue("sJPsi")*paramSPsiP;
  double npsip3s = signalPsiP->Integral(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);
  double nerrpsip3s = (fitTotal->GetParError(11)/fitTotal->GetParameter(11))*signalPsiP->IntegralError(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);

  Set("NofPsiP3s",npsip3s,nerrpsip3s);
  //_____________________________


  //_____Computation of bin significance and signal over background
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
  //__________________________

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2VWG2()
{
  /// Fit using 2 extended crystal balls (signal) + variable width gaussian (background)

  fHisto->GetListOfFunctions()->Delete();

  TString histoName = fHisto->GetTitle();
  TString sfitOption= histoName.Contains("Corrected") ? "0SERL" : "NO0SRLM"; //SERLM
  const char* fitOption = sfitOption.Data(); //We can add NO to avoid plotting
  const char* fitOptionBg = "0LSR"; //We can add NO to avoid plotting//SER


  //__________ Get tails parameters, fitting range and SigmaPsiP
  Double_t alphaLow     = GetValue("alJPsi");
  Double_t nLow         = GetValue("nlJPsi");
  Double_t alphaUp      = GetValue("auJPsi");
  Double_t nUp          = GetValue("nuJPsi");
  Double_t fitRangeLow  = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  Double_t paramSPsiP   = GetValue("FSigmaPsiP");
  Double_t meanJPsi     = GetValue("meanJPsi");
  Double_t sigmaJPsi    = GetValue("sigmaJPsi");
  Double_t binNormJPsi  = GetValue("binNormJPsi");
  // Double_t binNormPsiP  = GetValue("binNormPsiP");

  Double_t mVWG2_init    = IsValidValue(GetValue("mVWG2_init"))  ? GetValue("mVWG2_init")  : 2.;
  Double_t s1VWG2_init   = IsValidValue(GetValue("s1VWG2_init")) ? GetValue("s1VWG2_init") : -0.5;
  Double_t s2VWG2_init   = IsValidValue(GetValue("s2VWG2_init")) ? GetValue("s2VWG2_init") : -0.9;
  Double_t gVWG2_init    = IsValidValue(GetValue("gVWG2_init"))  ? GetValue("gVWG2_init")  : 0.15;

  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.2;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 3.9;

  TString msg;

  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  if (IsValidValue(nLow)) msg     += TString::Format("nLow=%e ",nLow);
  if (IsValidValue(alphaUp)) msg  += TString::Format("alphaUp=%e ",alphaUp);
  if (IsValidValue(nUp)) msg      += TString::Format("nUp=%e ",nUp);
  msg += TString::Format("mVWG2init=%.2f ",mVWG2_init);
  msg += TString::Format("s1VWG2init=%.2f ",s1VWG2_init);
  msg += TString::Format("s2VWG2init=%.2f ",s2VWG2_init);
  msg += TString::Format("gVWG2init=%.2f ",gVWG2_init);
  //__________


  AliDebug(1,Form("Fit with jpsi + psiprime VWG %s",msg.Data()));


  //__________ Define the function to fit the spectrum, and the background just for plotting
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG2,fitRangeLow,fitRangeHigh,13,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG2");

  fitTotal->SetParNames("kVWG2","mVWG2","s1VWG2","s2VWG2","gVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi");
//                        0        1       2        3        4      5      6        7       8           9       10
  fitTotal->SetParName(11, "nuJPsi");
//                            11
  fitTotal->SetParName(12, "kPsiP");
//                            12

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG2,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG2");
  //___________

  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG2,1.8,6.,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG2");
  Int_t bin = fHisto->FindBin(0.26);


  bckInit->SetParameters(fHisto->GetBinContent(bin),mVWG2_init,s1VWG2_init,s2VWG2_init,gVWG2_init);
  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOptionBg);
  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundVWG2",fitOptionBg);// Further attempts to fit bkg if the first one fails
  SetFitRejectRange();
  //____________

  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 5; ++i ) fitTotal->SetParameter(i, bckInit->GetParameter(i));

  if(IsValidValue(binNormJPsi)) bin = fHisto->FindBin(binNormJPsi); //
  else bin = fHisto->FindBin(3.09);
  fitTotal->SetParameter(5, fHisto->GetBinContent(bin)); // nor


  if(IsValidValue(meanJPsi)) fitTotal->FixParameter(6,meanJPsi); // mean
  else {
    fitTotal->SetParameter(6, 3.11); // mean
    fitTotal->SetParLimits(6, 2.95, 3.2);
  }

  if(IsValidValue(sigmaJPsi))fitTotal->FixParameter(7,sigmaJPsi);
  else {
    fitTotal->SetParameter(7, 0.08); // sigma
    fitTotal->SetParLimits(7, 0.03, 0.2);
  }

  if ( IsValidValue(alphaLow) ) {
    fitTotal->FixParameter(8, alphaLow);
  } else {
    fitTotal->SetParameter(8,0.9);
    fitTotal->SetParLimits(8,0.1,10.0);
  }

  if ( IsValidValue(nLow) ) {
    fitTotal->FixParameter(9, nLow);
  } else {
    fitTotal->SetParameter(9,5.0);
    fitTotal->SetParLimits(9,0.0,10.0);
  }

  if ( IsValidValue(alphaUp) ) {
    fitTotal->FixParameter(10, alphaUp);
  } else {
    fitTotal->SetParameter(10, 2.0);
    fitTotal->SetParLimits(10,0.1,10.0);
  }

  if ( IsValidValue(nUp) ) {
    fitTotal->FixParameter(11, nUp);
  } else {
    fitTotal->SetParameter(11,3.0);
    fitTotal->SetParLimits(11,0.0,10.0);
  }

  bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(12, fHisto->GetBinContent(bin)*0.2); //kPsi'
  // fitTotal->SetParLimits(12, 0.,fHisto->GetBinContent(bin));
  //______________


  //_____________First fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");

  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  //___________


  //___________ Further attempts to fit if the first one fails
  if ( ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 ) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,12,4);
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  //___________

  delete bckInit; //Delete the initial background funtion

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);

  Set("mJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  Set("sJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));

  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalJPsi->SetParameters(fitTotal->GetParameter(5),
                     fitTotal->GetParameter(6),
                     fitTotal->GetParameter(7),
                     fitTotal->GetParameter(8),
                     fitTotal->GetParameter(9),
                     fitTotal->GetParameter(10),
                     fitTotal->GetParameter(11));

  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalPsiP->SetParameters(fitTotal->GetParameter(12),
                        fitTotal->GetParameter(6) + (3.68609 - 3.096916),
                        fitTotal->GetParameter(7)*paramSPsiP, // /3.096916*3.68609,
                        fitTotal->GetParameter(8),
                        fitTotal->GetParameter(9),
                        fitTotal->GetParameter(10),
                        fitTotal->GetParameter(11));

  bck->SetParameters(fitTotal->GetParameter(0),
                     fitTotal->GetParameter(1),
                     fitTotal->GetParameter(2),
                     fitTotal->GetParameter(3),
                     fitTotal->GetParameter(4));


  Set("kVWG2",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("mVWG2",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("s1VWG2",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("s2VWG2",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("gVWG2",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set("kJPsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));

  Set("alJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("nlJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set("auJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("nuJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));

  Set("kPsiP",fitTotal->GetParameter(12),fitTotal->GetParError(12));


  SetFitRejectRange();

  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);


  Double_t cbParameters[7];
  Double_t covarianceMatrix[7][7];

  for ( int ix = 0; ix < 7; ++ix )
  {
    cbParameters[ix] = fitTotal->GetParameter(ix+5);
  }

  for ( int iy = 0; iy < 7; ++iy )
  {
    for ( int ix = 0; ix < 7; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+5,iy+5);
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

  double npsip = signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerrpsip = (fitTotal->GetParError(12)/fitTotal->GetParameter(12))*signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);

  Set("NofPsiP",npsip,nerrpsip);

  double mpsip = GetValue("mJPsi")+ (3.68609 - 3.096916);
  double spsip = GetValue("sJPsi")*paramSPsiP;
  double npsip3s = signalPsiP->Integral(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);
  double nerrpsip3s = (fitTotal->GetParError(12)/fitTotal->GetParameter(12))*signalPsiP->IntegralError(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);

  Set("NofPsiP3s",npsip3s,nerrpsip3s);
  //_____________________________


  //_____Computation of bin significance and signal over background
  Double_t bkgParameters[5];
  Double_t bkgcovarianceMatrix[5][5];

  for ( int ix = 0; ix < 5; ++ix ){
    bkgParameters[ix] = fitTotal->GetParameter(ix);
  }

  for ( int iy = 0; iy < 5; ++iy ){
    for ( int ix = 0; ix < 5; ++ix ) {
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
  //__________________________

  // playground
    // TF1* fitTotalCopy = new TF1("signal+bckCopy",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2VWG2,fitRangeLow,fitRangeHigh,13,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2VWG2");
    // for (int i = 0; i < 13; ++i)
    // {
    //   if (i == 12 ) fitTotalCopy->SetParameter(i,fitTotal->GetParameter(5)*0.024); // kJpsi * N_PsiP(13TeV)/JPsi(13TeV)
    //   else fitTotalCopy->SetParameter(i,fitTotal->GetParameter(i));
    // }

    // new TCanvas;
    // fitTotalCopy->SetLineColor(8);
    // printf("kPsiP = %f \n",fitTotal->GetParameter(4)*0.024);
    // fHisto->Draw("");
    // fitTotalCopy->DrawCopy("same");


}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2POL1POL2()
{
  /// Fit using 2 extended crystal balls (signal) + POL1/POL2 (background)

  fHisto->GetListOfFunctions()->Delete();

  TString histoName = fHisto->GetTitle();
  TString sfitOption= histoName.Contains("Corrected") ? "0SERL" : "NO0SRLM";
  const char* fitOption = sfitOption.Data();
  const char* fitOptionBg = "SRL";


  //__________ Get tails parameters, fitting range and SigmaPsiP
  Double_t alphaLow     = GetValue("alJPsi");
  Double_t nLow         = GetValue("nlJPsi");
  Double_t alphaUp      = GetValue("auJPsi");
  Double_t nUp          = GetValue("nuJPsi");
  Double_t fitRangeLow  = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  Double_t paramSPsiP   = GetValue("FSigmaPsiP");
  Double_t meanJPsi     = GetValue("meanJPsi");
  Double_t sigmaJPsi    = GetValue("sigmaJPsi");
  Double_t binNormJPsi  = GetValue("binNormJPsi");
  Double_t binNormPsiP  = GetValue("binNormPsiP");

  TString msg;

  if (IsValidValue(alphaLow)) msg  += TString::Format("alphaLow=%e ",alphaLow);
  if (IsValidValue(nLow))     msg  += TString::Format("nLow=%e ",nLow);
  if (IsValidValue(alphaUp))  msg  += TString::Format("alphaUp=%e ",alphaUp);
  if (IsValidValue(nUp))      msg  += TString::Format("nUp=%e ",nUp);
  //__________


  AliDebug(1,Form("Fit with jpsi + psiprime VWG %s",msg.Data()));


  //__________ Define the function to fit the spectrum, and the background just for plotting
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol1Pol2,fitRangeLow,fitRangeHigh,13,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2Pol1Pol2");

  fitTotal->SetParNames("a","b","a'","b'","c'","kJPsi","mJPsi","sJPsi",
//                        0  1   2    3    4      5      6       7
                        "alJPsi","nlJPsi","auJPsi");
//                         8        9        10
  fitTotal->SetParName(11, "nuJPsi");
  //                          11
  fitTotal->SetParName(12, "kPsiP");
//                            12

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol1Pol2,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol1Pol2");

  //___________


  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol1Pol2,1.7,6.,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol1Pol2");

  Int_t bin = fHisto->FindBin(0.70);
  // Int_t bin = fHisto->FindBin(1.2); // just after a low mass pick
  // Int_t bin = fHisto->FindBin(0.40);

  bckInit->SetParameters(0.,bin,0.,1.,1.);
  bckInit->SetParLimits(0,-200,50);
  bckInit->SetParLimits(2,-50,2000);
  bckInit->SetParLimits(3,-300,100);
  bckInit->FixParameter(4, 1.);

  SetFitRejectRange(2.4,3.2);

  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOptionBg);
  // CheckRoots(fitResultInit,bckInit,2,bckInit->GetParameter(2),bckInit->GetParameter(3),bckInit->GetParameter(4),0.,fitOptionBg);

  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

  //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundPol1Pol2",fitOptionBg);
  //___________


  SetFitRejectRange();

  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 5; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i) + 0.2*bckInit->GetParameter(i));// avoid finding the good parameters
    if(i==4)fitTotal->FixParameter(i, 1.);
  }
  // fitTotal->SetParLimits(0,-200,50);
  // fitTotal->SetParLimits(2,-50,2000);
  // fitTotal->SetParLimits(3,-300,100);

  // bin = fHisto->FindBin(3.09);
  // fitTotal->SetParameter(5, fHisto->GetBinContent(bin)); // norm

  // fitTotal->SetParameter(6, 3.15); // mean
  // fitTotal->SetParLimits(6, 2.95, 3.2);

  // fitTotal->SetParameter(7, 0.08); // sigma
  // fitTotal->SetParLimits(7, 0.05, 0.09);


  if(IsValidValue(binNormJPsi)) bin = fHisto->FindBin(binNormJPsi); //
  else bin = fHisto->FindBin(3.09);
  fitTotal->SetParameter(5, fHisto->GetBinContent(bin)); // norm

  if(IsValidValue(meanJPsi))fitTotal->SetParameter(5,meanJPsi); // mean
  else fitTotal->SetParameter(6, 3.15); // mean

  fitTotal->SetParLimits(6, 2.95, 3.2);

  if(IsValidValue(sigmaJPsi))fitTotal->SetParameter(5,sigmaJPsi);
  else fitTotal->SetParameter(7, 0.08); // sigma
  fitTotal->SetParLimits(7, 0.03, 0.2);


  if ( IsValidValue(alphaLow) )
  {
    fitTotal->FixParameter(8, alphaLow);
  }
  else
  {
    fitTotal->SetParameter(8,0.9);
    fitTotal->SetParLimits(8,0.1,10.0);
  }

  if ( IsValidValue(nLow) )
  {
    fitTotal->FixParameter(9, nLow);
  }
  else
  {
    fitTotal->SetParameter(9,5.0);
    fitTotal->SetParLimits(9,0.0,10.0);
  }

  if ( IsValidValue(alphaUp) )
  {
    fitTotal->FixParameter(10, alphaUp);
  }
  else
  {
    fitTotal->SetParameter(10, 2.0);
    fitTotal->SetParLimits(10,0.1,10.0);
  }

  if ( IsValidValue(nUp) )
  {
    fitTotal->FixParameter(11, nUp);
  }
  else
  {
    fitTotal->SetParameter(11,3.0);
    fitTotal->SetParLimits(11,0.0,10.0);
  }

  // bin = fHisto->FindBin(3.68);
  // fitTotal->SetParameter(12, fHisto->GetBinContent(bin)*0.5); //kPsi'
  // fitTotal->SetParLimits(12, 0.,fHisto->GetBinContent(bin));

  if(IsValidValue(binNormJPsi)) bin = fHisto->FindBin(binNormPsiP); //
  else bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(12, fHisto->GetBinContent(bin)*0.5); //kPsi'
  // fitTotal->SetParLimits(12, 0.,fHisto->GetBinContent(bin));
  //______________


  //_____________First fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
  // CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(2),fitTotal->GetParameter(3),fitTotal->GetParameter(4),0.,fitOption);

  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  //___________


  //___________ Further attempts to fit if the first one fails
  if ( ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 ) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,12,4);
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  // CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(2),fitTotal->GetParameter(3),fitTotal->GetParameter(4),0.,fitOption);
  // if ( static_cast<int>(fitResult) ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,12,4);
  //___________

  delete bckInit; //Delete the initial background funtion

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);

  Set("mJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  Set("sJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));

  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalJPsi->SetParameters(fitTotal->GetParameter(5),
                     fitTotal->GetParameter(6),
                     fitTotal->GetParameter(7),
                     fitTotal->GetParameter(8),
                     fitTotal->GetParameter(9),
                     fitTotal->GetParameter(10),
                     fitTotal->GetParameter(11));

  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalPsiP->SetParameters(fitTotal->GetParameter(12),
                        fitTotal->GetParameter(6) + (3.68609 - 3.096916),
                        fitTotal->GetParameter(7)*paramSPsiP, // /3.096916*3.68609,
                        fitTotal->GetParameter(8),
                        fitTotal->GetParameter(9),
                        fitTotal->GetParameter(10),
                        fitTotal->GetParameter(11));

  bck->SetParameters(fitTotal->GetParameter(0),
                     fitTotal->GetParameter(1),
                     fitTotal->GetParameter(2),
                     fitTotal->GetParameter(3),
                     fitTotal->GetParameter(4));


  Set("a",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("b",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("a'",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("b'",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("c'",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set("kJPsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));

  Set("alJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("nlJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set("auJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("nuJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));

  Set("kPsiP",fitTotal->GetParameter(12),fitTotal->GetParError(12));


  SetFitRejectRange();

  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);


  Double_t cbParameters[7];
  Double_t covarianceMatrix[7][7];

  for ( int ix = 0; ix < 7; ++ix )
  {
    cbParameters[ix] = fitTotal->GetParameter(ix+5);
  }

  for ( int iy = 0; iy < 7; ++iy )
  {
    for ( int ix = 0; ix < 7; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+5,iy+5);
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

  double npsip = signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerrpsip = (fitTotal->GetParError(12)/fitTotal->GetParameter(12))*signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);

  Set("NofPsiP",npsip,nerrpsip);

  double mpsip = GetValue("mJPsi")+ (3.68609 - 3.096916);
  double spsip = GetValue("sJPsi")*paramSPsiP;
  double npsip3s = signalPsiP->Integral(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);
  double nerrpsip3s = (fitTotal->GetParError(12)/fitTotal->GetParameter(12))*signalPsiP->IntegralError(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);

  Set("NofPsiP3s",npsip3s,nerrpsip3s);
  //_____________________________


  //_____Computation of bin significance and signal over background
  Double_t bkgParameters[5];
  Double_t bkgcovarianceMatrix[5][5];

  for ( int ix = 0; ix < 5; ++ix )
  {
    bkgParameters[ix] = fitTotal->GetParameter(ix);
  }

  for ( int iy = 0; iy < 5; ++iy )
  {
    for ( int ix = 0; ix < 5; ++ix )
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
  //__________________________

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2POL2POL3()
{
  /// Fit using 2 extended crystal balls (signal) + Pol2/Pol3 (background)

  fHisto->GetListOfFunctions()->Delete();

  TString histoName = fHisto->GetTitle();
  TString sfitOption= histoName.Contains("Corrected") ? "0SERL" : "NO0SLER";
  const char* fitOption = sfitOption.Data(); //We can add NO to avoid plotting
  const char* fitOptionBg = "0SRL"; //We can add NO to avoid plotting

  //__________ Get tails parameters, fitting range and SigmaPsiP
  Double_t alphaLow     = GetValue("alJPsi");
  Double_t nLow         = GetValue("nlJPsi");
  Double_t alphaUp      = GetValue("auJPsi");
  Double_t nUp          = GetValue("nuJPsi");
  Double_t fitRangeLow  = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  Double_t paramSPsiP   = GetValue("FSigmaPsiP");
  // Double_t meanJPsi     = GetValue("meanJPsi");
  // Double_t sigmaJPsi    = GetValue("sigmaJPsi");

  Double_t a_init    = IsValidValue(GetValue("a_init"))  ? GetValue("a_init")  : -130.;
  Double_t b_init    = IsValidValue(GetValue("b_init"))  ? GetValue("b_init")  : 350.;
  Double_t ap_init   = IsValidValue(GetValue("ap_init")) ? GetValue("ap_init") : -0.05;
  Double_t bp_init   = IsValidValue(GetValue("bp_init")) ? GetValue("bp_init") : 0.5;
  Double_t cp_init   = IsValidValue(GetValue("cp_init")) ? GetValue("cp_init") : -1.;

  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.2;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 3.9;

  TString msg;

  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  if (IsValidValue(nLow))     msg += TString::Format("nLow=%e ",nLow);
  if (IsValidValue(alphaUp))  msg += TString::Format("alphaUp=%e ",alphaUp);
  if (IsValidValue(nUp))      msg += TString::Format("nUp=%e ",nUp);
  msg += TString::Format("ainit=%.2f ",a_init);
  msg += TString::Format("binit=%.2f ",b_init);
  msg += TString::Format("a'init=%.2f ",ap_init);
  msg += TString::Format("b'init=%.2f ",bp_init);
  msg += TString::Format("c'init=%.2f ",cp_init);
  //__________

  AliDebug(1,Form("Fit with jpsi + psiprime POL2/POL3 %s",msg.Data()));

  //__________ Define the function to fit the spectrum, and the background just for plotting
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol2Pol3,fitRangeLow,fitRangeHigh,15,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2Pol2Pol3");

  fitTotal->SetParNames("a","b","c","a'","b'","c'","d'","kJPsi","mJPsi","sJPsi","alJPsi");
//                        0  1   2   3    4    5    6    7         8      9       10
  fitTotal->SetParName(11,"nlJPsi");
  //                         11
  fitTotal->SetParName(12,"auJPsi");
  //                         12
  fitTotal->SetParName(13,"nuJPsi");
  //                         13
  fitTotal->SetParName(14,"kPsiP");
//                            14

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Pol3,fitRangeLow,fitRangeHigh,7,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Pol3");
  //___________


  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Pol3,1.9,6.,7,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Pol3");

  Int_t bin = fHisto->FindBin(0.7);

  bckInit->SetParameters(a_init,b_init,bin,ap_init,bp_init,cp_init);
  bckInit->FixParameter(6.,1);

  // bckInit->SetParLimits(0.,-30,30);
  // bckInit->SetParLimits(1.,-10,10);
  // bckInit->SetParLimits(2.,-10,100);
  bckInit->SetParLimits(3.,-10,10);
  bckInit->SetParLimits(4.,-10,10);
  bckInit->SetParLimits(5.,-30,30);

  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOptionBg);
  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundPol2Pol3",fitOptionBg); // Further attempts to fit bkg if the first one fails
  SetFitRejectRange();
  //____________
  // new TCanvas;
  // fHisto->DrawCopy();
  // return;

  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 7; ++i ) {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
    if(i==6) fitTotal->FixParameter(6, 1.);
  }

  bin = fHisto->FindBin(3.10);
  fitTotal->SetParameter(7, fHisto->GetBinContent(bin)); // norm

  fitTotal->SetParameter(8, 3.15); // mean
  fitTotal->SetParLimits(8, 2.95, 3.2);

  fitTotal->SetParameter(9, 0.08); // sigma
  fitTotal->SetParLimits(9, 0.05, 0.09);

  if ( IsValidValue(alphaLow) ) {
    fitTotal->FixParameter(10, alphaLow);
  } else {
    fitTotal->SetParameter(10,0.9);
    fitTotal->SetParLimits(10,0.1,10.0);
  }

  if ( IsValidValue(nLow) ) {
    fitTotal->FixParameter(11, nLow);
  } else {
    fitTotal->SetParameter(11,5.0);
    fitTotal->SetParLimits(11,0.0,10.0);
  }

  if ( IsValidValue(alphaUp) ) {
    fitTotal->FixParameter(12, alphaUp);
  } else {
    fitTotal->SetParameter(12, 2.0);
    fitTotal->SetParLimits(12,0.1,10.0);
  }

  if ( IsValidValue(nUp) ) {
    fitTotal->FixParameter(13, nUp);
  } else {
    fitTotal->SetParameter(13,3.0);
    fitTotal->SetParLimits(13,0.0,10.0);
  }

  bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(14, fHisto->GetBinContent(bin)*0.5); //kPsi'
  // fitTotal->SetParLimits(14, 0.,1.5*fHisto->GetBinContent(bin));
  //______________

  //_____________Fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  if ( ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 ) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,14,6); // Further attempts to fit if the first one fails
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  //___________

  delete bckInit; //Delete the initial background funtion

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);

  Set("mJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("sJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));

   TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalJPsi->SetParameters(fitTotal->GetParameter(7),
                     fitTotal->GetParameter(8),
                     fitTotal->GetParameter(9),
                     fitTotal->GetParameter(10),
                     fitTotal->GetParameter(11),
                     fitTotal->GetParameter(12),
                     fitTotal->GetParameter(13));

  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),7,"AliAnalysisMuMuJpsiResult","FitFunctionSignalCrystalBallExtended");
  signalPsiP->SetParameters(fitTotal->GetParameter(14),
                        3.68609+(fitTotal->GetParameter(8)-3.096916)/3.096916*3.68609,
                        fitTotal->GetParameter(9)*paramSPsiP, // /3.096916*3.68609,
                        fitTotal->GetParameter(10),
                        fitTotal->GetParameter(11),
                        fitTotal->GetParameter(12),
                        fitTotal->GetParameter(13));

  bck->SetParameters(fitTotal->GetParameter(0),
                     fitTotal->GetParameter(1),
                     fitTotal->GetParameter(2),
                     fitTotal->GetParameter(3),
                     fitTotal->GetParameter(4),
                     fitTotal->GetParameter(5),
                     fitTotal->GetParameter(6));


  Set("a",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("b",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("c",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("a'",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("b'",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set("c'",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set("d'",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  Set("kJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));

  Set("alJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("nlJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));
  Set("auJPsi",fitTotal->GetParameter(12),fitTotal->GetParError(12));
  Set("nuJPsi",fitTotal->GetParameter(13),fitTotal->GetParError(13));

  Set("kPsiP",fitTotal->GetParameter(14),fitTotal->GetParError(14));


  SetFitRejectRange();
  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);

  Double_t cbParameters[7];
  Double_t covarianceMatrix[7][7];

  for ( int ix = 0; ix < 7; ++ix ) cbParameters[ix] = fitTotal->GetParameter(ix+7);

  for ( int iy = 0; iy < 7; ++iy ) {
    for ( int ix = 0; ix < 7; ++ix ) {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+7,iy+7);
    }
  }

  Double_t a   = fHisto->GetXaxis()->GetXmin();
  Double_t b   = fHisto->GetXaxis()->GetXmax();
  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr  = signalJPsi->IntegralError(a,b,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);

  Set("NofJPsi",njpsi,nerr);

  double m       = GetValue("mJPsi");
  double s       = GetValue("sJPsi");
  double njpsi3s = signalJPsi->Integral(m-3*s,m+3*s)/fHisto->GetBinWidth(1);
  double nerr3s  = signalJPsi->IntegralError(m-3*s,m+3*s,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);

  Set("NofJPsi3s",njpsi3s,nerr3s);
  //_____________________________

  //_____Computation of bin significance and signal over background
  Double_t bkgParameters[7];
  Double_t bkgcovarianceMatrix[7][7];

  for ( int ix = 0; ix < 7; ++ix ) bkgParameters[ix] = fitTotal->GetParameter(ix);

  for ( int iy = 0; iy < 7; ++iy ) {
    for ( int ix = 0; ix < 7; ++ix ) {
      bkgcovarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix,iy);
    }
  }

  double nbck3s      = bck->Integral(m-3.*s,m+3.*s)/fHisto->GetBinWidth(1);
  double nbck3sErr   = bck->IntegralError(m-3.*s,m+3.*s,&bkgParameters[0],&bkgcovarianceMatrix[0][0])/fHisto->GetBinWidth(1);

  double sOverB3s    = njpsi3s / nbck3s;
  double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));

  Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);

  double sig    = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
  double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
                               TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );

  Set("Significance3s",sig,sigErr);

  //__________________________

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2POL2POL3V2()
{
  /// Fit using 2 extended crystal balls (signal) + Pol2/Pol3 (background)

  fHisto->GetListOfFunctions()->Delete();

  const char* fitOption = "SERLI"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SRL"; //We can add NO to avoid plotting

  //__________ Get tails parameters, fitting range and SigmaPsiP
  Double_t alphaLow     = GetValue("alJPsi");
  Double_t nLow         = GetValue("nlJPsi");
  Double_t alphaUp      = GetValue("auJPsi");
  Double_t nUp          = GetValue("nuJPsi");
  Double_t fitRangeLow  = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  Double_t paramSPsiP   = GetValue("FSigmaPsiP");

  TString msg;

  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  if (IsValidValue(nLow))     msg += TString::Format("nLow=%e ",nLow);
  if (IsValidValue(alphaUp))  msg += TString::Format("alphaUp=%e ",alphaUp);
  if (IsValidValue(nUp))      msg += TString::Format("nUp=%e ",nUp);
  //__________

  AliDebug(1,Form("Fit with jpsi + psiprime POL2/POL3 %s",msg.Data()));

  //__________ Define the function to fit the spectrum, and the background just for plotting
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol2Pol3V2,fitRangeLow,fitRangeHigh,14,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2Pol2Pol3V2");

  fitTotal->SetParNames("a","b","c","a'","b'","c'","kJPsi","mJPsi","sJPsi","alJPsi");
//                        0  1   2   3    4    5     6       7         8      9
  fitTotal->SetParName(10,"nlJPsi");
  //                         10
  fitTotal->SetParName(11,"auJPsi");
  //                         11
  fitTotal->SetParName(12,"nuJPsi");
  //                         12
  fitTotal->SetParName(13,"kPsiP");
//                            13

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Pol3V2,fitRangeLow,fitRangeHigh,6,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Pol3V2");
  //___________


  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Pol3V2,1.7,8.,6,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Pol3V2");

  Int_t bin = fHisto->FindBin(0.7);
  bckInit->SetParameters(0.,0.,bin,1.,0.,0.);

  SetFitRejectRange(2.2,3.7);
  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOptionBg);
  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundPol2Pol3V2",fitOptionBg); // Further attempts to fit bkg if the first one fails
  SetFitRejectRange();
  //____________

  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 6; ++i ) fitTotal->SetParameter(i, bckInit->GetParameter(i));

  bin = fHisto->FindBin(3.10);
  fitTotal->SetParameter(6, fHisto->GetBinContent(bin)); // norm

  fitTotal->SetParameter(7, 3.15); // mean
  fitTotal->SetParLimits(7, 2.95, 3.2);

  fitTotal->SetParameter(8, 0.088); // sigma
  fitTotal->SetParLimits(8, 0.05, 0.09);

  if ( IsValidValue(alphaLow) ) {
    fitTotal->FixParameter(9, alphaLow);
  } else {
    fitTotal->SetParameter(9,0.9);
    fitTotal->SetParLimits(9,0.1,10.0);
  }

  if ( IsValidValue(nLow) ) {
    fitTotal->FixParameter(10, nLow);
  } else {
    fitTotal->SetParameter(10,5.0);
    fitTotal->SetParLimits(10,0.0,10.0);
  }

  if ( IsValidValue(alphaUp) ) {
    fitTotal->FixParameter(11, alphaUp);
  } else {
    fitTotal->SetParameter(11, 2.0);
    fitTotal->SetParLimits(11,0.1,10.0);
  }

  if ( IsValidValue(nUp) ) {
    fitTotal->FixParameter(12, nUp);
  } else {
    fitTotal->SetParameter(12,3.0);
    fitTotal->SetParLimits(12,0.0,10.0);
  }

  bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(13, fHisto->GetBinContent(bin)*0.5); //kPsi'
  // fitTotal->SetParLimits(13, 0.,1.5*fHisto->GetBinContent(bin));
  //______________

  //_____________Fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;


  if ( static_cast<int>(fitResult) /*||  static_cast<int>(fitResult->CovMatrixStatus())!=3*/ ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,13,5); // Further attempts to fit if the first one fails
  //___________

   new TCanvas;
  fHisto->DrawCopy();
  return;


  delete bckInit; //Delete the initial background funtion

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);

  Set("mJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("sJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));

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

  bck->SetParameters(fitTotal->GetParameter(0),
                     fitTotal->GetParameter(1),
                     fitTotal->GetParameter(2),
                     fitTotal->GetParameter(3),
                     fitTotal->GetParameter(4),
                     fitTotal->GetParameter(5));


  Set("a",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("b",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("c",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("a'",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("b'",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set("c'",fitTotal->GetParameter(5),fitTotal->GetParError(5));

  Set("kJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));

  Set("alJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set("nlJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("auJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));
  Set("nuJPsi",fitTotal->GetParameter(12),fitTotal->GetParError(12));

  Set("kPsiP",fitTotal->GetParameter(13),fitTotal->GetParError(13));


  SetFitRejectRange();
  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);

  Double_t cbParameters[7];
  Double_t covarianceMatrix[7][7];

  for ( int ix = 0; ix < 7; ++ix ) cbParameters[ix] = fitTotal->GetParameter(ix+6);

  for ( int iy = 0; iy < 7; ++iy ) {
    for ( int ix = 0; ix < 7; ++ix ) {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+6,iy+6);
    }
  }

  Double_t a   = fHisto->GetXaxis()->GetXmin();
  Double_t b   = fHisto->GetXaxis()->GetXmax();
  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr  = signalJPsi->IntegralError(a,b,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);

  Set("NofJPsi",njpsi,nerr);

  double m       = GetValue("mJPsi");
  double s       = GetValue("sJPsi");
  double njpsi3s = signalJPsi->Integral(m-3*s,m+3*s)/fHisto->GetBinWidth(1);
  double nerr3s  = signalJPsi->IntegralError(m-3*s,m+3*s,&cbParameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);

  Set("NofJPsi3s",njpsi3s,nerr3s);
  //_____________________________

  //_____Computation of bin significance and signal over background
  Double_t bkgParameters[6];
  Double_t bkgcovarianceMatrix[6][6];

  for ( int ix = 0; ix < 6; ++ix ) bkgParameters[ix] = fitTotal->GetParameter(ix);

  for ( int iy = 0; iy < 6; ++iy ) {
    for ( int ix = 0; ix < 6; ++ix ) {
      bkgcovarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix,iy);
    }
  }

  double nbck3s      = bck->Integral(m-3.*s,m+3.*s)/fHisto->GetBinWidth(1);
  double nbck3sErr   = bck->IntegralError(m-3.*s,m+3.*s,&bkgParameters[0],&bkgcovarianceMatrix[0][0])/fHisto->GetBinWidth(1);

  double sOverB3s    = njpsi3s / nbck3s;
  double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));

  Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);

  double sig    = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
  double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
                               TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );

  Set("Significance3s",sig,sigErr);
  //__________________________
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2VWGINDEPTAILS()
{
  /// Fit using 2 extended crystal balls with independent tails (signal) + variable width gaussian (background)
  // Was used as a test of the Psi' tails effect, since it was negligible is not used anymore

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

  bin = fHisto->FindBin(3.10);
  fitTotal->SetParameter(4, fHisto->GetBinContent(bin)); // norm(kJPsi)

  fitTotal->SetParameter(5, 3.12); // mean
  fitTotal->SetParLimits(5, 3.0, 3.2);

  fitTotal->SetParameter(6, 0.10); // sigma
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
  fitTotal->SetParameter(11, fHisto->GetBinContent(bin)*0.5); //kPsi'
  // fitTotal->SetParLimits(11, fHisto->GetBinContent(bin)*0.01,fHisto->GetBinContent(bin));

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

  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;

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
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  }


  if ( static_cast<int>(fitResult) )
  {
    if ( fitTotal->GetParameter(11) <= fitTotal->GetParError(11) ) //kPsi'
    {
      std::cout << "//-------Refitting again (setting Psi'norm=0)" << std::endl;

      fitTotal->FixParameter(11, 0.);
    }

    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  }

  if ( static_cast<int>(fitResult) )
  {
    std::cout << "//-------Refitting again (setting kVWG=kVWG/2)" << std::endl;

    fitTotal->SetParameter(0, fHisto->GetMaximum()/2.); // kVWG

    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  }

  if ( static_cast<int>(fitResult) )
  {
    std::cout << "//-------Cannot fit properly, try something else..." << std::endl;
  }

  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);

  Set("mJPsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set("sJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));

  Set("mPsiP",fitTotal->GetParameter(5) + (3.68609-3.096916),fitTotal->GetParError(5)/3.096916*3.68609);
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
                            fitTotal->GetParameter(5) + (3.68609-3.096916),
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

  double npsip = signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerrpsip = (fitTotal->GetParError(11)/fitTotal->GetParameter(11))*signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);

  Set("NofPsiP",npsip,nerrpsip);

  double mpsip = GetValue("mJPsi")+ (3.68609 - 3.096916);
  double spsip = GetValue("sJPsi")*paramSPsiP;
  double npsip3s = signalPsiP->Integral(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);
  double nerrpsip3s = (fitTotal->GetParError(11)/fitTotal->GetParameter(11))*signalPsiP->IntegralError(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);

  Set("NofPsiP3s",npsip3s,nerrpsip3s);
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2POL2EXP()
{
  /// Fit using 2 extended crystal balls (signal) + pol2 x exp (background)
  /// 13 parameters

  fHisto->GetListOfFunctions()->Delete();

  const char* fitOption = "SRM"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SR"; //We can add NO to avoid plotting

  //__________ Get tails parameters, fitting range and SigmaPsiP
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
  //__________

  AliDebug(1,Form("Fit with jpsi + psiprime Pol2 x Exp %s",msg.Data()));

  //__________ Define the function to fit the spectrum, and the background just for plotting
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoCB2Pol2Exp,fitRangeLow,fitRangeHigh,12,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoCB2Pol2Exp");

  fitTotal->SetParNames("pol0","pol1","pol2","exp","kJPsi","mJPsi","sJPsi","alJPsi",
  //                      0      1       2     3      4      5       6       7
                        "nlJPsi","auJPsi","nuJPsi");
  //                          8       9        10
  fitTotal->SetParName(11,"kPsiP");
  //                            11

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  //__________


  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,1.7,6.,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");

  Int_t bin = fHisto->FindBin(0.26);

  bckInit->SetParameters(fHisto->GetBinContent(bin),-fHisto->GetBinContent(bin)/3.,100.,0.05);//fHisto->GetBinContent(bin)

  SetFitRejectRange(2.7,4.0);

  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOptionBg);

  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

  //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundPol2Exp",fitOptionBg);
  //___________

  SetFitRejectRange();
  //____________


  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 4; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }

  bin = fHisto->FindBin(3.10);
  fitTotal->SetParameter(4, fHisto->GetBinContent(bin)); // norm(kJPsi)

  fitTotal->SetParameter(5, 3.12);
  fitTotal->SetParLimits(5, 3.07, 3.2);
  fitTotal->SetParameter(6, 0.10);
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

  bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(11, fHisto->GetBinContent(bin)*0.5); //kPsi'
  // fitTotal->SetParLimits(11, fHisto->GetBinContent(bin)*0.01,fHisto->GetBinContent(bin));


  //_____________First fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");

  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  //___________

  //___________ Further attempts to fit if the first one fails
  if ( ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 ) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,11,3);
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  //___________


  delete bckInit; //Delete the initial background funtion

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
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
                            fitTotal->GetParameter(5) + (3.68609-3.096916),
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

  double npsip = signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerrpsip = (fitTotal->GetParError(11)/fitTotal->GetParameter(11))*signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);

  Set("NofPsiP",npsip,nerrpsip);

  double mpsip = GetValue("mJPsi")+ (3.68609 - 3.096916);
  double spsip = GetValue("sJPsi")*paramSPsiP;
  double npsip3s = signalPsiP->Integral(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);
  double nerrpsip3s = (fitTotal->GetParError(11)/fitTotal->GetParameter(11))*signalPsiP->IntegralError(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);

  Set("NofPsiP3s",npsip3s,nerrpsip3s);
  //_____________________________


  //_____Computation of bin significance and signal over background
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
  //__________

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMECB2POL4EXP()
{
  /// Fit using 2 extended crystal balls (signal) + pol4 x exp (background)
  /// 15 parameters
  // Not used in pA Jpsi analysis: too many parameters for correct convergence (Error matrix not pos. def.)

  fHisto->GetListOfFunctions()->Delete();

  //__________ Get tails parameters, fitting range and SigmaPsiP
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
  //__________

  AliDebug(1,Form("Fit with jpsi + psiprime Pol4 x Exp %s",msg.Data()));

  //__________ Define the function to fit the spectrum, and the background just for plotting
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
  //__________


  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Exp,1.6,7.,6,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Exp");

  Int_t bin = fHisto->FindBin(1.6);

  bckInit->SetParameters(fHisto->GetBinContent(bin),-fHisto->GetBinContent(bin),fHisto->GetBinContent(bin)/2.,-fHisto->GetBinContent(bin)/10.,fHisto->GetBinContent(bin)/100.,-2.);

  SetFitRejectRange(2.6,4.0);

  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,"SRL");

  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

  //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundPol4Exp","SRL");
  //___________

  SetFitRejectRange();
  //____________


  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 6; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }

  bin = fHisto->FindBin(3.10);
  fitTotal->SetParameter(6, fHisto->GetBinContent(bin)); // norm(kJPsi)

  fitTotal->SetParameter(7, 3.15); // mean
  fitTotal->SetParLimits(7, 3.0, 3.2);

  fitTotal->SetParameter(8, 0.088); // sigma
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
  fitTotal->SetParameter(13, fHisto->GetBinContent(bin)*0.5); //kPsi'
  // fitTotal->SetParLimits(13, fHisto->GetBinContent(bin)*0.01,fHisto->GetBinContent(bin));

  const char* fitOption = "SER";

  //_____________First fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");

  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  //___________

  //___________ Further attempts to fit if the first one fails
  if ( ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 ) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,13,4);
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  //___________


  delete bckInit; //Delete the initial background funtion


  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
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

  double npsip = signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerrpsip = (fitTotal->GetParError(13)/fitTotal->GetParameter(13))*signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);

  Set("NofPsiP",npsip,nerrpsip);

  double mpsip = GetValue("mJPsi")+ (3.68609 - 3.096916);
  double spsip = GetValue("sJPsi")*paramSPsiP;
  double npsip3s = signalPsiP->Integral(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);
  double nerrpsip3s = (fitTotal->GetParError(13)/fitTotal->GetParameter(13))*signalPsiP->IntegralError(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);

  Set("NofPsiP3s",npsip3s,nerrpsip3s);
  //_____________________________


  //_____Computation of bin significance and signal over background
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
  //_____________________________

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMENA60NEWVWG()
{
  /// Fit using 2 NA60(new) (signal) + variable width gaussian (background)

  fHisto->GetListOfFunctions()->Delete();
  TString histoName = fHisto->GetTitle();
  TString sfitOption= histoName.Contains("Corrected") ? "0SERL" : "NO0SRLM";
  const char* fitOption = sfitOption.Data();
  const char* fitOptionBg = "SER";


  //__________ Get tails parameters, fitting range and SigmaPsiP
  Double_t p1Left       = GetValue("p1LJPsi");
  Double_t p2Left       = GetValue("p2LJPsi");
  Double_t p3Left       = GetValue("p3LJPsi");
  Double_t p1Right      = GetValue("p1RJPsi");
  Double_t p2Right      = GetValue("p2RJPsi");
  Double_t p3Right      = GetValue("p3RJPsi");

  Double_t alphaLeft    = GetValue("aLJPsi");
  Double_t alphaRight   = GetValue("aRJPsi");

  Double_t paramSPsiP   = GetValue("FSigmaPsiP");


  Double_t fitRangeLow  = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  Double_t meanJPsi     = GetValue("meanJPsi");
  Double_t sigmaJPsi    = GetValue("sigmaJPsi");
  Double_t binNormJPsi  = GetValue("binNormJPsi");
  Double_t binNormPsiP  = GetValue("binNormPsiP");

  TString msg;

  if (IsValidValue(p1Left)) msg += TString::Format("p1L=%e ",p1Left);
  if (IsValidValue(p2Left)) msg += TString::Format("p2L=%e ",p2Left);
  if (IsValidValue(p3Left)) msg += TString::Format("p3L=%e ",p3Left);
  if (IsValidValue(p1Right)) msg += TString::Format("p1R=%e ",p1Right);
  if (IsValidValue(p2Right)) msg += TString::Format("p2R=%e ",p2Right);
  if (IsValidValue(p3Right)) msg += TString::Format("p3R=%e ",p3Right);

  if (IsValidValue(alphaLeft)) msg += TString::Format("aL=%e ",alphaLeft);
  if (IsValidValue(alphaRight)) msg += TString::Format("aR=%e ",alphaRight);
  //__________

  AliDebug(1,Form("Fit with jpsi + psiprime NA60 new and VWG %s",msg.Data()));

  //__________ Define the function to fit the spectrum, and the background just for plotting
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

  //__________


  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG,1.7,6.,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG");

  Int_t bin = fHisto->FindBin(0.25);

  bckInit->SetParameters(fHisto->GetBinContent(bin),2.,0.5,0.3);

  SetFitRejectRange(2.2,3.4);

  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOptionBg);

  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

  //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundVWG",fitOptionBg);
  //___________

  SetFitRejectRange();
  //____________


  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 4; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }

  if(IsValidValue(binNormJPsi)) bin = fHisto->FindBin(binNormJPsi); //
  else bin = fHisto->FindBin(3.09);
  fitTotal->SetParameter(4, fHisto->GetBinContent(bin)); // norm

  if(IsValidValue(meanJPsi))fitTotal->SetParameter(5,meanJPsi); // mean
  else fitTotal->SetParameter(5, 3.1); // mean
  fitTotal->SetParLimits(5, 2.95, 3.2);

  if(IsValidValue(sigmaJPsi))fitTotal->SetParameter(6,sigmaJPsi);
  else fitTotal->SetParameter(6, 0.08); // sigma
  fitTotal->SetParLimits(6, 0.03, 0.2);

  fitTotal->FixParameter(7, p1Left);
  fitTotal->FixParameter(8, p2Left);
  fitTotal->FixParameter(9, p3Left);
  fitTotal->FixParameter(10, p1Right);
  fitTotal->FixParameter(11, p2Right);
  fitTotal->FixParameter(12, p3Right);

  fitTotal->FixParameter(13, alphaLeft);
  fitTotal->FixParameter(14, alphaRight);

  // bin = fHisto->FindBin(3.68);
  // fitTotal->SetParameter(15, fHisto->GetBinContent(bin)*0.5); //kPsi'
  // fitTotal->SetParLimits(15, 0.,fHisto->GetBinContent(bin));

  if(IsValidValue(binNormJPsi)) bin = fHisto->FindBin(binNormPsiP); //
  else bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(15, fHisto->GetBinContent(bin)*0.5); //kPsi'
  // fitTotal->SetParLimits(15, 0.,fHisto->GetBinContent(bin));

  //_____________First fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");

  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  //___________

  //___________ Further attempts to fit if the first one fails
  if ( static_cast<int>(fitResult) ||  (!static_cast<int>(fitResult)&&static_cast<int>(fitResult->CovMatrixStatus())!=3) ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,15,3);
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  //___________

  delete bckInit;//Delete the initial background funtion

  //___________Set parameters and fit functions to store in the result
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
                            fitTotal->GetParameter(5) + (3.68609-3.096916),
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

  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
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

  double npsip = signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerrpsip = (fitTotal->GetParError(15)/fitTotal->GetParameter(15))*signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);

  Set("NofPsiP",npsip,nerrpsip);

  double mpsip = GetValue("mJPsi")+ (3.68609-3.096916);
  double spsip = GetValue("sJPsi")*paramSPsiP;
  double npsip3s = signalPsiP->Integral(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);
  double nerrpsip3s = (fitTotal->GetParError(15)/fitTotal->GetParameter(15))*signalPsiP->IntegralError(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);

  Set("NofPsiP3s",npsip3s,nerrpsip3s);

  //_____________________________


  //_____Computation of bin significance and signal over background
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
  //___________________

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMENA60NEWVWG2()
{
  /// Fit using 2 NA60(new) (signal) + variable width gaussian (background)

  fHisto->GetListOfFunctions()->Delete();
  TString histoName = fHisto->GetTitle();
  TString sfitOption= histoName.Contains("Corrected") ? "0SERL" : "NO0SRLM";
  const char* fitOption = sfitOption.Data();
  const char* fitOptionBg = "0SR";


  //__________ Get tails parameters, fitting range and SigmaPsiP
  Double_t p1Left       = GetValue("p1LJPsi");
  Double_t p2Left       = GetValue("p2LJPsi");
  Double_t p3Left       = GetValue("p3LJPsi");
  Double_t p1Right      = GetValue("p1RJPsi");
  Double_t p2Right      = GetValue("p2RJPsi");
  Double_t p3Right      = GetValue("p3RJPsi");

  Double_t alphaLeft    = GetValue("aLJPsi");
  Double_t alphaRight   = GetValue("aRJPsi");

  Double_t paramSPsiP   = GetValue("FSigmaPsiP");
  Double_t meanJPsi     = GetValue("meanJPsi");
  Double_t sigmaJPsi    = GetValue("sigmaJPsi");
  Double_t fitRangeLow  = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  Double_t mVWG2_init    = IsValidValue(GetValue("mVWG2_init"))  ? GetValue("mVWG2_init")  : 2.;
  Double_t s1VWG2_init   = IsValidValue(GetValue("s1VWG2_init")) ? GetValue("s1VWG2_init") : -0.5;
  Double_t s2VWG2_init   = IsValidValue(GetValue("s2VWG2_init")) ? GetValue("s2VWG2_init") : -0.9;
  Double_t gVWG2_init    = IsValidValue(GetValue("gVWG2_init"))  ? GetValue("gVWG2_init")  : 0.15;

  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.2;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 3.9;

  TString msg;

  if (IsValidValue(p1Left)) msg     += TString::Format("p1L=%e ",p1Left);
  if (IsValidValue(p2Left)) msg     += TString::Format("p2L=%e ",p2Left);
  if (IsValidValue(p3Left)) msg     += TString::Format("p3L=%e ",p3Left);
  if (IsValidValue(p1Right)) msg    += TString::Format("p1R=%e ",p1Right);
  if (IsValidValue(p2Right)) msg    += TString::Format("p2R=%e ",p2Right);
  if (IsValidValue(p3Right)) msg    += TString::Format("p3R=%e ",p3Right);

  if (IsValidValue(alphaLeft)) msg  += TString::Format("aL=%e ",alphaLeft);
  if (IsValidValue(alphaRight)) msg += TString::Format("aR=%e ",alphaRight);

  msg += TString::Format("mVWG2init=%.2f ",mVWG2_init);
  msg += TString::Format("s1VWG2init=%.2f ",s1VWG2_init);
  msg += TString::Format("s2VWG2init=%.2f ",s2VWG2_init);
  msg += TString::Format("gVWG2init=%.2f ",gVWG2_init);
  //__________

  AliDebug(1,Form("Fit with jpsi + psiprime NA60 new and VWG %s",msg.Data()));

  //__________ Define the function to fit the spectrum, and the background just for plotting
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewVWG2,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoNA60NewVWG2");

  fitTotal->SetParNames("kVWG","mVWG","s1VWG2","s2VWG2","gVWG2","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi");
  //                       0      1      2       3       4       5       6      7         8         9        10
  fitTotal->SetParName(11, "p1RJPsi");
  //                           11
  fitTotal->SetParName(12, "p2RJPsi");
  //                           12
  fitTotal->SetParName(13, "p3RJPsi");
  //                           13
  fitTotal->SetParName(14, "aLJPsi");
  //                           14
  fitTotal->SetParName(15, "aRJPsi");
  //                           15
  fitTotal->SetParName(16, "kPsiP");
  //                           16

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG2,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG2");
  //__________


  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundVWG2,1.9,5.,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundVWG2");

  Int_t bin = fHisto->FindBin(0.26);

  bckInit->SetParameters(fHisto->GetBinContent(bin),mVWG2_init,s1VWG2_init,s2VWG2_init,gVWG2_init);

  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOptionBg);

  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

  //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundVWG2",fitOptionBg);
  //___________

  SetFitRejectRange();
  //____________

  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 5; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }

  if (IsValidValue(meanJPsi)) fitTotal->FixParameter(6,meanJPsi);
  else{
    fitTotal->SetParameter(6, 3.1); // mean
    fitTotal->SetParLimits(6, 3.0, 3.2);
  }

  if (IsValidValue(sigmaJPsi)) fitTotal->FixParameter(7, sigmaJPsi);
  else{
    fitTotal->SetParameter(7, 0.08); // sigma
    fitTotal->SetParLimits(7, 0.03, 0.2);
  }

  fitTotal->FixParameter(8, p1Left);
  fitTotal->FixParameter(9, p2Left);
  fitTotal->FixParameter(10, p3Left);
  fitTotal->FixParameter(11, p1Right);
  fitTotal->FixParameter(12, p2Right);
  fitTotal->FixParameter(13, p3Right);

  fitTotal->FixParameter(14, alphaLeft);
  fitTotal->FixParameter(15, alphaRight);

  bin = fHisto->FindBin(3.67);
  fitTotal->SetParameter(16, fHisto->GetBinContent(bin)*0.2); //kPsi'
  // fitTotal->SetParLimits(16, 0.,1.1*fHisto->GetBinContent(bin));

  //_____________First fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");

  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  //___________

  //___________ Further attempts to fit if the first one fails
  if ( ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 ) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,16,4);
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  //___________

  delete bckInit;//Delete the initial background funtion

  //___________Set parameters and fit functions to store in the result
  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","FitFunctionNA60New");

  signalJPsi->SetParameters(fitTotal->GetParameter(5),
                            fitTotal->GetParameter(6),
                            fitTotal->GetParameter(7),
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13),
                            fitTotal->GetParameter(14));

  signalJPsi->SetParameter(10,fitTotal->GetParameter(15));

  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","FitFunctionNA60New");

  signalPsiP->SetParameters(fitTotal->GetParameter(16),
                            fitTotal->GetParameter(6) + (3.68609 - 3.096916),
                            fitTotal->GetParameter(7)*paramSPsiP, // /3.096916*3.68609,
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13),
                            fitTotal->GetParameter(14));

  signalPsiP->SetParameter(10,fitTotal->GetParameter(15));

  bck->SetParameter(0,fitTotal->GetParameter(0));
  bck->SetParameter(1,fitTotal->GetParameter(1));
  bck->SetParameter(2,fitTotal->GetParameter(2));
  bck->SetParameter(3,fitTotal->GetParameter(3));
  bck->SetParameter(4,fitTotal->GetParameter(4));

  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);

  Set("kVWG2",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("mVWG2",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("s1VWG2",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("s2VWG2",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("gVWG2",fitTotal->GetParameter(4),fitTotal->GetParError(4));

  Set("kJPsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set("kPsiP",fitTotal->GetParameter(16),fitTotal->GetParError(16));

  Set("mJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  Set("sJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));

  Set("p1LJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("p2LJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set("p3LJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("p1RJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));
  Set("p2RJPsi",fitTotal->GetParameter(12),fitTotal->GetParError(12));
  Set("p3RJPsi",fitTotal->GetParameter(13),fitTotal->GetParError(13));

  Set("aLJPsi",fitTotal->GetParameter(14),fitTotal->GetParError(14));
  Set("aRJPsi",fitTotal->GetParameter(15),fitTotal->GetParError(15));

  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);

  Double_t na60Parameters[11];
  Double_t covarianceMatrix[11][11];

  for ( int ix = 0; ix < 11; ++ix )
  {
    na60Parameters[ix] = fitTotal->GetParameter(ix+5);
  }

  for ( int iy = 0; iy < 11; ++iy )
  {
    for ( int ix = 0; ix < 11; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+5,iy+5);
    }
  }

  Double_t a   = fHisto->GetXaxis()->GetXmin();
  Double_t b   = fHisto->GetXaxis()->GetXmax();
  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr  = signalJPsi->IntegralError(a,b,&na60Parameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);

  Set("NofJPsi",njpsi,nerr);

  double m       = GetValue("mJPsi");
  double s       = GetValue("sJPsi");
  double njpsi3s = signalJPsi->Integral(m-3*s,m+3*s)/fHisto->GetBinWidth(1);
  double nerr3s  = signalJPsi->IntegralError(m-3*s,m+3*s,&na60Parameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);

  Set("NofJPsi3s",njpsi3s,nerr3s);
  //_____________________________


  //_____Computation of bin significance and signal over background
  Double_t bkgParameters[5];
  Double_t bkgcovarianceMatrix[5][5];

  for ( int ix = 0; ix < 5; ++ix )
  {
    bkgParameters[ix] = fitTotal->GetParameter(ix);
  }

  for ( int iy = 0; iy < 5; ++iy )
  {
    for ( int ix = 0; ix < 5; ++ix )
    {
      bkgcovarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix,iy);
    }
  }

  double nbck3s      = bck->Integral(m-3.*s,m+3.*s)/fHisto->GetBinWidth(1);
  double nbck3sErr   = bck->IntegralError(m-3.*s,m+3.*s,&bkgParameters[0],&bkgcovarianceMatrix[0][0])/fHisto->GetBinWidth(1);

  double sOverB3s    = njpsi3s / nbck3s;
  double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));

  Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);

  double sig    = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
  double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
                              TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );

  Set("Significance3s",sig,sigErr);
  //___________________

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMENA60NEWPOL1POL2()
{
  /// Fit using 2 NA60(new) (signal) + POL1/POL2 (background)

  fHisto->GetListOfFunctions()->Delete();

  const char* fitOption = "SERL";
  const char* fitOptionBg = "SERLI";

  //__________ Get tails parameters, fitting range and SigmaPsiP
  Double_t p1Left       = GetValue("p1LJPsi");
  Double_t p2Left       = GetValue("p2LJPsi");
  Double_t p3Left       = GetValue("p3LJPsi");
  Double_t p1Right      = GetValue("p1RJPsi");
  Double_t p2Right      = GetValue("p2RJPsi");
  Double_t p3Right      = GetValue("p3RJPsi");

  Double_t alphaLeft    = GetValue("aLJPsi");
  Double_t alphaRight   = GetValue("aRJPsi");

  Double_t paramSPsiP   = GetValue("FSigmaPsiP");

  Double_t fitRangeLow  = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  Double_t meanJPsi     = GetValue("meanJPsi");
  Double_t sigmaJPsi    = GetValue("sigmaJPsi");
  Double_t binNormJPsi  = GetValue("binNormJPsi");
  Double_t binNormPsiP  = GetValue("binNormPsiP");

  TString msg;

  if (IsValidValue(p1Left)) msg += TString::Format("p1L=%e ",p1Left);
  if (IsValidValue(p2Left)) msg += TString::Format("p2L=%e ",p2Left);
  if (IsValidValue(p3Left)) msg += TString::Format("p3L=%e ",p3Left);
  if (IsValidValue(p1Right)) msg += TString::Format("p1R=%e ",p1Right);
  if (IsValidValue(p2Right)) msg += TString::Format("p2R=%e ",p2Right);
  if (IsValidValue(p3Right)) msg += TString::Format("p3R=%e ",p3Right);

  if (IsValidValue(alphaLeft)) msg += TString::Format("aL=%e ",alphaLeft);
  if (IsValidValue(alphaRight)) msg += TString::Format("aR=%e ",alphaRight);
  //__________

  AliDebug(1,Form("Fit with jpsi + psiprime NA60 new and pol1/pol2 %s",msg.Data()));

  //__________ Define the function to fit the spectrum, and the background just for plotting
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewPol1Pol2,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoNA60NewPol1Pol2");

  fitTotal->SetParNames("a","b","a'","b'","c'","kJPsi","mJPsi","sJPsi",
  //                     0   1   2    3    4      5       6       7
                        "p1LJPsi","p2LJPsi","p3LJPsi");
  //                        8         9         10
  fitTotal->SetParName(11, "p1RJPsi");
  //                           11
  fitTotal->SetParName(12, "p2RJPsi");
  //                           12
  fitTotal->SetParName(13, "p3RJPsi");
  //                           13
  fitTotal->SetParName(14, "aLJPsi");
  //                           14
  fitTotal->SetParName(15, "aRJPsi");
  //                           15
  fitTotal->SetParName(16, "kPsiP");
  //                           16


  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol1Pol2,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol1Pol2");
  //___________


  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol1Pol2,1.7,6.,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol1Pol2");
  Int_t bin = fHisto->FindBin(0.70);


  bckInit->SetParameters(0.,bin,0.,1.,1.);
  bckInit->SetParLimits(0,-200,50);
  bckInit->SetParLimits(2,-50,2000);
  bckInit->SetParLimits(3,-300,100);
  bckInit->FixParameter(4, 1.);

  SetFitRejectRange(2.45,3.4);

  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOptionBg);
  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

  //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundPol1Pol2",fitOptionBg);
  SetFitRejectRange();
  //____________

  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 5; ++i ) {
    fitTotal->SetParameter(i, bckInit->GetParameter(i) /*+0.1*bckInit->GetParameter(i)*/);
    if(i==4)fitTotal->FixParameter(i, 1.);
  }

  // fitTotal->SetParameter(6, 3.1); // mean
  // fitTotal->SetParLimits(6, 3.0, 3.2);

  // fitTotal->SetParameter(7, 0.08); // sigma
  // fitTotal->SetParLimits(7, 0.03, 0.2);

  if(IsValidValue(binNormJPsi)) bin = fHisto->FindBin(binNormJPsi); //
  else bin = fHisto->FindBin(3.09);
  fitTotal->SetParameter(5, fHisto->GetBinContent(bin)); // norm

  if(IsValidValue(meanJPsi))fitTotal->SetParameter(6,meanJPsi); // mean
  else fitTotal->SetParameter(6, 3.1); // mean

  fitTotal->SetParLimits(6, 2.95, 3.2);

  if(IsValidValue(sigmaJPsi))fitTotal->SetParameter(7,sigmaJPsi);
  else fitTotal->SetParameter(7, 0.08); // sigma
  fitTotal->SetParLimits(7, 0.03, 0.2);

  fitTotal->FixParameter(8, p1Left);
  fitTotal->FixParameter(9, p2Left);
  fitTotal->FixParameter(10, p3Left);
  fitTotal->FixParameter(11, p1Right);
  fitTotal->FixParameter(12, p2Right);
  fitTotal->FixParameter(13, p3Right);

  fitTotal->FixParameter(14, alphaLeft);
  fitTotal->FixParameter(15, alphaRight);

  // bin = fHisto->FindBin(3.68);
  // fitTotal->SetParameter(16, fHisto->GetBinContent(bin)*0.5); //kPsi'
  // fitTotal->SetParLimits(16, 0.,fHisto->GetBinContent(bin));

  if(IsValidValue(binNormJPsi)) bin = fHisto->FindBin(binNormPsiP); //
  else bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(16, fHisto->GetBinContent(bin)*0.6); //kPsi'
  // fitTotal->SetParLimits(16, 0.,fHisto->GetBinContent(bin));

  //_____________First fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");

  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  //___________

  //___________ Further attempts to fit if the first one fails
  if ( ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 ) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,16,4);
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  //___________

  delete bckInit;//Delete the initial background funtion

  //___________Set parameters and fit functions to store in the result
  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","FitFunctionNA60New");

  signalJPsi->SetParameters(fitTotal->GetParameter(5),
                            fitTotal->GetParameter(6),
                            fitTotal->GetParameter(7),
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13),
                            fitTotal->GetParameter(14));

  signalJPsi->SetParameter(10,fitTotal->GetParameter(15));

  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","FitFunctionNA60New");

  signalPsiP->SetParameters(fitTotal->GetParameter(16),
                            fitTotal->GetParameter(6) + (3.68609 - 3.096916),
                            fitTotal->GetParameter(7)*paramSPsiP, // /3.096916*3.68609,
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13),
                            fitTotal->GetParameter(14));

  signalPsiP->SetParameter(10,fitTotal->GetParameter(15));

  bck->SetParameter(0,fitTotal->GetParameter(0));
  bck->SetParameter(1,fitTotal->GetParameter(1));
  bck->SetParameter(2,fitTotal->GetParameter(2));
  bck->SetParameter(3,fitTotal->GetParameter(3));
  bck->SetParameter(4,fitTotal->GetParameter(4));

  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);

  Set("a",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("b",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("a'",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("b'",fitTotal->GetParameter(3),fitTotal->GetParError(3));

  Set("kJPsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set("kPsiP",fitTotal->GetParameter(16),fitTotal->GetParError(16));

  Set("mJPsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  Set("sJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));

  Set("p1LJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("p2LJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  Set("p3LJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("p1RJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));
  Set("p2RJPsi",fitTotal->GetParameter(12),fitTotal->GetParError(12));
  Set("p3RJPsi",fitTotal->GetParameter(13),fitTotal->GetParError(13));

  Set("aLJPsi",fitTotal->GetParameter(14),fitTotal->GetParError(14));
  Set("aRJPsi",fitTotal->GetParameter(15),fitTotal->GetParError(15));

  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);

  Double_t na60Parameters[11];
  Double_t covarianceMatrix[11][11];

  for ( int ix = 0; ix < 11; ++ix ) na60Parameters[ix] = fitTotal->GetParameter(ix+5);

  for ( int iy = 0; iy < 11; ++iy )
  {
    for ( int ix = 0; ix < 11; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+5,iy+5);
    }
  }


  Double_t a   = fHisto->GetXaxis()->GetXmin();
  Double_t b   = fHisto->GetXaxis()->GetXmax();
  double njpsi = signalJPsi->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerr  = signalJPsi->IntegralError(a,b,&na60Parameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);

  Set("NofJPsi",njpsi,nerr);

  double m = GetValue("mJPsi");
  double s = GetValue("sJPsi");
  double njpsi3s = signalJPsi->Integral(m-3*s,m+3*s)/fHisto->GetBinWidth(1);
  double nerr3s = signalJPsi->IntegralError(m-3*s,m+3*s,&na60Parameters[0],&covarianceMatrix[0][0])/fHisto->GetBinWidth(1);

  Set("NofJPsi3s",njpsi3s,nerr3s);
  //_____________________________

  double npsip = signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerrpsip  = (fitTotal->GetParError(16)/fitTotal->GetParameter(16))*signalPsiP->IntegralError(a,b)/fHisto->GetBinWidth(1);

  Set("NofPsiP",npsip,nerrpsip);

  double mpsip = GetValue("mJPsi")+ (3.68609 - 3.096916);
  double spsip = GetValue("sJPsi")*paramSPsiP;
  double npsip3s = signalPsiP->Integral(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);
  double nerrpsip3s = (fitTotal->GetParError(16)/fitTotal->GetParameter(16))*signalPsiP->IntegralError(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);

  Set("NofPsiP3s",npsip3s,nerrpsip3s);
  //_____________________________


  //_____Computation of bin significance and signal over background
  Double_t bkgParameters[5];
  Double_t bkgcovarianceMatrix[5][5];

  for ( int ix = 0; ix < 5; ++ix )bkgParameters[ix] = fitTotal->GetParameter(ix);

  for ( int iy = 0; iy < 5; ++iy )
  {
    for ( int ix = 0; ix < 5; ++ix )
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
  //___________________

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMENA60NEWPOL2POL3()
{
  /// Fit using 2 NA60(new) (signal) + POL2/POL3 (background)

  fHisto->GetListOfFunctions()->Delete();

  TString histoName = fHisto->GetTitle();
  TString sfitOption= histoName.Contains("Corrected") ? "0SERL" : "NO0SLER";
  const char* fitOption = sfitOption.Data();
  const char* fitOptionBg = "0SRL";


  //__________ Get tails parameters, fitting range and SigmaPsiP
  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");

  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");

  Double_t paramSPsiP = GetValue("FSigmaPsiP");
  Double_t meanJPsi     = GetValue("meanJPsi");
  Double_t sigmaJPsi    = GetValue("sigmaJPsi");
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  Double_t a_init    = IsValidValue(GetValue("a_init"))  ? GetValue("a_init")  : -130.;
  Double_t b_init    = IsValidValue(GetValue("b_init"))  ? GetValue("b_init")  : 350.;
  Double_t ap_init   = IsValidValue(GetValue("ap_init")) ? GetValue("ap_init") : -0.05;
  Double_t bp_init   = IsValidValue(GetValue("bp_init")) ? GetValue("bp_init") : 0.5;
  Double_t cp_init   = IsValidValue(GetValue("cp_init")) ? GetValue("cp_init") : -1.;

  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.2;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 3.9;

  TString msg;

  if (IsValidValue(p1Left)) msg += TString::Format("p1L=%e ",p1Left);
  if (IsValidValue(p2Left)) msg += TString::Format("p2L=%e ",p2Left);
  if (IsValidValue(p3Left)) msg += TString::Format("p3L=%e ",p3Left);
  if (IsValidValue(p1Right)) msg += TString::Format("p1R=%e ",p1Right);
  if (IsValidValue(p2Right)) msg += TString::Format("p2R=%e ",p2Right);
  if (IsValidValue(p3Right)) msg += TString::Format("p3R=%e ",p3Right);

  if (IsValidValue(alphaLeft)) msg += TString::Format("aL=%e ",alphaLeft);
  if (IsValidValue(alphaRight)) msg += TString::Format("aR=%e ",alphaRight);
  msg += TString::Format("ainit=%.2f",a_init);
  msg += TString::Format("binit=%.2f",b_init);
  msg += TString::Format("a'init=%.2f",ap_init);
  msg += TString::Format("b'init=%.2f",bp_init);
  msg += TString::Format("c'init=%.2f",cp_init);
  //__________

  AliDebug(1,Form("Fit with jpsi + psiprime NA60 new and pol2/pol3 %s",msg.Data()));

  //__________ Define the function to fit the spectrum, and the background just for plotting
  TF1* fitTotal = new TF1("signal+bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionTotalTwoNA60NewPol2Pol3,fitRangeLow,fitRangeHigh,19,"AliAnalysisMuMuJpsiResult","FitFunctionTotalTwoNA60NewPol2Pol3");

  fitTotal->SetParNames("a","b","c","a'","b'","c'","d'","kJPsi","mJPsi","sJPsi",
  //                     0   1   2   3    4    5    6      7      8       9
                        "p1LJPsi");
  //                        10
  fitTotal->SetParName(11, "p2LJPsi");
  //                           11
  fitTotal->SetParName(12, "p3LJPsi");
  //                           12
  fitTotal->SetParName(13, "p1RJPsi");
  //                           13
  fitTotal->SetParName(14, "p2RJPsi");
  //                           14
  fitTotal->SetParName(15, "p3RJPsi");
  //                           15
  fitTotal->SetParName(16, "aLJPsi");
  //                           16
  fitTotal->SetParName(17, "aRJPsi");
  //                           17
  fitTotal->SetParName(18, "kPsiP");
  //                           18

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Pol3,fitRangeLow,fitRangeHigh,7,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Pol3");


  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Pol3,1.9,6.,7,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Pol3");

  Int_t bin = fHisto->FindBin(0.7);

  bckInit->SetParameters(a_init,b_init,bin,ap_init,bp_init,cp_init);
  bckInit->FixParameter(6.,1);

  // bckInit->SetParLimits(0.,-30,30);
  // bckInit->SetParLimits(1.,-10,10);
  // bckInit->SetParLimits(2.,-10,100);
  bckInit->SetParLimits(3.,-10,10);
  bckInit->SetParLimits(4.,-10,10);
  bckInit->SetParLimits(5.,-30,30);

//  bckInit->SetParLimits(0,fHisto->GetBinContent(bin)*0.5,fHisto->GetBinContent(bin)*10);

  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOptionBg);

  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundPol2Pol3",fitOptionBg);

  SetFitRejectRange();
  //____________

  // new TCanvas,
  // fHisto->DrawCopy();
  // return;

  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 7; ++i ) {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
    if(i==6)fitTotal->FixParameter(i, 1.);
  }

  if (IsValidValue(meanJPsi)) fitTotal->FixParameter(8,meanJPsi);
  else{
    fitTotal->SetParameter(8, 3.1); // mean
    fitTotal->SetParLimits(8, 3.0, 3.2);
  }

  if (IsValidValue(sigmaJPsi)) fitTotal->FixParameter(9, sigmaJPsi);
  else{
    fitTotal->SetParameter(9, 0.08); // sigma
    fitTotal->SetParLimits(9, 0.03, 0.2);
  }

  fitTotal->FixParameter(10, p1Left);
  fitTotal->FixParameter(11, p2Left);
  fitTotal->FixParameter(12, p3Left);
  fitTotal->FixParameter(13, p1Right);
  fitTotal->FixParameter(14, p2Right);
  fitTotal->FixParameter(15, p3Right);

  fitTotal->FixParameter(16, alphaLeft);
  fitTotal->FixParameter(17, alphaRight);

  bin = fHisto->FindBin(3.68);
  fitTotal->SetParameter(18, fHisto->GetBinContent(bin)*0.5); //kPsi'
  // fitTotal->SetParLimits(18, fHisto->GetBinContent(bin)*0.1,fHisto->GetBinContent(bin)*1.5);

  //_____________First fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");
  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  //___________

  //___________ Further attempts to fit if the first one fails
  if ( ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 ) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,18,6);
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  delete bckInit;//Delete the initial background funtion

  //___________Set parameters and fit functions to store in the result
  TF1* signalJPsi = new TF1("signalJPsi",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","FitFunctionNA60New");

  signalJPsi->SetParameters(fitTotal->GetParameter(7),
                            fitTotal->GetParameter(8),
                            fitTotal->GetParameter(9),
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13),
                            fitTotal->GetParameter(14),
                            fitTotal->GetParameter(15),
                            fitTotal->GetParameter(16));

  signalJPsi->SetParameter(10,fitTotal->GetParameter(17));

  TF1* signalPsiP = new TF1("signalPsiP",this,&AliAnalysisMuMuJpsiResult::FitFunctionNA60New,fHisto->GetXaxis()->GetXmin(),fHisto->GetXaxis()->GetXmax(),11,"AliAnalysisMuMuJpsiResult","FitFunctionNA60New");

  signalPsiP->SetParameters(fitTotal->GetParameter(18),
                            3.68609+(fitTotal->GetParameter(8)-3.096916)/3.096916*3.68609,
                            fitTotal->GetParameter(9)*paramSPsiP, // /3.096916*3.68609,
                            fitTotal->GetParameter(10),
                            fitTotal->GetParameter(11),
                            fitTotal->GetParameter(12),
                            fitTotal->GetParameter(13),
                            fitTotal->GetParameter(14),
                            fitTotal->GetParameter(15),
                            fitTotal->GetParameter(16));

  signalPsiP->SetParameter(10,fitTotal->GetParameter(17));

  bck->SetParameter(0,fitTotal->GetParameter(0));
  bck->SetParameter(1,fitTotal->GetParameter(1));
  bck->SetParameter(2,fitTotal->GetParameter(2));
  bck->SetParameter(3,fitTotal->GetParameter(3));
  bck->SetParameter(4,fitTotal->GetParameter(4));
  bck->SetParameter(5,fitTotal->GetParameter(5));
  bck->SetParameter(6,fitTotal->GetParameter(6));

  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitTotal->GetChisquare()/fitTotal->GetNDF(),0.0);
  Set("FitNDF",fitTotal->GetNDF(),0.0);

  Set("a",fitTotal->GetParameter(0),fitTotal->GetParError(0));
  Set("b",fitTotal->GetParameter(1),fitTotal->GetParError(1));
  Set("c",fitTotal->GetParameter(2),fitTotal->GetParError(2));
  Set("a'",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("b'",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  Set("c'",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set("d'",fitTotal->GetParameter(6),fitTotal->GetParError(6));

  Set("kJPsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
  Set("kPsiP",fitTotal->GetParameter(18),fitTotal->GetParError(18));

  Set("mJPsi",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("sJPsi",fitTotal->GetParameter(9),fitTotal->GetParError(9));

  Set("p1LJPsi",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("p2LJPsi",fitTotal->GetParameter(11),fitTotal->GetParError(11));
  Set("p3LJPsi",fitTotal->GetParameter(12),fitTotal->GetParError(12));
  Set("p1RJPsi",fitTotal->GetParameter(13),fitTotal->GetParError(13));
  Set("p2RJPsi",fitTotal->GetParameter(14),fitTotal->GetParError(14));
  Set("p3RJPsi",fitTotal->GetParameter(15),fitTotal->GetParError(15));

  Set("aLJPsi",fitTotal->GetParameter(16),fitTotal->GetParError(16));
  Set("aRJPsi",fitTotal->GetParameter(17),fitTotal->GetParError(17));

  AttachFunctionsToHisto(signalJPsi,signalPsiP,bck,fitTotal,fitRangeLow,fitRangeHigh);

  Double_t na60Parameters[11];
  Double_t covarianceMatrix[11][11];

  for ( int ix = 0; ix < 11; ++ix )
  {
    na60Parameters[ix] = fitTotal->GetParameter(ix+7);
  }

  for ( int iy = 0; iy < 11; ++iy )
  {
    for ( int ix = 0; ix < 11; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+7,iy+7);
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

  double npsip = signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerrpsip = (fitTotal->GetParError(18)/fitTotal->GetParameter(18))*signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);

  Set("NofPsiP",npsip,nerrpsip);

  double mpsip = GetValue("mJPsi")+ (3.68609-3.096916);
  double spsip = GetValue("sJPsi")*paramSPsiP;
  double npsip3s = signalPsiP->Integral(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);
  double nerrpsip3s = (fitTotal->GetParError(18)/fitTotal->GetParameter(18))*signalPsiP->IntegralError(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);

  Set("NofPsiP3s",npsip3s,nerrpsip3s);

  //_____________________________


  //_____Computation of bin significance and signal over background
  Double_t bkgParameters[7];
  Double_t bkgcovarianceMatrix[7][7];

  for ( int ix = 0; ix < 7; ++ix )
  {
    bkgParameters[ix] = fitTotal->GetParameter(ix);
  }

  for ( int iy = 0; iy < 7; ++iy )
  {
    for ( int ix = 0; ix < 7; ++ix )
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
  //___________________

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMENA60NEWPOL2EXP()
{
  /// Fit using 2 NA60(new) (signal) + pol2 x exp (background)

  fHisto->GetListOfFunctions()->Delete();

  const char* fitOption = "SRM"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SR"; //We can add NO to avoid plotting

  //__________ Get tails parameters, fitting range and SigmaPsiP
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
  //__________

  AliDebug(1,Form("Fit with jpsi + psiprime NA60 new and pol4 x exp %s",msg.Data()));

  //__________ Define the function to fit the spectrum, and the background just for plotting
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

  //__________


  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,1.7,6.,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");

  Int_t bin = fHisto->FindBin(0.26);

  bckInit->SetParameters(fHisto->GetBinContent(bin),-fHisto->GetBinContent(bin)/3.,100.,0.05);//fHisto->GetBinContent(bin)

  SetFitRejectRange(2.7,4.0);

  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOptionBg);

  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

  //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundPol2Exp",fitOptionBg);
  //___________

  SetFitRejectRange();
  //____________


  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 4; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }

  bin = fHisto->FindBin(3.10);
  fitTotal->SetParameter(4, fHisto->GetBinContent(bin)); // kJPsi

  fitTotal->SetParameter(5, 3.12); // mean
  fitTotal->SetParLimits(5, 3.0, 3.2);

  fitTotal->SetParameter(6, 0.10); // sigma
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
  fitTotal->SetParameter(15, fHisto->GetBinContent(bin)*0.5); //kPsi'
  fitTotal->SetParLimits(15, fHisto->GetBinContent(bin)*0.01,fHisto->GetBinContent(bin));

  //_____________First fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");

  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  //___________

  //___________ Further attempts to fit if the first one fails
  if ( ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 ) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,15,3);
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  //___________

  delete bckInit;//Delete the initial background funtion


  //___________Set parameters and fit functions to store in the result
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
                            fitTotal->GetParameter(5) + (3.68609-3.096916),
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


  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
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

  double npsip = signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerrpsip = (fitTotal->GetParError(15)/fitTotal->GetParameter(15))*signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);

  Set("NofPsiP",npsip,nerrpsip);

  double mpsip = GetValue("mJPsi")+ (3.68609-3.096916);
  double spsip = GetValue("sJPsi")*paramSPsiP;
  double npsip3s = signalPsiP->Integral(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);
  double nerrpsip3s = (fitTotal->GetParError(15)/fitTotal->GetParameter(15))*signalPsiP->IntegralError(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);

  Set("NofPsiP3s",npsip3s,nerrpsip3s);

  //_____________________________


  //_____Computation of bin significance and signal over background
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
  //_____________________________

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitPSIPSIPRIMENA60NEWPOL4EXP()
{
  /// Fit using 2 NA60(new) (signal) + pol4 x exp (background)
  // Not used in pA Jpsi analysis: too many parameters for correct convergence (Error matrix not pos. def.)

  fHisto->GetListOfFunctions()->Delete();

  //__________ Get tails parameters, fitting range and SigmaPsiP
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
  //__________

  AliDebug(1,Form("Fit with jpsi + psiprime NA60 new and pol4 x exp %s",msg.Data()));

  //__________ Define the function to fit the spectrum, and the background just for plotting
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
  //________


  //__________ Fit background only for initial parameters
  TF1* bckInit = new TF1("bckInit",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Exp,1.6,7.,6,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Exp");

  Int_t bin = fHisto->FindBin(1.6);

  bckInit->SetParameters(-fHisto->GetBinContent(bin),fHisto->GetBinContent(bin),-fHisto->GetBinContent(bin)/2.,fHisto->GetBinContent(bin)/10.,fHisto->GetBinContent(bin)/100.,-2.);

  SetFitRejectRange(2.6,4.0);

  TFitResultPtr fitResultInit = fHisto->Fit(bckInit,"SRL");

  std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  //___________

  //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessBkgFit(fitResultInit,bckInit,"FitFunctionBackgroundPol4Exp","SRL");
  //___________

  SetFitRejectRange();
  //____________


  //__________ Set initial parameters in fitting function
  for ( Int_t i = 0; i < 6; ++i )
  {
    fitTotal->SetParameter(i, bckInit->GetParameter(i));
  }

  bin = fHisto->FindBin(3.10);
  fitTotal->SetParameter(6, fHisto->GetBinContent(bin)); // kJPsi

  fitTotal->SetParameter(7, 3.15); // mean
  fitTotal->SetParLimits(7, 3.0, 3.2);

  fitTotal->SetParameter(8, 0.088); // sigma
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
  // fitTotal->SetParameter(17, fHisto->GetBinContent(bin)*0.5); //kPsi'
  fitTotal->SetParLimits(17, fHisto->GetBinContent(bin)*0.01,fHisto->GetBinContent(bin));

  const char* fitOption = "SER";

  //_____________First fit attempt
  TFitResultPtr fitResult = fHisto->Fit(fitTotal,fitOption,"");

  std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  //___________

  //___________ Further attempts to fit if the first one fails
  if ( ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 ) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMinvFit(fitResult,fitTotal,bckInit,fitOption,17,4);
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  //___________

  delete bckInit;//Delete the initial background funtion


  //___________Set parameters and fit functions to store in the result
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


  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
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

  double npsip = signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);
  double nerrpsip = (fitTotal->GetParError(17)/fitTotal->GetParameter(17))*signalPsiP->Integral(a,b)/fHisto->GetBinWidth(1);

  Set("NofPsiP",npsip,nerrpsip);

  double mpsip = GetValue("mJPsi")+ (3.68609-3.096916);
  double spsip = GetValue("sJPsi")*paramSPsiP;
  double npsip3s = signalPsiP->Integral(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);
  double nerrpsip3s = (fitTotal->GetParError(17)/fitTotal->GetParameter(17))*signalPsiP->IntegralError(mpsip-3*spsip,mpsip+3*spsip)/fHisto->GetBinWidth(1);

  Set("NofPsiP3s",npsip3s,nerrpsip3s);

  //_____________________________


  //_____Computation of bin significance and signal over background
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
  //_____________________________

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2VWG_BKGMPTPOL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow       = GetValue("alJPsi");
  Double_t nLow           = GetValue("nlJPsi");
  Double_t alphaUp        = GetValue("auJPsi");
  Double_t nUp            = GetValue("nuJPsi");

  Double_t kVWG           = GetValue("kVWG");
  Double_t mVWG           = GetValue("mVWG");
  Double_t sVWG1          = GetValue("sVWG1");
  Double_t sVWG2          = GetValue("sVWG2");
  Double_t kJPsi          = GetValue("kJPsi");
  Double_t kPsiP          = GetValue("kPsiP");
  Double_t mJPsi          = GetValue("mJPsi");
  Double_t sJPsi          = GetValue("sJPsi");
  Double_t NofJPsi        = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate(); // Recalculates the error for low stat bins

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(3.,1.,0.);

  bck->SetParLimits(0, 0.,5.0);

  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);

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
  fitMeanpt->SetParLimits(12, 1.0,10.);

  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }


  Double_t psipPtLim = 10.;

  fitMeanpt->SetParameter(16, 3.);
  fitMeanpt->SetParLimits(16, 0.,psipPtLim);


  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(16) <= fitMeanpt->GetParError(16) || fitMeanpt->GetParError(16) >= 0.75*fitMeanpt->GetParameter(16) || (fitMeanpt->GetParameter(16)/psipPtLim > 0.9)) )
  {
    fitMeanpt->FixParameter(16, bck->Eval(3.68));
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
    if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(16) <= fitMeanpt->GetParError(16) || fitMeanpt->GetParError(16) >= 0.75*fitMeanpt->GetParameter(16) || (fitMeanpt->GetParameter(16)/psipPtLim > 0.9)) )
      printf("--------> Warning : problem with error estimation for back parameters <-------\n");
  }

  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(12),fitMeanpt->GetParError(12));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2POL1POL2_BKGMPTPOL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");

  Double_t a       = GetValue("a");
  Double_t b       = GetValue("b");
  Double_t aprime  = GetValue("a'");
  Double_t bprime  = GetValue("b'");
  Double_t cprime  = GetValue("c'");
  Double_t kJPsi   = GetValue("kJPsi");
  Double_t kPsiP   = GetValue("kPsiP");
  Double_t mJPsi   = GetValue("mJPsi");
  Double_t sJPsi   = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate(); // Recalculates the error for low stat bins

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(3.,1.,0.);

  bck->SetParLimits(0, 0.,5.0);

  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);

  SetFitRejectRange();


  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2POL1POL2POL2,fitRangeLow,fitRangeHigh,18,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2POL1POL2POL2");

  fitMeanpt->SetParNames("a","b","a'","b'","c'","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi");
  fitMeanpt->SetParName(11,"nuJPsi");
  fitMeanpt->SetParName(12,"kPsiP");
  fitMeanpt->SetParName(13,"<pt>JPsi");
  fitMeanpt->SetParName(14,"<pt>BG0");
  fitMeanpt->SetParName(15,"<pt>BG1");
  fitMeanpt->SetParName(16,"<pt>BG2");
  fitMeanpt->SetParName(17,"<pt>PsiP");

  fitMeanpt->FixParameter(0,a);
  fitMeanpt->FixParameter(1,b);
  fitMeanpt->FixParameter(2,aprime);
  fitMeanpt->FixParameter(3,bprime);
  fitMeanpt->FixParameter(4,cprime);
  fitMeanpt->FixParameter(5,kJPsi);
  fitMeanpt->FixParameter(6,mJPsi);
  fitMeanpt->FixParameter(7,sJPsi);
  fitMeanpt->FixParameter(8,alphaLow);
  fitMeanpt->FixParameter(9,nLow);
  fitMeanpt->FixParameter(10,alphaUp);
  fitMeanpt->FixParameter(11,nUp);
  fitMeanpt->FixParameter(12,kPsiP);

  fitMeanpt->SetParameter(13, 3.);
  fitMeanpt->SetParLimits(13, 1.0,10.);

  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanpt->SetParameter(i + 14, bck->GetParameter(i));
  }


  Double_t psipPtLim = 10.;

  fitMeanpt->SetParameter(17, 3.);
  fitMeanpt->SetParLimits(17, 0.,psipPtLim);


  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 14, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(17) <= fitMeanpt->GetParError(17) || fitMeanpt->GetParError(17) >= 0.75*fitMeanpt->GetParameter(17) || (fitMeanpt->GetParameter(17)/psipPtLim > 0.9)) )
  {
    fitMeanpt->FixParameter(17, bck->Eval(3.68));
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 14, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  bck->SetParameters(fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15),fitMeanpt->GetParameter(16));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(13),fitMeanpt->GetParError(13));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(17),fitMeanpt->GetParError(17));
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2VWG_BKGMPTPOL2EXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow       = GetValue("alJPsi");
  Double_t nLow           = GetValue("nlJPsi");
  Double_t alphaUp        = GetValue("auJPsi");
  Double_t nUp            = GetValue("nuJPsi");

  Double_t kVWG           = GetValue("kVWG");
  Double_t mVWG           = GetValue("mVWG");
  Double_t sVWG1          = GetValue("sVWG1");
  Double_t sVWG2          = GetValue("sVWG2");
  Double_t kJPsi          = GetValue("kJPsi");
  Double_t kPsiP          = GetValue("kPsiP");
  Double_t mJPsi          = GetValue("mJPsi");
  Double_t sJPsi          = GetValue("sJPsi");
  Double_t NofJPsi        = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");

  Double_t fitRangeLow    = GetValue(kFitRangeLow);
  Double_t fitRangeHigh   = GetValue(kFitRangeHigh);

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();

  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");

  bck->SetParameters(3.,-2.,0.4,-0.0);

  bck->SetParLimits(0, 0.,5.0);

  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);

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
  fitMeanpt->SetParLimits(12, 1.0,10.);

  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }

  fitMeanpt->SetParameter(17, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(17, 0.,psipPtLim);

  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult= " << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(17) <= fitMeanpt->GetParError(17) || fitMeanpt->GetParError(17) >= 0.75*fitMeanpt->GetParameter(17) || (fitMeanpt->GetParameter(17)/psipPtLim > 0.9) )
  {
    fitMeanpt->FixParameter(17, bck->Eval(3.68));
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(13) <= fitMeanpt->GetParError(13) || fitMeanpt->GetParError(13) >= 0.75*fitMeanpt->GetParameter(13)) )
  {
    fitMeanpt->SetParameter(13, 2.);

    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
    if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(13) <= fitMeanpt->GetParError(13) || fitMeanpt->GetParError(13) >= 0.75*fitMeanpt->GetParameter(13)) )
      printf("--------> Warning : problem with error estimation for back parameters <-------\n");
  }

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15),fitMeanpt->GetParameter(16));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(12),fitMeanpt->GetParError(12));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(17),fitMeanpt->GetParError(17));

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2POL1POL2_BKGMPTPOL2EXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow       = GetValue("alJPsi");
  Double_t nLow           = GetValue("nlJPsi");
  Double_t alphaUp        = GetValue("auJPsi");
  Double_t nUp            = GetValue("nuJPsi");

  Double_t a              = GetValue("a");
  Double_t b              = GetValue("b");
  Double_t aprime         = GetValue("a'");
  Double_t bprime         = GetValue("b'");
  Double_t cprime         = GetValue("c'");
  Double_t kJPsi          = GetValue("kJPsi");
  Double_t kPsiP          = GetValue("kPsiP");
  Double_t mJPsi          = GetValue("mJPsi");
  Double_t sJPsi          = GetValue("sJPsi");
  Double_t NofJPsi        = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");

  Double_t fitRangeLow    = GetValue(kFitRangeLow);
  Double_t fitRangeHigh   = GetValue(kFitRangeHigh);

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();

  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");

  bck->SetParameters(3.,-2.,0.4,-0.0);

  bck->SetParLimits(0, 0.,5.0);

  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SR","",fitRangeLow,fitRangeHigh);

  SetFitRejectRange();

  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2POL1POL2POL2EXP,fitRangeLow,fitRangeHigh,19,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2POL1POL2POL2EXP");

  fitMeanpt->SetParNames("a","b","a'","b'","c'","kJPsi","mJPsi","sJPsi","alJPsi","auJPsi","nlJPsi");
  fitMeanpt->SetParName(11,"nuJPsi");
  fitMeanpt->SetParName(12,"kPsiP");
  fitMeanpt->SetParName(13,"<pt>JPsi");
  fitMeanpt->SetParName(14,"<pt>BG0");
  fitMeanpt->SetParName(15,"<pt>BG1");
  fitMeanpt->SetParName(16,"<pt>BG2");
  fitMeanpt->SetParName(17,"<pt>BGEXP");
  fitMeanpt->SetParName(18,"<pt>PsiP");

  fitMeanpt->FixParameter(0,a);
  fitMeanpt->FixParameter(1,b);
  fitMeanpt->FixParameter(2,aprime);
  fitMeanpt->FixParameter(3,bprime);
  fitMeanpt->FixParameter(4,cprime);
  fitMeanpt->FixParameter(5,kJPsi);
  fitMeanpt->FixParameter(6,mJPsi);
  fitMeanpt->FixParameter(7,sJPsi);
  fitMeanpt->FixParameter(8,alphaLow);
  fitMeanpt->FixParameter(9,nLow);
  fitMeanpt->FixParameter(10,alphaUp);
  fitMeanpt->FixParameter(11,nUp);
  fitMeanpt->FixParameter(12,kPsiP);

  fitMeanpt->SetParameter(13, 3.);
  fitMeanpt->SetParLimits(13, 1.0,10.);

  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanpt->SetParameter(i + 14, bck->GetParameter(i));
  }

  fitMeanpt->SetParameter(18, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(18, 0.,psipPtLim);

  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult= " << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 14, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(18) <= fitMeanpt->GetParError(18) || fitMeanpt->GetParError(18) >= 0.75*fitMeanpt->GetParameter(18) || (fitMeanpt->GetParameter(18)/psipPtLim > 0.9) )
  {
    fitMeanpt->FixParameter(18, bck->Eval(3.68));
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 14, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(14) <= fitMeanpt->GetParError(14) || fitMeanpt->GetParError(14) >= 0.75*fitMeanpt->GetParameter(14)) )
  {
    fitMeanpt->SetParameter(14, 2.);

    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15),fitMeanpt->GetParameter(16),fitMeanpt->GetParameter(17));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(13),fitMeanpt->GetParError(13));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(18),fitMeanpt->GetParError(18));

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


  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(3.,1.,0.);

  bck->SetParLimits(0, 0.,5.0);


  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);

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
  fitMeanpt->SetParLimits(12, 1.0,10.);

  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }

  fitMeanpt->SetParameter(16, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(16, 0.,psipPtLim);

  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(16) <= fitMeanpt->GetParError(16) || fitMeanpt->GetParError(16) >= 0.75*fitMeanpt->GetParameter(16) || (fitMeanpt->GetParameter(16)/psipPtLim > 0.9)) )
  {
    fitMeanpt->FixParameter(16, bck->Eval(3.68));
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
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


  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");

  bck->SetParameters(3.,-1.,0.4,-0.1);
  bck->SetParLimits(0, 0.,5.0);


  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);

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
  fitMeanpt->SetParLimits(12, 1.0,10.);

  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }

  fitMeanpt->SetParameter(17, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(17, 0.,psipPtLim);


  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(17) <= fitMeanpt->GetParError(17) || fitMeanpt->GetParError(17) >= 0.75*fitMeanpt->GetParameter(17) || (fitMeanpt->GetParameter(17)/psipPtLim > 0.9)) )
  {
    fitMeanpt->FixParameter(17, bck->Eval(3.68));
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(13) <= fitMeanpt->GetParError(13) || fitMeanpt->GetParError(13) >= 0.75*fitMeanpt->GetParameter(13)) )
  {
    fitMeanpt->SetParameter(13, 2.);

    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15),fitMeanpt->GetParameter(16));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(12),fitMeanpt->GetParError(12));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(17),fitMeanpt->GetParError(17));

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();

  Double_t p1Left         = GetValue("p1LJPsi");
  Double_t p2Left         = GetValue("p2LJPsi");
  Double_t p3Left         = GetValue("p3LJPsi");
  Double_t p1Right        = GetValue("p1RJPsi");
  Double_t p2Right        = GetValue("p2RJPsi");
  Double_t p3Right        = GetValue("p3RJPsi");

  Double_t alphaLeft      = GetValue("aLJPsi");
  Double_t alphaRight     = GetValue("aRJPsi");

  Double_t kVWG           = GetValue("kVWG");
  Double_t mVWG           = GetValue("mVWG");
  Double_t sVWG1          = GetValue("sVWG1");
  Double_t sVWG2          = GetValue("sVWG2");
  Double_t kJPsi          = GetValue("kJPsi");
  Double_t kPsiP          = GetValue("kPsiP");
  Double_t mJPsi          = GetValue("mJPsi");
  Double_t sJPsi          = GetValue("sJPsi");
  Double_t NofJPsi        = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");

  Double_t fitRangeLow    = GetValue(kFitRangeLow);
  Double_t fitRangeHigh   = GetValue(kFitRangeHigh);


  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(3.,1.,0.);
  bck->SetParLimits(0, 0.,5.0);


  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);

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
  fitMeanpt->SetParLimits(16, 1.0,10.);


  for ( Int_t i = 0; i < 3; ++i )
  {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i));
  }

  fitMeanpt->SetParameter(20, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(20, 0.,psipPtLim);


  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( static_cast<int>(fitResult)&& (fitMeanpt->GetParameter(20) <= fitMeanpt->GetParError(20) || fitMeanpt->GetParError(20) >= 0.75*fitMeanpt->GetParameter(20) || (fitMeanpt->GetParameter(20)/psipPtLim > 0.9)) )
  {
    fitMeanpt->FixParameter(20, bck->Eval(3.68));
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
    if ( static_cast<int>(fitResult)&& (fitMeanpt->GetParameter(20) <= fitMeanpt->GetParError(20) || fitMeanpt->GetParError(20) >= 0.75*fitMeanpt->GetParameter(20) || (fitMeanpt->GetParameter(20)/psipPtLim > 0.9)) )
      printf("--------> Warning : problem with <pt>PsiP error estimation <-------\n");
  }

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanpt->GetParameter(17),fitMeanpt->GetParameter(18),fitMeanpt->GetParameter(19));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(20),fitMeanpt->GetParError(20));

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMENA60NEWPOL1POL2_BKGMPTPOL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();

  Double_t p1Left         = GetValue("p1LJPsi");
  Double_t p2Left         = GetValue("p2LJPsi");
  Double_t p3Left         = GetValue("p3LJPsi");
  Double_t p1Right        = GetValue("p1RJPsi");
  Double_t p2Right        = GetValue("p2RJPsi");
  Double_t p3Right        = GetValue("p3RJPsi");

  Double_t alphaLeft      = GetValue("aLJPsi");
  Double_t alphaRight     = GetValue("aRJPsi");

  Double_t a              = GetValue("a");
  Double_t b              = GetValue("b");
  Double_t aprime         = GetValue("a'");
  Double_t bprime         = GetValue("b'");
  Double_t cprime         = GetValue("c'");
  Double_t kJPsi          = GetValue("kJPsi");
  Double_t kPsiP          = GetValue("kPsiP");
  Double_t mJPsi          = GetValue("mJPsi");
  Double_t sJPsi          = GetValue("sJPsi");
  Double_t NofJPsi        = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");

  Double_t fitRangeLow    = GetValue(kFitRangeLow);
  Double_t fitRangeHigh   = GetValue(kFitRangeHigh);


  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();

  // --- fit the background ---

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(2.,1.,0.);
  bck->SetParLimits(0, 0.,5.0);

  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);

  SetFitRejectRange();

  //---



  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2NA60NEWPOL1POL2POL2,fitRangeLow,fitRangeHigh,22,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2NA60NEWPOL1POL2POL2");

  fitMeanpt->SetParNames("a","b","a'","b'","c'","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi");
  fitMeanpt->SetParName(11,"p1RJPsi");
  fitMeanpt->SetParName(12,"p2RJPsi");
  fitMeanpt->SetParName(13,"p3RJPsi");
  fitMeanpt->SetParName(14,"aLJPsi");
  fitMeanpt->SetParName(15,"aRJPsi");
  fitMeanpt->SetParName(16,"kPsiP");
  fitMeanpt->SetParName(17,"<pt>JPsi");
  fitMeanpt->SetParName(18,"<pt>BG0");
  fitMeanpt->SetParName(19,"<pt>BG1");
  fitMeanpt->SetParName(20,"<pt>BG2");
  fitMeanpt->SetParName(21,"<pt>PsiP");

  fitMeanpt->FixParameter(0,a);
  fitMeanpt->FixParameter(1,b);
  fitMeanpt->FixParameter(2,aprime);
  fitMeanpt->FixParameter(3,bprime);
  fitMeanpt->FixParameter(4,cprime);
  fitMeanpt->FixParameter(5,kJPsi);
  fitMeanpt->FixParameter(6,mJPsi);
  fitMeanpt->FixParameter(7,sJPsi);
  fitMeanpt->FixParameter(8,p1Left);
  fitMeanpt->FixParameter(9,p2Left);
  fitMeanpt->FixParameter(10,p3Left);
  fitMeanpt->FixParameter(11,p1Right);
  fitMeanpt->FixParameter(12,p2Right);
  fitMeanpt->FixParameter(13,p3Right);
  fitMeanpt->FixParameter(14,alphaLeft);
  fitMeanpt->FixParameter(15,alphaRight);
  fitMeanpt->FixParameter(16,kPsiP);

  fitMeanpt->SetParameter(17, 3.);
  fitMeanpt->SetParLimits(17, 1.0,10.);


  for ( Int_t i = 0; i < 3; ++i )
  {
      fitMeanpt->SetParameter(i + 18, bck->GetParameter(i));
  }

  fitMeanpt->SetParameter(21, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(21, 0.,psipPtLim);


  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 18, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( static_cast<int>(fitResult)&& (fitMeanpt->GetParameter(21) <= fitMeanpt->GetParError(21) || fitMeanpt->GetParError(21) >= 0.75*fitMeanpt->GetParameter(21) || (fitMeanpt->GetParameter(21)/psipPtLim > 0.9)) )
  {
    fitMeanpt->FixParameter(21, bck->Eval(3.68));
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 18, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanpt->GetParameter(18),fitMeanpt->GetParameter(19),fitMeanpt->GetParameter(20));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(17),fitMeanpt->GetParError(17));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(21),fitMeanpt->GetParError(21));

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


  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");

  bck->SetParameters(2.5,-1.,0.2,-0.05);
  bck->SetParLimits(0, 0.,5.0);


  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);

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
  fitMeanpt->SetParLimits(16, 1.0,10.);

  for ( Int_t i = 0; i < 4; ++i )
  {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i));
  }

  fitMeanpt->SetParameter(21, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(21, 0.,psipPtLim);


  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( static_cast<int>(fitResult)&& (fitMeanpt->GetParameter(21) <= fitMeanpt->GetParError(21) || fitMeanpt->GetParError(21) >= 0.75*fitMeanpt->GetParameter(21) || (fitMeanpt->GetParameter(21)/psipPtLim > 0.9)) )
  {
    fitMeanpt->FixParameter(21, bck->Eval(3.68));
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(17) <= fitMeanpt->GetParError(17) || fitMeanpt->GetParError(17) >= 0.75*fitMeanpt->GetParameter(13)) )
  {
    fitMeanpt->SetParameter(17, 2.);

    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
    if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(17) <= fitMeanpt->GetParError(17) || fitMeanpt->GetParError(17) >= 0.75*fitMeanpt->GetParameter(13)) )
      printf("--------> Warning : problem with error estimation for back parameters <-------\n");

  }


  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanpt->GetParameter(17),fitMeanpt->GetParameter(18),fitMeanpt->GetParameter(19),fitMeanpt->GetParameter(20));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(21),fitMeanpt->GetParError(21));

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMENA60NEWPOL1POL2_BKGMPTPOL2EXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();

  Double_t p1Left         = GetValue("p1LJPsi");
  Double_t p2Left         = GetValue("p2LJPsi");
  Double_t p3Left         = GetValue("p3LJPsi");
  Double_t p1Right        = GetValue("p1RJPsi");
  Double_t p2Right        = GetValue("p2RJPsi");
  Double_t p3Right        = GetValue("p3RJPsi");

  Double_t alphaLeft      = GetValue("aLJPsi");
  Double_t alphaRight     = GetValue("aRJPsi");

  Double_t a              = GetValue("a");
  Double_t b              = GetValue("b");
  Double_t aprime         = GetValue("a'");
  Double_t bprime         = GetValue("b'");
  Double_t cprime         = GetValue("c'");
  Double_t kJPsi          = GetValue("kJPsi");
  Double_t kPsiP          = GetValue("kPsiP");
  Double_t mJPsi          = GetValue("mJPsi");
  Double_t sJPsi          = GetValue("sJPsi");
  Double_t NofJPsi        = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetErrorStat("NofJPsi");

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);


  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");

  bck->SetParameters(2.5,-1.,0.2,-0.05);
  bck->SetParLimits(0, 0.,5.0);


  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SR","",fitRangeLow,fitRangeHigh);

  SetFitRejectRange();

  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2NA60NEWPOL1POL2POL2EXP,fitRangeLow,fitRangeHigh,23,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2NA60NEWPOL1POL2POL2EXP");

  fitMeanpt->SetParNames("a","b","a'","b'","c'","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi");
  fitMeanpt->SetParName(11,"p1RJPsi");
  fitMeanpt->SetParName(12,"p2RJPsi");
  fitMeanpt->SetParName(13,"p3RJPsi");
  fitMeanpt->SetParName(14,"aLJPsi");
  fitMeanpt->SetParName(15,"aRJPsi");
  fitMeanpt->SetParName(16,"kPsiP");
  fitMeanpt->SetParName(17,"<pt>JPsi");
  fitMeanpt->SetParName(18,"<pt>BG0");
  fitMeanpt->SetParName(19,"<pt>BG1");
  fitMeanpt->SetParName(20,"<pt>BG2");
  fitMeanpt->SetParName(21,"<pt>BGEXP");
  fitMeanpt->SetParName(22,"<pt>PsiP");

  fitMeanpt->FixParameter(0,a);
  fitMeanpt->FixParameter(1,b);
  fitMeanpt->FixParameter(2,aprime);
  fitMeanpt->FixParameter(3,bprime);
  fitMeanpt->FixParameter(4,cprime);
  fitMeanpt->FixParameter(5,kJPsi);
  fitMeanpt->FixParameter(6,mJPsi);
  fitMeanpt->FixParameter(7,sJPsi);
  fitMeanpt->FixParameter(8,p1Left);
  fitMeanpt->FixParameter(9,p2Left);
  fitMeanpt->FixParameter(10,p3Left);
  fitMeanpt->FixParameter(11,p1Right);
  fitMeanpt->FixParameter(12,p2Right);
  fitMeanpt->FixParameter(13,p3Right);
  fitMeanpt->FixParameter(14,alphaLeft);
  fitMeanpt->FixParameter(15,alphaRight);
  fitMeanpt->FixParameter(16,kPsiP);

  fitMeanpt->SetParameter(17, 3.);
  fitMeanpt->SetParLimits(17, 1.0,10.);

  for ( Int_t i = 0; i < 4; ++i )
  {
      fitMeanpt->SetParameter(i + 18, bck->GetParameter(i));
  }

  fitMeanpt->SetParameter(22, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(22, 0.,psipPtLim);


  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 18, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( static_cast<int>(fitResult)&& (fitMeanpt->GetParameter(22) <= fitMeanpt->GetParError(22) || fitMeanpt->GetParError(22) >= 0.75*fitMeanpt->GetParameter(22) || (fitMeanpt->GetParameter(22)/psipPtLim > 0.9)) )
  {
    fitMeanpt->FixParameter(22, bck->Eval(3.68));
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 18, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(18) <= fitMeanpt->GetParError(18) || fitMeanpt->GetParError(18) >= 0.75*fitMeanpt->GetParameter(14)) )
  {
    fitMeanpt->SetParameter(17, 2.);

    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }


  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanpt->GetParameter(18),fitMeanpt->GetParameter(19),fitMeanpt->GetParameter(20),fitMeanpt->GetParameter(21));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(17),fitMeanpt->GetParError(17));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(22),fitMeanpt->GetParError(22));

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


  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(3.,1.,0.);
  bck->SetParLimits(0, 0.,5.0);


  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);

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
  fitMeanpt->SetParLimits(16, 1.0,10.);

  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanpt->SetParameter(i + 17, bck->GetParameter(i));
  }

  fitMeanpt->SetParameter(20, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(20, 0.,psipPtLim);

  //  TProfile::Approximate();

  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(20) <= fitMeanpt->GetParError(20) || fitMeanpt->GetParError(20) >= 0.75*fitMeanpt->GetParameter(20) || (fitMeanpt->GetParameter(20)/psipPtLim > 0.9)) )
  {
    fitMeanpt->FixParameter(20, bck->Eval(3.68));
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanpt->GetParameter(17),fitMeanpt->GetParameter(18),fitMeanpt->GetParameter(19));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
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

  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");

  bck->SetParameters(2.8,-0.5,0.2,0.05);
  bck->SetParLimits(0, 0.,5.0);


  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);

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
  fitMeanpt->SetParLimits(16, 1.0,10.);

  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanpt->SetParameter(i + 17, bck->GetParameter(i));
  }

  fitMeanpt->SetParameter(21, 3.);
  Double_t psipPtLim = 10.;
  fitMeanpt->SetParLimits(21, 0.,psipPtLim);


  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(21) <= fitMeanpt->GetParError(21) || fitMeanpt->GetParError(21) >= 0.75*fitMeanpt->GetParameter(21) || (fitMeanpt->GetParameter(21)/psipPtLim > 0.9)) )
  {
    fitMeanpt->FixParameter(21, bck->Eval(3.68));
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 17, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(17) <= fitMeanpt->GetParError(17) || fitMeanpt->GetParError(17) >= 0.75*fitMeanpt->GetParameter(13)) )
  {
    fitMeanpt->SetParameter(17, 2.);

    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanpt->GetParameter(17),fitMeanpt->GetParameter(18),fitMeanpt->GetParameter(19),fitMeanpt->GetParameter(20));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
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


  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol3,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol3");

  bck->SetParameters(3.,1.,-0.4,0.05);
  bck->SetParLimits(0, 0.,5.0);


  SetFitRejectRange(2.6,4.0);

  p->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);

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
  fitMeanpt->SetParLimits(12, 1.0,10.);

  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }

  fitMeanpt->SetParameter(17, 3.);
  fitMeanpt->SetParLimits(17, 0.,10.);


  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( static_cast<int>(fitResult) && (fitMeanpt->GetParameter(17) <= fitMeanpt->GetParError(17) || fitMeanpt->GetParError(17) >= 0.75*fitMeanpt->GetParameter(17)) )
  {
    fitMeanpt->FixParameter(17, 0.);
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14),fitMeanpt->GetParameter(15),fitMeanpt->GetParameter(16));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
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

  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundLin,fitRangeLow,fitRangeHigh,2,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundLin");

  bck->SetParameters(3.,0.);
  bck->SetParLimits(0, 0.,5.0);
  SetFitRejectRange(2.0,4.0);

  fHisto->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);

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

  fitMeanpt->SetParameter(15, 3.);
  fitMeanpt->SetParLimits(15, 0.5,6.);


  //  TProfile::Approximate();

  const char* fitOption = "SER"; //+";

  TFitResultPtr fitResult = fHisto->Fit(fitMeanpt,fitOption,"");

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanpt->GetParameter(13),fitMeanpt->GetParameter(14));

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
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
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(20),fitMeanpt->GetParError(20));
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSIPSIPRIMECB2VWG_BKGMPTPOL4()
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


  //  TString resultName(Form("MEANPTFIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4");

  bck->SetParameters(3.,-0.5,-0.5,1.5,-0.01);
  bck->SetParLimits(0, 0.,5.0);
  //  bck->SetParLimits(0, 1.,8.0);

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

  const char* fitOption = "SER";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanpt->GetParameter(18) <= fitMeanpt->GetParError(18) || fitMeanpt->GetParError(18) >= 0.75*fitMeanpt->GetParameter(18) )
  {
    fitMeanpt->FixParameter(18, 0.);
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanpt->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  for ( Int_t i = 0; i < 5; ++i )
  {
    bck->SetParameter(i, fitMeanpt->GetParameter(i+13));
  }

  AttachFunctionsToHisto(0x0,bck,fitMeanpt,fitRangeLow,fitRangeHigh);//
  printf("Final fit status : %d (Cov. Mat. : %d)\n",static_cast<int>(fitResult),fitResult->CovMatrixStatus());
  Set("FitResult",static_cast<int>(fitResult),0);
  Set("CovMatrixStatus",static_cast<int>(fitResult->CovMatrixStatus()),0);
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("MeanPtJPsi",fitMeanpt->GetParameter(12),fitMeanpt->GetParError(12));
  Set("MeanPtPsiP",fitMeanpt->GetParameter(16),fitMeanpt->GetParError(16));
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2VWG_BKGMV2POL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2

  const char* fitOption = "SERI"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SERI"; //We can add NO to avoid plotting

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

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  // bck->SetParameters(.3,-0.5,0.5,1.5,-0.01);
  // bck->SetParameters(-5.,-90.,90.,-30,3.);
  // bck->SetParameters(4.,-6.,3.);
  bck->SetParameters(.3,-.17,.02);

  // bck->SetParLimits(0, -1.,1.0);
 // bck->SetParLimits(0, 1.,8.0);
  // SetFitRejectRange(2.3,4.0);
  SetFitRejectRange(2.5,3.5);
  // SetFitRejectRange(2.6,4.0);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol2",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL2,fitRangeLow,fitRangeHigh,17,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2VWGPOL2");

  fitMeanv2->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanv2->SetParName(11,"kPsiP");
  fitMeanv2->SetParName(12,"<v2>JPsi");
  fitMeanv2->SetParName(13,"<v2>BG0");
  fitMeanv2->SetParName(14,"<v2>BG1");
  fitMeanv2->SetParName(15,"<v2>BG2");
  fitMeanv2->SetParName(16,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG);
  fitMeanv2->FixParameter(1,mVWG);
  fitMeanv2->FixParameter(2,sVWG1);
  fitMeanv2->FixParameter(3,sVWG2);
  fitMeanv2->FixParameter(4,kJPsi);
  fitMeanv2->FixParameter(5,mJPsi);
  fitMeanv2->FixParameter(6,sJPsi);
  fitMeanv2->FixParameter(7,alphaLow);
  fitMeanv2->FixParameter(8,nLow);
  fitMeanv2->FixParameter(9,alphaUp);
  fitMeanv2->FixParameter(10,nUp);
  fitMeanv2->FixParameter(11,kPsiP);

  fitMeanv2->SetParameter(12, 0.01);
  fitMeanv2->SetParLimits(12, -1.,1.);


  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanv2->SetParameter(i + 13, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(16, 0.01);
  fitMeanv2->SetParLimits(16, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  // else if ( fitMeanv2->GetParameter(16) <= fitMeanv2->GetParError(16) || fitMeanv2->GetParError(16) >= 0.75*fitMeanv2->GetParameter(16) )
  // {
  //   fitMeanv2->FixParameter(16, 0.);
  //   // for ( Int_t i = 0; i < 3; ++i )
  //   // {
  //   //   fitMeanv2->SetParameter(i + 13, bck->GetParameter(i)*0.6);
  //   // }
  //   fitResult = p->Fit(fitMeanv2,fitOption,"");

  //   std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  // }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

    //___________ Further attempts to fit if the first one fails
  if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,11,3); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(5),fitMeanv2->GetParError(5));
  Set("sJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  // Set("<v2>JPsi",fitMeanv2->GetParameter(12),fitMeanv2->GetParError(12));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 3; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+13));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(12),fitMeanv2->GetParError(12));
  Set("<v2>PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2VWG_BKGMV2POL2EXP()
{
  const char* fitOption = "SERI"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SERI"; //We can add NO to avoid plotting

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

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");

  bck->SetParameters(.1,1.,0.4,-0.0);
  // bck->SetParameters(.1,1.,0.);

  bck->SetParLimits(0, -2.,2.);

  // SetFitRejectRange(2.7,4.); //2.6,3.4
  SetFitRejectRange(2.5,3.5);

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  //
  //  bck->SetParameters(1,0,0,0);


  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"",fitRangeLow,fitRangeHigh);

  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol2Exp",fitOptionBg);
  //___________

  SetFitRejectRange();
  std::cout << "range : " << fitRangeLow << " - " << fitRangeHigh << std::endl;
  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL2EXP,fitRangeLow,fitRangeHigh,18,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2VWGPOL2EXP");

  fitMeanv2->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  //                        0      1      2       3       4       5       6       7         8        9        10
  fitMeanv2->SetParName(11,"kPsiP");
  fitMeanv2->SetParName(12,"<v2>JPsi");
  fitMeanv2->SetParName(13,"<v2>BG0");
  fitMeanv2->SetParName(14,"<v2>BG1");
  fitMeanv2->SetParName(15,"<v2>BG2");
  fitMeanv2->SetParName(16,"<v2>BGEXP");
  fitMeanv2->SetParName(17,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG);
  fitMeanv2->FixParameter(1,mVWG);
  fitMeanv2->FixParameter(2,sVWG1);
  fitMeanv2->FixParameter(3,sVWG2);
  fitMeanv2->FixParameter(4,kJPsi);
  fitMeanv2->FixParameter(5,mJPsi);
  fitMeanv2->FixParameter(6,sJPsi);
  fitMeanv2->FixParameter(7,alphaLow);
  fitMeanv2->FixParameter(8,nLow);
  fitMeanv2->FixParameter(9,alphaUp);
  fitMeanv2->FixParameter(10,nUp);
  fitMeanv2->FixParameter(11,kPsiP);

  fitMeanv2->SetParameter(12, 0.01);
  fitMeanv2->SetParLimits(12, -1.0,1.);


  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanv2->SetParameter(i + 13, bck->GetParameter(i));
  }

  //  fitMeanv2->SetParameter(13, 3.);
  //  fitMeanv2->SetParLimits(13, 0.5,10.);
  //
  //  fitMeanv2->SetParameter(14, 0.2);
  //  //  fitMeanv2->SetParLimits(14, 0.1,0.2);
  //
  //  fitMeanv2->SetParameter(15, 0.1);
  //  //  fitMeanv2->SetParLimits(15, 0.,1.);

  fitMeanv2->SetParameter(17, 0.01);
  // Double_t psipPtLim = 0.1;
  fitMeanv2->SetParLimits(17, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

    //___________ Further attempts to fit if the first one fails
  if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,11,3); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(5),fitMeanv2->GetParError(5));
  Set("sJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  // Set("<v2>JPsi",fitMeanv2->GetParameter(12),fitMeanv2->GetParError(12));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanv2->GetParameter(13),fitMeanv2->GetParameter(14),fitMeanv2->GetParameter(15),fitMeanv2->GetParameter(16));

  AttachFunctionsToHisto(0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(12),fitMeanv2->GetParError(12));
  Set("<v2>PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));
}


//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2VWG_BKGMV2POL3()
{

  const char* fitOption = "SERL"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SER"; //We can add NO to avoid plotting
  //Fit mean dimuon mean v2 to get Jpsi mean v2 using the CB2 signal parameters, VWG for the Bkg, Jpsi mv2 = cte and Bkg mv2 = pol2
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

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean V2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();

    TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol3,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol3");

  bck->SetParameters(3.,1.,0.,0.05);

  bck->SetParLimits(0, 0.,5.0);
  // bck->SetParLimits(0, 1.,8.0);

  SetFitRejectRange(2.6,4.0);
  // SetFitRejectRange(2.7,4.0);


  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"",fitRangeLow,fitRangeHigh);
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol2",fitOptionBg);
  //___________

  SetFitRejectRange();


  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL3,fitRangeLow,fitRangeHigh,18,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2VWGPOL3");

  fitMeanv2->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanv2->SetParName(11,"kPsiP");
  fitMeanv2->SetParName(12,"<v2>JPsi");
  fitMeanv2->SetParName(13,"<v2>BG0");
  fitMeanv2->SetParName(14,"<v2>BG1");
  fitMeanv2->SetParName(15,"<v2>BG2");
  fitMeanv2->SetParName(16,"<v2>BG3");
  fitMeanv2->SetParName(17,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG);
  fitMeanv2->FixParameter(1,mVWG);
  fitMeanv2->FixParameter(2,sVWG1);
  fitMeanv2->FixParameter(3,sVWG2);
  fitMeanv2->FixParameter(4,kJPsi);
  fitMeanv2->FixParameter(5,mJPsi);
  fitMeanv2->FixParameter(6,sJPsi);
  fitMeanv2->FixParameter(7,alphaLow);
  fitMeanv2->FixParameter(8,nLow);
  fitMeanv2->FixParameter(9,alphaUp);
  fitMeanv2->FixParameter(10,nUp);
  fitMeanv2->FixParameter(11,kPsiP);

  fitMeanv2->SetParameter(12, 0.01);
  fitMeanv2->SetParLimits(12, -1.0,1.);

  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanv2->SetParameter(i + 13, bck->GetParameter(i));
  }
    fitMeanv2->SetParameter(17, 0.01);
  fitMeanv2->SetParLimits(17, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  if ( static_cast<int>(fitResult) && (fitMeanv2->GetParameter(17) <= fitMeanv2->GetParError(17) || fitMeanv2->GetParError(17) >= 0.75*fitMeanv2->GetParameter(17)) )
  {
    fitMeanv2->FixParameter(17, 0.);
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

    //___________ Further attempts to fit if the first one fails
  if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,11,3); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(5),fitMeanv2->GetParError(5));
  Set("sJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  // Set("<v2>JPsi",fitMeanv2->GetParameter(12),fitMeanv2->GetParError(12));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  //  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(fitMeanv2->GetParameter(13),fitMeanv2->GetParameter(14),fitMeanv2->GetParameter(15),fitMeanv2->GetParameter(16));

  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(12),fitMeanv2->GetParError(12));
  Set("<v2>PsiP",fitMeanv2->GetParameter(17),fitMeanv2->GetParError(17));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2VWG_BKGMV2POL4()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2

  const char* fitOption = "SERL"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SER"; //We can add NO to avoid plotting

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

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4");

  // bck->SetParameters(.3,-0.5,0.5,1.5,-0.01);
  // bck->SetParameters(-5.,-90.,90.,-30,3.);
  bck->SetParameters(4.,-6.,3.,-.7,.02);
  // bck->SetParLimits(0, -1.,1.0);
 // bck->SetParLimits(0, 1.,8.0);
  // SetFitRejectRange(2.3,4.0);
  SetFitRejectRange(2.6,3.4);
  // SetFitRejectRange(2.6,4.0);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol4",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtS2CB2VWGPOL4,fitRangeLow,fitRangeHigh,19,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtS2CB2VWGPOL4");

  fitMeanv2->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanv2->SetParName(11,"kPsiP");
  fitMeanv2->SetParName(12,"<v2>JPsi");
  fitMeanv2->SetParName(13,"<v2>BG0");
  fitMeanv2->SetParName(14,"<v2>BG1");
  fitMeanv2->SetParName(15,"<v2>BG2");
  fitMeanv2->SetParName(16,"<v2>BG3");
  fitMeanv2->SetParName(17,"<v2>BG4");
  fitMeanv2->SetParName(18,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG);
  fitMeanv2->FixParameter(1,mVWG);
  fitMeanv2->FixParameter(2,sVWG1);
  fitMeanv2->FixParameter(3,sVWG2);
  fitMeanv2->FixParameter(4,kJPsi);
  fitMeanv2->FixParameter(5,mJPsi);
  fitMeanv2->FixParameter(6,sJPsi);
  fitMeanv2->FixParameter(7,alphaLow);
  fitMeanv2->FixParameter(8,nLow);
  fitMeanv2->FixParameter(9,alphaUp);
  fitMeanv2->FixParameter(10,nUp);
  fitMeanv2->FixParameter(11,kPsiP);

  fitMeanv2->SetParameter(12, 0.01);
  fitMeanv2->SetParLimits(12, -1.,1.);

  for ( Int_t i = 0; i < 5; ++i )
  {
    fitMeanv2->SetParameter(i + 13, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(18, 0.01);
  fitMeanv2->SetParLimits(18, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 13, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  else if ( fitMeanv2->GetParameter(18) <= fitMeanv2->GetParError(18) || fitMeanv2->GetParError(18) >= 0.75*fitMeanv2->GetParameter(18) )
  {
    fitMeanv2->FixParameter(18, 0.);
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 13, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

    //___________ Further attempts to fit if the first one fails
  if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,11,3); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(5),fitMeanv2->GetParError(5));
  Set("sJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  // Set("<v2>JPsi",fitMeanv2->GetParameter(12),fitMeanv2->GetParError(12));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));


  for ( Int_t i = 0; i < 5; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+13));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(12),fitMeanv2->GetParError(12));
  Set("<v2>PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMPTPSI_HFUNCTION()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2
  fHisto->GetListOfFunctions()->Delete();

  Double_t ah = GetValue("ah");
  if(!IsValidValue(ah))  ah=0; // Defalut value if nothing in the config string
  Double_t bh = GetValue("bh");
  if(!IsValidValue(bh))  bh=0; // Defalut value if nothing in the config string
  Double_t ch = GetValue("ch");
  if(!IsValidValue(ch))  ch=0; // Defalut value if nothing in the config string
  Double_t dh = GetValue("dh");
  if(!IsValidValue(dh))  dh=0; // Defalut value if nothing in the config string
  Double_t eh = GetValue("eh");
  if(!IsValidValue(eh))  eh=0; // Defalut value if nothing in the config string
  Double_t fh = GetValue("fh");
  if(!IsValidValue(fh))  fh=0; // Defalut value if nothing in the config string
  Double_t gh = GetValue("gh");
  if(!IsValidValue(gh))  gh=0; // Defalut value if nothing in the config string

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean pt histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::hFunction,2.5,3.29,6,"AliAnalysisMuMuJpsiResult","hFunction");
  bck->SetParameter(0,12864.7);
  bck->SetParameter(1,ah);
  bck->SetParameter(2,bh);
  bck->SetParameter(3,ch);
  bck->SetParameter(4,dh);
  bck->SetParameter(5,eh);

  p->Fit(bck,"SRL","");

  TF1* fitMeanpt = new TF1("fitMeanpt",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtHFunction,fitRangeLow,fitRangeHigh,8,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtHFunction");

  fitMeanpt->SetParNames("<pt>JPsi","ah","bh","ch","dh","eh","fh","gh");

  fitMeanpt->SetParameter(0,bck->GetParameter(0));
  fitMeanpt->SetParameter(1,bck->GetParameter(1));
  fitMeanpt->SetParameter(2,bck->GetParameter(2));
  fitMeanpt->SetParameter(3,bck->GetParameter(3));
  fitMeanpt->SetParameter(4,bck->GetParameter(4));
  fitMeanpt->SetParameter(5,bck->GetParameter(5));
  fitMeanpt->SetParameter(6,bck->GetParameter(6));
  fitMeanpt->SetParameter(7,bck->GetParameter(7));

  const char* fitOption = "SERM";

  TFitResultPtr fitResult = p->Fit(fitMeanpt,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 )
  {
    for ( Int_t i = 1; i < 6; ++i )
    {
      fitMeanpt->SetParameter(i,fitMeanpt->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanpt,fitOption,"");

    std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  }

  AttachFunctionsToHisto(0x0,0x0,fitMeanpt,fitRangeLow,fitRangeHigh);//

  Set("MeanPtJPsi",fitMeanpt->GetParameter(0),fitMeanpt->GetParError(0));
}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2VWG2_BKGMV2POLEXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SER"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SER"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  Double_t kVWG2 = GetValue("kVWG2");
  Double_t mVWG2 = GetValue("mVWG2");
  Double_t s1VWG2 = GetValue("s1VWG2");
  Double_t s2VWG2 = GetValue("s2VWG2");
  Double_t gVWG2 = GetValue("gVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");



  Double_t BG0_init    = IsValidValue(GetValue("BG0_init"))  ? GetValue("BG0_init")  : .4;
  Double_t BG1_init    = IsValidValue(GetValue("BG1_init"))  ? GetValue("BG1_init")  : -.2;
  Double_t BGexp_init    = IsValidValue(GetValue("BGexp_init"))  ? GetValue("BGexp_init")  : -.2;
  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.6;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 4.;

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : alJPsi = %f nlJPsi = %f auJPsi = %f nuJPsi = %f",alphaLow,nLow,alphaUp,nUp));
  AliDebug(1,Form("values : kVWG2 = %f mVG2 = %f s1VWG2 = %f s2VWG2 = %f gVWG2 = %f",kVWG2, mVWG2, s1VWG2,s2VWG2,gVWG2));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));
  if(IsValidValue(GetValue("BG0_init"))) AliInfo(Form("Bckg initialised using parameters : %f,%f,%f",BG0_init,BG1_init,BGexp_init));
  if(IsValidValue(GetValue("rejRl"))) AliInfo(Form("Bckg initialised using reject range: %f - %f",rejRl,rejRh));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPolExp,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPolExp");
   bck->SetParameters(BG0_init,BG1_init,BGexp_init);
  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPolExp",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2VWG2POLEXP,fitRangeLow,fitRangeHigh,18,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSCB2VWG2POLEXP");

  fitMeanv2->SetParNames("kVWG2","mVWG2","s1VWG2","s2VWG2","gVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi");
  //                        0      1      2       3       4       5       6       7         8        9        10
  fitMeanv2->SetParName(11,"nuJPsi");
  fitMeanv2->SetParName(12,"kPsiP");
  fitMeanv2->SetParName(13,"<v2>JPsi");
  fitMeanv2->SetParName(14,"<v2>BG0");
  fitMeanv2->SetParName(15,"<v2>BG1");
  fitMeanv2->SetParName(16,"<v2>BGEXP");
  fitMeanv2->SetParName(17,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG2);
  fitMeanv2->FixParameter(1,mVWG2);
  fitMeanv2->FixParameter(2,s1VWG2);
  fitMeanv2->FixParameter(3,s2VWG2);
  fitMeanv2->FixParameter(4,gVWG2);
  fitMeanv2->FixParameter(5,kJPsi);
  fitMeanv2->FixParameter(6,mJPsi);
  fitMeanv2->FixParameter(7,sJPsi);
  fitMeanv2->FixParameter(8,alphaLow);
  fitMeanv2->FixParameter(9,nLow);
  fitMeanv2->FixParameter(10,alphaUp);
  fitMeanv2->FixParameter(11,nUp);
  fitMeanv2->FixParameter(12,kPsiP);

  fitMeanv2->SetParameter(13, 0.01);
  fitMeanv2->SetParLimits(13, -1.,1.);

  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanv2->SetParameter(i + 14, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(17, 0.01);
  fitMeanv2->SetParLimits(17, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 14, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(17) <= fitMeanv2->GetParError(17) || fitMeanv2->GetParError(17) >= 0.75*fitMeanv2->GetParameter(17) )
  {
    fitMeanv2->FixParameter(17, 0.);
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 14, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  Set("sJPsi",fitMeanv2->GetParameter(7),fitMeanv2->GetParError(7));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 3; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+14));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(13),fitMeanv2->GetParError(13));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(17),fitMeanv2->GetParError(17));

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2VWG2_BKGMV2POL2EXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2

  const char* fitOption = "SERIM"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SER"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  Double_t kVWG2 = GetValue("kVWG2");
  Double_t mVWG2 = GetValue("mVWG2");
  Double_t s1VWG2 = GetValue("s1VWG2");
  Double_t s2VWG2 = GetValue("s2VWG2");
  Double_t gVWG2 = GetValue("gVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : alJPsi = %f nlJPsi = %f auJPsi = %f nuJPsi = %f",alphaLow,nLow,alphaUp,nUp));
  AliDebug(1,Form("values : kVWG2 = %f mVG2 = %f s1VWG2 = %f s2VWG2 = %f gVWG2 = %f",kVWG2, mVWG2, s1VWG2,s2VWG2,gVWG2));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));

   TString msg;

   if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
   if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
   if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
   if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);

   AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol2 for Bkg <pt> %s",msg.Data()));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();
  //_____________
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");
  bck->SetParameters(1.7,-1.3,0.2,-0.7);
  SetFitRejectRange(2.4,3.6);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol2Exp",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2VWG2POL2EXP,fitRangeLow,fitRangeHigh,19,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSCB2VWG2POL2EXP");

  fitMeanv2->SetParNames("kVWG2","mVWG2","s1VWG2","s2VWG2","gVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi");
  //                        0      1      2       3       4       5       6       7         8        9        10
  fitMeanv2->SetParName(11,"nuJPsi");
  fitMeanv2->SetParName(12,"kPsiP");
  fitMeanv2->SetParName(13,"<v2>JPsi");
  fitMeanv2->SetParName(14,"<v2>BG0");
  fitMeanv2->SetParName(15,"<v2>BG1");
  fitMeanv2->SetParName(16,"<v2>BG2");
  fitMeanv2->SetParName(17,"<v2>BGEXP");
  fitMeanv2->SetParName(18,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG2);
  fitMeanv2->FixParameter(1,mVWG2);
  fitMeanv2->FixParameter(2,s1VWG2);
  fitMeanv2->FixParameter(3,s2VWG2);
  fitMeanv2->FixParameter(4,gVWG2);
  fitMeanv2->FixParameter(5,kJPsi);
  fitMeanv2->FixParameter(6,mJPsi);
  fitMeanv2->FixParameter(7,sJPsi);
  fitMeanv2->FixParameter(8,alphaLow);
  fitMeanv2->FixParameter(9,nLow);
  fitMeanv2->FixParameter(10,alphaUp);
  fitMeanv2->FixParameter(11,nUp);
  fitMeanv2->FixParameter(12,kPsiP);

  fitMeanv2->SetParameter(13, 0.01);
  fitMeanv2->SetParLimits(13, -1.,1.);

  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanv2->SetParameter(i + 14, bck->GetParameter(i));
  }
  // fitMeanv2->FixParameter(17,0.);
  fitMeanv2->SetParameter(18, 0.01);
  fitMeanv2->SetParLimits(18, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    std::cout << " ========= Refitting with bck param * 0.9 =========" << std::endl;
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 14, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  }

  if ( fitMeanv2->GetParameter(18) <= fitMeanv2->GetParError(18) || fitMeanv2->GetParError(18) >= 0.75*fitMeanv2->GetParameter(18) )
  {
    std::cout << " ========= Refitting with bck param * 0.9 and v2(psi2S)=0 =========" << std::endl;
    fitMeanv2->FixParameter(18, 0.);
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 14, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  }
  if ( fitMeanv2->GetParameter(18) <= fitMeanv2->GetParError(18) || fitMeanv2->GetParError(18) >= 0.75*fitMeanv2->GetParameter(18) )
  {
    std::cout << " ========= Refitting with bck param * 0.75 and v2(psi2S)=0 =========" << std::endl;
    fitMeanv2->FixParameter(18, 0.);
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 14, bck->GetParameter(i)*0.75);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  }
  if (  static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3  )
  {
    std::cout << " ========= Refitting with bck param * 0.6 and v2(psi2S)=0 =========" << std::endl;
    fitMeanv2->FixParameter(18, 0.);
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 14, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  }
  if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3  )
  {
    std::cout << " ========= Refitting with bck param * 0.5 and v2(psi2S)=0 =========" << std::endl;
    fitMeanv2->FixParameter(18, 0.);
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 14, bck->GetParameter(i)*0.5);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  }
  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;

    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  Set("sJPsi",fitMeanv2->GetParameter(7),fitMeanv2->GetParError(7));
  Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 4; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+14));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(13),fitMeanv2->GetParError(13));
  Set("<v2>PsiP",fitMeanv2->GetParameter(18),fitMeanv2->GetParError(18));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2VWG2_BKGMV2POL4()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2

  const char* fitOption = "SERI"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SERI"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  Double_t kVWG2 = GetValue("kVWG2");
  Double_t mVWG2 = GetValue("mVWG2");
  Double_t s1VWG2 = GetValue("s1VWG2");
  Double_t s2VWG2 = GetValue("s2VWG2");
  Double_t gVWG2 = GetValue("gVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : alJPsi = %f nlJPsi = %f auJPsi = %f nuJPsi = %f",alphaLow,nLow,alphaUp,nUp));
  AliDebug(1,Form("values : kVWG2 = %f mVG2 = %f s1VWG2 = %f s2VWG2 = %f gVWG2 = %f",kVWG2, mVWG2, s1VWG2,s2VWG2,gVWG2));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));

   TString msg;

   if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
   if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
   if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
   if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);

   AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol4 for Bkg <pt> %s",msg.Data()));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();
  //_____________
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4");
  bck->SetParameters(.1,-.1,.01,0.,0.);
  SetFitRejectRange(2.5,3.5);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol4",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2VWG2POL4,fitRangeLow,fitRangeHigh,19,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSCB2VWG2POL4");

  fitMeanv2->SetParNames("kVWG2","mVWG2","s1VWG2","s2VWG2","gVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi");
  //                        0      1      2       3       4       5       6       7         8        9        10
  fitMeanv2->SetParName(11,"nuJPsi");
  fitMeanv2->SetParName(12,"kPsiP");
  fitMeanv2->SetParName(13,"<v2>JPsi");
  fitMeanv2->SetParName(14,"<v2>BG0");
  fitMeanv2->SetParName(15,"<v2>BG1");
  fitMeanv2->SetParName(16,"<v2>BG2");
  fitMeanv2->SetParName(17,"<v2>BG3");
  fitMeanv2->SetParName(18,"<v2>BG4");
  fitMeanv2->SetParName(19,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG2);
  fitMeanv2->FixParameter(1,mVWG2);
  fitMeanv2->FixParameter(2,s1VWG2);
  fitMeanv2->FixParameter(3,s2VWG2);
  fitMeanv2->FixParameter(4,gVWG2);
  fitMeanv2->FixParameter(5,kJPsi);
  fitMeanv2->FixParameter(6,mJPsi);
  fitMeanv2->FixParameter(7,sJPsi);
  fitMeanv2->FixParameter(8,alphaLow);
  fitMeanv2->FixParameter(9,nLow);
  fitMeanv2->FixParameter(10,alphaUp);
  fitMeanv2->FixParameter(11,nUp);
  fitMeanv2->FixParameter(12,kPsiP);

  fitMeanv2->SetParameter(13, 0.01);
  fitMeanv2->SetParLimits(13, -1.,1.);

  for ( Int_t i = 0; i < 5; ++i )
  {
    fitMeanv2->SetParameter(i + 14, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(19, 0.01);
  fitMeanv2->SetParLimits(19, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 14, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }
  else if (  static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3  )
  {
    fitMeanv2->FixParameter(19, 0.);
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 14, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");

    std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;

    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  Set("sJPsi",fitMeanv2->GetParameter(7),fitMeanv2->GetParError(7));
  Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 5; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+14));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(13),fitMeanv2->GetParError(13));
  Set("<v2>PsiP",fitMeanv2->GetParameter(19),fitMeanv2->GetParError(19));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2VWG2_BKGMV2POL4Cheb()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol2

  const char* fitOption = "SER"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SER"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");
  Double_t kVWG2 = GetValue("kVWG2");
  Double_t mVWG2 = GetValue("mVWG2");
  Double_t s1VWG2 = GetValue("s1VWG2");
  Double_t s2VWG2 = GetValue("s2VWG2");
  Double_t gVWG2 = GetValue("gVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t BG0_init    = IsValidValue(GetValue("BG0_init"))  ? GetValue("BG0_init")  : .1;
  Double_t BG1_init    = IsValidValue(GetValue("BG1_init"))  ? GetValue("BG1_init")  : -.1;
  Double_t BG2_init    = IsValidValue(GetValue("BG2_init"))  ? GetValue("BG2_init")  : -.05;
  Double_t BG3_init    = IsValidValue(GetValue("BG3_init"))  ? GetValue("BG3_init")  : 0.;
  Double_t BG4_init    = IsValidValue(GetValue("BG4_init"))  ? GetValue("BG4_init")  : 0.;
  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.5;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 4.;

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : alJPsi = %f nlJPsi = %f auJPsi = %f nuJPsi = %f",alphaLow,nLow,alphaUp,nUp));
  AliDebug(1,Form("values : kVWG2 = %f mVG2 = %f s1VWG2 = %f s2VWG2 = %f gVWG2 = %f",kVWG2, mVWG2, s1VWG2,s2VWG2,gVWG2));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));

  TString msg;

  if (IsValidValue(alphaLow)) msg += TString::Format("alphaLow=%e ",alphaLow);
  if (IsValidValue(nLow)) msg += TString::Format("nLow=%e ",nLow);
  if (IsValidValue(alphaUp)) msg += TString::Format("alphaUp=%e ",alphaUp);
  if (IsValidValue(nUp)) msg += TString::Format("nUp=%e ",nUp);
  if(IsValidValue(GetValue("BG0_init"))) AliInfo(Form("Bckg initialised using parameters : %f,%f,%f,%f",BG0_init,BG1_init,BG3_init,BG4_init));
  if(IsValidValue(GetValue("rejRl"))) AliInfo(Form("Bckg initialised using reject range: %f - %f",rejRl,rejRh));

   AliDebug(1,Form("Mean pt fit with jpsi + psiprime (CB2),Bkg VWG and Pol4Cheb for Bkg <pt> %s",msg.Data()));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();
  //_____________
  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Cheb,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Cheb");


  bck->SetParameters(BG0_init,BG1_init,BG2_init,BG3_init,BG4_init);
  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol4Cheb",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2VWG2POL4Cheb,fitRangeLow,fitRangeHigh,19,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSCB2VWG2POL4Cheb");

  fitMeanv2->SetParNames("kVWG2","mVWG2","s1VWG2","s2VWG2","gVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi");
  //                        0      1      2       3       4       5       6       7         8        9        10
  fitMeanv2->SetParName(11,"nuJPsi");
  fitMeanv2->SetParName(12,"kPsiP");
  fitMeanv2->SetParName(13,"<v2>JPsi");
  fitMeanv2->SetParName(14,"<v2>BGC0");
  fitMeanv2->SetParName(15,"<v2>BGC1");
  fitMeanv2->SetParName(16,"<v2>BGC2");
  fitMeanv2->SetParName(17,"<v2>BGC3");
  fitMeanv2->SetParName(18,"<v2>BGC4");
  fitMeanv2->SetParName(19,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG2);
  fitMeanv2->FixParameter(1,mVWG2);
  fitMeanv2->FixParameter(2,s1VWG2);
  fitMeanv2->FixParameter(3,s2VWG2);
  fitMeanv2->FixParameter(4,gVWG2);
  fitMeanv2->FixParameter(5,kJPsi);
  fitMeanv2->FixParameter(6,mJPsi);
  fitMeanv2->FixParameter(7,sJPsi);
  fitMeanv2->FixParameter(8,alphaLow);
  fitMeanv2->FixParameter(9,nLow);
  fitMeanv2->FixParameter(10,alphaUp);
  fitMeanv2->FixParameter(11,nUp);
  fitMeanv2->FixParameter(12,kPsiP);

  fitMeanv2->SetParameter(13, 0.01);
  fitMeanv2->SetParLimits(13, -1.,1.);

  for ( Int_t i = 0; i < 5; ++i )
  {
    fitMeanv2->SetParameter(i + 14, bck->GetParameter(i));
  }
  // fitMeanv2->FixParameter(18,0.);
  // fitMeanv2->FixParameter(17,0.);

  fitMeanv2->SetParameter(19, 0.01);
  fitMeanv2->SetParLimits(19, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 14, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(19) <= fitMeanv2->GetParError(19) || fitMeanv2->GetParError(19) >= 0.75*fitMeanv2->GetParameter(19) )
  {
    fitMeanv2->FixParameter(19, 0.);
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 14, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;

    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  Set("sJPsi",fitMeanv2->GetParameter(7),fitMeanv2->GetParError(7));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 5; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+14));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(13),fitMeanv2->GetParError(13));
  Set("<v2>PsiP",fitMeanv2->GetParameter(19),fitMeanv2->GetParError(19));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2POL2POL3_BKGMV2POL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SER"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SER"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");

  Double_t a  = GetValue("a");
  Double_t b  = GetValue("b");
  Double_t c  = GetValue("c");
  Double_t ap = GetValue("a'");
  Double_t bp = GetValue("b'");
  Double_t cp = GetValue("c'");
  Double_t dp = GetValue("d'");

  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t BG0_init    = IsValidValue(GetValue("BG0_init"))  ? GetValue("BG0_init")  : -0.15;
  Double_t BG1_init    = IsValidValue(GetValue("BG1_init"))  ? GetValue("BG1_init")  : .2;
  Double_t BG2_init    = IsValidValue(GetValue("BG2_init"))  ? GetValue("BG2_init")  : -.004;
  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.6;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 3.6;

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : alJPsi = %f nlJPsi = %f auJPsi = %f nuJPsi = %f",alphaLow,nLow,alphaUp,nUp));
  AliDebug(1,Form("values : a = %f b = %f c = %f a' = %f b' = %f c' = %f d' = %f ",a,b,c,ap,bp,cp,dp));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));
  if(IsValidValue(GetValue("BG0_init"))) AliInfo(Form("Bckg initialised using parameters : %f,%f,%f",BG0_init,BG1_init,BG2_init));
  if(IsValidValue(GetValue("rejRl"))) AliInfo(Form("Bckg initialised using reject range: %f - %f",rejRl,rejRh));

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(BG0_init,BG1_init,BG2_init);

  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol2",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2POL2POL3POL2,fitRangeLow,fitRangeHigh,20,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSCB2POL2POL3POL2");

  fitMeanv2->SetParNames("a","b","c","ap","bp","cp","dp","kJPsi","mJPsi","sJPsi","alJPsi");
  //                      0   1   2    3    4    5   6       7       8       9      10
  fitMeanv2->SetParName(11,"nlJPsi");
  fitMeanv2->SetParName(12,"auJPsi");
  fitMeanv2->SetParName(13,"nuJPsi");
  fitMeanv2->SetParName(14,"kPsiP");
  fitMeanv2->SetParName(15,"<v2>JPsi");
  fitMeanv2->SetParName(16,"<v2>BG0");
  fitMeanv2->SetParName(17,"<v2>BG1");
  fitMeanv2->SetParName(18,"<v2>BG2");
  fitMeanv2->SetParName(19,"<v2>PsiP");

  fitMeanv2->FixParameter(0,a);
  fitMeanv2->FixParameter(1,b);
  fitMeanv2->FixParameter(2,c);
  fitMeanv2->FixParameter(3,ap);
  fitMeanv2->FixParameter(4,bp);
  fitMeanv2->FixParameter(5,cp);
  fitMeanv2->FixParameter(6,dp);
  fitMeanv2->FixParameter(7,kJPsi);
  fitMeanv2->FixParameter(8,mJPsi);
  fitMeanv2->FixParameter(9,sJPsi);
  fitMeanv2->FixParameter(10,alphaLow);
  fitMeanv2->FixParameter(11,nLow);
  fitMeanv2->FixParameter(12,alphaUp);
  fitMeanv2->FixParameter(13,nUp);
  fitMeanv2->FixParameter(14,kPsiP);

  fitMeanv2->SetParameter(15, 0.01);
  fitMeanv2->SetParLimits(15, -1.,1.);

  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanv2->SetParameter(i + 16, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(19, 0.01);
  fitMeanv2->SetParLimits(19, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 16, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(19) <= fitMeanv2->GetParError(19) || fitMeanv2->GetParError(19) >= 0.75*fitMeanv2->GetParameter(19) )
  {
    fitMeanv2->FixParameter(19, 0.);
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 16, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(8),fitMeanv2->GetParError(8));
  Set("sJPsi",fitMeanv2->GetParameter(9),fitMeanv2->GetParError(9));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 3; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+16));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(15),fitMeanv2->GetParError(15));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(19),fitMeanv2->GetParError(19));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2POL2POL3_BKGMV2POLEXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SER"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SER"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");

  Double_t a  = GetValue("a");
  Double_t b  = GetValue("b");
  Double_t c  = GetValue("c");
  Double_t ap = GetValue("a'");
  Double_t bp = GetValue("b'");
  Double_t cp = GetValue("c'");
  Double_t dp = GetValue("d'");

  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");


  Double_t BG0_init    = IsValidValue(GetValue("BG0_init"))  ? GetValue("BG0_init")  : .3;
  Double_t BG1_init    = IsValidValue(GetValue("BG1_init"))  ? GetValue("BG1_init")  : -.08;
  Double_t BGexp_init    = IsValidValue(GetValue("BGexp_init"))  ? GetValue("BGexp_init")  : 0.;
  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.6;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 4.;

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : alJPsi = %f nlJPsi = %f auJPsi = %f nuJPsi = %f",alphaLow,nLow,alphaUp,nUp));
  AliDebug(1,Form("values : a = %f b = %f c = %f a' = %f b' = %f c' = %f d' = %f ",a,b,c,ap,bp,cp,dp));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));
  if(IsValidValue(GetValue("BG0_init"))) AliInfo(Form("Bckg initialised using parameters : %f,%f,%f",BG0_init,BG1_init,BGexp_init));
  if(IsValidValue(GetValue("rejRl"))) AliInfo(Form("Bckg initialised using reject range: %f - %f",rejRl,rejRh));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPolExp,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPolExp");
  bck->SetParameters(BG0_init,BG1_init,BGexp_init);

  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPolExp",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2POL2POL3_POLEXP,fitRangeLow,fitRangeHigh,20,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSCB2POL2POL3_POLEXP");

  fitMeanv2->SetParNames("a","b","c","ap","bp","cp","dp","kJPsi","mJPsi","sJPsi","alJPsi");
  //                      0   1   2    3    4    5   6       7       8       9      10
  fitMeanv2->SetParName(11,"nlJPsi");
  fitMeanv2->SetParName(12,"auJPsi");
  fitMeanv2->SetParName(13,"nuJPsi");
  fitMeanv2->SetParName(14,"kPsiP");
  fitMeanv2->SetParName(15,"<v2>JPsi");
  fitMeanv2->SetParName(16,"<v2>BG0");
  fitMeanv2->SetParName(17,"<v2>BG1");
  fitMeanv2->SetParName(18,"<v2>BGEXP");
  fitMeanv2->SetParName(19,"<v2>PsiP");

  fitMeanv2->FixParameter(0,a);
  fitMeanv2->FixParameter(1,b);
  fitMeanv2->FixParameter(2,c);
  fitMeanv2->FixParameter(3,ap);
  fitMeanv2->FixParameter(4,bp);
  fitMeanv2->FixParameter(5,cp);
  fitMeanv2->FixParameter(6,dp);
  fitMeanv2->FixParameter(7,kJPsi);
  fitMeanv2->FixParameter(8,mJPsi);
  fitMeanv2->FixParameter(9,sJPsi);
  fitMeanv2->FixParameter(10,alphaLow);
  fitMeanv2->FixParameter(11,nLow);
  fitMeanv2->FixParameter(12,alphaUp);
  fitMeanv2->FixParameter(13,nUp);
  fitMeanv2->FixParameter(14,kPsiP);

  fitMeanv2->SetParameter(15, 0.01);
  fitMeanv2->SetParLimits(15, -1.,1.);

  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanv2->SetParameter(i + 16, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(19, 0.01);
  fitMeanv2->SetParLimits(19, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 16, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(19) <= fitMeanv2->GetParError(19) || fitMeanv2->GetParError(19) >= 0.75*fitMeanv2->GetParameter(19) )
  {
    fitMeanv2->FixParameter(19, 0.);
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 16, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(8),fitMeanv2->GetParError(8));
  Set("sJPsi",fitMeanv2->GetParameter(9),fitMeanv2->GetParError(9));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 3; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+16));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(15),fitMeanv2->GetParError(15));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(19),fitMeanv2->GetParError(19));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2POL2POL3_BKGMV2POL2EXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SERI"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SERI"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");

  Double_t a  = GetValue("a");
  Double_t b  = GetValue("b");
  Double_t c  = GetValue("c");
  Double_t ap = GetValue("a'");
  Double_t bp = GetValue("b'");
  Double_t cp = GetValue("c'");
  Double_t dp = GetValue("d'");

  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : alJPsi = %f nlJPsi = %f auJPsi = %f nuJPsi = %f",alphaLow,nLow,alphaUp,nUp));
  AliDebug(1,Form("values : a = %f b = %f c = %f a' = %f b' = %f c' = %f d' = %f ",a,b,c,ap,bp,cp,dp));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");

  bck->SetParameters(.4,-.2,.02);
  SetFitRejectRange(2.6,3.4);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol2",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2POL2POL3POL2EXP,fitRangeLow,fitRangeHigh,21,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSCB2POL2POL3POL2EXP");

  fitMeanv2->SetParNames("a","b","c","ap","bp","cp","dp","kJPsi","mJPsi","sJPsi","alJPsi");
  //                      0   1   2    3    4    5   6       7       8       9      10
  fitMeanv2->SetParName(11,"nlJPsi");
  fitMeanv2->SetParName(12,"auJPsi");
  fitMeanv2->SetParName(13,"nuJPsi");
  fitMeanv2->SetParName(14,"kPsiP");
  fitMeanv2->SetParName(15,"<v2>JPsi");
  fitMeanv2->SetParName(16,"<v2>BG0");
  fitMeanv2->SetParName(17,"<v2>BG1");
  fitMeanv2->SetParName(18,"<v2>BG2");
  fitMeanv2->SetParName(19,"<v2>BGEXP");
  fitMeanv2->SetParName(20,"<v2>PsiP");

  fitMeanv2->FixParameter(0,a);
  fitMeanv2->FixParameter(1,b);
  fitMeanv2->FixParameter(2,c);
  fitMeanv2->FixParameter(3,ap);
  fitMeanv2->FixParameter(4,bp);
  fitMeanv2->FixParameter(5,cp);
  fitMeanv2->FixParameter(6,dp);
  fitMeanv2->FixParameter(7,kJPsi);
  fitMeanv2->FixParameter(8,mJPsi);
  fitMeanv2->FixParameter(9,sJPsi);
  fitMeanv2->FixParameter(10,alphaLow);
  fitMeanv2->FixParameter(11,nLow);
  fitMeanv2->FixParameter(12,alphaUp);
  fitMeanv2->FixParameter(13,nUp);
  fitMeanv2->FixParameter(14,kPsiP);

  fitMeanv2->SetParameter(15, 0.01);
  fitMeanv2->SetParLimits(15, -1.,1.);

  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanv2->SetParameter(i + 16, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(20, 0.01);
  fitMeanv2->SetParLimits(20, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 16, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(20) <= fitMeanv2->GetParError(20) || fitMeanv2->GetParError(20) >= 0.75*fitMeanv2->GetParameter(20) )
  {
    fitMeanv2->FixParameter(20, 0.);
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 16, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(8),fitMeanv2->GetParError(8));
  Set("sJPsi",fitMeanv2->GetParameter(9),fitMeanv2->GetParError(9));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 4; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+16));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(15),fitMeanv2->GetParError(15));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(20),fitMeanv2->GetParError(20));

}

//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2POL2POL3_BKGMV2POL4Cheb()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SER"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SER"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");

  Double_t a  = GetValue("a");
  Double_t b  = GetValue("b");
  Double_t c  = GetValue("c");
  Double_t ap = GetValue("a'");
  Double_t bp = GetValue("b'");
  Double_t cp = GetValue("c'");
  Double_t dp = GetValue("d'");

  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");


  Double_t BG0_init    = IsValidValue(GetValue("BG0_init"))  ? GetValue("BG0_init")  : .3;
  Double_t BG1_init    = IsValidValue(GetValue("BG1_init"))  ? GetValue("BG1_init")  : -.08;
  Double_t BG2_init    = IsValidValue(GetValue("BG2_init"))  ? GetValue("BG2_init")  : 0.;
  Double_t BG3_init    = IsValidValue(GetValue("BG3_init"))  ? GetValue("BG3_init")  : 0.;
  Double_t BG4_init    = IsValidValue(GetValue("BG4_init"))  ? GetValue("BG4_init")  : 0.;
  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.6;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 4.;

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : alJPsi = %f nlJPsi = %f auJPsi = %f nuJPsi = %f",alphaLow,nLow,alphaUp,nUp));
  AliDebug(1,Form("values : a = %f b = %f c = %f a' = %f b' = %f c' = %f d' = %f ",a,b,c,ap,bp,cp,dp));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));
  if(IsValidValue(GetValue("BG0_init"))) AliInfo(Form("Bckg initialised using parameters : %f,%f,%f,%f",BG0_init,BG1_init,BG3_init,BG4_init));
  if(IsValidValue(GetValue("rejRl"))) AliInfo(Form("Bckg initialised using reject range: %f - %f",rejRl,rejRh));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Cheb,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Cheb");

  bck->SetParameters(BG0_init,BG1_init,BG2_init,BG3_init,BG4_init);

  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol4Cheb",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2POL2POL3POL4Cheb,fitRangeLow,fitRangeHigh,22,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSCB2POL2POL3POL4Cheb");

  fitMeanv2->SetParNames("a","b","c","ap","bp","cp","dp","kJPsi","mJPsi","sJPsi","alJPsi");
  //                      0   1   2    3    4    5   6       7       8       9      10
  fitMeanv2->SetParName(11,"nlJPsi");
  fitMeanv2->SetParName(12,"auJPsi");
  fitMeanv2->SetParName(13,"nuJPsi");
  fitMeanv2->SetParName(14,"kPsiP");
  fitMeanv2->SetParName(15,"<v2>JPsi");
  fitMeanv2->SetParName(16,"<v2>BGC0");
  fitMeanv2->SetParName(17,"<v2>BGC1");
  fitMeanv2->SetParName(18,"<v2>BGC2");
  fitMeanv2->SetParName(19,"<v2>BGC3");
  fitMeanv2->SetParName(20,"<v2>BGC4");
  fitMeanv2->SetParName(21,"<v2>PsiP");

  fitMeanv2->FixParameter(0,a);
  fitMeanv2->FixParameter(1,b);
  fitMeanv2->FixParameter(2,c);
  fitMeanv2->FixParameter(3,ap);
  fitMeanv2->FixParameter(4,bp);
  fitMeanv2->FixParameter(5,cp);
  fitMeanv2->FixParameter(6,dp);
  fitMeanv2->FixParameter(7,kJPsi);
  fitMeanv2->FixParameter(8,mJPsi);
  fitMeanv2->FixParameter(9,sJPsi);
  fitMeanv2->FixParameter(10,alphaLow);
  fitMeanv2->FixParameter(11,nLow);
  fitMeanv2->FixParameter(12,alphaUp);
  fitMeanv2->FixParameter(13,nUp);
  fitMeanv2->FixParameter(14,kPsiP);

  fitMeanv2->SetParameter(15, 0.01);
  fitMeanv2->SetParLimits(15, -1.,1.);

  for ( Int_t i = 0; i < 5; ++i )
  {
    fitMeanv2->SetParameter(i + 16, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(21, 0.01);
  fitMeanv2->SetParLimits(21, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 16, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(21) <= fitMeanv2->GetParError(21) || fitMeanv2->GetParError(21) >= 0.75*fitMeanv2->GetParameter(21) )
  {
    fitMeanv2->FixParameter(21, 0.);
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 16, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(8),fitMeanv2->GetParError(8));
  Set("sJPsi",fitMeanv2->GetParameter(9),fitMeanv2->GetParError(9));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 5; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+16));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(15),fitMeanv2->GetParError(15));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(21),fitMeanv2->GetParError(21));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMECB2POL2POL3_BKGMV2POL4()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SERI"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SERI"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t alphaLow = GetValue("alJPsi");
  Double_t nLow = GetValue("nlJPsi");
  Double_t alphaUp = GetValue("auJPsi");
  Double_t nUp = GetValue("nuJPsi");

  Double_t a  = GetValue("a");
  Double_t b  = GetValue("b");
  Double_t c  = GetValue("c");
  Double_t ap = GetValue("a'");
  Double_t bp = GetValue("b'");
  Double_t cp = GetValue("c'");
  Double_t dp = GetValue("d'");

  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : alJPsi = %f nlJPsi = %f auJPsi = %f nuJPsi = %f",alphaLow,nLow,alphaUp,nUp));
  AliDebug(1,Form("values : a = %f b = %f c = %f a' = %f b' = %f c' = %f d' = %f ",a,b,c,ap,bp,cp,dp));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4");

  bck->SetParameters(.4,-.2,.02);
  SetFitRejectRange(2.6,3.4);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol2",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSCB2POL2POL3POL4,fitRangeLow,fitRangeHigh,22,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSCB2POL2POL3POL4");

  fitMeanv2->SetParNames("a","b","c","ap","bp","cp","dp","kJPsi","mJPsi","sJPsi","alJPsi");
  //                      0   1   2    3    4    5   6       7       8       9      10
  fitMeanv2->SetParName(11,"nlJPsi");
  fitMeanv2->SetParName(12,"auJPsi");
  fitMeanv2->SetParName(13,"nuJPsi");
  fitMeanv2->SetParName(14,"kPsiP");
  fitMeanv2->SetParName(15,"<v2>JPsi");
  fitMeanv2->SetParName(16,"<v2>BG0");
  fitMeanv2->SetParName(17,"<v2>BG1");
  fitMeanv2->SetParName(18,"<v2>BG2");
  fitMeanv2->SetParName(19,"<v2>BG3");
  fitMeanv2->SetParName(20,"<v2>BG4");
  fitMeanv2->SetParName(21,"<v2>PsiP");

  fitMeanv2->FixParameter(0,a);
  fitMeanv2->FixParameter(1,b);
  fitMeanv2->FixParameter(2,c);
  fitMeanv2->FixParameter(3,ap);
  fitMeanv2->FixParameter(4,bp);
  fitMeanv2->FixParameter(5,cp);
  fitMeanv2->FixParameter(6,dp);
  fitMeanv2->FixParameter(7,kJPsi);
  fitMeanv2->FixParameter(8,mJPsi);
  fitMeanv2->FixParameter(9,sJPsi);
  fitMeanv2->FixParameter(10,alphaLow);
  fitMeanv2->FixParameter(11,nLow);
  fitMeanv2->FixParameter(12,alphaUp);
  fitMeanv2->FixParameter(13,nUp);
  fitMeanv2->FixParameter(14,kPsiP);

  fitMeanv2->SetParameter(15, 0.01);
  fitMeanv2->SetParLimits(15, -1.,1.);

  for ( Int_t i = 0; i < 5; ++i )
  {
    fitMeanv2->SetParameter(i + 16, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(21, 0.01);
  fitMeanv2->SetParLimits(21, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 16, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(21) <= fitMeanv2->GetParError(21) || fitMeanv2->GetParError(21) >= 0.75*fitMeanv2->GetParameter(21) )
  {
    fitMeanv2->FixParameter(21, 0.);
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 16, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(8),fitMeanv2->GetParError(8));
  Set("sJPsi",fitMeanv2->GetParameter(9),fitMeanv2->GetParError(9));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 5; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+16));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(15),fitMeanv2->GetParError(15));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(21),fitMeanv2->GetParError(21));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMENA60NEWVWG2_BKGMV2POL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SERI"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SERI"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");

  Double_t kVWG2 = GetValue("kVWG2");
  Double_t mVWG2 = GetValue("mVWG2");
  Double_t s1VWG2 = GetValue("s1VWG2");
  Double_t s2VWG2 = GetValue("s2VWG2");
  Double_t gVWG2 = GetValue("gVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");


  Double_t BG0_init    = IsValidValue(GetValue("BG0_init"))  ? GetValue("BG0_init")  : 0.4;
  Double_t BG1_init    = IsValidValue(GetValue("BG1_init"))  ? GetValue("BG1_init")  : -.2;
  Double_t BG2_init    = IsValidValue(GetValue("BG2_init"))  ? GetValue("BG2_init")  : .02;
  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.6;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 3.6;

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);



  AliDebug(1,Form("values : p1JPsi = %f p2LJPsi = %f p3LJPsi = %f p1RJPsi = %f p2RJPsi = %f p3RJPsi = %f  aLJPsi= %f  aRJPsi= %f",p1Left,p2Left,p3Left,p1Right,p2Right,p3Right,alphaLeft,alphaRight));
  AliDebug(1,Form("values : kVWG2 = %f mVG2 = %f s1VWG2 = %f s2VWG2 = %f gVWG2 = %f",kVWG2, mVWG2, s1VWG2,s2VWG2,gVWG2));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));
  if(IsValidValue(GetValue("BG0_init"))) AliInfo(Form("Bckg initialised using parameters : %f,%f,%f",BG0_init,BG1_init,BG2_init));
  if(IsValidValue(GetValue("rejRl"))) AliInfo(Form("Bckg initialised using reject range: %f - %f",rejRl,rejRh));

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(BG0_init,BG1_init,BG2_init);
  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol2",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWVWG2POL2,fitRangeLow,fitRangeHigh,22,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSNA60NEWVWG2POL2");

  fitMeanv2->SetParNames("kVWG2","mVWG2","s1VWG2","s2VWG2","gVWG2","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi");
  //                        0      1      2       3       4       5       6       7         8        9        10
  fitMeanv2->SetParName(11,"p1RJPsi");
  fitMeanv2->SetParName(12,"p2RJPsi");
  fitMeanv2->SetParName(13,"p3RJPsi");
  fitMeanv2->SetParName(14,"aLJPsi");
  fitMeanv2->SetParName(15,"aRJPsi");
  fitMeanv2->SetParName(16,"kPsiP");
  fitMeanv2->SetParName(17,"<v2>JPsi");
  fitMeanv2->SetParName(18,"<v2>BG0");
  fitMeanv2->SetParName(19,"<v2>BG1");
  fitMeanv2->SetParName(20,"<v2>BG2");
  fitMeanv2->SetParName(21,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG2);
  fitMeanv2->FixParameter(1,mVWG2);
  fitMeanv2->FixParameter(2,s1VWG2);
  fitMeanv2->FixParameter(3,s2VWG2);
  fitMeanv2->FixParameter(4,gVWG2);
  fitMeanv2->FixParameter(5,kJPsi);
  fitMeanv2->FixParameter(6,mJPsi);
  fitMeanv2->FixParameter(7,sJPsi);

  fitMeanv2->FixParameter(8,p1Left);
  fitMeanv2->FixParameter(9,p2Left);
  fitMeanv2->FixParameter(10,p3Left);
  fitMeanv2->FixParameter(11,p1Right);
  fitMeanv2->FixParameter(12,p2Right);
  fitMeanv2->FixParameter(13,p3Right);
  fitMeanv2->FixParameter(14,alphaLeft);
  fitMeanv2->FixParameter(15,alphaRight);

  fitMeanv2->FixParameter(16,kPsiP);

  fitMeanv2->SetParameter(17, 0.01);
  fitMeanv2->SetParLimits(17, -1.,1.);

  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanv2->SetParameter(i + 18, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(21, 0.01);
  fitMeanv2->SetParLimits(21, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 18, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(21) <= fitMeanv2->GetParError(21) || fitMeanv2->GetParError(21) >= 0.75*fitMeanv2->GetParameter(21) )
  {
    fitMeanv2->FixParameter(21, bck->Eval(3.68));
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 18, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  Set("sJPsi",fitMeanv2->GetParameter(7),fitMeanv2->GetParError(7));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 3; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+18));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(17),fitMeanv2->GetParError(17));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(17),fitMeanv2->GetParError(17));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMENA60NEWVWG2_BKGMV2POLEXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SER"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SER"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");

  Double_t kVWG2 = GetValue("kVWG2");
  Double_t mVWG2 = GetValue("mVWG2");
  Double_t s1VWG2 = GetValue("s1VWG2");
  Double_t s2VWG2 = GetValue("s2VWG2");
  Double_t gVWG2 = GetValue("gVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t BG0_init    = IsValidValue(GetValue("BG0_init"))  ? GetValue("BG0_init")  : .2;
  Double_t BG1_init    = IsValidValue(GetValue("BG1_init"))  ? GetValue("BG1_init")  : -.05;
  Double_t BGexp_init    = IsValidValue(GetValue("BGexp_init"))  ? GetValue("BGexp_init")  : -0.4;
  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.6;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 4.;


  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);



  AliDebug(1,Form("values : p1JPsi = %f p2LJPsi = %f p3LJPsi = %f p1RJPsi = %f p2RJPsi = %f p3RJPsi = %f  aLJPsi= %f  aRJPsi= %f",p1Left,p2Left,p3Left,p1Right,p2Right,p3Right,alphaLeft,alphaRight));
  AliDebug(1,Form("values : kVWG2 = %f mVG2 = %f s1VWG2 = %f s2VWG2 = %f gVWG2 = %f",kVWG2, mVWG2, s1VWG2,s2VWG2,gVWG2));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));
  if(IsValidValue(GetValue("BG0_init"))) AliInfo(Form("Bckg initialised using parameters : %f,%f,%f",BG0_init,BG1_init,BGexp_init));
  if(IsValidValue(GetValue("rejRl"))) AliInfo(Form("Bckg initialised using reject range: %f - %f",rejRl,rejRh));
  //  TString resultName(Form("MEANV2FIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPolExp,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPolExp");

  bck->SetParameters(BG0_init,BG1_init,BGexp_init);
  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol2",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWVWG2POLEXP,fitRangeLow,fitRangeHigh,22,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSNA60NEWVWG2POLEXP");

  fitMeanv2->SetParNames("kVWG2","mVWG2","s1VWG2","s2VWG2","gVWG2","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi");
  //                        0      1      2       3       4       5       6       7         8        9        10
  fitMeanv2->SetParName(11,"p1RJPsi");
  fitMeanv2->SetParName(12,"p2RJPsi");
  fitMeanv2->SetParName(13,"p3RJPsi");
  fitMeanv2->SetParName(14,"aLJPsi");
  fitMeanv2->SetParName(15,"aRJPsi");
  fitMeanv2->SetParName(16,"kPsiP");
  fitMeanv2->SetParName(17,"<v2>JPsi");
  fitMeanv2->SetParName(18,"<v2>BG0");
  fitMeanv2->SetParName(19,"<v2>BG1");
  fitMeanv2->SetParName(20,"<v2>BGEXP");
  fitMeanv2->SetParName(21,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG2);
  fitMeanv2->FixParameter(1,mVWG2);
  fitMeanv2->FixParameter(2,s1VWG2);
  fitMeanv2->FixParameter(3,s2VWG2);
  fitMeanv2->FixParameter(4,gVWG2);
  fitMeanv2->FixParameter(5,kJPsi);
  fitMeanv2->FixParameter(6,mJPsi);
  fitMeanv2->FixParameter(7,sJPsi);

  fitMeanv2->FixParameter(8,p1Left);
  fitMeanv2->FixParameter(9,p2Left);
  fitMeanv2->FixParameter(10,p3Left);
  fitMeanv2->FixParameter(11,p1Right);
  fitMeanv2->FixParameter(12,p2Right);
  fitMeanv2->FixParameter(13,p3Right);
  fitMeanv2->FixParameter(14,alphaLeft);
  fitMeanv2->FixParameter(15,alphaRight);

  fitMeanv2->FixParameter(16,kPsiP);

  fitMeanv2->SetParameter(17, 0.01);
  fitMeanv2->SetParLimits(17, -1.,1.);

  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanv2->SetParameter(i + 18, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(21, 0.01);
  fitMeanv2->SetParLimits(21, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 18, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(21) <= fitMeanv2->GetParError(21) || fitMeanv2->GetParError(21) >= 0.75*fitMeanv2->GetParameter(21) )
  {
    fitMeanv2->FixParameter(21, bck->Eval(3.68));
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 18, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  Set("sJPsi",fitMeanv2->GetParameter(7),fitMeanv2->GetParError(7));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 3; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+18));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//
  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(17),fitMeanv2->GetParError(17));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(17),fitMeanv2->GetParError(17));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMENA60NEWVWG2_BKGMV2POL2EXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SERI"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SERI"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");

  Double_t kVWG2 = GetValue("kVWG2");
  Double_t mVWG2 = GetValue("mVWG2");
  Double_t s1VWG2 = GetValue("s1VWG2");
  Double_t s2VWG2 = GetValue("s2VWG2");
  Double_t gVWG2 = GetValue("gVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);


  //  TString resultName(Form("MEANV2FIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");

  bck->SetParameters(1.8,-1.,.1,-.7);
  bck->SetParLimits(0,0.,10.);
  SetFitRejectRange(2.4,3.6);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol2Exp",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWVWG2POL2EXP,fitRangeLow,fitRangeHigh,23,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSNA60NEWVWG2POL2EXP");

  fitMeanv2->SetParNames("kVWG2","mVWG2","s1VWG2","s2VWG2","gVWG2","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi");
  //                        0      1      2       3       4       5       6       7         8        9        10
  fitMeanv2->SetParName(11,"p1RJPsi");
  fitMeanv2->SetParName(12,"p2RJPsi");
  fitMeanv2->SetParName(13,"p3RJPsi");
  fitMeanv2->SetParName(14,"aLJPsi");
  fitMeanv2->SetParName(15,"aRJPsi");
  fitMeanv2->SetParName(16,"kPsiP");
  fitMeanv2->SetParName(17,"<v2>JPsi");
  fitMeanv2->SetParName(18,"<v2>BG0");
  fitMeanv2->SetParName(19,"<v2>BG1");
  fitMeanv2->SetParName(20,"<v2>BG2");
  fitMeanv2->SetParName(21,"<v2>BGEXP");
  fitMeanv2->SetParName(22,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG2);
  fitMeanv2->FixParameter(1,mVWG2);
  fitMeanv2->FixParameter(2,s1VWG2);
  fitMeanv2->FixParameter(3,s2VWG2);
  fitMeanv2->FixParameter(4,gVWG2);
  fitMeanv2->FixParameter(5,kJPsi);
  fitMeanv2->FixParameter(6,mJPsi);
  fitMeanv2->FixParameter(7,sJPsi);

  fitMeanv2->FixParameter(8,p1Left);
  fitMeanv2->FixParameter(9,p2Left);
  fitMeanv2->FixParameter(10,p3Left);
  fitMeanv2->FixParameter(11,p1Right);
  fitMeanv2->FixParameter(12,p2Right);
  fitMeanv2->FixParameter(13,p3Right);
  fitMeanv2->FixParameter(14,alphaLeft);
  fitMeanv2->FixParameter(15,alphaRight);

  fitMeanv2->FixParameter(16,kPsiP);

  fitMeanv2->SetParameter(17, 0.01);
  fitMeanv2->SetParLimits(17, -1.,1.);

  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanv2->SetParameter(i + 18, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(22, 0.01);
  fitMeanv2->SetParLimits(22, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 18, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(22) <= fitMeanv2->GetParError(22) || fitMeanv2->GetParError(22) >= 0.75*fitMeanv2->GetParameter(22) )
  {
    fitMeanv2->FixParameter(22, 0.);
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 18, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  Set("sJPsi",fitMeanv2->GetParameter(7),fitMeanv2->GetParError(7));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 4; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+18));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(17),fitMeanv2->GetParError(17));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(17),fitMeanv2->GetParError(17));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMENA60NEWVWG2_BKGMV2POL4()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SERI"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SERI"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");

  Double_t kVWG2 = GetValue("kVWG2");
  Double_t mVWG2 = GetValue("mVWG2");
  Double_t s1VWG2 = GetValue("s1VWG2");
  Double_t s2VWG2 = GetValue("s2VWG2");
  Double_t gVWG2 = GetValue("gVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);


  //  TString resultName(Form("MEANV2FIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4");

  bck->SetParameters(.6,-.4,.03,0.015,-0.003);
  bck->SetParLimits(0,0.,2.);
  SetFitRejectRange(2.5,3.5);
  // SetFitRejectRange(2.6,4.0);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol4",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWVWG2POL4,fitRangeLow,fitRangeHigh,24,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSNA60NEWVWG2POL4");

  fitMeanv2->SetParNames("kVWG2","mVWG2","s1VWG2","s2VWG2","gVWG2","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi");
  //                        0      1      2       3       4       5       6       7         8        9        10
  fitMeanv2->SetParName(11,"p1RJPsi");
  fitMeanv2->SetParName(12,"p2RJPsi");
  fitMeanv2->SetParName(13,"p3RJPsi");
  fitMeanv2->SetParName(14,"aLJPsi");
  fitMeanv2->SetParName(15,"aRJPsi");
  fitMeanv2->SetParName(16,"kPsiP");
  fitMeanv2->SetParName(17,"<v2>JPsi");
  fitMeanv2->SetParName(18,"<v2>BG0");
  fitMeanv2->SetParName(19,"<v2>BG1");
  fitMeanv2->SetParName(20,"<v2>BG2");
  fitMeanv2->SetParName(21,"<v2>BG3");
  fitMeanv2->SetParName(22,"<v2>BG4");
  fitMeanv2->SetParName(23,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG2);
  fitMeanv2->FixParameter(1,mVWG2);
  fitMeanv2->FixParameter(2,s1VWG2);
  fitMeanv2->FixParameter(3,s2VWG2);
  fitMeanv2->FixParameter(4,gVWG2);
  fitMeanv2->FixParameter(5,kJPsi);
  fitMeanv2->FixParameter(6,mJPsi);
  fitMeanv2->FixParameter(7,sJPsi);

  fitMeanv2->FixParameter(8,p1Left);
  fitMeanv2->FixParameter(9,p2Left);
  fitMeanv2->FixParameter(10,p3Left);
  fitMeanv2->FixParameter(11,p1Right);
  fitMeanv2->FixParameter(12,p2Right);
  fitMeanv2->FixParameter(13,p3Right);
  fitMeanv2->FixParameter(14,alphaLeft);
  fitMeanv2->FixParameter(15,alphaRight);

  fitMeanv2->FixParameter(16,kPsiP);

  fitMeanv2->SetParameter(17, 0.01);
  fitMeanv2->SetParLimits(17, -1.,1.);

  for ( Int_t i = 0; i < 5; ++i )
  {
    fitMeanv2->SetParameter(i + 18, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(23, 0.01);
  fitMeanv2->SetParLimits(23, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 18, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(23) <= fitMeanv2->GetParError(23) || fitMeanv2->GetParError(23) >= 0.75*fitMeanv2->GetParameter(23) )
  {
    fitMeanv2->FixParameter(23, bck->Eval(3.68));
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 18, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  Set("sJPsi",fitMeanv2->GetParameter(7),fitMeanv2->GetParError(7));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 5; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+18));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(17),fitMeanv2->GetParError(17));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(17),fitMeanv2->GetParError(17));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMENA60NEWVWG2_BKGMV2POL4Cheb()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SER"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SER"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();

  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");

  Double_t kVWG2 = GetValue("kVWG2");
  Double_t mVWG2 = GetValue("mVWG2");
  Double_t s1VWG2 = GetValue("s1VWG2");
  Double_t s2VWG2 = GetValue("s2VWG2");
  Double_t gVWG2 = GetValue("gVWG2");
  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t BG0_init    = IsValidValue(GetValue("BG0_init"))  ? GetValue("BG0_init")  : 0.4;
  Double_t BG1_init    = IsValidValue(GetValue("BG1_init"))  ? GetValue("BG1_init")  : -.2;
  Double_t BG2_init    = IsValidValue(GetValue("BG2_init"))  ? GetValue("BG2_init")  : .02;
  Double_t BG3_init    = IsValidValue(GetValue("BG3_init"))  ? GetValue("BG3_init")  : -.2;
  Double_t BG4_init    = IsValidValue(GetValue("BG4_init"))  ? GetValue("BG4_init")  : .02;
  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.6;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 3.6;
  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : p1JPsi = %f p2LJPsi = %f p3LJPsi = %f p1RJPsi = %f p2RJPsi = %f p3RJPsi = %f  aLJPsi= %f  aRJPsi= %f",p1Left,p2Left,p3Left,p1Right,p2Right,p3Right,alphaLeft,alphaRight));
  AliDebug(1,Form("values : kVWG2 = %f mVG2 = %f s1VWG2 = %f s2VWG2 = %f gVWG2 = %f",kVWG2, mVWG2, s1VWG2,s2VWG2,gVWG2));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));
  if(IsValidValue(GetValue("BG0_init"))) AliInfo(Form("Bckg initialised using parameters : %f,%f,%f,%f",BG0_init,BG1_init,BG3_init,BG4_init));
  if(IsValidValue(GetValue("rejRl"))) AliInfo(Form("Bckg initialised using reject range: %f - %f",rejRl,rejRh));
  //  TString resultName(Form("MEANV2FIT_%salphalow=%5.2fnlow=%5.2falphaup=%5.2fnup=%5.2f",fitName.Data(),par[7],par[8],par[9],par[10]));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Cheb,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Cheb");
  bck->SetParameters(BG0_init,BG1_init,BG2_init,BG3_init,BG4_init);
  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol4Cheb",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWVWG2POL4Cheb,fitRangeLow,fitRangeHigh,24,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSNA60NEWVWG2POL4Cheb");

  fitMeanv2->SetParNames("kVWG2","mVWG2","s1VWG2","s2VWG2","gVWG2","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi");
  //                        0      1      2       3       4       5       6       7         8        9        10
  fitMeanv2->SetParName(11,"p1RJPsi");
  fitMeanv2->SetParName(12,"p2RJPsi");
  fitMeanv2->SetParName(13,"p3RJPsi");
  fitMeanv2->SetParName(14,"aLJPsi");
  fitMeanv2->SetParName(15,"aRJPsi");
  fitMeanv2->SetParName(16,"kPsiP");
  fitMeanv2->SetParName(17,"<v2>JPsi");
  fitMeanv2->SetParName(18,"<v2>BGC0");
  fitMeanv2->SetParName(19,"<v2>BGC1");
  fitMeanv2->SetParName(20,"<v2>BGC2");
  fitMeanv2->SetParName(21,"<v2>BGC3");
  fitMeanv2->SetParName(22,"<v2>BGC4");
  fitMeanv2->SetParName(23,"<v2>PsiP");

  fitMeanv2->FixParameter(0,kVWG2);
  fitMeanv2->FixParameter(1,mVWG2);
  fitMeanv2->FixParameter(2,s1VWG2);
  fitMeanv2->FixParameter(3,s2VWG2);
  fitMeanv2->FixParameter(4,gVWG2);
  fitMeanv2->FixParameter(5,kJPsi);
  fitMeanv2->FixParameter(6,mJPsi);
  fitMeanv2->FixParameter(7,sJPsi);

  fitMeanv2->FixParameter(8,p1Left);
  fitMeanv2->FixParameter(9,p2Left);
  fitMeanv2->FixParameter(10,p3Left);
  fitMeanv2->FixParameter(11,p1Right);
  fitMeanv2->FixParameter(12,p2Right);
  fitMeanv2->FixParameter(13,p3Right);
  fitMeanv2->FixParameter(14,alphaLeft);
  fitMeanv2->FixParameter(15,alphaRight);

  fitMeanv2->FixParameter(16,kPsiP);

  fitMeanv2->SetParameter(17, 0.01);
  fitMeanv2->SetParLimits(17, -1.,1.);

  for ( Int_t i = 0; i < 5; ++i )
  {
    fitMeanv2->SetParameter(i + 18, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(23, 0.01);
  fitMeanv2->SetParLimits(23, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 18, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(23) <= fitMeanv2->GetParError(23) || fitMeanv2->GetParError(23) >= 0.75*fitMeanv2->GetParameter(23) )
  {
    fitMeanv2->FixParameter(23, bck->Eval(3.68));
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 18, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(6),fitMeanv2->GetParError(6));
  Set("sJPsi",fitMeanv2->GetParameter(7),fitMeanv2->GetParError(7));
  // Set("mv2PsiP",fitMeanv2->GetParameter(16),fitMeanv2->GetParError(16));

  for ( Int_t i = 0; i < 5; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+18));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(17),fitMeanv2->GetParError(17));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(17),fitMeanv2->GetParError(17));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMENA60NEWPOL2POL3_BKGMV2POL2()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SER"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SERI"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();


  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");

  Double_t a  = GetValue("a");
  Double_t b  = GetValue("b");
  Double_t c  = GetValue("c");
  Double_t ap = GetValue("a'");
  Double_t bp = GetValue("b'");
  Double_t cp = GetValue("c'");
  Double_t dp = GetValue("d'");

  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t BG0_init    = IsValidValue(GetValue("BG0_init"))  ? GetValue("BG0_init")  : 0.4;
  Double_t BG1_init    = IsValidValue(GetValue("BG1_init"))  ? GetValue("BG1_init")  : -.2;
  Double_t BG2_init    = IsValidValue(GetValue("BG2_init"))  ? GetValue("BG2_init")  : .03;
  Double_t rejRl       = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.6;
  Double_t rejRh       = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 3.6;

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : p1JPsi = %f p2LJPsi = %f p3LJPsi = %f p1RJPsi = %f p2RJPsi = %f p3RJPsi = %f  aLJPsi= %f  aRJPsi= %f",p1Left,p2Left,p3Left,p1Right,p2Right,p3Right,alphaLeft,alphaRight));
  AliDebug(1,Form("values : a = %f b = %f c = %f a' = %f b' = %f c' = %f d' = %f ",a,b,c,ap,bp,cp,dp));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));
  if(IsValidValue(GetValue("BG0_init"))) AliInfo(Form("Bckg initialised using parameters : %f,%f,%f",BG0_init,BG1_init,BG2_init));
  if(IsValidValue(GetValue("rejRl"))) AliInfo(Form("Bckg initialised using reject range: %f - %f",rejRl,rejRh));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2");

  bck->SetParameters(BG0_init,BG1_init,BG2_init);
  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol2",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWPOL2POL3POL2,fitRangeLow,fitRangeHigh,24,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSNA60NEWPOL2POL3POL2");

  fitMeanv2->SetParNames("a","b","c","ap","bp","cp","dp","kJPsi","mJPsi","sJPsi","p1LJPsi");
  //                      0   1   2    3    4    5   6       7       8       9      10

  fitMeanv2->SetParName(11,"p2LJPsi");
  fitMeanv2->SetParName(12,"p3LJPsi");
  fitMeanv2->SetParName(13,"p1RJPsi");
  fitMeanv2->SetParName(14,"p2RJPsi");
  fitMeanv2->SetParName(15,"p3RJPsi");
  fitMeanv2->SetParName(16,"aLJPsi");
  fitMeanv2->SetParName(17,"aRJPsi");

  fitMeanv2->SetParName(18,"kPsiP");
  fitMeanv2->SetParName(19,"<v2>JPsi");
  fitMeanv2->SetParName(20,"<v2>BG0");
  fitMeanv2->SetParName(21,"<v2>BG1");
  fitMeanv2->SetParName(22,"<v2>BG2");
  fitMeanv2->SetParName(23,"<v2>PsiP");

  fitMeanv2->FixParameter(0,a);
  fitMeanv2->FixParameter(1,b);
  fitMeanv2->FixParameter(2,c);
  fitMeanv2->FixParameter(3,ap);
  fitMeanv2->FixParameter(4,bp);
  fitMeanv2->FixParameter(5,cp);
  fitMeanv2->FixParameter(6,dp);
  fitMeanv2->FixParameter(7,kJPsi);
  fitMeanv2->FixParameter(8,mJPsi);
  fitMeanv2->FixParameter(9,sJPsi);
  fitMeanv2->FixParameter(10,p1Left);
  fitMeanv2->FixParameter(11,p2Left);
  fitMeanv2->FixParameter(12,p3Left);
  fitMeanv2->FixParameter(13,p1Right);
  fitMeanv2->FixParameter(14,p2Right);
  fitMeanv2->FixParameter(15,p3Right);
  fitMeanv2->FixParameter(16,alphaLeft);
  fitMeanv2->FixParameter(17,alphaRight);

  fitMeanv2->FixParameter(18,kPsiP);

  fitMeanv2->SetParameter(19, 0.01);
  fitMeanv2->SetParLimits(19, -1.,1.);

  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanv2->SetParameter(i + 20, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(23, 0.01);
  fitMeanv2->SetParLimits(23, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 20, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(23) <= fitMeanv2->GetParError(23) || fitMeanv2->GetParError(23) >= 0.75*fitMeanv2->GetParameter(23) )
  {
    fitMeanv2->FixParameter(23, 0.);
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 20, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(8),fitMeanv2->GetParError(8));
  Set("sJPsi",fitMeanv2->GetParameter(9),fitMeanv2->GetParError(9));
  // Set("mv2PsiP",fitMeanv2->GetParameter(23),fitMeanv2->GetParError(23));

  for ( Int_t i = 0; i < 3; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+20));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(19),fitMeanv2->GetParError(19));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(23),fitMeanv2->GetParError(23));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMENA60NEWPOL2POL3_BKGMV2POLEXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SER"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SER"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();


  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");

  Double_t a  = GetValue("a");
  Double_t b  = GetValue("b");
  Double_t c  = GetValue("c");
  Double_t ap = GetValue("a'");
  Double_t bp = GetValue("b'");
  Double_t cp = GetValue("c'");
  Double_t dp = GetValue("d'");

  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t BG0_init    = IsValidValue(GetValue("BG0_init"))  ? GetValue("BG0_init")  : 1.;
  Double_t BG1_init    = IsValidValue(GetValue("BG1_init"))  ? GetValue("BG1_init")  : -3.;
  Double_t BGexp_init    = IsValidValue(GetValue("BGexp_init"))  ? GetValue("BGexp_init")  : -1.7;
  Double_t rejRl    = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.6;
  Double_t rejRh    = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 4.;

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : p1JPsi = %f p2LJPsi = %f p3LJPsi = %f p1RJPsi = %f p2RJPsi = %f p3RJPsi = %f  aLJPsi= %f  aRJPsi= %f",p1Left,p2Left,p3Left,p1Right,p2Right,p3Right,alphaLeft,alphaRight));
  AliDebug(1,Form("values : a = %f b = %f c = %f a' = %f b' = %f c' = %f d' = %f ",a,b,c,ap,bp,cp,dp));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));
  if(IsValidValue(GetValue("BG0_init"))) AliInfo(Form("Bckg initialised using parameters : %f,%f,%f",BG0_init,BG1_init,BGexp_init));
  if(IsValidValue(GetValue("rejRl"))) AliInfo(Form("Bckg initialised using reject range: %f - %f",rejRl,rejRh));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPolExp,fitRangeLow,fitRangeHigh,3,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPolExp");

  // bck->SetParameters(1.2,-.3,-1.);
  bck->SetParameters(BG0_init,BG1_init,BGexp_init);

   // bck->SetParameters(-.15,.2,-.004);

  SetFitRejectRange(rejRl,rejRh);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPolExp",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWPOL2POL3_POLEXP,fitRangeLow,fitRangeHigh,24,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSNA60NEWPOL2POL3_POLEXP");

  fitMeanv2->SetParNames("a","b","c","ap","bp","cp","dp","kJPsi","mJPsi","sJPsi","p1LJPsi");
  //                      0   1   2    3    4    5   6       7       8       9      10

  fitMeanv2->SetParName(11,"p2LJPsi");
  fitMeanv2->SetParName(12,"p3LJPsi");
  fitMeanv2->SetParName(13,"p1RJPsi");
  fitMeanv2->SetParName(14,"p2RJPsi");
  fitMeanv2->SetParName(15,"p3RJPsi");
  fitMeanv2->SetParName(16,"aLJPsi");
  fitMeanv2->SetParName(17,"aRJPsi");

  fitMeanv2->SetParName(18,"kPsiP");
  fitMeanv2->SetParName(19,"<v2>JPsi");
  fitMeanv2->SetParName(20,"<v2>BG0");
  fitMeanv2->SetParName(21,"<v2>BG1");
  fitMeanv2->SetParName(22,"<v2>BGEXP");
  fitMeanv2->SetParName(23,"<v2>PsiP");

  fitMeanv2->FixParameter(0,a);
  fitMeanv2->FixParameter(1,b);
  fitMeanv2->FixParameter(2,c);
  fitMeanv2->FixParameter(3,ap);
  fitMeanv2->FixParameter(4,bp);
  fitMeanv2->FixParameter(5,cp);
  fitMeanv2->FixParameter(6,dp);
  fitMeanv2->FixParameter(7,kJPsi);
  fitMeanv2->FixParameter(8,mJPsi);
  fitMeanv2->FixParameter(9,sJPsi);
  fitMeanv2->FixParameter(10,p1Left);
  fitMeanv2->FixParameter(11,p2Left);
  fitMeanv2->FixParameter(12,p3Left);
  fitMeanv2->FixParameter(13,p1Right);
  fitMeanv2->FixParameter(14,p2Right);
  fitMeanv2->FixParameter(15,p3Right);
  fitMeanv2->FixParameter(16,alphaLeft);
  fitMeanv2->FixParameter(17,alphaRight);

  fitMeanv2->FixParameter(18,kPsiP);

  fitMeanv2->SetParameter(19, 0.01);
  fitMeanv2->SetParLimits(19, -1.,1.);

  for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanv2->SetParameter(i + 20, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(23, 0.01);
  fitMeanv2->SetParLimits(23, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 20, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(23) <= fitMeanv2->GetParError(23) || fitMeanv2->GetParError(23) >= 0.75*fitMeanv2->GetParameter(23) )
  {
    fitMeanv2->FixParameter(23, 0.);
    for ( Int_t i = 0; i < 3; ++i )
    {
      fitMeanv2->SetParameter(i + 20, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(8),fitMeanv2->GetParError(8));
  Set("sJPsi",fitMeanv2->GetParameter(9),fitMeanv2->GetParError(9));
  // Set("mv2PsiP",fitMeanv2->GetParameter(23),fitMeanv2->GetParError(23));

  for ( Int_t i = 0; i < 3; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+20));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(19),fitMeanv2->GetParError(19));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(23),fitMeanv2->GetParError(23));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMENA60NEWPOL2POL3_BKGMV2POL2EXP()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SERI"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SERI"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();


  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");

  Double_t a  = GetValue("a");
  Double_t b  = GetValue("b");
  Double_t c  = GetValue("c");
  Double_t ap = GetValue("a'");
  Double_t bp = GetValue("b'");
  Double_t cp = GetValue("c'");
  Double_t dp = GetValue("d'");

  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : p1JPsi = %f p2LJPsi = %f p3LJPsi = %f p1RJPsi = %f p2RJPsi = %f p3RJPsi = %f  aLJPsi= %f  aRJPsi= %f",p1Left,p2Left,p3Left,p1Right,p2Right,p3Right,alphaLeft,alphaRight));
  AliDebug(1,Form("values : a = %f b = %f c = %f a' = %f b' = %f c' = %f d' = %f ",a,b,c,ap,bp,cp,dp));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol2Exp,fitRangeLow,fitRangeHigh,4,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol2Exp");

  bck->SetParameters(1.8,-1.,.1,-.7);

  SetFitRejectRange(2.6,3.4);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol2Exp",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWPOL2POL3POL2EXP,fitRangeLow,fitRangeHigh,25,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSNA60NEWPOL2POL3POL2EXP");

  fitMeanv2->SetParNames("a","b","c","ap","bp","cp","dp","kJPsi","mJPsi","sJPsi","p1LJPsi");
  //                      0   1   2    3    4    5   6       7       8       9      10

  fitMeanv2->SetParName(11,"p2LJPsi");
  fitMeanv2->SetParName(12,"p3LJPsi");
  fitMeanv2->SetParName(13,"p1RJPsi");
  fitMeanv2->SetParName(14,"p2RJPsi");
  fitMeanv2->SetParName(15,"p3RJPsi");
  fitMeanv2->SetParName(16,"aLJPsi");
  fitMeanv2->SetParName(17,"aRJPsi");

  fitMeanv2->SetParName(18,"kPsiP");
  fitMeanv2->SetParName(19,"<v2>JPsi");
  fitMeanv2->SetParName(20,"<v2>BG0");
  fitMeanv2->SetParName(21,"<v2>BG1");
  fitMeanv2->SetParName(22,"<v2>BG2");
  fitMeanv2->SetParName(23,"<v2>BGEXP");
  fitMeanv2->SetParName(24,"<v2>PsiP");

  fitMeanv2->FixParameter(0,a);
  fitMeanv2->FixParameter(1,b);
  fitMeanv2->FixParameter(2,c);
  fitMeanv2->FixParameter(3,ap);
  fitMeanv2->FixParameter(4,bp);
  fitMeanv2->FixParameter(5,cp);
  fitMeanv2->FixParameter(6,dp);
  fitMeanv2->FixParameter(7,kJPsi);
  fitMeanv2->FixParameter(8,mJPsi);
  fitMeanv2->FixParameter(9,sJPsi);
  fitMeanv2->FixParameter(10,p1Left);
  fitMeanv2->FixParameter(11,p2Left);
  fitMeanv2->FixParameter(12,p3Left);
  fitMeanv2->FixParameter(13,p1Right);
  fitMeanv2->FixParameter(14,p2Right);
  fitMeanv2->FixParameter(15,p3Right);
  fitMeanv2->FixParameter(16,alphaLeft);
  fitMeanv2->FixParameter(17,alphaRight);

  fitMeanv2->FixParameter(18,kPsiP);

  fitMeanv2->SetParameter(19, 0.01);
  fitMeanv2->SetParLimits(19, -1.,1.);

  for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanv2->SetParameter(i + 20, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(24, 0.01);
  fitMeanv2->SetParLimits(24, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 20, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(24) <= fitMeanv2->GetParError(24) || fitMeanv2->GetParError(24) >= 0.75*fitMeanv2->GetParameter(24) )
  {
    fitMeanv2->FixParameter(24, 0.);
    for ( Int_t i = 0; i < 4; ++i )
    {
      fitMeanv2->SetParameter(i + 20, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(8),fitMeanv2->GetParError(8));
  Set("sJPsi",fitMeanv2->GetParameter(9),fitMeanv2->GetParError(9));
  // Set("mv2PsiP",fitMeanv2->GetParameter(23),fitMeanv2->GetParError(23));

  for ( Int_t i = 0; i < 4; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+20));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(19),fitMeanv2->GetParError(19));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(24),fitMeanv2->GetParError(24));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMENA60NEWPOL2POL3_BKGMV2POL4()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SERI"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SERI"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();


  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");

  Double_t a  = GetValue("a");
  Double_t b  = GetValue("b");
  Double_t c  = GetValue("c");
  Double_t ap = GetValue("a'");
  Double_t bp = GetValue("b'");
  Double_t cp = GetValue("c'");
  Double_t dp = GetValue("d'");

  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : p1JPsi = %f p2LJPsi = %f p3LJPsi = %f p1RJPsi = %f p2RJPsi = %f p3RJPsi = %f  aLJPsi= %f  aRJPsi= %f",p1Left,p2Left,p3Left,p1Right,p2Right,p3Right,alphaLeft,alphaRight));
  AliDebug(1,Form("values : a = %f b = %f c = %f a' = %f b' = %f c' = %f d' = %f ",a,b,c,ap,bp,cp,dp));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));

  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4");

  bck->SetParameters(.4,-.2,.02,0.,0.);
  SetFitRejectRange(2.6,3.4);

  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol4",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWPOL2POL3POL2,fitRangeLow,fitRangeHigh,26,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSNA60NEWPOL2POL3POL2");

  fitMeanv2->SetParNames("a","b","c","ap","bp","cp","dp","kJPsi","mJPsi","sJPsi","p1LJPsi");
  //                      0   1   2    3    4    5   6       7       8       9      10

  fitMeanv2->SetParName(11,"p2LJPsi");
  fitMeanv2->SetParName(12,"p3LJPsi");
  fitMeanv2->SetParName(13,"p1RJPsi");
  fitMeanv2->SetParName(14,"p2RJPsi");
  fitMeanv2->SetParName(15,"p3RJPsi");
  fitMeanv2->SetParName(16,"aLJPsi");
  fitMeanv2->SetParName(17,"aRJPsi");

  fitMeanv2->SetParName(18,"kPsiP");
  fitMeanv2->SetParName(19,"<v2>JPsi");
  fitMeanv2->SetParName(20,"<v2>BG0");
  fitMeanv2->SetParName(21,"<v2>BG1");
  fitMeanv2->SetParName(22,"<v2>BG2");
  fitMeanv2->SetParName(23,"<v2>BG3");
  fitMeanv2->SetParName(24,"<v2>BG4");
  fitMeanv2->SetParName(25,"<v2>PsiP");

  fitMeanv2->FixParameter(0,a);
  fitMeanv2->FixParameter(1,b);
  fitMeanv2->FixParameter(2,c);
  fitMeanv2->FixParameter(3,ap);
  fitMeanv2->FixParameter(4,bp);
  fitMeanv2->FixParameter(5,cp);
  fitMeanv2->FixParameter(6,dp);
  fitMeanv2->FixParameter(7,kJPsi);
  fitMeanv2->FixParameter(8,mJPsi);
  fitMeanv2->FixParameter(9,sJPsi);

  fitMeanv2->FixParameter(10,p1Left);
  fitMeanv2->FixParameter(11,p2Left);
  fitMeanv2->FixParameter(12,p3Left);
  fitMeanv2->FixParameter(13,p1Right);
  fitMeanv2->FixParameter(14,p2Right);
  fitMeanv2->FixParameter(15,p3Right);
  fitMeanv2->FixParameter(16,alphaLeft);
  fitMeanv2->FixParameter(17,alphaRight);

  fitMeanv2->FixParameter(18,kPsiP);

  fitMeanv2->SetParameter(19, 0.01);
  fitMeanv2->SetParLimits(19, -1.,1.);

  for ( Int_t i = 0; i < 5; ++i )
  {
    fitMeanv2->SetParameter(i + 20, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(25, 0.01);
  fitMeanv2->SetParLimits(25, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 20, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(25) <= fitMeanv2->GetParError(25) || fitMeanv2->GetParError(25) >= 0.75*fitMeanv2->GetParameter(25) )
  {
    fitMeanv2->FixParameter(25, 0.);
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 20, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(8),fitMeanv2->GetParError(8));
  Set("sJPsi",fitMeanv2->GetParameter(9),fitMeanv2->GetParError(9));
  // Set("mv2PsiP",fitMeanv2->GetParameter(25),fitMeanv2->GetParError(25));

  for ( Int_t i = 0; i < 5; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+20));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(19),fitMeanv2->GetParError(19));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(23),fitMeanv2->GetParError(23));

}
//_____________________________________________________________________________
void AliAnalysisMuMuJpsiResult::FitMV2PSIPSIPRIMENA60NEWPOL2POL3_BKGMV2POL4Cheb()
{
  //Fit mean dimuon mean pt to get Jpsi mean pt using the CB2 signal parameters, VWG2 for the Bkg, Jpsi mpt = cte and Bkg mpt = pol4

  const char* fitOption = "SER"; //We can add NO to avoid plotting
  const char* fitOptionBg = "SER"; //We can add NO to avoid plotting

  fHisto->GetListOfFunctions()->Delete();


  Double_t p1Left = GetValue("p1LJPsi");
  Double_t p2Left = GetValue("p2LJPsi");
  Double_t p3Left = GetValue("p3LJPsi");
  Double_t p1Right = GetValue("p1RJPsi");
  Double_t p2Right = GetValue("p2RJPsi");
  Double_t p3Right = GetValue("p3RJPsi");
  Double_t alphaLeft = GetValue("aLJPsi");
  Double_t alphaRight = GetValue("aRJPsi");

  Double_t a  = GetValue("a");
  Double_t b  = GetValue("b");
  Double_t c  = GetValue("c");
  Double_t ap = GetValue("a'");
  Double_t bp = GetValue("b'");
  Double_t cp = GetValue("c'");
  Double_t dp = GetValue("d'");

  Double_t kJPsi = GetValue("kJPsi");
  Double_t kPsiP = GetValue("kPsiP");
  Double_t mJPsi = GetValue("mJPsi");
  Double_t sJPsi = GetValue("sJPsi");
  Double_t NofJPsi = GetValue("NofJPsi");
  Double_t ErrStatNofJPsi = GetValue("ErrStatNofJPsi");


  Double_t BG0_init    = IsValidValue(GetValue("BG0_init"))  ? GetValue("BG0_init")  : 0.4;
  Double_t BG1_init    = IsValidValue(GetValue("BG1_init"))  ? GetValue("BG1_init")  : -.2;
  Double_t BG2_init    = IsValidValue(GetValue("BG2_init"))  ? GetValue("BG2_init")  : .03;
  Double_t BG3_init    = IsValidValue(GetValue("BG3_init"))  ? GetValue("BG3_init")  : -.2;
  Double_t BG4_init    = IsValidValue(GetValue("BG4_init"))  ? GetValue("BG4_init")  : .03;
  Double_t rejRl       = IsValidValue(GetValue("rejRl"))  ? GetValue("rejRl")  : 2.6;
  Double_t rejRh       = IsValidValue(GetValue("rejRh"))  ? GetValue("rejRh")  : 3.6;

  Double_t fitRangeLow = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);

  AliDebug(1,Form("values : p1JPsi = %f p2LJPsi = %f p3LJPsi = %f p1RJPsi = %f p2RJPsi = %f p3RJPsi = %f  aLJPsi= %f  aRJPsi= %f",p1Left,p2Left,p3Left,p1Right,p2Right,p3Right,alphaLeft,alphaRight));
  AliDebug(1,Form("values : a = %f b = %f c = %f a' = %f b' = %f c' = %f d' = %f ",a,b,c,ap,bp,cp,dp));
  AliDebug(1,Form("values : kJPsi = %f mJPsi = %f sJPsi = %f kPsiP = %f NofJPsi= %f ErrStatNofJPsi = %f",kJPsi,mJPsi,sJPsi,kPsiP,NofJPsi,ErrStatNofJPsi));
  if(IsValidValue(GetValue("BG0_init"))) AliInfo(Form("Bckg initialised using parameters : %f,%f,%f,%f",BG0_init,BG1_init,BG3_init,BG4_init));
  if(IsValidValue(GetValue("rejRl"))) AliInfo(Form("Bckg initialised using reject range: %f - %f",rejRl,rejRh));
  TProfile* p(0x0);
  if ( fHisto->IsA() == TProfile::Class() ) p = static_cast<TProfile*>(fHisto);
  else
  {
    AliError("Mean v2 histo has to be a TProfile");
    return;
  }

  TProfile::Approximate();


  //_____________

  TF1* bck = new TF1("bck",this,&AliAnalysisMuMuJpsiResult::FitFunctionBackgroundPol4Cheb,fitRangeLow,fitRangeHigh,5,"AliAnalysisMuMuJpsiResult","FitFunctionBackgroundPol4Cheb");

  bck->SetParameters(.4,-.2,.02,0.,0.);
  SetFitRejectRange(2.6,3.4);
  bck->SetParameters(BG0_init,BG1_init,BG2_init,BG3_init,BG4_init);
  // bck->SetParLimits(0, -1.,1.0);
 // bck->SetParLimits(0, 1.,8.0);
  // SetFitRejectRange(2.3,4.0);
  SetFitRejectRange(rejRl,rejRh);
  TFitResultPtr fitResultInit = p->Fit(bck,fitOptionBg,"");
  std::cout << "FitResultmBkgInit=" << static_cast<int>(fitResultInit) << std::endl;

    //___________ Further attempts to fit bkg if the first one fails
  if ( static_cast<int>(fitResultInit) ) ProcessmBkgFit(fitResultInit,bck,"FitFunctionBackgroundPol4Cheb",fitOptionBg);
  //___________

  SetFitRejectRange();

  TF1* fitMeanv2 = new TF1("fitMeanv2",this,&AliAnalysisMuMuJpsiResult::FitFunctionMeanPtSNA60NEWPOL2POL3POL4Cheb,fitRangeLow,fitRangeHigh,26,"AliAnalysisMuMuJpsiResult","FitFunctionMeanPtSNA60NEWPOL2POL3POL4Cheb");

  fitMeanv2->SetParNames("a","b","c","ap","bp","cp","dp","kJPsi","mJPsi","sJPsi","p1LJPsi");
  //                      0   1   2    3    4    5   6       7       8       9      10

  fitMeanv2->SetParName(11,"p2LJPsi");
  fitMeanv2->SetParName(12,"p3LJPsi");
  fitMeanv2->SetParName(13,"p1RJPsi");
  fitMeanv2->SetParName(14,"p2RJPsi");
  fitMeanv2->SetParName(15,"p3RJPsi");
  fitMeanv2->SetParName(16,"aLJPsi");
  fitMeanv2->SetParName(17,"aRJPsi");

  fitMeanv2->SetParName(18,"kPsiP");
  fitMeanv2->SetParName(19,"<v2>JPsi");
  fitMeanv2->SetParName(20,"<v2>BGC0");
  fitMeanv2->SetParName(21,"<v2>BGC1");
  fitMeanv2->SetParName(22,"<v2>BGC2");
  fitMeanv2->SetParName(23,"<v2>BGC3");
  fitMeanv2->SetParName(24,"<v2>BGC4");
  fitMeanv2->SetParName(25,"<v2>PsiP");

  fitMeanv2->FixParameter(0,a);
  fitMeanv2->FixParameter(1,b);
  fitMeanv2->FixParameter(2,c);
  fitMeanv2->FixParameter(3,ap);
  fitMeanv2->FixParameter(4,bp);
  fitMeanv2->FixParameter(5,cp);
  fitMeanv2->FixParameter(6,dp);
  fitMeanv2->FixParameter(7,kJPsi);
  fitMeanv2->FixParameter(8,mJPsi);
  fitMeanv2->FixParameter(9,sJPsi);

  fitMeanv2->FixParameter(10,p1Left);
  fitMeanv2->FixParameter(11,p2Left);
  fitMeanv2->FixParameter(12,p3Left);
  fitMeanv2->FixParameter(13,p1Right);
  fitMeanv2->FixParameter(14,p2Right);
  fitMeanv2->FixParameter(15,p3Right);
  fitMeanv2->FixParameter(16,alphaLeft);
  fitMeanv2->FixParameter(17,alphaRight);

  fitMeanv2->FixParameter(18,kPsiP);

  fitMeanv2->SetParameter(19, 0.01);
  fitMeanv2->SetParLimits(19, -1.,1.);

  for ( Int_t i = 0; i < 5; ++i )
  {
    fitMeanv2->SetParameter(i + 20, bck->GetParameter(i));
  }

  fitMeanv2->SetParameter(25, 0.01);
  fitMeanv2->SetParLimits(25, -1.,1.);

  TFitResultPtr fitResult = p->Fit(fitMeanv2,fitOption,"");

  std::cout << "FitResult=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  if ( static_cast<int>(fitResult) )
  {
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 20, bck->GetParameter(i)*0.9);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }
  else if ( fitMeanv2->GetParameter(25) <= fitMeanv2->GetParError(25) || fitMeanv2->GetParError(25) >= 0.75*fitMeanv2->GetParameter(25) )
  {
    fitMeanv2->FixParameter(25, 0.);
    for ( Int_t i = 0; i < 5; ++i )
    {
      fitMeanv2->SetParameter(i + 20, bck->GetParameter(i)*0.6);
    }
    fitResult = p->Fit(fitMeanv2,fitOption,"");
  }

  std::cout << "FitResultt=" << static_cast<int>(fitResult) << std::endl;
  std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  std::cout << "FitChi2PerNDF = " << fitMeanv2->GetChisquare()<<"/"<<fitMeanv2->GetNDF() << std::endl;
    //___________ Further attempts to fit if the first one fails
  // if ( static_cast<int>(fitResult) ||  static_cast<int>(fitResult->CovMatrixStatus())!=3 ) ProcessMv2Fit(fitResult,fitMeanv2,bck,fitOption,12,4); // Int_t iParKPsip, Int_t iLastParBkg);
  //___________

  //___________Set parameters and fit functions to store in the result
  Int_t fitStatus = 0;
  if(!CheckFitStatus(fitResult)) fitStatus =-1;
  else Set("FitStatus",fitStatus,0.);
  Set("FitChi2PerNDF",fitMeanv2->GetChisquare()/fitMeanv2->GetNDF(),0.0);
  Set("FitNDF",fitMeanv2->GetNDF(),0.0);
  Set("mJPsi",fitMeanv2->GetParameter(8),fitMeanv2->GetParError(8));
  Set("sJPsi",fitMeanv2->GetParameter(9),fitMeanv2->GetParError(9));
  // Set("mv2PsiP",fitMeanv2->GetParameter(25),fitMeanv2->GetParError(25));

  for ( Int_t i = 0; i < 5; ++i )
  {
    bck->SetParameter(i, fitMeanv2->GetParameter(i+20));
  }
  AttachFunctionsToHisto(0x0,0x0,bck,fitMeanv2,fitRangeLow,fitRangeHigh);//

  Set("NofJPsi",NofJPsi,ErrStatNofJPsi);
  Set("<v2>JPsi",fitMeanv2->GetParameter(19),fitMeanv2->GetParError(19));
  // Set("<v2>PsiP",fitMeanv2->GetParameter(23),fitMeanv2->GetParError(23));

}
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
    callEnv.Execute(r);// here fit Method ("fit<SOMETHING>") is called and the fit is proceed.
  }
  else
  {
    AliError(Form("Could not get the method %s",fittingMethod.Data()));
    delete r;
    return kFALSE;
  }

  if ( r->IsValid() )
  {
    StdoutToAliDebug(1,r->Print(););
    r->SetBin(Bin());
    r->SetNofTriggers(NofTriggers());
    r->SetNofRuns(NofRuns());

    Bool_t adoptOK = AdoptSubResult(r);
    if ( adoptOK ) {

      std::cout << "Subresult " << r->GetName() << " adopted in " << GetName() <<  std::endl;
      if(IsValidValue(r->Weight()))  SetWeight(Weight()+r->Weight());
      else SetWeight(Weight()+1);
    }
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
  Int_t rebin         =1;
  Double_t fitMinvMin =2.0;
  Double_t fitMinvMax =5.0;
  Double_t paramSPsiP = 1.05 ;/*3.68609/3.096916;*/
  Double_t Weight = 1.;
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
    AliDebug(1,Form("key = %s, value = %s",key.Data(),value.Data()));

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
    else if ( key.CompareTo(kKeyRebin,TString::kIgnoreCase) == 0 )rebin = value.Atoi();

    else if ( key.CompareTo(kKeyHistoType,TString::kIgnoreCase) == 0 ) //FIXME::Is really necesary to save the histoType? I think I dont use it
    {
      histoType = value;

      if ( histoType.CompareTo("minv",TString::kIgnoreCase) == 0 ) Set(kKeyHistoType,0.,0.0); //histoType=0 means minv histo
      else if ( histoType.CompareTo("mpt",TString::kIgnoreCase) == 0 ) Set(kKeyHistoType,1.,0.0); //histoType=1 means mpt histo
      else if ( histoType.CompareTo("mV2",TString::kIgnoreCase) == 0 ) Set(kKeyHistoType,1.,0.0); //histoType=1 means mpt histo
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
    else if ( key.CompareTo(kKeySPsiP,TString::kIgnoreCase) == 0 )paramSPsiP = value.Atof();
    else if ( key.CompareTo(kKeyWeight,TString::kIgnoreCase) == 0 )Weight = value.Atof();
    else if ( key.CompareTo(kKeyMinvRS,TString::kIgnoreCase) == 0 )fMinvRS = value.Data();

    else Set(key.Data(),value.Atof(),0.0);

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
    Set(kKeyWeight,Weight,0.0);
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

//________________________
void AliAnalysisMuMuJpsiResult::ProcessMinvFit(TFitResultPtr& fitResult, TF1* fitTotal, TF1* bckInit, const char* fitOption, Int_t iParKPsip, Int_t iLastParBkg)
{
  // If a Minv fit fails this algorithm changes some initial parameters to get the fit converged

  Int_t bin(0);

  if ( (static_cast<int>(fitResult) && (static_cast<int>(fitResult)!=4000||static_cast<int>(fitResult)!=0)  ) || static_cast<int>(fitResult->CovMatrixStatus())!=3 /*|| static_cast<int>(fitResult->CovMatrixStatus())!=2*/)
  {
    if ( (0.5*fitTotal->GetParameter(iParKPsip) <= fitTotal->GetParError(iParKPsip))) { //kPsi'
      std::cout << "-----------------------------------------------" << std::endl;
      std::cout << "------- Setting Psi'norm= Psi' norm*0.8 -------" << std::endl;
      std::cout << "-----------------------------------------------" << std::endl;
      bin = fHisto->FindBin(3.68);
      // fitTotal->SetParLimits(iParKPsip, 0.,fHisto->GetBinContent(bin)*1.5); // we further restrict the range of psi' norm
      fitTotal->SetParameter(iParKPsip, fHisto->GetBinContent(bin)*0.8);
    }

    if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) {
      std::cout << "-----------------------------------------------" << std::endl;
      std::cout << "-------       Setting p0=MAX/2          -------" << std::endl;
      std::cout << "-----------------------------------------------" << std::endl;
    }

    std::cout << "================================\\" << std::endl;
    std::cout << "======== Refitting again =======\\" << std::endl;
    std::cout << "================================\\" << std::endl;
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

    //Check if there are poles somewhere/
    if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption);
    if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
   }

  if ( (static_cast<int>(fitResult) && (static_cast<int>(fitResult)!=4000||static_cast<int>(fitResult)!=0)) || static_cast<int>(fitResult->CovMatrixStatus())!=3 /*|| static_cast<int>(fitResult->CovMatrixStatus())!=2*/)
  {
    if ( (0.5*fitTotal->GetParameter(iParKPsip) <= fitTotal->GetParError(iParKPsip))  ){ //kPsi'
      std::cout << "------------------------------------------------" << std::endl;
      std::cout << "------- Setting Psi'norm= Psi' norm*0.5) -------" << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      bin = fHisto->FindBin(3.68);
      // fitTotal->SetParLimits(iParKPsip, 0.,fHisto->GetBinContent(bin)*0.9); // we further restrict the range of psi' norm
      fitTotal->SetParameter(iParKPsip, fHisto->GetBinContent(bin)*0.5);
    }
    if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) {
      std::cout << "------------------------------------------------" << std::endl;
      std::cout << "-------         Setting p0=MAX/2)        -------" << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      fitTotal->SetParameter(0, fHisto->GetMaximum()*2.); // kVWG
    }

    std::cout << "================================\\" << std::endl;
    std::cout << "======== Refitting again =======\\" << std::endl;
    std::cout << "================================\\" << std::endl;
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

     //Check if there are poles somewhere/
    if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption) ;
    if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  }

  if ( (static_cast<int>(fitResult) && (static_cast<int>(fitResult)!=4000||static_cast<int>(fitResult)!=0)) || static_cast<int>(fitResult->CovMatrixStatus())!=3 /*|| static_cast<int>(fitResult->CovMatrixStatus())!=2*/) {
    std::cout << "============================================================================================\\" << std::endl;
    std::cout << "======== Refitting bkg again (setting range rejected 2.5-3.7, and fit range 1.7-4.5) =======\\" << std::endl;
    std::cout << "============================================================================================\\" << std::endl;
    SetFitRejectRange(2.5,3.7);
    TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOption,"",1.7,4.5);
    SetFitRejectRange();

    for ( Int_t i = 0; i < iLastParBkg+1 ; ++i ) fitTotal->SetParameter(i, bckInit->GetParameter(i)); //set initial background parameters

    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

     //Check if there are poles somewhere/
    if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption) ;
    if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  }

  if ( (static_cast<int>(fitResult) && (static_cast<int>(fitResult)!=4000||static_cast<int>(fitResult)!=0)) || static_cast<int>(fitResult->CovMatrixStatus())!=3 /*|| static_cast<int>(fitResult->CovMatrixStatus())!=2*/) {
    std::cout << "============================================================================================\\" << std::endl;
    std::cout << "======== Refitting bkg again (setting range rejected 2.5-3.5, and fit range 1.5-5)  ========\\" << std::endl;
    std::cout << "============================================================================================\\" << std::endl;

    SetFitRejectRange(2.5,3.5);
    TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOption,"",1.5,5.);
    SetFitRejectRange();

    for ( Int_t i = 0; i < iLastParBkg+1 ; ++i ) fitTotal->SetParameter(i, bckInit->GetParameter(i)); //set initial background parameters

    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

    if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption) ;
    if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  }

  if ( (static_cast<int>(fitResult) && (static_cast<int>(fitResult)!=4000||static_cast<int>(fitResult)!=0)) || static_cast<int>(fitResult->CovMatrixStatus())!=3.)
  {

    for ( Int_t i = 0; i < iLastParBkg+1 ; ++i ) fitTotal->SetParameter(i, bckInit->GetParameter(i));

    if ( (0.5*fitTotal->GetParameter(iParKPsip) <= fitTotal->GetParError(iParKPsip))  ) { //kPsi'
      std::cout << "------------------------------------------------" << std::endl;
      std::cout << "------- Setting Psi'norm= Psi' norm*0.3) -------" << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      bin = fHisto->FindBin(3.68);
      // fitTotal->SetParLimits(iParKPsip, 0.,fHisto->GetBinContent(bin)*0.7); // we further restrict the range of psi' norm
      fitTotal->SetParameter(iParKPsip, fHisto->GetBinContent(bin)*0.3);
    }

    if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ){
      std::cout << "------------------------------------------------" << std::endl;
      std::cout << "-------         Setting p0=MAX*0.        -------" << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      // fitTotal->SetParLimits(0,bckInit->GetParameter(0)*0.1,bckInit->GetParameter(0)*1.5);
      fitTotal->SetParameter(0, fHisto->GetMaximum()*0.6); // kVWG
    }

    std::cout << "================================\\" << std::endl;
    std::cout << "======== Refitting again =======\\" << std::endl;
    std::cout << "================================\\" << std::endl;    fitResult = fHisto->Fit(fitTotal,fitOption,"");

    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

    if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption) ;
    if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  }

  // if ( !static_cast<int>(fitResult) && static_cast<int>(fitResult->CovMatrixStatus())!=3.){
  //   // change fit option if the only problem is the error calculation
  //   std::cout << "//------- Erro estimation problem, changing fit option" << std::endl;
  //   std::cout << "//======== Refitting again =======\\" << std::endl;

  //   fitResult = fHisto->Fit(fitTotal,Form("%sM",fitOption),"");
  //   std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  //   std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  // }

  if ( (static_cast<int>(fitResult) && (static_cast<int>(fitResult)!=4000||static_cast<int>(fitResult)!=0)) || static_cast<int>(fitResult->CovMatrixStatus())!=3.){
    std::cout << "===========================================================\\" << std::endl;
    std::cout << "======== Cannot fit properly, try something else... =======\\" << std::endl;
    std::cout << "===========================================================\\" << std::endl;
  }
}

//________________________
void AliAnalysisMuMuJpsiResult::ProcessBkgFit(TFitResultPtr &fitResultInit, TF1* bckInit, const char* bkgFuncName, const char* fitOption)
{
  // If a Bkg fit fails this algorithm changes some initial parameters or fitting range to get the fit converged. Can be refined

  TString sbkgFuncName(bkgFuncName);
  Int_t nParsBkg(0);
  TArrayD* initPars(0x0);

  Int_t bin = fHisto->FindBin(0.82); // We change the bin from where we get the initial value for 1st parameter
  if ( !sbkgFuncName.CompareTo("FitFunctionBackgroundVWG")|| !sbkgFuncName.CompareTo("FitFunctionBackgroundVWG2") ) {
    nParsBkg = 4;
    initPars = new TArrayD(nParsBkg);
    initPars->AddAt(fHisto->GetBinContent(bin),0);
    initPars->AddAt(2.,1);
    initPars->AddAt(0.5,2);
    initPars->AddAt(0.3,3);
  } else if ( !sbkgFuncName.CompareTo("FitFunctionBackgroundPol2Exp") ) {
    nParsBkg = 4;
    initPars = new TArrayD(nParsBkg);
    initPars->AddAt(fHisto->GetBinContent(bin),0);
    initPars->AddAt(-fHisto->GetBinContent(bin)/3.,1);
    initPars->AddAt(100.,2);
    initPars->AddAt(0.05,3);
  } else if ( !sbkgFuncName.CompareTo("FitFunctionBackgroundPol4Exp") ) {
    nParsBkg = 6;
    bin = fHisto->FindBin(1.6);
    bckInit->SetParameters(0.,1.,1.,1.,2.,0.5);
    initPars = new TArrayD(nParsBkg);
    initPars->AddAt(-fHisto->GetBinContent(bin),0);
    initPars->AddAt(fHisto->GetBinContent(bin),1);
    initPars->AddAt(-fHisto->GetBinContent(bin)/2.,2);
    initPars->AddAt(fHisto->GetBinContent(bin)/10.,3);
    initPars->AddAt(fHisto->GetBinContent(bin)/100.,4);
    initPars->AddAt(-2.,5);
  //    initPars->AddAt(0.,0);
  //    initPars->AddAt(1.,1);
  //    initPars->AddAt(1.,2);
  //    initPars->AddAt(1.,3);
  //    initPars->AddAt(2.,4);
  //    initPars->AddAt(0.5,5);
  } else if ( !sbkgFuncName.CompareTo("FitFunctionBackgroundPol1Pol2") ) {
    nParsBkg = 4;
    initPars = new TArrayD(nParsBkg);

    initPars->AddAt(0,0);
    initPars->AddAt(fHisto->GetBinContent(bin),1);
    initPars->AddAt(0.,2);
    initPars->AddAt(1.,3);
    // initPars->AddAt(1,4);
  } else if ( !sbkgFuncName.CompareTo("FitFunctionBackgroundPol2Pol3V2") ) {
    bin = fHisto->FindBin(1.5);
    nParsBkg = 6;
    initPars = new TArrayD(nParsBkg);
    initPars->AddAt(0.9,0);
    initPars->AddAt(-9,1);
    initPars->AddAt(bin,2);
    initPars->AddAt(1,3);
    initPars->AddAt(-6,4);
    initPars->AddAt(9,5);
  } else if ( !sbkgFuncName.CompareTo("FitFunctionBackgroundPol2Pol3") ) {
    bin = fHisto->FindBin(1.5);
    nParsBkg = 6;
    initPars = new TArrayD(nParsBkg);
    initPars->AddAt(0.9,0);
    initPars->AddAt(-9,1);
    initPars->AddAt(bin,2);
    initPars->AddAt(1,3);
    initPars->AddAt(-6,4);
    initPars->AddAt(9,5);
  } else {
    AliError("Unrecognized Background function");
    return;
  }

  if(!sbkgFuncName.CompareTo("FitFunctionBackgroundPol2Pol3")){
    for (int i = 0; i < 5; ++i) {
      if(CheckRoots(fitResultInit,bckInit,3,bckInit->GetParameter(3),bckInit->GetParameter(4),bckInit->GetParameter(5),bckInit->GetParameter(6),fitOption)) continue;
      printf("attempt to remove pole n° %d\n", i+1 );
    }
  }

  if(!sbkgFuncName.CompareTo("FitFunctionBackgroundPol1Pol2") ){
    for (int i = 0; i < 5; ++i){
      if(CheckRoots(fitResultInit,bckInit,2,bckInit->GetParameter(3),bckInit->GetParameter(4),bckInit->GetParameter(5),0.,fitOption)) continue;
      printf("attempt to remove pole n° %d\n", i+1 );
    }
  }

  if ( static_cast<int>(fitResultInit) ) {

    std::cout << "//------- Resetting default initial parameters" << std::endl;

    for ( Int_t i = 0 ; i < nParsBkg ; i++ ) {
      Double_t par = initPars->At(i);
      if ( i == 0 ) par =2*par + 1.;
      bckInit->SetParameter(i,par);
    }

    std::cout << "//======== Fitting background again =======\\" << std::endl;
    fitResultInit = fHisto->Fit(bckInit,fitOption);
    std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  }

  if ( static_cast<int>(fitResultInit) )
  {
    // Same initial parameters but change the fitting range

    std::cout << "//------- Changing fitting range to (1.5,5.0)" << std::endl;

    for ( Int_t i = 0 ; i < nParsBkg ; i++ ) bckInit->SetParameter(i,initPars->At(i));
    SetFitRejectRange(2.6,3.5);

    std::cout << "//======== Fitting background again =======\\" << std::endl;
    fitResultInit = fHisto->Fit(bckInit,fitOption,"",1.5,5.);

    std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  }

  if ( static_cast<int>(fitResultInit) ) {
    // Chage first initial parameter and change the fitting range
    std::cout << "Fitting background again" << std::endl;

    for ( Int_t i = 0 ; i < nParsBkg ; i++ ) {
      Double_t par = initPars->At(i);
      if ( i == 0 ) par =2*par + 1.;
      bckInit->SetParameter(i,par);
    }
    SetFitRejectRange(2.7,4.0);

    std::cout << "//======== Fitting background again =======\\" << std::endl;
    fitResultInit = fHisto->Fit(bckInit,fitOption,"",1.5,5.);
    std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  }

  if ( static_cast<int>(fitResultInit) ) {
    // Chage first initial parameter and change the fitting range
    std::cout << "Fitting background again with range à la Momo" << std::endl;

    for ( Int_t i = 0 ; i < nParsBkg ; i++ ) {
      Double_t par = initPars->At(i);
      if ( i == 0 ) par =2*par + 1.;
      bckInit->SetParameter(i,par);
    }
    // SetFitRejectRange(2.7,4.0);

    std::cout << "//======== Fitting background again =======\\" << std::endl;
    fitResultInit = fHisto->Fit(bckInit,fitOption,"",2.2,2.8);
    std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  }
  if ( static_cast<int>(fitResultInit) ) {
    std::cout << std::endl;
    std::cout << "Cannot fit background properly, try something else" << std::endl;
    std::cout << std::endl;
  }

  delete initPars;
}
//________________________
void AliAnalysisMuMuJpsiResult::ProcessMv2Fit(TFitResultPtr& fitResult, TF1* fitTotal, TF1* bckInit, const char* fitOption, Int_t iParKPsip, Int_t iLastParBkg)
{
  //

  Int_t bin(0);
  // TString minuitStatus = gMinuit->fCstatu;
  // if ( static_cast<int>(fitResult) || static_cast<int>(fitResult->CovMatrixStatus())!=3 /*|| static_cast<int>(fitResult->CovMatrixStatus())!=2*/)
  // {
  //   if ( (0.5*fitTotal->GetParameter(iParKPsip) <= fitTotal->GetParError(iParKPsip))) { //kPsi'
  //     std::cout << "-----------------------------------------------" << std::endl;
  //     std::cout << "------- Setting Psi'norm= Psi' norm*0.8 -------" << std::endl;
  //     std::cout << "-----------------------------------------------" << std::endl;
  //     bin = fHisto->FindBin(3.68);
  //     // fitTotal->SetParLimits(iParKPsip, 0.,fHisto->GetBinContent(bin)*1.5); // we further restrict the range of psi' norm
  //     fitTotal->SetParameter(iParKPsip, fHisto->GetBinContent(bin)*0.8);
  //   }

  //   if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) {
  //     std::cout << "-----------------------------------------------" << std::endl;
  //     std::cout << "-------       Setting p0=MAX/2          -------" << std::endl;
  //     std::cout << "-----------------------------------------------" << std::endl;
  //   }

  //   std::cout << "================================\\" << std::endl;
  //   std::cout << "======== Refitting again =======\\" << std::endl;
  //   std::cout << "================================\\" << std::endl;
  //   fitResult = fHisto->Fit(fitTotal,fitOption,"");
  //   std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  //   std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  //   //Check if there are poles somewhere/
  //   if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption);
  //   if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  // }

  //  //

  // if ( static_cast<int>(fitResult) || static_cast<int>(fitResult->CovMatrixStatus())!=3 /*|| static_cast<int>(fitResult->CovMatrixStatus())!=2*/)
  // {
  //   if ( (0.5*fitTotal->GetParameter(iParKPsip) <= fitTotal->GetParError(iParKPsip))  ){ //kPsi'
  //     std::cout << "------------------------------------------------" << std::endl;
  //     std::cout << "------- Setting Psi'norm= Psi' norm*0.5) -------" << std::endl;
  //     std::cout << "------------------------------------------------" << std::endl;
  //     bin = fHisto->FindBin(3.68);
  //     // fitTotal->SetParLimits(iParKPsip, 0.,fHisto->GetBinContent(bin)*0.9); // we further restrict the range of psi' norm
  //     fitTotal->SetParameter(iParKPsip, fHisto->GetBinContent(bin)*0.5);
  //   }
  //   if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) {
  //     std::cout << "------------------------------------------------" << std::endl;
  //     std::cout << "-------         Setting p0=MAX/2)        -------" << std::endl;
  //     std::cout << "------------------------------------------------" << std::endl;
  //     fitTotal->SetParameter(0, fHisto->GetMaximum()*2.); // kVWG
  //   }

  //   std::cout << "================================\\" << std::endl;
  //   std::cout << "======== Refitting again =======\\" << std::endl;
  //   std::cout << "================================\\" << std::endl;
  //   fitResult = fHisto->Fit(fitTotal,fitOption,"");
  //   std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  //   std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

  //    //Check if there are poles somewhere/
  //   if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption) ;
  //   if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  // }

  if ( static_cast<int>(fitResult) ||static_cast<int>(fitResult->CovMatrixStatus())!=3 /*|| static_cast<int>(fitResult->CovMatrixStatus())!=2*/) {
    std::cout << "============================================================================================\\" << std::endl;
    std::cout << "======== Refitting bkg again (setting range rejected 2.5-3.7, and fit range 1.7-4.5) =======\\" << std::endl;
    std::cout << "============================================================================================\\" << std::endl;
    SetFitRejectRange(2.5,3.7);
    TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOption,"",1.7,4.5);
    SetFitRejectRange();

    for ( Int_t i = 0; i < iLastParBkg+1 ; ++i ) fitTotal->SetParameter(i, bckInit->GetParameter(i)); //set initial background parameters

    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

     //Check if there are poles somewhere/
    // if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption) ;
    // if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  }

  if ( static_cast<int>(fitResult) ||static_cast<int>(fitResult->CovMatrixStatus())!=3 /*|| static_cast<int>(fitResult->CovMatrixStatus())!=2*/) {
    std::cout << "============================================================================================\\" << std::endl;
    std::cout << "======== Refitting bkg again (setting range rejected 2.5-3.5, and fit range 1.5-5)  ========\\" << std::endl;
    std::cout << "============================================================================================\\" << std::endl;

    SetFitRejectRange(2.5,3.5);
    TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOption,"",1.5,5.);
    SetFitRejectRange();

    for ( Int_t i = 0; i < iLastParBkg+1 ; ++i ) fitTotal->SetParameter(i, bckInit->GetParameter(i)); //set initial background parameters

    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

    if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption) ;
    if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  }

  if ( static_cast<int>(fitResult) ||static_cast<int>(fitResult->CovMatrixStatus())!=3.)
  {

    for ( Int_t i = 0; i < iLastParBkg+1 ; ++i ) fitTotal->SetParameter(i, bckInit->GetParameter(i));

    if ( (0.5*fitTotal->GetParameter(iParKPsip) <= fitTotal->GetParError(iParKPsip))  ) { //kPsi'
      std::cout << "------------------------------------------------" << std::endl;
      std::cout << "------- Setting Psi'norm= Psi' norm*0.3) -------" << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      bin = fHisto->FindBin(3.68);
      // fitTotal->SetParLimits(iParKPsip, 0.,fHisto->GetBinContent(bin)*0.7); // we further restrict the range of psi' norm
      fitTotal->SetParameter(iParKPsip, fHisto->GetBinContent(bin)*0.3);
    }

    if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ){
      std::cout << "------------------------------------------------" << std::endl;
      std::cout << "-------         Setting p0=MAX*0.        -------" << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      // fitTotal->SetParLimits(0,bckInit->GetParameter(0)*0.1,bckInit->GetParameter(0)*1.5);
      fitTotal->SetParameter(0, fHisto->GetMaximum()*0.6); // kVWG
    }

    std::cout << "================================\\" << std::endl;
    std::cout << "======== Refitting again =======\\" << std::endl;
    std::cout << "================================\\" << std::endl;    fitResult = fHisto->Fit(fitTotal,fitOption,"");

    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

    if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption) ;
    if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  }

  // if ( !static_cast<int>(fitResult) && static_cast<int>(fitResult->CovMatrixStatus())!=3.){
  //   // change fit option if the only problem is the error calculation
  //   std::cout << "//------- Erro estimation problem, changing fit option" << std::endl;
  //   std::cout << "//======== Refitting again =======\\" << std::endl;

  //   fitResult = fHisto->Fit(fitTotal,Form("%sM",fitOption),"");
  //   std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  //   std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  // }

  if ( static_cast<int>(fitResult) ||static_cast<int>(fitResult->CovMatrixStatus())!=3.){
    std::cout << "===========================================================\\" << std::endl;
    std::cout << "======== Cannot fit properly, try something else... =======\\" << std::endl;
    std::cout << "===========================================================\\" << std::endl;
  }
}
//________________________
void AliAnalysisMuMuJpsiResult::ProcessmBkgFit(TFitResultPtr &fitResultInit, TF1* bckInit, const char* bkgFuncName, const char* fitOption)
{
  // If a Bkg fit fails this algorithm changes some initial parameters or fitting range to get the fit converged. Can be refined

  TString sbkgFuncName(bkgFuncName);
  Int_t nParsBkg(0);
  TArrayD* initPars(0x0);

  Int_t bin = fHisto->FindBin(0.82); // We change the bin from where we get the initial value for 1st parameter
  if ( !sbkgFuncName.CompareTo("FitFunctionBackgroundPol2")) {
    nParsBkg = 3;
    initPars = new TArrayD(nParsBkg);
    initPars->AddAt(fHisto->GetBinContent(bin),0);
    initPars->AddAt(2.,1);
    initPars->AddAt(0.5,2);
    initPars->AddAt(0.3,3);
  } else if ( !sbkgFuncName.CompareTo("FitFunctionBackgroundPol2Exp") ) {
    nParsBkg = 4;
    initPars = new TArrayD(nParsBkg);
    initPars->AddAt(fHisto->GetBinContent(bin),0);
    initPars->AddAt(-fHisto->GetBinContent(bin)/3.,1);
    initPars->AddAt(100.,2);
    initPars->AddAt(0.05,3);
  } else if ( !sbkgFuncName.CompareTo("FitFunctionBackgroundPol3") ) {
    nParsBkg = 4;
    bin = fHisto->FindBin(1.6);
    bckInit->SetParameters(0.,1.,1.,1.,2.,0.5);
    initPars = new TArrayD(nParsBkg);
    initPars->AddAt(-fHisto->GetBinContent(bin),0);
    initPars->AddAt(fHisto->GetBinContent(bin),1);
    initPars->AddAt(-fHisto->GetBinContent(bin)/2.,2);
    initPars->AddAt(fHisto->GetBinContent(bin)/10.,3);
    initPars->AddAt(fHisto->GetBinContent(bin)/100.,4);
    initPars->AddAt(-2.,5);
  } else if ( !sbkgFuncName.CompareTo("FitFunctionBackgroundPol4") ) {
    nParsBkg = 5;
    initPars = new TArrayD(nParsBkg);

    initPars->AddAt(0,0);
    initPars->AddAt(fHisto->GetBinContent(bin),1);
    initPars->AddAt(0.,2);
    initPars->AddAt(1.,3);
  } else {
    AliError("Unrecognized Background function");
    return;
  }
  if ( static_cast<int>(fitResultInit) ) {

    std::cout << "//------- Resetting default initial parameters" << std::endl;

    for ( Int_t i = 0 ; i < nParsBkg ; i++ ) {
      Double_t par = initPars->At(i);
      if ( i == 0 ) par =2*par + 1.;
      bckInit->SetParameter(i,par);
    }

    std::cout << "//======== Fitting background again =======\\" << std::endl;
    fitResultInit = fHisto->Fit(bckInit,fitOption);
    std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  }

  if ( static_cast<int>(fitResultInit) ){
    // Same initial parameters but change the fitting range

    std::cout << "//------- Changing fitting range to (1.5,5.0)" << std::endl;

    for ( Int_t i = 0 ; i < nParsBkg ; i++ ) bckInit->SetParameter(i,initPars->At(i));
    SetFitRejectRange(2.6,3.5);

    std::cout << "//======== Fitting background again =======\\" << std::endl;
    fitResultInit = fHisto->Fit(bckInit,fitOption,"",1.5,5.);

    std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  }

  if ( static_cast<int>(fitResultInit) ) {
    // Chage first initial parameter and change the fitting range
    std::cout << "Fitting background again" << std::endl;

    for ( Int_t i = 0 ; i < nParsBkg ; i++ ) {
      Double_t par = initPars->At(i);
      if ( i == 0 ) par =2*par + 1.;
      bckInit->SetParameter(i,par);
    }
    SetFitRejectRange(2.7,4.0);

    std::cout << "//======== Fitting background again =======\\" << std::endl;
    fitResultInit = fHisto->Fit(bckInit,fitOption,"",1.5,5.);
    std::cout << "FitResultBkgInit=" << static_cast<int>(fitResultInit) << std::endl;
  }

  if ( static_cast<int>(fitResultInit) ) {
    std::cout << std::endl;
    std::cout << "Cannot fit background properly, try something else" << std::endl;
    std::cout << std::endl;
  }

  delete initPars;
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

    AliDebug(1,Form("i=%d particle %s n %e sigma %.2f",i,particleNames[i],n,sigma));

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

  Double_t npart     = GetValue(Form("Nof%s",particle));
  Double_t ninput    = GetValue(Form("NofInput%s",particle));

  Set(Form("AccEff%s",particle),
      npart/ninput,TMath::Max(1./ninput,TMath::Sqrt(npart/ninput*TMath::Abs(1.-npart/ninput)/ninput)));

  TIter next(SubResults());
  AliAnalysisMuMuJpsiResult* r;

  while ( ( r = static_cast<AliAnalysisMuMuJpsiResult*>(next())) )
  {
    r->Set(Form("NofInput%s",particle),n,TMath::Sqrt(n));

    npart = r->GetValue(Form("Nof%s",particle));

    Double_t value = (ninput>0.) ? npart/ninput : 0 ;
    Double_t error = (ninput>0.) ? TMath::Max(1./ninput, TMath::Sqrt(value*TMath::Abs(1.-value)/ninput)) : 1.;

    r->Set(Form("AccEff%s",particle),value,error);
    printf("AccEff%s : %f +/- %f \n",particle,value,error);

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
Bool_t AliAnalysisMuMuJpsiResult::CheckRoots(TFitResultPtr &fitResult, TF1* fitFunction, Int_t deg, Double_t a, Double_t b, Double_t c, Double_t d,const char* fitOption  )
{
  /// Check if there is no roots founs in the x axis inside the fit range. If there are some, refit again.
  /// Meant to work for polynomial ratio background function.

  if(!fitFunction){
    AliError("Cannot check roots without functons ...");
    return kFALSE;
  }

  Double_t fitRangeLow  = GetValue(kFitRangeLow);
  Double_t fitRangeHigh = GetValue(kFitRangeHigh);
  Double_t epsilon      = 2.;

  if(deg == 2) {
    Double_t delta = b*b -4*a*c;
    if(delta <= 0.) return kTRUE;
    Double_t x1 = (-b - TMath::Sqrt(delta ))/(2.*a);
    Double_t x2 = (-b + TMath::Sqrt(delta ))/(2.*a);
    if ((fitRangeLow < x1 && x1 < fitRangeHigh )){
       printf(" !!!!!! Found roots at %f !!!!!! \n",x1);
       fitFunction->SetParameters(2,a+a*epsilon);
       fitResult = fHisto->Fit(fitFunction,fitOption,"");
       printf("fit Result       = %d\n",static_cast<int>(fitResult) );
       printf("Covariant matrix = %d\n",static_cast<int>(fitResult->CovMatrixStatus()) );
       return kFALSE;
    }
    if ((fitRangeLow < x2 && x2 < fitRangeHigh )){
       printf(" !!!!!! Found roots at %f !!!!!! \n",x2);
       fitFunction->SetParameters(2,a+a*epsilon);
       fitResult = fHisto->Fit(fitFunction,fitOption,"");
       printf("fit Result       = %d\n",static_cast<int>(fitResult) );
       printf("Covariant matrix = %d\n",static_cast<int>(fitResult->CovMatrixStatus()) );
       return kFALSE;
    }
  }

  if(deg == 3 && a != 0.) {
    // Cardan methods, see https://fr.wikipedia.org/wiki/M%C3%A9thode_de_Cardan

    Double_t p = -b*b/(3.*a*a) + c/a;
    // printf("p =%f\n",p );
    Double_t q =b/(27.*a) * (2.*b*b/a/a -9*c/a ) +d/a;
    // printf("q =%f\n",q );
    Double_t delta = -(4*p*p*p + 27*q*q);
    // printf("delta =%f\n",delta );

    Double_t z1 = 0.5*(-q - TMath::Sqrt(-0.037*delta));
    // printf("z1 =%f\n",z1 );
    Double_t z2 = 0.5*(-q + TMath::Sqrt(-0.037*delta));
    // printf("z2 =%f\n",z2 );

    Double_t r1 = z1/TMath::Abs(z1)*TMath::Power(TMath::Abs(z1),0.33) ;
    // printf("r1 =%f\n",r1 );
    Double_t r2 = z2/TMath::Abs(z2)*TMath::Power(TMath::Abs(z2),0.33) ;
    // printf("r2 =%f\n",r2 );

    Double_t R1 =0.;
    Double_t R2 =0.;
    Double_t R3 =0.;

    if( delta < 0) R1 = r1 + r2 -0.333*b/a;
    else if(delta == 0){
      R1 = 3*q/p -0.333*b/a;
      R2 = -1.5*q*p -0.333*b/a;
      R3 = -1.5*q*p -0.333*b/a;
    }
    else if(delta > 0.){
      R1 = 2*TMath::Sqrt(-0.333*p)*TMath::Cos(0.333*TMath::ACos(-0.5*q*TMath::Sqrt(-27/p/p/p) + 0.*0.666*3.14)) -0.333*b/a;
      // R2 = 2*TMath::Sqrt(-0.333*p)*TMath::Cos(0.333*TMath::ACos(-0.5*q*TMath::Sqrt(-27/p/p/p) + 1.*0.666*3.14)) -0.333*b/a;
      // R3 = 2*TMath::Sqrt(-0.333*p)*TMath::Cos(0.333*TMath::ACos(-0.5*q*TMath::Sqrt(-27/p/p/p) + 2.*0.666*3.14)) -0.333*b/a;
    }

    if (fitRangeLow < R1  && R1 < fitRangeHigh){
      printf(" !!!!!! Found roots at %f !!!!!! \n", R1);
      fitFunction->SetParameter(3,a+a*epsilon);
      printf("Parameter 3 = %f \n",fitFunction->GetParameter(3) );
      fitResult = fHisto->Fit(fitFunction,fitOption,"");
      printf("fit Result       = %d\n",static_cast<int>(fitResult) );
      printf("Covariant matrix = %d\n",static_cast<int>(fitResult->CovMatrixStatus()) );
      // if ( static_cast<int>(fitResult) ||static_cast<int>(fitResult->CovMatrixStatus())!=3.) return kFALSE;
      return kFALSE;
    }
    if (fitRangeLow < R2  && R2 < fitRangeHigh){
      printf(" !!!!!! Found roots at %f !!!!!! \n", R2);
      fitFunction->SetParameter(3,a+a*epsilon);
      printf("Parameter 3 = %f \n",fitFunction->GetParameter(3) );
      fitResult = fHisto->Fit(fitFunction,fitOption,"");
      printf("fit Result = %d\n",static_cast<int>(fitResult) );
      printf("Covariant matrix = %d\n",static_cast<int>(fitResult->CovMatrixStatus()) );
      // if ( static_cast<int>(fitResult) ||static_cast<int>(fitResult->CovMatrixStatus())!=3.) return kFALSE;
      return kFALSE;
    }
    if (fitRangeLow < R3  && R3 < fitRangeHigh){
      printf(" !!!!!! Found roots at %f !!!!!! \n", R3);
      fitFunction->SetParameter(3,a+a*epsilon);
      printf("Parameter 3 = %f \n",fitFunction->GetParameter(3) );
      fitFunction->SetParameter(4,b+b*epsilon);
      printf("Parameter 4 = %f \n",fitFunction->GetParameter(4) );
      fitResult = fHisto->Fit(fitFunction,fitOption,"");
      printf("fit Result = %d\n",static_cast<int>(fitResult) );
      printf("Covariant matrix = %d\n",static_cast<int>(fitResult->CovMatrixStatus()) );
      // if ( static_cast<int>(fitResult) || static_cast<int>(fitResult->CovMatrixStatus())!=3.) return kFALSE;
      return kFALSE;
    }
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisMuMuJpsiResult::CheckFitStatus(TFitResultPtr &fitResult)
{
  Bool_t isok = kTRUE;
  if ( static_cast<int>(fitResult) && static_cast<int>(fitResult)!=4000 ){
    isok =kFALSE;
    AliDebug(1,Form("Fit rejected because of fitresult : %d",static_cast<int>(fitResult)));
  }

  if(static_cast<int>(fitResult->CovMatrixStatus())!=3 ){
    isok =kFALSE;
    AliDebug(1,Form("Fit rejected because of covariant matrix : %d",fitResult->CovMatrixStatus()));
  }
  TString minuitStatus = gMinuit->fCstatu;
  if(!minuitStatus.Contains("SUCCESSFUL") && !minuitStatus.Contains("OK") && !minuitStatus.Contains("PROBLEMS")){
    isok =kFALSE;
    AliDebug(1,"Minuit status is not ok !");
  }
  return isok;
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
  /// return kTRUE if npar-th parameter of fit function has a big error,
  /// and in that case fix the parameter's value to fixValueIfWrong

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
