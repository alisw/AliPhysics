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
/// author: Laurent Aphecetche (Subatech)
///

#include "AliAnalysisMuMuJpsiResult.h"

ClassImp(AliAnalysisMuMuJpsiResult)

#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TH2.h"
#include "THashList.h"
#include "TLine.h"
#include "TList.h"
#include "TMap.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TParameter.h"
#include "AliAnalysisMuMuBinning.h"
#include "AliLog.h"
#include <map>

namespace {
  
  const std::map<std::string,Double_t>& MassMap()
  {
    /// a simple map of masses...
    static std::map<std::string,Double_t> massMap;
    // not super elegant, but that way we do not depend on TDatabasePDG and thus
    // can decide on our particle naming
    if (massMap.empty())
    {
      massMap["Jpsi"] = 3.096916e+00;
      massMap["PsiPrime"] = 3.68609e+00;
      massMap["Upsilon"] = 9.46030e+00;
      massMap["UpsilonPrime"] = 1.00233e+01;
    }
    return massMap;
  }
  
  
  Double_t funcCB(Double_t* xx, Double_t* par)
  {
    /// Crystal ball
    
    Double_t norm = par[0];
    Double_t alpha = par[1];
    Double_t n = par[2];
    Double_t mean = par[3];
    Double_t sigma = par[4];
    
    Double_t x = xx[0];
    
    Double_t a = TMath::Power(n/TMath::Abs(alpha),n)*TMath::Exp(-0.5*alpha*alpha);
    Double_t b = n/TMath::Abs(alpha) - TMath::Abs(alpha);
    
    Double_t y = ( TMath::Abs(sigma) > 1E-12 ? (x-mean)/sigma : 0 );
    
    if ( y > alpha*-1.0 ) 
    {
      return norm*TMath::Exp(-0.5*y*y);
    }
    else 
    {
      return norm*a*TMath::Power(b-y,-n);
    }
  }
  
  Double_t funcJpsiGCBE(Double_t* xx, Double_t* par)
  {
    /// crystal ball + expo + gaussian
    Double_t x = xx[0];
    
    Double_t g = par[0]*TMath::Gaus(x,par[1],par[2]);
    
    Double_t jpsi = funcCB(xx,par+3);
    
    Double_t expo = par[8]*TMath::Exp(par[9]*x);
    
    return g+expo+jpsi;
  }
  
  Double_t funcCB2(Double_t* xx, Double_t* par)
  {
    /// CB2 = extended crystal ball
    
    Double_t norm = par[0];
    Double_t alpha = par[1];
    Double_t n = par[2];
    Double_t mean = par[3];
    Double_t sigma = par[4];
    Double_t alphaprime = par[5];
    Double_t nprime = par[6];
    
    Double_t x = xx[0];
    
    Double_t a = TMath::Power(n/TMath::Abs(alpha),n)*TMath::Exp(-0.5*alpha*alpha);
    Double_t b = n/TMath::Abs(alpha) - TMath::Abs(alpha);
    Double_t c = TMath::Power(nprime/TMath::Abs(alphaprime),nprime)*TMath::Exp(-0.5*alphaprime*alphaprime);
    Double_t d = nprime/TMath::Abs(alphaprime) - TMath::Abs(alphaprime);
    
    Double_t y = ( TMath::Abs(sigma) > 1E-12 ? (x-mean)/sigma : 0 );
    
    if ( y > alphaprime )
    {
      return norm*c*TMath::Power(d+y,-nprime);
    }
    else if ( y > alpha*-1.0 ) 
    {
      return norm*TMath::Exp(-0.5*y*y);
    }
    else 
    {
      return norm*a*TMath::Power(b-y,-n);
    }
  }
  
  
  Double_t funcJpsiNA48(Double_t*x, Double_t* par)
  {
    /// Fit function from e.q. 4.8 of Ruben's PhD.
    Double_t c1 = par[0];
    Double_t c2 = par[1];
    Double_t d1 = par[2];
    Double_t d2 = par[3];
    Double_t g1 = par[4];
    Double_t g2 = par[5];
    Double_t m0 = par[6];
    Double_t sigma1 = par[7];
    Double_t sigma2 = par[8];
    Double_t b1 = par[9];
    Double_t b2 = par[10];
    Double_t norm = par[11];
    
    Double_t m = x[0];
    
    Double_t rv(0);
    
    if ( m <= c1*m0 )
    {
      Double_t e = d1-g1*TMath::Sqrt(c1*m0-m);
      rv = TMath::Power(b1*(c1*m0-m),e);
      rv += sigma1;
    }
    else if( m >= c1*m0 && m <= m0 )
    {
      rv = sigma1;
    }
    else if ( m >= m0 && m < c2*m0 )
    {
      rv = sigma2;
    }
    else if( m >= c2*m0 )
    {
      Double_t e = d2-g2*TMath::Sqrt(m-c2*m0);
      rv = TMath::Power(b2*(m-c2*m0),e);
      rv += sigma2;
    }
    
    return norm*TMath::Exp(-(m-m0)*(m-m0)/(2.0*rv*rv));
  }
  
  //------------------------------------------------------------------------------
  Double_t BackgroundVWG(Double_t *x, Double_t *par)
  {
    // gaussian variable width
    Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
    return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
    
  }
  
  //------------------------------------------------------------------------------
  Double_t CrystalBallExtended(Double_t *x,Double_t *par)
  {
    //par[0] = Normalization
    //par[1] = mean
    //par[2] = sigma
    //par[3] = alpha
    //par[4] = n
    //par[5] = alpha'
    //par[6] = n'
    
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
  
  
  //---------------------------------------------------------------------------
//  Double_t fitFunctionVWG(Double_t *x, Double_t *par)
//  {
//    if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
//    return BackgroundVWG(x, par);
//  }
  
  //---------------------------------------------------------------------------
  Double_t fitFunctionCB2VWG(Double_t *x, Double_t *par)
  {
    return BackgroundVWG(x, par) + CrystalBallExtended(x, &par[4]);
  }
  
  //---------------------------------------------------------------------------
  Double_t func2CB2VWG(Double_t *x, Double_t *par)
  {
    /// 2 extended crystal balls + variable width gaussian
    /// width of the second CB related to the first (free) one.
    
    Double_t par2[7] = {
      par[11],
      3.68609+(par[5]-3.096916)/3.096916*3.68609,
      par[6]/3.096916*3.68609,
      par[7],
      par[8],
      par[9],
      par[10]
    };
    return BackgroundVWG(x, par) + CrystalBallExtended(x, &par[4]) + CrystalBallExtended(x, par2);
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult::AliAnalysisMuMuJpsiResult(TRootIOCtor* /*io*/) :
AliAnalysisMuMuResult("",""),
fNofRuns(),
fNofTriggers(-1),
fMinv(0x0),
fBin(),
fRebin(0),
fTriggerClass(),
fEventSelection(),
fPairSelection(),
fCentralitySelection()
{
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult::AliAnalysisMuMuJpsiResult(const TH1& hminv) :
AliAnalysisMuMuResult("",""),
fNofRuns(1),
fNofTriggers(-1),
fMinv(0x0),
fBin(),
fRebin(0),
fTriggerClass(),
fEventSelection(),
fPairSelection(),
fCentralitySelection()
{
  SetMinv(hminv);
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult::AliAnalysisMuMuJpsiResult(const TH1& hminv,
                                             const char* fitType,
                                             Int_t nrebin)
:
AliAnalysisMuMuResult(Form("%s:%d",fitType,nrebin),""),
fNofRuns(1),
fNofTriggers(-1),
fMinv(0x0),
fBin(),
fRebin(nrebin),
fTriggerClass(),
fEventSelection(),
fPairSelection(),
fCentralitySelection()
{
  SetMinv(hminv);
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult::AliAnalysisMuMuJpsiResult(const TH1& hminv,
                                             const char* triggerName,
                                             const char* eventSelection,
                                             const char* pairSelection,
                                             const char* centSelection,
                                             const AliAnalysisMuMuBinning::Range& bin)
:
AliAnalysisMuMuResult(Form("%s-%s-%s-%s",triggerName,eventSelection,pairSelection,centSelection),""),
fNofRuns(1),
fNofTriggers(-1),
fMinv(0x0),
fBin(bin),
fRebin(1),
fTriggerClass(triggerName),
fEventSelection(eventSelection),
fPairSelection(pairSelection),
fCentralitySelection(centSelection)
{
  SetMinv(hminv);
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult::AliAnalysisMuMuJpsiResult(const AliAnalysisMuMuJpsiResult& rhs)
:
AliAnalysisMuMuResult(rhs),
fNofRuns(rhs.NofRuns()),
fNofTriggers(rhs.NofTriggers()),
fMinv(0x0),
fBin(rhs.Bin()),
fRebin(rhs.fRebin),
fTriggerClass(rhs.fTriggerClass),
fEventSelection(rhs.fEventSelection),
fPairSelection(rhs.fPairSelection),
fCentralitySelection(rhs.fCentralitySelection)
{
  /// copy ctor
  /// Note that the mother is lost
  /// fKeys remains 0x0 so it will be recomputed if need be

  if ( rhs.fMinv )
  {
    fMinv = static_cast<TH1*>(rhs.fMinv->Clone());
  }  
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult& AliAnalysisMuMuJpsiResult::operator=(const AliAnalysisMuMuJpsiResult& rhs)
{
  /// Assignment operator
  
  if (this!=&rhs)
  {
    static_cast<AliAnalysisMuMuResult&>(*this) = static_cast<const AliAnalysisMuMuResult&>(rhs);
    delete fMinv;

    if ( rhs.fMinv )
    {
      fMinv = static_cast<TH1*>(rhs.fMinv->Clone());
    }
    
    fNofRuns = rhs.NofRuns();
    fNofTriggers = rhs.NofTriggers();
    fBin = rhs.Bin();
    fRebin = rhs.fRebin;    
  }
  
  return *this;
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult::~AliAnalysisMuMuJpsiResult()
{
  // dtor
  delete fMinv;
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

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult* AliAnalysisMuMuJpsiResult::FitJpsiGCBE(TH1& h)
{
  /// Fit Jpsi spectra with crystal ball + gaussian + exponential
  
  std::cout << "Fit with jpsi alone (gaus + CB + expo)" << std::endl;
  
  Int_t nrebin = fMinv->GetXaxis()->GetNbins() / h.GetXaxis()->GetNbins();
  
  AliAnalysisMuMuJpsiResult* r = new AliAnalysisMuMuJpsiResult(h,"JPSIGCBE",nrebin);
  
  TH1* hfit = r->Minv();
  
  const Double_t xmin(1.0);
  const Double_t xmax(8.0);
  
  TF1* fitTotal = new TF1("fitTotal",funcJpsiGCBE,xmin,xmax,10);
  fitTotal->SetParNames("cste","x0","sigma0","N","alpha","n","mean","sigma","expocste","exposlope");
  
  fitTotal->SetParLimits(3,0,h.GetMaximum()*2); // N
  
  const Double_t cbalpha(0.98);
  const Double_t cbn(5.2);
  
  fitTotal->FixParameter(4,cbalpha);
  fitTotal->FixParameter(5,cbn);
  
  fitTotal->SetParLimits(6,2.8,3.2); // mean
  fitTotal->SetParLimits(7,0.02,0.3); // sigma
  
  TF1* fg = new TF1("fg","gaus",xmin,xmax);
  
  hfit->Fit(fg,"","",0.75,3.0);
  
  fitTotal->SetParameter(0,fg->GetParameter(0));
  fitTotal->SetParameter(1,fg->GetParameter(1));
  fitTotal->SetParameter(2,fg->GetParameter(2));
  
  TF1* fexpo = new TF1("expo","expo",xmin,xmax);
  
  hfit->Fit(fexpo,"","",3.5,5);
  
  fitTotal->SetParameter(8,fexpo->GetParameter(0));
  fitTotal->SetParameter(9,fexpo->GetParameter(1));
  
  fitTotal->SetParameter(3,h.GetMaximum()),
  fitTotal->SetParameter(4,cbalpha);
  fitTotal->SetParameter(5,cbn);
  fitTotal->SetParameter(6,3.15);
  fitTotal->SetParameter(7,0.1);
  
  const char* fitOption = "QSI+";
  
  TFitResultPtr fitResult = hfit->Fit(fitTotal,fitOption,"",2,5);
  
  r->Set("MeanJpsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  r->Set("SigmaJpsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
  
  double m = r->GetValue("MeanJpsi");
  double s = r->GetValue("SigmaJpsi");
  double n = 3.0;
  
  TF1* fcb = new TF1("fcb",funcCB,xmin,xmax,5);
  fcb->SetParameters(fitTotal->GetParameter(3),
                     fitTotal->GetParameter(4),
                     fitTotal->GetParameter(5),
                     fitTotal->GetParameter(6),
                     fitTotal->GetParameter(7));
  
  fcb->SetLineColor(6);
  fcb->SetNpx(100);
  TLine* l1 = new TLine(m-n*s,0,m-n*s,fitTotal->GetParameter(3));
  TLine* l2 = new TLine(m+n*s,0,m+n*s,fitTotal->GetParameter(3));
  l1->SetLineColor(6);
  l2->SetLineColor(6);
  h.GetListOfFunctions()->Add(fcb);
  h.GetListOfFunctions()->Add(l1);
  h.GetListOfFunctions()->Add(l2);
  
  
  Double_t cbParameters[5];
  Double_t covarianceMatrix[5][5];
  
  cbParameters[0] = fitTotal->GetParameter(3);
  cbParameters[1] = fitTotal->GetParameter(4);
  cbParameters[2] = fitTotal->GetParameter(5);
  cbParameters[3] = fitTotal->GetParameter(6);
  cbParameters[4] = fitTotal->GetParameter(7);
  
  for ( int iy = 0; iy < 5; ++iy )
  {
    for ( int ix = 0; ix < 5; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+3,iy+3);
    }
  }
  
  double njpsi = fcb->Integral(m-n*s,m+n*s)/h.GetBinWidth(1);
  
  double nerr = fcb->IntegralError(m-n*s,m+n*s,&cbParameters[0],&covarianceMatrix[0][0])/h.GetBinWidth(1);
  
  r->Set("NofJpsi",njpsi,nerr);
  
  return r;
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult* AliAnalysisMuMuJpsiResult::FitJpsi(TH1& h)
{
  /// Fit Jpsi spectra using extended crystall ball (CB2) with free tails
  
  StdoutToAliDebug(1,std::cout << "Fit with jpsi alone" << std::endl;);

  Int_t nrebin = fMinv->GetXaxis()->GetNbins() / h.GetXaxis()->GetNbins();
  
  AliAnalysisMuMuJpsiResult* r = new AliAnalysisMuMuJpsiResult(h,"JPSI",nrebin);
  
  TH1* hfit = r->Minv();

  const Double_t xmin(1.5);
  const Double_t xmax(8.0);

  TF1* fitTotal = new TF1("fitTotal",funcCB2,xmin,xmax,7);
  fitTotal->SetParNames("N","alphaLow","nLow","mean","sigma","alphaUp","nUp");
  fitTotal->SetParameters(h.GetMaximum(),1,5,3.1,0.07,1.5,3);
  fitTotal->SetParLimits(0,0,h.GetMaximum()*2); // N
  fitTotal->SetParLimits(1,0,10); // alpha
  fitTotal->SetParLimits(2,0.1,10); // n
  fitTotal->SetParLimits(3,3,3.15); // mean
  fitTotal->SetParLimits(4,0.01,1); // sigma
  fitTotal->SetParLimits(5,0,10); // alpha
  fitTotal->SetParLimits(6,0.1,10); // n
  
  hfit->Fit(fitTotal,"QSER+","",2,5);
  
  
  r->Set("MeanJpsi",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  r->Set("SigmaJpsi",fitTotal->GetParameter(4),fitTotal->GetParError(4));

  double m = r->GetValue("MeanJpsi");
  double s = r->GetValue("SigmaJpsi");
  double n = 10.0;

  r->Set("NofJpsi",fitTotal->Integral(m-n*s,m+n*s)/h.GetBinWidth(1),fitTotal->IntegralError(m-n*s,m+n*s)/h.GetBinWidth(1));

  return r;
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult* AliAnalysisMuMuJpsiResult::FitJpsiCB2VWG(const TH1& h)
{
  /// Fit Jpsi spectra using extended crystal ball (CB2) + variable width gaussian (VWG)
  
  StdoutToAliDebug(1,std::cout << "Fit with jpsi VWG" << std::endl;);
  
  Int_t nrebin = fMinv->GetXaxis()->GetNbins() / h.GetXaxis()->GetNbins();
  
  AliAnalysisMuMuJpsiResult* r = new AliAnalysisMuMuJpsiResult(h,"JPSICB2VWG",nrebin);
  
  
  TH1* hfit = r->Minv();
  
  const Double_t xmin(2.0);
  const Double_t xmax(5.0);
  
//  // gaussian variable width
//  Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
//  return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
//  Double_t CrystalBallExtended(Double_t *x,Double_t *par)
//  //par[0] = Normalization
//  //par[1] = mean
//  //par[2] = sigma
//  //par[3] = alpha
//  //par[4] = n
//  //par[5] = alpha'
//  //par[6] = n'

  TF1* fitTotal = new TF1("fitTotal",fitFunctionCB2VWG,xmin,xmax,11);
  fitTotal->SetParNames("kVWG","mVWG","sVWG1","sVWG2","norm","mean","sigma","alpha","n","alpha'","n'");
  
  fitTotal->SetParameter(0, 10000.); // kVWG
  fitTotal->SetParameter(1, 1.9); // mVWG
  
  fitTotal->SetParameter(2, 0.5); // sVWG1
  fitTotal->SetParLimits(2, 0., 100.);
  
  fitTotal->SetParameter(3, 0.3); // sVWG2
  fitTotal->SetParLimits(3, 0., 100.);
  
  fitTotal->SetParameter(4, h.GetMaximum()); // norm
  
  fitTotal->SetParameter(5, 3.1); // mean
  fitTotal->SetParLimits(5, 3.0, 3.2);
  
  fitTotal->SetParameter(6, 0.08); // sigma
  fitTotal->SetParLimits(6, 0.04, 0.20);
  
  fitTotal->SetParameter(7,1.0); // alpha
  fitTotal->SetParameter(8,5); // n
  fitTotal->SetParameter(9,2.0); // alpha'
  fitTotal->SetParameter(10,4); // n'
  
//  fitTotal->FixParameter(7, 0.93);
//  fitTotal->FixParameter(8, 5.59);
//  fitTotal->FixParameter(9, 2.32);
//  fitTotal->FixParameter(10, 3.39);
//  fitTotal->SetParameter(11, 10.);
  
  const char* fitOption = "QSIER"; //+";
  
  TFitResultPtr fitResult = hfit->Fit(fitTotal,fitOption,"");
  
  r->Set("MeanJpsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  r->Set("SigmaJpsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  
  double m = r->GetValue("MeanJpsi");
  double s = r->GetValue("SigmaJpsi");
  double n = 3.0;
  
  TF1* fcb = new TF1("fcb",CrystalBallExtended,xmin,xmax,7);
  fcb->SetParameters(fitTotal->GetParameter(4),
                     fitTotal->GetParameter(5),
                     fitTotal->GetParameter(6),
                     fitTotal->GetParameter(7),
                     fitTotal->GetParameter(8),
                     fitTotal->GetParameter(9),
                     fitTotal->GetParameter(10));
  
  
  fcb->SetLineColor(1);
  fcb->SetNpx(1000);
  TLine* l1 = new TLine(m-n*s,0,m-n*s,fitTotal->GetParameter(4));
  TLine* l2 = new TLine(m+n*s,0,m+n*s,fitTotal->GetParameter(4));
  l1->SetLineColor(6);
  l2->SetLineColor(6);
  hfit->GetListOfFunctions()->Add(fcb);
  hfit->GetListOfFunctions()->Add(l1);
  hfit->GetListOfFunctions()->Add(l2);
  
  Double_t cbParameters[7];
  Double_t covarianceMatrix[7][7];
  
  for ( int ix = 0; ix < 7; ++ix )
  {
    cbParameters[ix] = fitTotal->GetParameter(ix+4);
  }
  
  for ( int iy = 0; iy < 5; ++iy )
  {
    for ( int ix = 0; ix < 5; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+4,iy+4);
    }
  }
  
  double njpsi = fcb->Integral(m-n*s,m+n*s)/h.GetBinWidth(1);
  double nerr = fcb->IntegralError(m-n*s,m+n*s,&cbParameters[0],&covarianceMatrix[0][0])/h.GetBinWidth(1);
  
  r->Set("NofJpsi",njpsi,nerr);
  
  return r;
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult* AliAnalysisMuMuJpsiResult::FitJpsi2CB2VWG(const TH1& h,
                                                             Double_t alphaLow,
                                                             Double_t nLow,
                                                             Double_t alphaUp,
                                                             Double_t nUp)
{
  /// Fit using extended crystal ball + variable width gaussian
  
  StdoutToAliDebug(1,std::cout << Form("Fit with jpsi + psiprime VWG alphaLow=%5.2f nLow=%5.2f alphaUp=%5.2f nUp=%5.2f",
                                       alphaLow,nLow,alphaUp,nUp) << std::endl;);
  
  Int_t nrebin = fMinv->GetXaxis()->GetNbins() / h.GetXaxis()->GetNbins();
  
  TString resultName("JPSI2CB2VWG");
  if ( alphaLow > 0 )
  {
    resultName += TString::Format("alphaLow=%5.2f",alphaLow);
  }
  if ( nLow > 0 )
  {
    resultName += TString::Format("nLow=%5.2f",nLow);
  }
  if ( alphaUp > 0 )
  {
    resultName += TString::Format("alphaUp=%5.2f",alphaUp);
  }
  if ( nUp > 0 )
  {
    resultName += TString::Format("nUp=%5.2f",nUp);
  }
  resultName.ReplaceAll(" ","");
  
  AliAnalysisMuMuJpsiResult* r = new AliAnalysisMuMuJpsiResult(h,resultName.Data(),nrebin);
  
  TH1* hfit = r->Minv();
  
  const Double_t xmin(2.2);
  const Double_t xmax(5.0);
  
  TF1* fitTotal = new TF1("fitTotal",func2CB2VWG,xmin,xmax,12);
  fitTotal->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");
  fitTotal->SetParName(11, "kPsi'");

  fitTotal->SetParameter(0, 10000.);
  fitTotal->SetParameter(1, 1.9);
  fitTotal->SetParameter(2, 0.5);
  fitTotal->SetParLimits(2, 0., 100.);
  fitTotal->SetParameter(3, 0.3);
  fitTotal->SetParLimits(3, 0., 100.);
  fitTotal->SetParameter(4, 100.);
  fitTotal->SetParameter(5, 3.1);
  fitTotal->SetParLimits(5, 3.08, 3.2);
  fitTotal->SetParameter(6, 0.08);
  fitTotal->SetParLimits(6, 0.05, 0.15);
  
//  r = FitJpsi2CB2VWG(*hminv,0.93,5.59,2.32,3.39);

  if ( alphaLow > 0 )
  {
    fitTotal->FixParameter(7, alphaLow);
  }
  else
  {
    fitTotal->SetParameter(7,0.9);
    fitTotal->SetParLimits(7,0.1,10.0);
  }
  
  if ( nLow > 0 )
  {
    fitTotal->FixParameter(8, nLow);
  }
  else
  {
    fitTotal->SetParameter(8,5.0);
    fitTotal->SetParLimits(8,0.0,10.0);
  }
  
  if ( alphaUp > 0 )
  {
    fitTotal->FixParameter(9, alphaUp);
  }
  else
  {
    fitTotal->SetParameter(9, 2.0);
    fitTotal->SetParLimits(9,0.1,10.0);
  }
  
  if ( nUp > 0 )
  {
    fitTotal->FixParameter(10, nUp);    
  }
  else
  {
    fitTotal->SetParameter(10,3.0);
    fitTotal->SetParLimits(10,0.0,10.0);
  }
  
  fitTotal->SetParameter(11, 10.);

  const char* fitOption = "QSER"; //+";
  
  TFitResultPtr fitResult = hfit->Fit(fitTotal,fitOption,"");
  
  r->Set("MeanJpsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  r->Set("SigmaJpsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  
  double m = r->GetValue("MeanJpsi");
  double s = r->GetValue("SigmaJpsi");
  double n = 3.0;
    
  TF1* fcb = new TF1("fcb",CrystalBallExtended,xmin,xmax,7);
  fcb->SetParameters(fitTotal->GetParameter(4),
                     fitTotal->GetParameter(5),
                     fitTotal->GetParameter(6),
                     fitTotal->GetParameter(7),
                     fitTotal->GetParameter(8),
                     fitTotal->GetParameter(9),
                     fitTotal->GetParameter(10));
                     
  
  fcb->SetLineColor(1);
  fcb->SetNpx(1000);
  TLine* l1 = new TLine(m-n*s,0,m-n*s,fitTotal->GetParameter(4));
  TLine* l2 = new TLine(m+n*s,0,m+n*s,fitTotal->GetParameter(4));
  l1->SetLineColor(6);
  l2->SetLineColor(6);
  hfit->GetListOfFunctions()->Add(fcb);
  hfit->GetListOfFunctions()->Add(l1);
  hfit->GetListOfFunctions()->Add(l2);
  
  Double_t cbParameters[7];
  Double_t covarianceMatrix[7][7];
  
  for ( int ix = 0; ix < 7; ++ix )
  {
    cbParameters[ix] = fitTotal->GetParameter(ix+4);
  }
  
  for ( int iy = 0; iy < 5; ++iy )
  {
    for ( int ix = 0; ix < 5; ++ix )
    {
      covarianceMatrix[ix][iy] = (fitResult->GetCovarianceMatrix())(ix+4,iy+4);
    }
  }
  
  double njpsi = fcb->Integral(m-n*s,m+n*s)/h.GetBinWidth(1);
  double nerr = fcb->IntegralError(m-n*s,m+n*s,&cbParameters[0],&covarianceMatrix[0][0])/h.GetBinWidth(1);
  
  r->Set("NofJpsi",njpsi,nerr);
  
  return r;
}

//_____________________________________________________________________________
AliAnalysisMuMuJpsiResult* AliAnalysisMuMuJpsiResult::FitJpsiNA48(const TH1& h)
{
  /// fit using functional form from Ruben Shahoyan's thesis (2001) (eq. 4.8.)
  
  StdoutToAliDebug(1,std::cout << "Fit with jpsi NA50 Ruben eq. 4.8" << std::endl;);
  
  Int_t nrebin = fMinv->GetXaxis()->GetNbins() / h.GetXaxis()->GetNbins();
  
  AliAnalysisMuMuJpsiResult* r = new AliAnalysisMuMuJpsiResult(h,"JPSINA",nrebin);
  
  TH1* hfit = r->Minv();
  
  const Double_t xmin(2.0);
  const Double_t xmax(5.0);
  
  TF1* fitTotal = new TF1("fitTotal",funcJpsiNA48,xmin,xmax,12);
  
  fitTotal->SetParName( 0, "c1");
  fitTotal->FixParameter(0,0.97);
  
  fitTotal->SetParName( 1, "c2");
  fitTotal->FixParameter(1,1.05);
  
  fitTotal->SetParName( 2, "d1");
  fitTotal->SetParameter(2,0.0);
  fitTotal->SetParLimits(2,0,1);
  
  fitTotal->SetParName( 3, "d2");
  fitTotal->SetParameter(3,0.0);
  fitTotal->SetParLimits(3,0,1);
  
  fitTotal->SetParName( 4, "g1");
  fitTotal->SetParameter(4,0.0);
  fitTotal->SetParLimits(4,0.0,1);
  
  fitTotal->SetParName( 5, "g2");
  fitTotal->SetParameter(5,0.0);
  fitTotal->SetParLimits(5,0.0,1);
  
  fitTotal->SetParName( 6, "m0");
  fitTotal->SetParameter(6,3.1);
  fitTotal->SetParLimits(6,2.8,3.4);

  fitTotal->SetParName( 7, "sigma1");
  fitTotal->SetParameter(7,0.05);
  fitTotal->SetParLimits(7,0.01,0.2);
  
  fitTotal->SetParName( 8, "sigma2");
  fitTotal->SetParameter(8,0.05);
  fitTotal->SetParLimits(8,0.01,0.2);

  fitTotal->SetParName( 9, "b1");
  fitTotal->SetParameter(9,1.0);
  fitTotal->SetParLimits(9,0,1);
  
  fitTotal->SetParName(10, "b2");
  fitTotal->SetParameter(10,1.0);
  fitTotal->SetParLimits(10,0,1);
  
  fitTotal->SetParName(11, "norm");
  fitTotal->SetParameter(11,h.GetMaximum());
  
  const char* fitOption = "QSIER"; //+";
  
  TFitResultPtr fitResult = hfit->Fit(fitTotal,fitOption,"");
  
  r->Set("MeanJpsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  r->Set("SigmaJpsi",
         fitTotal->GetParameter(7)+fitTotal->GetParameter(8),
         0.0);

  double m = r->GetValue("MeanJpsi");
  double s = r->GetValue("SigmaJpsi");
  double n = 3.0;
  
  TLine* l1 = new TLine(m-n*s,0,m-n*s,fitTotal->GetParameter(11));
  TLine* l2 = new TLine(m+n*s,0,m+n*s,fitTotal->GetParameter(11));
  l1->SetLineColor(6);
  l2->SetLineColor(6);
  hfit->GetListOfFunctions()->Add(l1);
  hfit->GetListOfFunctions()->Add(l2);
  
  double njpsi = fitTotal->Integral(m-n*s,m+n*s)/h.GetBinWidth(1);
  double nerr = fitTotal->IntegralError(m-n*s,m+n*s)/h.GetBinWidth(1);
  
  r->Set("NofJpsi",njpsi,nerr);
  
  return r;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuJpsiResult::AddFit(const char* fitType, Int_t npar, Double_t* par)
{
  // Add a fit to this result
  
  TString msg(Form("fitType=%s npar=%d par[]=",fitType,npar));
  
  for ( Int_t i = 0; i < npar; ++i )
  {
    msg += TString::Format("%e,",par[i]);
  }
  
  msg += TString::Format(" minv=%p %d",fMinv,fMinv?TMath::Nint(fMinv->GetEntries()):0);
  
  if ( !fMinv ) return kFALSE;
  
  TString sFitType(fitType);
  sFitType.ToUpper();
  Int_t nrebin(1);
  
  if (sFitType.CountChar(':'))
  {
    Int_t index = sFitType.Index(":");
    nrebin = TString(sFitType(index+1,sFitType.Length()-index-1)).Atoi();
    sFitType = sFitType(0,index);
  }
  
  msg += TString::Format(" nrebin=%d",nrebin);
  
  AliDebug(1,msg.Data());
  

  if ( fMinv->GetEntries()<100 && !sFitType.Contains("COUNT")) return kFALSE;
  
  TH1* hminv = static_cast<TH1*>(fMinv->Clone());
  
  hminv->Rebin(nrebin);
  hminv->SetDirectory(0);

  AliAnalysisMuMuJpsiResult* r(0x0);
  
  if ( sFitType=="PSI1")
  {
    r = FitJpsi(*hminv);
  }
  else if ( sFitType == "PSILOW")
  {
    r = FitJpsi2CB2VWG(*hminv,-1,-1,-1,-1); // free tails
  }
  else if ( sFitType == "PSILOWMCTAILS" )
  {
    if ( npar!= 4 )
    {
      AliError("Cannot use PSILOWMCTAILS without being given the MC tails !");
      delete hminv;
      return kFALSE;
    }
    r = FitJpsi2CB2VWG(*hminv,par[0],par[1],par[2],par[3]);
    if (r)
    {
      r->SetAlias("MCTAILS");
    }
  }
  else if ( sFitType.BeginsWith("PSILOWALPHA") )
  {
    Float_t lpar[] = { 0.0, 0.0, 0.0, 0.0 };
    
    AliDebug(1,Form("sFitType=%s",sFitType.Data()));
    
    sscanf(sFitType.Data(),"PSILOWALPHALOW%fNLOW%fALPHAUP%fNUP%f",
           &lpar[0],&lpar[1],&lpar[2],&lpar[3]);
    
    AliDebug(1,Form("PSILOW ALPHALOW=%f NLOW=%f ALPHAUP=%f NUP=%f",lpar[0],lpar[1],lpar[2],lpar[3]));
    
    if ( lpar[0] == 0.0 && lpar[1] == 0.0 && lpar[0] == 0.0 && lpar[1] == 0.0 )
    {
      AliError("Cannot work with zero tails !");
    }
    else
    {
      r = FitJpsi2CB2VWG(*hminv,lpar[0],lpar[1],lpar[2],lpar[3]);      
    }
  }
  else if ( sFitType == "COUNTJPSI" )
  {
    r = new AliAnalysisMuMuJpsiResult(*hminv,"COUNTJPSI",1);
    Double_t n = CountParticle(*hminv,"Jpsi");
    r->Set("NofJpsi",n,TMath::Sqrt(n));
  }
  
  
  if ( r )
  {
    StdoutToAliDebug(1,r->Print(););
    r->SetBin(Bin());
    r->SetNofTriggers(NofTriggers());
    r->SetNofRuns(NofRuns());

    AdoptSubResult(r);
  }
  
  delete hminv;
  
  return (r!=0x0);
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
void AliAnalysisMuMuJpsiResult::SetNofInputParticles(const TH1& hminv)
{
  /// Set the number of input particle from the invariant mass spectra
  
  static const char* particleNames[] = { "Jpsi" , "PsiPrime", "Upsilon","UpsilonPrime" };

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
void AliAnalysisMuMuJpsiResult::SetMinv(const TH1& hminv)
{
    /// Set the inv. mass spectrum to be fitted.
  static UInt_t n(0);
  
  delete fMinv;
  fMinv = static_cast<TH1*>(hminv.Clone(Form("Minv%u",n++)));
  fMinv->SetDirectory(0);
}
