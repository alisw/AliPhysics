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


// $Id$

#include "AliAnalysisMuMuResult.h"

ClassImp(AliAnalysisMuMuResult)

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
#include "AliHistogramCollection.h"
#include "AliLog.h"
#include <map>

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

Double_t funcJpsiPsiPrimeCustom(Double_t* xx, Double_t* par)
{
  // custom fit for jpsi + psi prime
  
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
  
  Double_t cb(0);
  
  if ( y > alphaprime )
  {
    cb = norm*c*TMath::Power(d+y,-nprime);
  }
  else if ( y > alpha*-1.0 ) 
  {
    cb = norm*TMath::Exp(-0.5*y*y);
  }
  else 
  {
    cb = norm*a*TMath::Power(b-y,-n);
  }
  
  if ( x < mean )
  {
    return cb + par[7] + par[8]*x; // gaus + pol1
  }
  else
  {
    Double_t yprime = (x-par[10])/par[11];
    return cb + par[9]*TMath::Exp(-0.5*yprime*yprime)+par[12]*TMath::Exp(-par[13]*x);
    // gaus (j/psi) + gaus (psi') + expo
  }
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

Double_t funcJpsiPsiPrime(Double_t* xx, Double_t* par)
{
  /// CB + CB2
  
  Double_t jpsi = funcCB(xx,par);
  Double_t psiprime = funcCB2(xx,par+5);
  
  int n = 10;
  Double_t x = xx[0];
    
  Double_t e1 = par[n]*TMath::Exp(par[n+1]*x);
  Double_t e2 = par[n+2]*TMath::Exp(par[n+3]*x);    
  
  Double_t e = e1;
  
  if ( x > par[3] ) e=e2;
  
  return jpsi+psiprime+e;
}

Double_t funcJpsiCBE(Double_t* xx, Double_t* par)
{
  // CB + expo
  
  Double_t jpsi = funcCB(xx,par);
  
  Double_t x = xx[0];
  
  Double_t e1 = par[5]*TMath::Exp(par[6]*x);
  
  return jpsi+e1;
}


Double_t funcJpsiPCBE(Double_t* xx, Double_t* par)
{
  // CB + expo + pol2
  
  Double_t x = xx[0];

  Double_t pol2 = par[0] + par[1]*x + par[2]*x*x;

  Double_t jpsi = funcCB(xx,par+3);
  
  Double_t expo = par[8]*TMath::Exp(par[9]*x);
  
  return pol2+jpsi+expo;
}

Double_t funcJpsiECBE(Double_t* xx, Double_t* par)
{
  // CB + expo
  
  Double_t jpsi = funcCB(xx,par+2);
  
  Double_t x = xx[0];
  
  Double_t e1 = par[0]*TMath::Exp(par[1]*x);
  
  Double_t e2 = par[7]*TMath::Exp(par[8]*x);
  
  return e1+e2+jpsi;
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

const char* NormalizeName(const char* name, const char* suffix)
{
  /// Remove - and / from the name, and adds _suffix
  
  TString str(Form("%s_%s",name,suffix));
  
  str.ReplaceAll("-","_");
  str.ReplaceAll("/","%");
  
  return str.Data();
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

//------------------------------------------------------------------------------
Double_t Gaus(Double_t *x, Double_t *par)
{
  // gaussian
  return par[0]/TMath::Sqrt(2.*TMath::Pi())/par[2]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*par[2]*par[2]));
  
}

//------------------------------------------------------------------------------
Double_t Exp(Double_t *x, Double_t *par)
{
  // exponential
  return par[0]*TMath::Exp(par[1]*x[0]);
  
}

//------------------------------------------------------------------------------
Double_t Pow(Double_t *x, Double_t *par)
{
  // power law
  return par[0]*TMath::Power(x[0],par[1]);
  
}

//------------------------------------------------------------------------------
Double_t fitFunctionVWG(Double_t *x, Double_t *par)
{
  if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  return BackgroundVWG(x, par);
}

//------------------------------------------------------------------------------
Double_t fitFunctionCB2VWG(Double_t *x, Double_t *par)
{
  return BackgroundVWG(x, par) + CrystalBallExtended(x, &par[4]);
}

//------------------------------------------------------------------------------
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

//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________

//_____________________________________________________________________________
AliAnalysisMuMuResult::AliAnalysisMuMuResult(TRootIOCtor* /*io*/) :
TNamed("",""),
fNofRuns(),
fNofTriggers(-1),
fMinv(0x0),
fBin(),
fSubResults(0x0),
fMap(0x0),
fMother(0x0),
fKeys(0x0),
fWeight(0.0),
fRebin(0),
fTriggerClass(),
fEventSelection(),
fPairSelection(),
fCentralitySelection(),
fAlias()
{
}

//_____________________________________________________________________________
AliAnalysisMuMuResult::AliAnalysisMuMuResult(const TH1& hminv)
:
  TNamed("",""),
  fNofRuns(1),
  fNofTriggers(-1),
  fMinv(0x0),
  fBin(),
  fSubResults(0x0),
  fMap(0x0),
  fMother(0x0),
fKeys(0x0),
fWeight(0.0),
fRebin(0),
fTriggerClass(),
fEventSelection(),
fPairSelection(),
fCentralitySelection(),
fAlias()
{
  SetMinv(hminv);
}

//_____________________________________________________________________________
AliAnalysisMuMuResult::AliAnalysisMuMuResult(const TH1& hminv,
                                             const char* fitType,
                                             Int_t nrebin)
:
TNamed(Form("%s:%d",fitType,nrebin),""),
fNofRuns(1),
fNofTriggers(-1),
fMinv(0x0),
fBin(),
fSubResults(0x0),
fMap(0x0),
fMother(0x0),
fKeys(0x0),
fWeight(0.0),
fRebin(nrebin),
fTriggerClass(),
fEventSelection(),
fPairSelection(),
fCentralitySelection(),
fAlias()
{
  SetMinv(hminv);
}

//_____________________________________________________________________________
AliAnalysisMuMuResult::AliAnalysisMuMuResult(const TH1& hminv,
                                             const char* triggerName,
                                             const char* eventSelection,
                                             const char* pairSelection,
                                             const char* centSelection,
                                             const AliAnalysisMuMuBinning::Range& bin)
:
TNamed(Form("%s-%s-%s-%s",triggerName,eventSelection,pairSelection,centSelection),""),
fNofRuns(1),
fNofTriggers(-1),
fMinv(0x0),
fBin(bin),
fSubResults(0x0),
fMap(0x0),
fMother(0x0),
fKeys(0x0),
fWeight(0.0),
fRebin(1),
fTriggerClass(triggerName),
fEventSelection(eventSelection),
fPairSelection(pairSelection),
fCentralitySelection(centSelection),
fAlias()
{
  SetMinv(hminv);
}

//_____________________________________________________________________________
AliAnalysisMuMuResult::AliAnalysisMuMuResult(const AliAnalysisMuMuResult& rhs)
:
TNamed(rhs),
fNofRuns(rhs.NofRuns()),
fNofTriggers(rhs.NofTriggers()),
fMinv(0x0),
fBin(rhs.Bin()),
fSubResults(0x0),
fMap(0x0),
fMother(0x0),
fKeys(0x0),
fWeight(rhs.fWeight),
fRebin(rhs.fRebin),
fAlias()
{
  /// copy ctor
  /// Note that the mother is lost
  /// fKeys remains 0x0 so it will be recomputed if need be

  if ( rhs.fMinv )
  {
    fMinv = static_cast<TH1*>(rhs.fMinv->Clone());
  }
  
  if (rhs.fSubResults)
  {
    fSubResults = static_cast<TObjArray*>(rhs.fSubResults->Clone());
  }
  
  if ( rhs.fMap )
  {
    fMap = static_cast<TMap*>(rhs.fMap->Clone());
  }
  
  if ( rhs.fAlias.Length() > 0 )
  {
    fAlias = rhs.fAlias;
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuResult& AliAnalysisMuMuResult::operator=(const AliAnalysisMuMuResult& rhs)
{
  /// Assignment operator
  
  if (this!=&rhs)
  {
    delete fMinv;
    delete fMap;
    delete fSubResults;
    
    fMinv = 0x0;
    fMap = 0x0;
    fSubResults = 0x0;
    fKeys = 0x0;
    
    if ( rhs.fMinv )
    {
      fMinv = static_cast<TH1*>(rhs.fMinv->Clone());
    }
    
    if (rhs.fSubResults)
    {
      fSubResults = static_cast<TObjArray*>(rhs.fSubResults->Clone());
    }
    
    if ( rhs.fMap )
    {
      fMap = static_cast<TMap*>(rhs.fMap->Clone());
    }

    static_cast<TNamed&>(*this)=rhs;
    
    fNofRuns = rhs.NofRuns();
    fNofTriggers = rhs.NofTriggers();
    fBin = rhs.Bin();
    fWeight = rhs.fWeight;
    fRebin = rhs.fRebin;
    
    fAlias="";
    
    if ( rhs.fAlias.Length() > 0 )
    {
      fAlias = rhs.fAlias;
    }
  }
  
  return *this;
}

//_____________________________________________________________________________
AliAnalysisMuMuResult::~AliAnalysisMuMuResult()
{
  // dtor
  delete fMap;
  delete fMinv;
  delete fSubResults;
  delete fKeys;
}

//_____________________________________________________________________________
const AliAnalysisMuMuBinning::Range& AliAnalysisMuMuResult::Bin() const
{
  /// Get the bin of this result
  if ( !Mother() ) return fBin;
  else return Mother()->Bin();
}

//_____________________________________________________________________________
TObject* AliAnalysisMuMuResult::Clone(const char* /*newname*/) const
{
  /// Clone this result
  return new AliAnalysisMuMuResult(*this);
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuResult::Correct(const AliAnalysisMuMuResult& other, const char* particle, const char* subResultName)
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
Double_t AliAnalysisMuMuResult::CountParticle(const TH1& hminv, const char* particle, Double_t sigma)
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
  
  TAxis* x = hminv.GetXaxis();

  Int_t b1 = x->FindBin(mass-sigma);
  Int_t b2 = x->FindBin(mass+sigma);
  
  AliDebugClass(1,Form("hminv getentries %e integral %e",hminv.GetEntries(),hminv.Integral(b1,b2)));
  
  return hminv.Integral(b1,b2);
}

/*
//_____________________________________________________________________________
void AliAnalysisMuMuResult::FitJpsiPsiPrimeCustom(TH1& h)
{
  std::cout << "Fit with jpsi + psiprime (custom)" << std::endl;
  
  const Double_t xmin(1.5);
  const Double_t xmax(8.0);
  
  fitTotal = new TF1("fitTotal",funcJpsiPsiPrimeCustom,xmin,xmax,14);
  fitTotal->SetLineColor(4);
  
  fitTotal->SetParName(0,"cstecb");
  fitTotal->SetParName(1,"alpha");
  fitTotal->SetParName(2,"n");
  fitTotal->SetParName(3,"meanjpsi");
  fitTotal->SetParName(4,"sigmajpsi");
  fitTotal->SetParName(5,"alphaprime");
  fitTotal->SetParName(6,"nprime");
  fitTotal->SetParName(7,"cstepol1");
  fitTotal->SetParName(8,"slopepol1");
  fitTotal->SetParName(9,"cstegaus");
  fitTotal->SetParName(10,"meanpsiprime");
  fitTotal->SetParName(11,"sigmapsiprime");
  fitTotal->SetParName(12,"csteexpo");
  fitTotal->SetParName(13,"slopeexpo");
  
  fitTotal->SetParameter( 0,1);
    
  const char* fitOption = "SQBR+";
  const Double_t alphaMC = 0.936;
  const Double_t nMC = 4.44;
  const Double_t alphaprimeMC = 1.60;
  const Double_t nprimeMC = 3.23;
  
  TF1* fcb = new TF1("cb",funcCB2,2.9,3.3,7);
  fcb->SetParameters(1,1.0,4.0,3.1,0.1,1.5,3);

  fcb->SetParLimits(3,3,4); 
  fcb->SetParLimits(4,0,1); 

  fcb->FixParameter(1,alphaMC);
  fcb->FixParameter(2,nMC);
  fcb->FixParameter(5,alphaprimeMC);
  fcb->FixParameter(6,nprimeMC);
  
  TFitResultPtr rcb = h.Fit(fcb,fitOption,"",2.9,3.3);

  if (!rcb.Get())
  {
    return;
  }
  
  fitTotal->SetParameter(0,rcb->Parameter(0));
  fitTotal->SetParameter(1,rcb->Parameter(1)); fitTotal->SetParLimits(1,0,10); // alpha
  fitTotal->SetParameter(2,rcb->Parameter(2)); fitTotal->SetParLimits(2,1,10); // n
  fitTotal->SetParameter(3,rcb->Parameter(3)); fitTotal->SetParLimits(3,3.0,3.5); // mean
  fitTotal->SetParameter(4,rcb->Parameter(4)); fitTotal->SetParLimits(4,0,1); // sigma
  fitTotal->SetParameter(5,rcb->Parameter(5)); fitTotal->SetParLimits(1,0,10); // alphaprime
  fitTotal->SetParameter(6,rcb->Parameter(6)); fitTotal->SetParLimits(2,1,10); // nprime

  fitTotal->FixParameter(1,alphaMC);
  fitTotal->FixParameter(2,nMC);
  fitTotal->FixParameter(5,alphaprimeMC);
  fitTotal->FixParameter(6,nprimeMC);
  
  TF1* fge = new TF1("fge","gaus(0)+expo(3)",3.5,4.4);
  fge->SetParameters(1,3.6,0.25,1,1);
  TFitResultPtr rpsiprime = h.Fit(fge,fitOption,"",3.5,4.4);
  
  if (static_cast<int>(rpsiprime))
  {
    AliInfo("Will fix psiprime parameters");
    fitTotal->FixParameter(9,0);
    fitTotal->FixParameter(10,3.7);
    fitTotal->FixParameter(11,0.1);
  }
  else
  {
    fitTotal->SetParameter(10,rpsiprime->Parameter(1)); fitTotal->SetParLimits(10,3.5,3.8); // mean'
    fitTotal->SetParameter(11,rpsiprime->Parameter(2)); fitTotal->SetParLimits(11,0.05,0.7); // sigma'
  }
  
  TFitResultPtr rpol1 = h.Fit("pol1",fitOption,"",1.5,2.5);
  fitTotal->SetParameter( 7,rpol1->Parameter(0));
  fitTotal->SetParameter( 8,rpol1->Parameter(1));
  
  TFitResultPtr rexpo = h.Fit("expo",fitOption,"",4.5,7.0);
  fitTotal->SetParameter(12,rexpo->Parameter(0));
  fitTotal->SetParameter(13,rexpo->Parameter(1));
  
  
  TFitResultPtr r = h.Fit(fitTotal,fitOption,"",1.5,7);
  
  TF1* signal = new TF1("signal","gaus",2,6);  
  signal->SetParameters(fitTotal->GetParameter(0),
                        fitTotal->GetParameter(3),
                        fitTotal->GetParameter(4));

  TF1* signalPrime = new TF1("signalPrime","gaus",2,6);  
  signalPrime->SetParameters(fitTotal->GetParameter(9),
                             fitTotal->GetParameter(10),
                             fitTotal->GetParameter(11));
  
  Double_t gausParameters[3];
  Double_t covarianceMatrix[3][3];
  Double_t gausParametersPrime[3];
  Double_t covarianceMatrixPrime[3][3];
  
  covarianceMatrix[0][0] = (r->GetCovarianceMatrix())(0,0);
  covarianceMatrix[1][0] = (r->GetCovarianceMatrix())(3,0);
  covarianceMatrix[2][0] = (r->GetCovarianceMatrix())(4,0);
  covarianceMatrix[0][1] = (r->GetCovarianceMatrix())(0,3);
  covarianceMatrix[0][2] = (r->GetCovarianceMatrix())(0,4);  
  
  for ( int iy = 1; iy < 3; ++iy )
  {
    for ( int ix = 1; ix < 3; ++ix )
    {
      covarianceMatrix[ix][iy] = (r->GetCovarianceMatrix())(ix+2,iy+2);
    }
  }
  
  gausParameters[0] = fitTotal->GetParameter(0);
  gausParameters[1] = fitTotal->GetParameter(3);
  gausParameters[2] = fitTotal->GetParameter(4);

  gausParametersPrime[0] = fitTotal->GetParameter(9);
  gausParametersPrime[1] = fitTotal->GetParameter(10);
  gausParametersPrime[2] = fitTotal->GetParameter(11);
  
  covarianceMatrixPrime[0][0] = (r->GetCovarianceMatrix())(9,9);
  covarianceMatrixPrime[1][0] = (r->GetCovarianceMatrix())(10,9);
  covarianceMatrixPrime[2][0] = (r->GetCovarianceMatrix())(11,9);
  covarianceMatrixPrime[0][1] = (r->GetCovarianceMatrix())(9,10);
  covarianceMatrixPrime[0][2] = (r->GetCovarianceMatrix())(9,11);  
  
  for ( int iy = 1; iy < 3; ++iy )
  {
    for ( int ix = 1; ix < 3; ++ix )
    {
      covarianceMatrixPrime[ix][iy] = (r->GetCovarianceMatrix())(ix+2,iy+2);
    }
  }
  
  double n = signal->Integral(2,6)/h.GetBinWidth(10);
  double nerr = signal->IntegralError(2,6,&gausParameters[0],&covarianceMatrix[0][0])/h.GetBinWidth(10);
  Set("NofJpsi",n,nerr);      
  Set("MeanJpsi",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("SigmaJpsi",fitTotal->GetParameter(4),fitTotal->GetParError(4));

  double nprime = signalPrime->Integral(2,6)/h.GetBinWidth(10);
  double nerrprime = signalPrime->IntegralError(2,6,&gausParametersPrime[0],&covarianceMatrixPrime[0][0])/h.GetBinWidth(10);
  Set("NofPsiPrime",nprime,nerrprime);
  Set("MeanPsiPrime",fitTotal->GetParameter(10),fitTotal->GetParError(10));
  Set("SigmaPsiPrime",fitTotal->GetParameter(11),fitTotal->GetParError(11));
}
*/

/*
//_____________________________________________________________________________
AliAnalysisMuMuResult::SubResult AliAnalysisMuMuResult::FitJpsiPsiPrimeCB(TH1& h)
{
  std::cout << "Fit with jpsi + psiprime (CB) " << std::endl;

  const Double_t xmin(1.5);
  const Double_t xmax(8.0);
  
  fitTotal = new TF1("fitTotal",funcJpsiPsiPrime,xmin,xmax,14);

//  Double_t N = par[0];
//  Double_t alpha = par[1];
//  Double_t n = par[2];
//  Double_t mean = par[3];
//  Double_t sigma = par[4];
  
  fitTotal->SetParameter( 0,1); // N
  fitTotal->FixParameter( 1,0.936); // alpha
  fitTotal->FixParameter( 2,4.44); // n
  fitTotal->SetParameter( 3,3.1); fitTotal->SetParLimits(3,3.0,3.2); // mean
  fitTotal->SetParameter( 4,0.07); fitTotal->SetParLimits(4,0.02,1); // sigma

  fitTotal->SetParameter( 5,0.01); // N'
  fitTotal->FixParameter( 6,0.936); // alpha'
  fitTotal->FixParameter( 7,4.44); // n'
  fitTotal->SetParameter( 8,3.7); fitTotal->SetParLimits(8,3.5,3.8); // mean'
  fitTotal->SetParameter( 9,0.1); fitTotal->SetParLimits(9,0.02,1.0); // sigma'
  
  fitTotal->SetParameter(10,h.GetMaximum());
  fitTotal->SetParameter(11,-1);

  fitTotal->SetParameter(12,h.GetMaximum()/100);
  fitTotal->SetParameter(13,-1);

  TFitResultPtr r = h.Fit(fitTotal,"SQBI","",1.5,6);
  
//  for ( int ix = 0; ix < fitTotal->GetNpar(); ++ix )
//  {
//    for ( int iy = 0; iy < fitTotal->GetNpar(); ++iy )
//    {      
//      std::cout << Form("COV(%d,%d)=%e ",ix,iy,r->GetCovarianceMatrix()(ix,iy));        
//    }
//    std::cout << std::endl;
//  }
  
  
  TF1* signal = new TF1("signal","gaus",2,8);
  
  signal->SetParameters(fitTotal->GetParameter(0),
                        fitTotal->GetParameter(3),
                        fitTotal->GetParameter(4));

  TF1* signalPrime = new TF1("signalPrime","gaus",2,8);
  
  signalPrime->SetParameters(fitTotal->GetParameter(0),
                             fitTotal->GetParameter(8),
                             fitTotal->GetParameter(9));
  
  Double_t gausParameters[3];
  Double_t gausParametersPrime[3];
  Double_t covarianceMatrix[3][3];
  Double_t covarianceMatrixPrime[3][3];
  
  gausParameters[0] = fitTotal->GetParameter(0);
  gausParameters[1] = fitTotal->GetParameter(3);
  gausParameters[2] = fitTotal->GetParameter(4);

  covarianceMatrix[0][0] = (r->GetCovarianceMatrix())(0,0);
  covarianceMatrix[1][0] = (r->GetCovarianceMatrix())(3,0);
  covarianceMatrix[2][0] = (r->GetCovarianceMatrix())(4,0);
  covarianceMatrix[0][1] = (r->GetCovarianceMatrix())(0,3);
  covarianceMatrix[0][2] = (r->GetCovarianceMatrix())(0,4);
  
  for ( int iy = 1; iy < 3; ++iy )
  {
    for ( int ix = 1; ix < 3; ++ix )
    {
      covarianceMatrix[ix][iy] = (r->GetCovarianceMatrix())(ix+2,iy+2);
    }
  }

  gausParametersPrime[0] = fitTotal->GetParameter(0);
  gausParametersPrime[1] = fitTotal->GetParameter(8);
  gausParametersPrime[2] = fitTotal->GetParameter(9);

  covarianceMatrixPrime[0][0] = (r->GetCovarianceMatrix())(0,0);
  covarianceMatrixPrime[1][0] = (r->GetCovarianceMatrix())(8,0);
  covarianceMatrixPrime[2][0] = (r->GetCovarianceMatrix())(9,0);
  covarianceMatrixPrime[0][1] = (r->GetCovarianceMatrix())(0,8);
  covarianceMatrixPrime[0][2] = (r->GetCovarianceMatrix())(0,9);
  
  for ( int iy = 1; iy < 3; ++iy )
  {
    for ( int ix = 1; ix < 3; ++ix )
    {
      covarianceMatrixPrime[ix][iy] = (r->GetCovarianceMatrix())(ix+2,iy+2);
    }
  }
  
  double n = signal->Integral(2,6)/h.GetBinWidth(10);
  double nerr = signal->IntegralError(2,6,&gausParameters[0],&covarianceMatrix[0][0])/h.GetBinWidth(10);

  Set("NofJpsi",n,nerr);      
  Set("MeanJpsi",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("SigmaJpsi",fitTotal->GetParameter(4),fitTotal->GetParError(4));

  double nprime = signalPrime->Integral(2,6)/h.GetBinWidth(10);
  double nerrprime = signalPrime->IntegralError(2,6,&gausParametersPrime[0],&covarianceMatrixPrime[0][0])/h.GetBinWidth(10);
  
  Set("NofPsiPrime",nprime,nerrprime);
  Set("MeanPsiPrime",fitTotal->GetParameter(8),fitTotal->GetParError(8));
  Set("SigmaPsiPrime",fitTotal->GetParameter(9),fitTotal->GetParError(9));
  
}
  */

//_____________________________________________________________________________
AliAnalysisMuMuResult* AliAnalysisMuMuResult::FitJpsiGCBE(TH1& h)
{
  /// Fit Jpsi spectra with crystal ball + gaussian + exponential
  
  std::cout << "Fit with jpsi alone (gaus + CB + expo)" << std::endl;
  
  Int_t nrebin = fMinv->GetXaxis()->GetNbins() / h.GetXaxis()->GetNbins();
  
  AliAnalysisMuMuResult* r = new AliAnalysisMuMuResult(h,"JPSIGCBE",nrebin);
  
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
  
  const char* fitOption = "SI+";
  
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

/*
//_____________________________________________________________________________
void AliAnalysisMuMuResult::FitJpsiPCBE(TH1& h)
{
  std::cout << "Fit with jpsi alone (pol2 + CB + expo)" << std::endl;
  
  const Double_t xmin(2.0);
  const Double_t xmax(5.0);
  
  fitTotal = new TF1("fitTotal",funcJpsiPCBE,xmin,xmax,10);
  fitTotal->SetParNames("p0","p1","p2","N","alpha","n","mean","sigma","expocste","exposlope");
  
  fitTotal->SetParLimits(3,0,h.GetMaximum()*2); // N

  const Double_t cbalpha(0.98);
  const Double_t cbn(5.2);
  
  fitTotal->FixParameter(4,cbalpha);
  fitTotal->FixParameter(5,cbn);
  
  fitTotal->SetParLimits(6,2,4); // mean
  fitTotal->SetParLimits(7,0.05,0.2); // sigma
  
  TF1* fpol2 = new TF1("pol2","pol2",xmin,xmax);
                       
  h.Fit(fpol2,"+","",2,2.8);
  
  fitTotal->SetParameter(0,fpol2->GetParameter(0));
  fitTotal->SetParameter(1,fpol2->GetParameter(1));
  fitTotal->SetParameter(2,fpol2->GetParameter(2));

  TF1* fexpo = new TF1("expo","expo",xmin,xmax);
  
  h.Fit(fexpo,"+","",3.5,4.5);
  
  fitTotal->SetParameter(8,fexpo->GetParameter(0));
  fitTotal->SetParameter(9,fexpo->GetParameter(1));
    
  fitTotal->SetParameter(3,h.GetMaximum()),
  fitTotal->SetParameter(4,cbalpha);
  fitTotal->SetParameter(5,cbn);
  fitTotal->SetParameter(6,3.15);
  fitTotal->SetParameter(7,0.1);
  
  h.Fit(fitTotal,"+","",2.5,5);
    
  Set("MeanJpsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  Set("SigmaJpsi",fitTotal->GetParameter(7),fitTotal->GetParError(7));
  
  double m = GetValue("MeanJpsi");
  double s = GetValue("SigmaJpsi");
  double n = 2.0;
  
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
  
  
  Set("NofJpsi",fitTotal->Integral(m-n*s,m+n*s)/h.GetBinWidth(1),fitTotal->IntegralError(m-n*s,m+n*s)/h.GetBinWidth(1));
  
  //  Set("NofJpsi",fitTotal->Integral(0,10)/h.GetBinWidth(1),fitTotal->IntegralError(0,10)/h.GetBinWidth(1));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::FitJpsiCBE(TH1& h)
{
  std::cout << "Fit with jpsi alone" << std::endl;
  
  const Double_t xmin(1.5);
  const Double_t xmax(8.0);
  
  fitTotal = new TF1("fitTotal",funcJpsiCBE,xmin,xmax,7);
  fitTotal->SetParNames("N","alpha","n","mean","sigma","expocste","exposlope");
  
//  fitTotal->SetParameters(h.GetMaximum(),1,5,3.0,0.07,1.5,3,1,0);

  fitTotal->SetParameters(1,1.15,3.6,3.0,0.07,1,-1);

  fitTotal->SetParLimits(0,0,h.GetMaximum()); // N
//  fitTotal->SetParLimits(1,0.1,2); // alpha
  fitTotal->FixParameter(1,0.98);
//  fitTotal->SetParLimits(2,0.01,5); // n
  fitTotal->FixParameter(2,5.2);
  fitTotal->SetParLimits(3,2.8,3.5); // mean
  fitTotal->SetParLimits(4,0.05,0.2); // sigma
  
  TF1* fexpo = new TF1("expo","expo",xmin,xmax);
  
  h.Fit(fexpo,"+","",2,3);
  
  fitTotal->SetParameter(5,fexpo->GetParameter(0));
  fitTotal->SetParameter(6,fexpo->GetParameter(1));
  
  h.Fit(fitTotal,"+","",2,5);
  
  
  Set("MeanJpsi",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("SigmaJpsi",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  
  double m = GetValue("MeanJpsi");
  double s = GetValue("SigmaJpsi");
  double n = 3.0;
  
  TF1* fcb = new TF1("fcb",funcCB,xmin,xmax,5);
  fcb->SetParameters(fitTotal->GetParameter(0),
                     fitTotal->GetParameter(1),
                     fitTotal->GetParameter(2),
                     fitTotal->GetParameter(3),
                     fitTotal->GetParameter(4));

  fcb->SetLineColor(6);
  fcb->SetNpx(1000);
  TLine* l1 = new TLine(m-n*s,0,m-n*s,fitTotal->GetParameter(0));
  TLine* l2 = new TLine(m+n*s,0,m+n*s,fitTotal->GetParameter(0));
  l1->SetLineColor(6);
  l2->SetLineColor(6);
  h.GetListOfFunctions()->Add(fcb);
  h.GetListOfFunctions()->Add(l1);
  h.GetListOfFunctions()->Add(l2);
  
  
  Set("NofJpsi",fitTotal->Integral(m-n*s,m+n*s)/h.GetBinWidth(1),fitTotal->IntegralError(m-n*s,m+n*s)/h.GetBinWidth(1));
  
  //  Set("NofJpsi",fitTotal->Integral(0,10)/h.GetBinWidth(1),fitTotal->IntegralError(0,10)/h.GetBinWidth(1));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::FitJpsiECBE(TH1& h)
{
  std::cout << "Fit with jpsi alone (expo + CB + expo)" << std::endl;
  
  const Double_t xmin(1.5);
  const Double_t xmax(8.0);
  
  fitTotal = new TF1("fitTotal",funcJpsiECBE,xmin,xmax,9);
  fitTotal->SetParNames("e0","s0","N","alpha","n","mean","sigma","e1","s1");
  
  fitTotal->SetParameters(1,-1,1,1.15,3.6,3.2,0.06,-1);

  fitTotal->SetParLimits(0,0,h.GetMaximum()*2);
  
  fitTotal->FixParameter(3,0.98); // alpha
  fitTotal->FixParameter(4,5.2); // n
  fitTotal->SetParLimits(5,2.8,3.5); // mean
  fitTotal->SetParLimits(6,0.05,0.2); // sigma
  
  TF1* fexpo1 = new TF1("expo1","expo",xmin,xmax);
  TF1* fexpo2 = new TF1("expo2","expo",xmin,xmax);
  
  h.Fit(fexpo1,"","",1.5,3);
  
  fitTotal->SetParameter(0,fexpo1->GetParameter(0));
  fitTotal->SetParameter(1,fexpo1->GetParameter(1));

  h.Fit(fexpo2,"","",3.5,5.0);

  fitTotal->SetParameter(7,fexpo2->GetParameter(0));
  fitTotal->SetParameter(8,fexpo2->GetParameter(1));

  const char* fitOption = "SI+";
  
  TFitResultPtr r = h.Fit(fitTotal,fitOption,"",2,5);
  
  Set("MeanJpsi",fitTotal->GetParameter(5),fitTotal->GetParError(5));
  Set("SigmaJpsi",fitTotal->GetParameter(6),fitTotal->GetParError(6));
  
  double m = GetValue("MeanJpsi");
  double s = GetValue("SigmaJpsi");
  double n = 3.0;
  
  TF1* fcb = new TF1("fcb",funcCB,xmin,xmax,5);
  fcb->SetParameters(fitTotal->GetParameter(2),
                     fitTotal->GetParameter(3),
                     fitTotal->GetParameter(4),
                     fitTotal->GetParameter(5),
                     fitTotal->GetParameter(6));

  fcb->SetParError(0,fitTotal->GetParError(2));
  fcb->SetParError(1,fitTotal->GetParError(3));
  fcb->SetParError(2,fitTotal->GetParError(4));
  fcb->SetParError(3,fitTotal->GetParError(5));
  fcb->SetParError(4,fitTotal->GetParError(6));
  
  fcb->SetLineColor(6);
  fcb->SetNpx(1000);
  TLine* l1 = new TLine(m-n*s,0,m-n*s,fitTotal->GetParameter(2));
  TLine* l2 = new TLine(m+n*s,0,m+n*s,fitTotal->GetParameter(2));
  l1->SetLineColor(6);
  l2->SetLineColor(6);
  h.GetListOfFunctions()->Add(fcb);
  h.GetListOfFunctions()->Add(l1);
  h.GetListOfFunctions()->Add(l2);
  
  
  Double_t cbParameters[5];
  Double_t covarianceMatrix[5][5];
  
  cbParameters[0] = fitTotal->GetParameter(2);
  cbParameters[1] = fitTotal->GetParameter(3);
  cbParameters[2] = fitTotal->GetParameter(4);
  cbParameters[3] = fitTotal->GetParameter(5);
  cbParameters[4] = fitTotal->GetParameter(6);
  
  for ( int iy = 0; iy < 5; ++iy )
  {
    for ( int ix = 0; ix < 5; ++ix )
    {
      covarianceMatrix[ix][iy] = (r->GetCovarianceMatrix())(ix+2,iy+2);
    }
  }
  

  double njpsi = fcb->Integral(m-n*s,m+n*s)/h.GetBinWidth(1);
  
  double nerr = fcb->IntegralError(m-n*s,m+n*s,&cbParameters[0],&covarianceMatrix[0][0])/h.GetBinWidth(1);
  

  Set("NofJpsi",njpsi,nerr);
}

 */

//_____________________________________________________________________________
AliAnalysisMuMuResult* AliAnalysisMuMuResult::FitJpsi(TH1& h)
{
  /// Fit Jpsi spectra using extended crystall ball (CB2) with free tails
  
  std::cout << "Fit with jpsi alone" << std::endl;

  Int_t nrebin = fMinv->GetXaxis()->GetNbins() / h.GetXaxis()->GetNbins();
  
  AliAnalysisMuMuResult* r = new AliAnalysisMuMuResult(h,"JPSI",nrebin);
  
  TH1* hfit = r->Minv();

  const Double_t xmin(1.5);
  const Double_t xmax(8.0);

  TF1* fitTotal = new TF1("fitTotal",funcCB2,xmin,xmax,7);
  fitTotal->SetParNames("N","alphaLow","nLow","mean","sigma","alphaUp","nUp");
  fitTotal->SetParameters(h.GetMaximum(),1,5,3.1,0.07,1.5,3);
  fitTotal->SetParLimits(0,0,h.GetMaximum()*2); // N
  fitTotal->SetParLimits(1,0,10); // alpha
  fitTotal->SetParLimits(2,0.1,10); // n
  fitTotal->SetParLimits(3,3,3.1); // mean
  fitTotal->SetParLimits(4,0.01,1); // sigma
  fitTotal->SetParLimits(5,0,10); // alpha
  fitTotal->SetParLimits(6,0.1,10); // n
  
  hfit->Fit(fitTotal,"+","",2,5);
  
  
  r->Set("MeanJpsi",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  r->Set("SigmaJpsi",fitTotal->GetParameter(4),fitTotal->GetParError(4));

  double m = r->GetValue("MeanJpsi");
  double s = r->GetValue("SigmaJpsi");
  double n = 10.0;

  r->Set("NofJpsi",fitTotal->Integral(m-n*s,m+n*s)/h.GetBinWidth(1),fitTotal->IntegralError(m-n*s,m+n*s)/h.GetBinWidth(1));

  return r;
}

//_____________________________________________________________________________
AliAnalysisMuMuResult* AliAnalysisMuMuResult::FitJpsiCB2VWG(const TH1& h)
{
  /// Fit Jpsi spectra using extended crystal ball (CB2) + variable width gaussian (VWG)
  
  std::cout << "Fit with jpsi VWG" << std::endl;
  
  Int_t nrebin = fMinv->GetXaxis()->GetNbins() / h.GetXaxis()->GetNbins();
  
  AliAnalysisMuMuResult* r = new AliAnalysisMuMuResult(h,"JPSICB2VWG",nrebin);
  
  
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
  
  const char* fitOption = "SIER"; //+";
  
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
AliAnalysisMuMuResult* AliAnalysisMuMuResult::FitJpsi2CB2VWG(const TH1& h,
                                                             Double_t alphaLow,
                                                             Double_t nLow,
                                                             Double_t alphaUp,
                                                             Double_t nUp)
{
  /// Fit using extended crystal ball + variable width gaussian
  
  std::cout << Form("Fit with jpsi + psiprime VWG alphaLow=%5.2f nLow=%5.2f alphaUp=%5.2f nUp=%5.2f",
                    alphaLow,nLow,alphaUp,nUp) << std::endl;
  
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
  
  AliAnalysisMuMuResult* r = new AliAnalysisMuMuResult(h,resultName.Data(),nrebin);
  
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

  const char* fitOption = "SER"; //+";
  
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
AliAnalysisMuMuResult* AliAnalysisMuMuResult::FitJpsiNA48(const TH1& h)
{
  /// fit using functional form from Ruben Shahoyan's thesis (2001) (eq. 4.8.)
  
  std::cout << "Fit with jpsi NA50 Ruben eq. 4.8" << std::endl;
  
  Int_t nrebin = fMinv->GetXaxis()->GetNbins() / h.GetXaxis()->GetNbins();
  
  AliAnalysisMuMuResult* r = new AliAnalysisMuMuResult(h,"JPSINA",nrebin);
  
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
  
  const char* fitOption = "SIER"; //+";
  
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

/*
//_____________________________________________________________________________
void AliAnalysisMuMuResult::FitUpsilon(TH1& h)
{
  std::cout << "Fit with upsilon alone" << std::endl;
  
  const Double_t xmin(6.0);
  const Double_t xmax(12.0);
  
  fitTotal = new TF1("fitTotal",funcCB2,xmin,xmax,7);
  fitTotal->SetParNames("N","alpha","n","mean","sigma","alphaprime","nprime");
  fitTotal->SetParameters(h.GetMaximum(),1,5,9.46,0.2,1.5,3);
  fitTotal->SetParLimits(0,0,h.GetMaximum()*2); // N
  fitTotal->SetParLimits(1,0,10); // alpha
  fitTotal->SetParLimits(2,1,10); // n
  fitTotal->SetParLimits(3,8,12); // mean
  fitTotal->SetParLimits(4,0.01,1); // sigma
  fitTotal->SetParLimits(5,0,10); // alpha
  fitTotal->SetParLimits(6,1,10); // n
  
  h.Fit(fitTotal,"+","",6,12);
  
  
  Set("MeanUpsilon",fitTotal->GetParameter(3),fitTotal->GetParError(3));
  Set("SigmaUpsilon",fitTotal->GetParameter(4),fitTotal->GetParError(4));
  
  double m = GetValue("MeanUpsilon");
  double s = GetValue("SigmaUpsilon");
  double n = 3.0;
  
  Set("NofUpsilon",fitTotal->Integral(m-n*s,m+n*s)/h.GetBinWidth(1),fitTotal->IntegralError(m-n*s,m+n*s)/h.GetBinWidth(1));  
}
*/

//_____________________________________________________________________________
Double_t AliAnalysisMuMuResult::ErrorAB(Double_t a, Double_t aerr, Double_t b, Double_t berr)
{
  /// Compute the quadratic sum of 2 errors

  Double_t e(0.0);
  
  if ( TMath::Abs(a) > 1E-12 )
  {
    e += (aerr*aerr)/(a*a);
  }
  
  if ( TMath::Abs(b) > 1E-12 )
  {
    e += (berr*berr)/(b*b);
  }
  
  return TMath::Sqrt(e);
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuResult::ErrorABC(Double_t a, Double_t aerr, Double_t b, Double_t berr, Double_t c, Double_t cerror)
{
  /// Compute the quadratic sum of 3 errors
  
  Double_t e(0.0);
  
  if ( TMath::Abs(a) > 1E-12 )
  {
    e += (aerr*aerr)/(a*a);
  }
  
  if ( TMath::Abs(b) > 1E-12 )
  {
    e += (berr*berr)/(b*b);
  }
  
  if ( TMath::Abs(b) > 1E-12 )
  {
    e += (cerror*cerror)/(c*c);
  }
  
  return TMath::Sqrt(e);
}


//_____________________________________________________________________________
Bool_t AliAnalysisMuMuResult::AddFit(const char* fitType, Int_t npar, Double_t* par)
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

  AliAnalysisMuMuResult* r(0x0);
  
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
    r = new AliAnalysisMuMuResult(*hminv,"COUNTJPSI",1);
    Double_t n = CountParticle(*hminv,"Jpsi");
    r->Set("NofJpsi",n,TMath::Sqrt(n));
  }
  
  
  if ( r )
  {
    r->Print();
    r->SetBin(Bin());
    r->SetNofTriggers(NofTriggers());
    r->SetNofRuns(NofRuns());

    if (!fSubResults)
    {
      fSubResults = new TObjArray;
      fSubResults->SetOwner(kTRUE);
    }

    fSubResults->Add(r);
  }
  
  delete hminv;
  
  return (r!=0x0);
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuResult::GetErrorStat(const char* name, const char* subResultName) const
{
  // compute the mean value from all subresults
  
  if ( strlen(subResultName) > 0 )
  {
    if ( !fSubResults)
    {
      AliError(Form("No subresult from which I could get the %s one...",subResultName));
      return TMath::Limits<Double_t>::Max();
    }
    AliAnalysisMuMuResult* sub = static_cast<AliAnalysisMuMuResult*>(fSubResults->FindObject(subResultName));
    if (!sub)
    {
      AliError(Form("Could not get subresult named %s",subResultName));
      return TMath::Limits<Double_t>::Max();
    }
    return sub->GetErrorStat(name);
  }
  
  if ( fMap )
  {
    TObjArray* p = static_cast<TObjArray*>(fMap->GetValue(name));
    if (p)
    {
      TParameter<double>* val = static_cast<TParameter<double>*>(p->At(kErrorStat));
      return val->GetVal();
    }
  }
  
  TIter next(fSubResults);
  AliAnalysisMuMuResult* r;
  Int_t n(0);
  Double_t mean(0);
  
  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
  {
    if ( r->HasValue(name) )
    {
      mean += r->GetErrorStat(name);
      ++n;
    }
  }
  return ( n ? mean/n : 0.0 );
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuResult::GetValue(const char* name, const char* subResultName) const
{
  // get a value (either directly or by computing the mean of the subresults)
  
  if ( strlen(subResultName) > 0 )
  {
    if ( !fSubResults)
    {
      AliError(Form("No subresult from which I could get the %s one...",subResultName));
      return TMath::Limits<Double_t>::Max();
    }
    AliAnalysisMuMuResult* sub = static_cast<AliAnalysisMuMuResult*>(fSubResults->FindObject(subResultName));
    if (!sub)
    {
      AliError(Form("Could not get subresult named %s",subResultName));
      return TMath::Limits<Double_t>::Max();
    }
    return sub->GetValue(name);
  }
  
  if (fMap)
  {
    TObjArray* p = static_cast<TObjArray*>(fMap->GetValue(name));
    if (p)
    {
      TParameter<double>* val = static_cast<TParameter<double>*>(p->At(kValue));
      return val->GetVal();
    }
  }
  
  // compute the mean value from all subresults
  TIter next(fSubResults);
  AliAnalysisMuMuResult* r;
  Int_t n(0);
  Double_t mean(0);
  
  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
  {
    if ( r->HasValue(name) )
    {
      mean += r->GetValue(name);
      ++n;
    }
  }
  return ( n ? mean/n : 0.0 );
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuResult::HasValue(const char* name, const char* subResultName) const
{
  /// Whether this result (or subresult if subResultName is provided) has a property
  /// named "name"
  
  if ( strlen(subResultName) > 0 )
  {
    if ( !fSubResults)
    {
      AliError(Form("No subresult from which I could get the %s one...",subResultName));
      return kFALSE;
    }
    AliAnalysisMuMuResult* sub = static_cast<AliAnalysisMuMuResult*>(fSubResults->FindObject(subResultName));
    if (!sub)
    {
      AliError(Form("Could not get subresult named %s",subResultName));
      return kFALSE;
    }
    return sub->HasValue(name);
  }

  if ( fMap && ( fMap->GetValue(name) != 0x0 ) )
  {
    return kTRUE;
  }
  
  TIter next(fSubResults);
  AliAnalysisMuMuResult* r;

  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
  {
    if ( r->HasValue(name) ) return kTRUE;
  }
  
  return kFALSE;
}

//_____________________________________________________________________________
THashList* AliAnalysisMuMuResult::Keys() const
{
  /// Return the complete list of keys we're using
  if (!fKeys)
  {
    fKeys = new THashList;
    fKeys->SetOwner(kTRUE);
    TIter next(fMap);
    TObjString* key;
    
    while ( ( key = static_cast<TObjString*>(next()) ) )
    {
      if ( !fKeys->FindObject(key->String()) )
      {
        fKeys->Add(new TObjString(key->String()));
      }
    }
    
    AliAnalysisMuMuResult* r;
    TIter nextResult(fSubResults);
    
    while ( ( r = static_cast<AliAnalysisMuMuResult*>(nextResult())) )
    {
      TIter nextHL(r->Keys());
      TObjString* s;
      
      while ( ( s = static_cast<TObjString*>(nextHL())) )
      {
        if ( !fKeys->FindObject(s->String()) )
        {
          fKeys->Add(new TObjString(s->String()));
        }
      }
    }

  }
  return fKeys;
}

//_____________________________________________________________________________
Long64_t AliAnalysisMuMuResult::Merge(TCollection* list)
{
  /// Merge method
  ///
  /// Merge a list of AliAnalysisMuMuResult objects with this
  /// Returns the number of merged objects (including this).
  ///
  /// Note that the merging is to be understood here as an average operation
  
  AliInfo(Form("this=%p",this));
  if (!list) return 0;
  
  if (list->IsEmpty()) return 1;
  
  TIter next(list);
  TObject* currObj;
  TList hminvList;
  hminvList.SetOwner(kTRUE);
  
  Double_t thisWeight = Weight();
  Double_t sumOfWeights = thisWeight;

  Double_t nofTriggers = fNofTriggers*thisWeight;
  Double_t nofRuns = fNofRuns*thisWeight;
  
  fNofRuns *= thisWeight;
  
  while ( ( currObj = next() ) )
  {
    AliAnalysisMuMuResult* result = dynamic_cast<AliAnalysisMuMuResult*>(currObj);
    if (!result)
    {
      AliFatal(Form("object named \"%s\" is a %s instead of an AliAnalysisMuMuResult!", currObj->GetName(), currObj->ClassName()));
      continue;
    }
    
    Double_t w = result->Weight();
    
    nofRuns += result->NofRuns()*w;
    nofTriggers += result->NofTriggers()*w;
    fWeight += result->fWeight;
    sumOfWeights += w;
  }
  
  thisWeight/= sumOfWeights;
  fNofRuns = nofRuns/sumOfWeights;
  fNofTriggers = nofTriggers/sumOfWeights;
  fWeight /= sumOfWeights;
  
  AliInfo(Form("thisWeight=%e sumOfWeight=%8.2f noftriggers=%e weight=%e",thisWeight,sumOfWeights,1.0*fNofTriggers,fWeight));
  
  TIter nextKey(fMap);
  TObjString* key;
  
  while ( ( key = static_cast<TObjString*>(nextKey())) )
  {
    AliInfo(key->String().Data());

    Double_t value = GetValue(key->String())*thisWeight;
    Double_t estat = GetErrorStat(key->String())*GetErrorStat(key->String())*thisWeight*thisWeight;

    Double_t test(thisWeight);

    next.Reset();
    
    while ( ( currObj = next() ) )
    {
      AliAnalysisMuMuResult* result = dynamic_cast<AliAnalysisMuMuResult*>(currObj);
    
      if (!result)
      { 
        continue;
      }
      
      if (!result->HasValue(key->String()))
      {
        AliError(Form("Did not find key %s in of the result to merge",key->String().Data()));
        continue;
      }
      
      // can only merge under the condition we have the same bin
    
      if ( fBin != result->Bin() )
      {
        AliError("Cannot merge results of different bin");
        continue;
      }
    
      Double_t w = result->Weight()/sumOfWeights;
      
      Double_t w2 = w*w;
      
      test += w;
      
      value += result->GetValue(key->String())*w;
      estat += result->GetErrorStat(key->String())*result->GetErrorStat(key->String())*w2;
      
    }
    
    Set(key->String(),
        value,
        TMath::Sqrt(estat));
    
    AliInfo(Form("test=%e",test));
  }

  if ( fMinv )
  {
    fMinv->Scale(thisWeight);    
  }
  
  next.Reset();
  
  while ( ( currObj = next() ) )
  {
    AliAnalysisMuMuResult* result = dynamic_cast<AliAnalysisMuMuResult*>(currObj);

    if ( result->Minv() )
    {
      TH1* h = static_cast<TH1*>(result->Minv()->Clone());
      AliInfo(Form("Nbins %d xmin %e xmax %e",h->GetXaxis()->GetNbins(),h->GetXaxis()->GetXmin(),
                   h->GetXaxis()->GetXmax()));
      h->Scale(result->Weight());
      hminvList.Add(h);
    }
  }
  
  if ( fMinv )
  {
    fMinv->Merge(&hminvList);
    fMinv->Scale(1.0/sumOfWeights);
  }
  
  TIter nextSubresult(fSubResults);
  AliAnalysisMuMuResult* r;
  
  while ( ( r = static_cast<AliAnalysisMuMuResult*>(nextSubresult())) )
  {
    TList sublist;

    next.Reset();
    
    while ( ( currObj = next() ) )
    {
      sublist.Add(currObj);
    }

    r->Merge(&sublist);
  }
  
  return list->GetEntries()+1;
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuResult::NofRuns() const
{
  /// Get the number of runs
  if ( !Mother() ) return fNofRuns;
  else return Mother()->NofRuns();
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuResult::NofTriggers() const
{
  /// Get the number of triggers
  
  if ( !Mother() ) return fNofTriggers;
  else return Mother()->NofTriggers();
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::Print(Option_t* opt) const
{
  /// printout
  
  TString sopt(opt);
  sopt.ToUpper();
  
  for ( Int_t i = 0; i < 9; ++i )
  {
    sopt.ReplaceAll(Form("%d",i),"");
  }
  
  TString pot(sopt);
  pot.ReplaceAll("ALL","");
  pot.ReplaceAll("FULL","");

  std::cout << pot.Data();
  
  if ( fAlias.Length() > 0 )
  {
    std::cout << Form("%s - ",fAlias.Data());
  }
  
  std::cout << Form("%50s - NRUNS %d - NTRIGGER %10d - %s%s",
               GetName(),
                    NofRuns(),
                    NofTriggers(),
                                  fWeight > 0.0 ? Form(" WEIGHT %e -",fWeight) : "",
                           fBin.AsString().Data());
  
    if ( fSubResults && fSubResults->GetEntries()>1 )
    {
      std::cout << " (" << fSubResults->GetEntries()-1 << " subresults)";
    }
  
  if (!fSubResults )
  {
    std::cout << Form(" - REBIN %d",fRebin) << std::endl;
  }
  
  std::cout << std::endl;

  if ( sopt.Contains("DUMP") )
  {
    TIter next(fMap);
    TObjString* str;
    while ( ( str = static_cast<TObjString*>(next()) ) )
    {
      TObjArray* a = static_cast<TObjArray*>(fMap->GetValue(str->String()));

      std::cout << Form("%s %e %e",
                        str->String().Data(),
                        static_cast<TParameter<Double_t>*> (a->At(kValue))->GetVal(),
                        static_cast<TParameter<Double_t>*> (a->At(kErrorStat))->GetVal()) << std::endl;
    }
  }
  else
  {
    
    PrintParticle("Jpsi",pot.Data());
    PrintParticle("PsiPrime",pot.Data());
    PrintParticle("Upsilon",pot.Data());
  }
  
  if ( HasValue("MBR"))
  {
    std::cout << Form("\t\tMBR %e +- %e",GetValue("MBR"),GetErrorStat("MBR")) << std::endl;
  }
  
  if ( fSubResults /* && fSubResults->GetEntries() > 1 */ && ( sopt.Contains("ALL") || sopt.Contains("FULL") ) )
  {
    std::cout << pot.Data() << "\t===== sub results ===== " << std::endl;
    
    sopt += "\t\t";
    
    TIter next(fSubResults);
    AliAnalysisMuMuResult* r;
    
    while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
    {
      r->Print(sopt.Data());
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::PrintParticle(const char* particle, const char* opt) const
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
    
    PrintValue(key->String(),opt,GetValue(key->String()),GetErrorStat(key->String()));
  } 
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::PrintValue(const char* key, const char* opt, Double_t value, Double_t errorStat)
{
  // print one value and its associated error
  
  if ( TString(key).Contains("AccEff") )
  {
    std::cout << opt << Form("\t\t%20s %9.2f +- %5.2f %%",key,value*100,errorStat*100) << std::endl;
  }
  else if ( TString(key).BeginsWith("Sigma") ||  TString(key).BeginsWith("Mean") )
  {
    std::cout << opt << Form("\t\t%20s %9.2f +- %5.2f MeV/c^2",key,value*1E3,1E3*errorStat) << std::endl;
  }
  else if ( TString(key).Contains("Nof") )
  {
    std::cout << opt << Form("\t\t%20s %9.2f +- %5.2f",key,value,errorStat) << std::endl;
  }
  else if ( value > 1E-3 && value < 1E3 )
  {
    std::cout << opt << Form("\t\t%20s %9.2f +- %5.2f",key,value,errorStat) << std::endl;
  }
  else
  {
    std::cout << opt << Form("\t\t%20s %9.2e +- %9.2e",key,value,errorStat) << std::endl;    
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::Set(const char* name, Double_t value, Double_t errorStat)
{
  /// Set a (value,error) pair with a given name
  
  if ( !fMap )
  {
    fMap = new TMap;
    fMap->SetOwnerKeyValue(kTRUE,kTRUE);
  }
  
  TObjArray* p = static_cast<TObjArray*>(fMap->GetValue(name));
  if (!p)
  {
    p = new TObjArray(4);
    
    p->SetOwner(kTRUE);
    
    p->AddAt(new TParameter<Double_t>(name,value),kValue);
    p->AddAt(new TParameter<Double_t>(name,errorStat),kErrorStat);
    
    fMap->Add(new TObjString(name),p);
  }
  else
  {
    static_cast<TParameter<double>*>(p->At(kValue))->SetVal(value);
    static_cast<TParameter<double>*>(p->At(kErrorStat))->SetVal(errorStat);
  }
  
  if ( TString(name)=="NofJpsi" )
  {
    if ( NofTriggers() > 0 )
    {
      Double_t rate = value/NofTriggers();
      Double_t rateError = rate*ErrorAB(value,errorStat,NofTriggers(),TMath::Sqrt(NofTriggers()));
      Set("RateJpsi",rate,rateError);
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::SetBin(const AliAnalysisMuMuBinning::Range& bin)
{
  /// Set the bin
  
  if (!Mother()) fBin = bin;
  else Mother()->SetBin(bin);
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::SetNofInputParticles(const TH1& hminv)
{
  /// Set the number of input particle from the invariant mass spectra
  
  static const char* particleNames[] = { "Jpsi" , "PsiPrime", "Upsilon","UpsilonPrime" };

  for ( Int_t i = 0; i < 4; ++i )
  {
    Double_t n = CountParticle(hminv,particleNames[i]);

    AliDebug(1,Form("i=%d particle %s n %e",i,particleNames[i],n));
    
    if ( n > 0 )
    {
      SetNofInputParticles(particleNames[i],n);
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::SetNofInputParticles(const char* particle, int n)
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
  
  TIter next(fSubResults);
  AliAnalysisMuMuResult* r;
  
  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next())) )
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
void AliAnalysisMuMuResult::SetNofRuns(Int_t n)
{
  if ( !Mother() ) fNofRuns=n;
  else Mother()->SetNofRuns(n);
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::SetNofTriggers(Int_t n)
{
  if ( !Mother() ) fNofTriggers=n;
  else Mother()->SetNofTriggers(n);
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::SetMinv(const TH1& hminv)
{
    /// Set the inv. mass spectrum to be fitted.
  static UInt_t n(0);
  
  delete fMinv;
  fMinv = static_cast<TH1*>(hminv.Clone(Form("Minv%u",n++)));
  fMinv->SetDirectory(0);
}

//_____________________________________________________________________________
AliAnalysisMuMuResult*
AliAnalysisMuMuResult::SubResult(const char* subResultName) const
{
  /// get a given subresult
  if (!fSubResults)
  {
    return 0x0;
  }
  TIter next(fSubResults);
  AliAnalysisMuMuResult* r;
  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next())) )
  {
    if ( r->Alias() == subResultName )
    {
      return r;
    }
  }
  return 0x0;
}

