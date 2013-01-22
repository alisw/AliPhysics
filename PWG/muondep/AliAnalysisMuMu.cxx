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

#include "AliAnalysisMuMu.h"

#include "AliAnalysisTriggerScalers.h"
#include "AliCounterCollection.h"
#include "AliHistogramCollection.h"
#include "AliLog.h"
#include "Riostream.h"
#include "TArrayL64.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGrid.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TList.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include <map>
#include <set>
#include <string>
#include "TParameter.h"
#include "TMap.h"
#include "TFitResult.h"
#include "TLine.h"
#include "TASImage.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TColor.h"

ClassImp(AliAnalysisMuMu)

void Add(TObjArray* a, AliAnalysisMuMu::Result* r)
{
  if ( r ) a->Add(r);
}

void ALICEseal(Double_t xPad, Double_t yPad, TList* extraLines)
{
  TVirtualPad* currPad = gPad;
  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.76,0.7,0.90,0.87);
  myPadLogo->SetFillColor(kWhite);
  myPadLogo->SetFillStyle(1001);
  myPadLogo->SetBorderMode(0);
  myPadLogo->SetBorderSize(2);
  myPadLogo->SetFrameBorderMode(0);
  myPadLogo->SetLeftMargin(0.0);
  myPadLogo->SetTopMargin(0.0);
  myPadLogo->SetBottomMargin(0.0);
  myPadLogo->SetRightMargin(0.0);
  myPadLogo->Draw();
  myPadLogo->cd();
  TASImage *myAliceLogo = new TASImage("$HOME/Pictures/2011-Nov-24-ALICE_PERFORMANCE_logo_BLACK_small_usage_design.gif");
  myAliceLogo->Draw();
  currPad->cd();
  Double_t x1 = xPad - 0.07, y1 = yPad - 0.06;
  Double_t x2 = x1 + 0.25, y2 = y1 + 0.08;
  TPaveText* t2=new TPaveText(x1+0.06,y1-0.06,x2-0.06,y2-0.06,"NDC");
  //  t2->SetFillStyle(0);
  t2->SetFillColor(kWhite);
  t2->SetBorderSize(0);
  t2->SetTextColor(kBlack);
  t2->SetTextFont(52);
  t2->SetTextSize(0.035);
  TDatime dt;
  TString today = Form("%02i/%02i/%4i", dt.GetDay(), dt.GetMonth(), dt.GetYear());
  t2->AddText(0.,0.,today.Data());
  t2->Draw();
  
  if ( extraLines ) 
  {
    int n = extraLines->GetSize();
    TPaveText* t3 = new TPaveText(xPad-0.07,yPad-0.06*2,xPad-0.07+0.25,yPad-0.06*(n+3),"NDC");
    t3->SetFillColor(kWhite);
    t3->SetBorderSize(0);
    t3->SetTextColor(kBlack);
    t3->SetTextFont(42);
    t3->SetTextSize(0.035);
    
    TIter next(extraLines);
    TObjString* str;
    
    while ( ( str = static_cast<TObjString*>(next()) ) )
    {
      t3->AddText(str->String().Data());
    }
    t3->SetFillColor(kWhite);
    t3->Draw();
  }
}

Double_t funcCB(Double_t* xx, Double_t* par)
{ 
  Double_t N = par[0];
  Double_t alpha = par[1];
  Double_t n = par[2];
  Double_t mean = par[3];
  Double_t sigma = par[4];
  
  Double_t x = xx[0];
  
  Double_t A = TMath::Power(n/TMath::Abs(alpha),n)*TMath::Exp(-0.5*alpha*alpha);
  Double_t B = n/TMath::Abs(alpha) - TMath::Abs(alpha);
  
  Double_t y = ( TMath::Abs(sigma) > 1E-12 ? (x-mean)/sigma : 0 );
  
  if ( y > alpha*-1.0 ) 
  {
    return N*TMath::Exp(-0.5*y*y);
  }
  else 
  {
    return N*A*TMath::Power(B-y,-n);
  }
}

Double_t funcJpsiGCBE(Double_t* xx, Double_t* par)
{
  Double_t x = xx[0];
  
  Double_t g = par[0]*TMath::Gaus(x,par[1],par[2]);
  
  Double_t jpsi = funcCB(xx,par+3);
  
  Double_t expo = par[8]*TMath::Exp(par[9]*x);
  
  return g+expo+jpsi;
}

Double_t funcJpsiJpsiPrimeCustom(Double_t* xx, Double_t* par)
{ 
  Double_t N = par[0];
  Double_t alpha = par[1];
  Double_t n = par[2];
  Double_t mean = par[3];
  Double_t sigma = par[4];
  Double_t alphaprime = par[5];
  Double_t nprime = par[6];
  
  Double_t x = xx[0];
  
  Double_t A = TMath::Power(n/TMath::Abs(alpha),n)*TMath::Exp(-0.5*alpha*alpha);
  Double_t B = n/TMath::Abs(alpha) - TMath::Abs(alpha);
  Double_t C = TMath::Power(nprime/TMath::Abs(alphaprime),nprime)*TMath::Exp(-0.5*alphaprime*alphaprime);
  Double_t D = nprime/TMath::Abs(alphaprime) - TMath::Abs(alphaprime);
  
  Double_t y = ( TMath::Abs(sigma) > 1E-12 ? (x-mean)/sigma : 0 );
  
  Double_t cb(0);
  
  if ( y > alphaprime )
  {
    cb = N*C*TMath::Power(D+y,-nprime);
  }
  else if ( y > alpha*-1.0 ) 
  {
    cb = N*TMath::Exp(-0.5*y*y);
  }
  else 
  {
    cb = N*A*TMath::Power(B-y,-n);
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
  Double_t N = par[0];
  Double_t alpha = par[1];
  Double_t n = par[2];
  Double_t mean = par[3];
  Double_t sigma = par[4];
  Double_t alphaprime = par[5];
  Double_t nprime = par[6];
  
  Double_t x = xx[0];
  
  Double_t A = TMath::Power(n/TMath::Abs(alpha),n)*TMath::Exp(-0.5*alpha*alpha);
  Double_t B = n/TMath::Abs(alpha) - TMath::Abs(alpha);
  Double_t C = TMath::Power(nprime/TMath::Abs(alphaprime),nprime)*TMath::Exp(-0.5*alphaprime*alphaprime);
  Double_t D = nprime/TMath::Abs(alphaprime) - TMath::Abs(alphaprime);
  
  Double_t y = ( TMath::Abs(sigma) > 1E-12 ? (x-mean)/sigma : 0 );
  
  if ( y > alphaprime )
  {
    return N*C*TMath::Power(D+y,-nprime);
  }
  else if ( y > alpha*-1.0 ) 
  {
    return N*TMath::Exp(-0.5*y*y);
  }
  else 
  {
    return N*A*TMath::Power(B-y,-n);
  }
}

Double_t funcJpsiJpsiPrime(Double_t* xx, Double_t* par)
{
  Double_t jpsi = funcCB(xx,par);
  Double_t jpsiprime = funcCB2(xx,par+5);
  
  int n = 10;
  Double_t x = xx[0];
    
  Double_t e1 = par[n]*TMath::Exp(par[n+1]*x);
  Double_t e2 = par[n+2]*TMath::Exp(par[n+3]*x);    
  
  Double_t e = e1;
  
  if ( x > par[3] ) e=e2;
  
  return jpsi+jpsiprime+e;
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

const char* NormalizeName(const char* name, const char* suffix)
{
  TString str(Form("%s_%s",name,suffix));
  
  str.ReplaceAll("-","_");
  str.ReplaceAll("/","%");
  
  return str.Data();
}

//_____________________________________________________________________________
AliAnalysisMuMu::Result::~Result()
{
  delete fHC;
  delete fMap;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Result::FitJpsiJpsiPrimeCustom(TH1& h)
{
  std::cout << "Fit with jpsi + jpsiprime' (custom)" << std::endl;
  
  const Double_t xmin(1.5);
  const Double_t xmax(8.0);
  
  fFitTotal = new TF1("fFitTotal",funcJpsiJpsiPrimeCustom,xmin,xmax,14);
  fFitTotal->SetLineColor(4);
  
  fFitTotal->SetParName(0,"cstecb");
  fFitTotal->SetParName(1,"alpha");
  fFitTotal->SetParName(2,"n");
  fFitTotal->SetParName(3,"meanjpsi");
  fFitTotal->SetParName(4,"sigmajpsi");
  fFitTotal->SetParName(5,"alphaprime");
  fFitTotal->SetParName(6,"nprime");
  fFitTotal->SetParName(7,"cstepol1");
  fFitTotal->SetParName(8,"slopepol1");
  fFitTotal->SetParName(9,"cstegaus");
  fFitTotal->SetParName(10,"meanpsiprime");
  fFitTotal->SetParName(11,"sigmapsiprime");
  fFitTotal->SetParName(12,"csteexpo");
  fFitTotal->SetParName(13,"slopeexpo");
  
  fFitTotal->SetParameter( 0,1);
    
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
  
  fFitTotal->SetParameter(0,rcb->Parameter(0));
  fFitTotal->SetParameter(1,rcb->Parameter(1)); fFitTotal->SetParLimits(1,0,10); // alpha
  fFitTotal->SetParameter(2,rcb->Parameter(2)); fFitTotal->SetParLimits(2,1,10); // n
  fFitTotal->SetParameter(3,rcb->Parameter(3)); fFitTotal->SetParLimits(3,3.0,3.5); // mean
  fFitTotal->SetParameter(4,rcb->Parameter(4)); fFitTotal->SetParLimits(4,0,1); // sigma
  fFitTotal->SetParameter(5,rcb->Parameter(5)); fFitTotal->SetParLimits(1,0,10); // alphaprime
  fFitTotal->SetParameter(6,rcb->Parameter(6)); fFitTotal->SetParLimits(2,1,10); // nprime

  fFitTotal->FixParameter(1,alphaMC);
  fFitTotal->FixParameter(2,nMC);
  fFitTotal->FixParameter(5,alphaprimeMC);
  fFitTotal->FixParameter(6,nprimeMC);
  
  TF1* fge = new TF1("fge","gaus(0)+expo(3)",3.5,4.4);
  fge->SetParameters(1,3.6,0.25,1,1);
  TFitResultPtr rpsiprime = h.Fit(fge,fitOption,"",3.5,4.4);
  
  if (static_cast<int>(rpsiprime))
  {
    AliInfo("Will fix psiprime parameters");
    fFitTotal->FixParameter(9,0);
    fFitTotal->FixParameter(10,3.7);
    fFitTotal->FixParameter(11,0.1);
  }
  else
  {
    fFitTotal->SetParameter(10,rpsiprime->Parameter(1)); fFitTotal->SetParLimits(10,3.5,3.8); // mean'
    fFitTotal->SetParameter(11,rpsiprime->Parameter(2)); fFitTotal->SetParLimits(11,0.05,0.7); // sigma'
  }
  
  TFitResultPtr rpol1 = h.Fit("pol1",fitOption,"",1.5,2.5);
  fFitTotal->SetParameter( 7,rpol1->Parameter(0));
  fFitTotal->SetParameter( 8,rpol1->Parameter(1));
  
  TFitResultPtr rexpo = h.Fit("expo",fitOption,"",4.5,7.0);
  fFitTotal->SetParameter(12,rexpo->Parameter(0));
  fFitTotal->SetParameter(13,rexpo->Parameter(1));
  
  
  TFitResultPtr r = h.Fit(fFitTotal,fitOption,"",1.5,7);
  
  TF1* signal = new TF1("signal","gaus",2,6);  
  signal->SetParameters(fFitTotal->GetParameter(0),
                        fFitTotal->GetParameter(3),
                        fFitTotal->GetParameter(4));

  TF1* signalPrime = new TF1("signalPrime","gaus",2,6);  
  signalPrime->SetParameters(fFitTotal->GetParameter(9),
                             fFitTotal->GetParameter(10),
                             fFitTotal->GetParameter(11));
  
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
  
  gausParameters[0] = fFitTotal->GetParameter(0);
  gausParameters[1] = fFitTotal->GetParameter(3);
  gausParameters[2] = fFitTotal->GetParameter(4);

  gausParametersPrime[0] = fFitTotal->GetParameter(9);
  gausParametersPrime[1] = fFitTotal->GetParameter(10);
  gausParametersPrime[2] = fFitTotal->GetParameter(11);
  
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
  Set("MeanJpsi",fFitTotal->GetParameter(3),fFitTotal->GetParError(3));
  Set("SigmaJpsi",fFitTotal->GetParameter(4),fFitTotal->GetParError(4));

  double nprime = signalPrime->Integral(2,6)/h.GetBinWidth(10);
  double nerrprime = signalPrime->IntegralError(2,6,&gausParametersPrime[0],&covarianceMatrixPrime[0][0])/h.GetBinWidth(10);
  Set("NofJpsiPrime",nprime,nerrprime);      
  Set("MeanJpsiPrime",fFitTotal->GetParameter(10),fFitTotal->GetParError(10));
  Set("SigmaJpsiPrime",fFitTotal->GetParameter(11),fFitTotal->GetParError(11));
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Result::FitJpsiJpsiPrimeCB(TH1& h)
{
  std::cout << "Fit with jpsi + jpsiprime' (CB) " << std::endl;

  const Double_t xmin(1.5);
  const Double_t xmax(8.0);
  
  fFitTotal = new TF1("fFitTotal",funcJpsiJpsiPrime,xmin,xmax,14);

//  Double_t N = par[0];
//  Double_t alpha = par[1];
//  Double_t n = par[2];
//  Double_t mean = par[3];
//  Double_t sigma = par[4];
  
  fFitTotal->SetParameter( 0,1); // N
  fFitTotal->FixParameter( 1,0.936); // alpha
  fFitTotal->FixParameter( 2,4.44); // n
  fFitTotal->SetParameter( 3,3.1); fFitTotal->SetParLimits(3,3.0,3.2); // mean
  fFitTotal->SetParameter( 4,0.07); fFitTotal->SetParLimits(4,0.02,1); // sigma

  fFitTotal->SetParameter( 5,0.01); // N'
  fFitTotal->FixParameter( 6,0.936); // alpha'
  fFitTotal->FixParameter( 7,4.44); // n'
  fFitTotal->SetParameter( 8,3.7); fFitTotal->SetParLimits(8,3.5,3.8); // mean'
  fFitTotal->SetParameter( 9,0.1); fFitTotal->SetParLimits(9,0.02,1.0); // sigma'
  
  fFitTotal->SetParameter(10,h.GetMaximum());
  fFitTotal->SetParameter(11,-1);

  fFitTotal->SetParameter(12,h.GetMaximum()/100);
  fFitTotal->SetParameter(13,-1);

  TFitResultPtr r = h.Fit(fFitTotal,"SQBI","",1.5,6);
  
//  for ( int ix = 0; ix < fFitTotal->GetNpar(); ++ix )
//  {
//    for ( int iy = 0; iy < fFitTotal->GetNpar(); ++iy )
//    {      
//      std::cout << Form("COV(%d,%d)=%e ",ix,iy,r->GetCovarianceMatrix()(ix,iy));        
//    }
//    std::cout << std::endl;
//  }
  
  
  TF1* signal = new TF1("signal","gaus",2,8);
  
  signal->SetParameters(fFitTotal->GetParameter(0),
                        fFitTotal->GetParameter(3),
                        fFitTotal->GetParameter(4));

  TF1* signalPrime = new TF1("signalPrime","gaus",2,8);
  
  signalPrime->SetParameters(fFitTotal->GetParameter(0),
                             fFitTotal->GetParameter(8),
                             fFitTotal->GetParameter(9));
  
  Double_t gausParameters[3];
  Double_t gausParametersPrime[3];
  Double_t covarianceMatrix[3][3];
  Double_t covarianceMatrixPrime[3][3];
  
  gausParameters[0] = fFitTotal->GetParameter(0);
  gausParameters[1] = fFitTotal->GetParameter(3);
  gausParameters[2] = fFitTotal->GetParameter(4);

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

  gausParametersPrime[0] = fFitTotal->GetParameter(0);
  gausParametersPrime[1] = fFitTotal->GetParameter(8);
  gausParametersPrime[2] = fFitTotal->GetParameter(9);

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
  Set("MeanJpsi",fFitTotal->GetParameter(3),fFitTotal->GetParError(3));
  Set("SigmaJpsi",fFitTotal->GetParameter(4),fFitTotal->GetParError(4));

  double nprime = signalPrime->Integral(2,6)/h.GetBinWidth(10);
  double nerrprime = signalPrime->IntegralError(2,6,&gausParametersPrime[0],&covarianceMatrixPrime[0][0])/h.GetBinWidth(10);
  
  Set("NofJpsiPrime",nprime,nerrprime);
  Set("MeanJpsiPrime",fFitTotal->GetParameter(8),fFitTotal->GetParError(8));
  Set("SigmaJpsiPrime",fFitTotal->GetParameter(9),fFitTotal->GetParError(9));
  
}
  
//_____________________________________________________________________________
void AliAnalysisMuMu::Result::FitJpsiGCBE(TH1& h)
{
  std::cout << "Fit with jpsi alone (gaus + CB + expo)" << std::endl;
  
  const Double_t xmin(1.0);
  const Double_t xmax(8.0);
  
  fFitTotal = new TF1("fFitTotal",funcJpsiGCBE,xmin,xmax,10);
  fFitTotal->SetParNames("cste","x0","sigma0","N","alpha","n","mean","sigma","expocste","exposlope");
  
  fFitTotal->SetParLimits(3,0,h.GetMaximum()*2); // N
  
  const Double_t cbalpha(0.98);
  const Double_t cbn(5.2);
  
  fFitTotal->FixParameter(4,cbalpha);
  fFitTotal->FixParameter(5,cbn);
  
  fFitTotal->SetParLimits(6,2.8,3.2); // mean
  fFitTotal->SetParLimits(7,0.02,0.3); // sigma
  
  TF1* fg = new TF1("fg","gaus",xmin,xmax);
  
  h.Fit(fg,"","",0.75,3.0);
  
  fFitTotal->SetParameter(0,fg->GetParameter(0));
  fFitTotal->SetParameter(1,fg->GetParameter(1));
  fFitTotal->SetParameter(2,fg->GetParameter(2));
  
  TF1* fexpo = new TF1("expo","expo",xmin,xmax);
  
  h.Fit(fexpo,"","",3.5,5);
  
  fFitTotal->SetParameter(8,fexpo->GetParameter(0));
  fFitTotal->SetParameter(9,fexpo->GetParameter(1));
  
  fFitTotal->SetParameter(3,h.GetMaximum()),
  fFitTotal->SetParameter(4,cbalpha);
  fFitTotal->SetParameter(5,cbn);
  fFitTotal->SetParameter(6,3.15);
  fFitTotal->SetParameter(7,0.1);
  
  const char* fitOption = "SI+";
  
  TFitResultPtr r = h.Fit(fFitTotal,fitOption,"",2,5);
  
  Set("MeanJpsi",fFitTotal->GetParameter(6),fFitTotal->GetParError(6));
  Set("SigmaJpsi",fFitTotal->GetParameter(7),fFitTotal->GetParError(7));
  
  double m = MeanJpsi();
  double s = SigmaJpsi();
  double n = 3.0;
  
  TF1* fcb = new TF1("fcb",funcCB,xmin,xmax,5);
  fcb->SetParameters(fFitTotal->GetParameter(3),
                     fFitTotal->GetParameter(4),
                     fFitTotal->GetParameter(5),
                     fFitTotal->GetParameter(6),
                     fFitTotal->GetParameter(7));
  
  fcb->SetLineColor(6);
  fcb->SetNpx(100);
  TLine* l1 = new TLine(m-n*s,0,m-n*s,fFitTotal->GetParameter(3));
  TLine* l2 = new TLine(m+n*s,0,m+n*s,fFitTotal->GetParameter(3));
  l1->SetLineColor(6);
  l2->SetLineColor(6);
  h.GetListOfFunctions()->Add(fcb);
  h.GetListOfFunctions()->Add(l1);
  h.GetListOfFunctions()->Add(l2);
  
  
  Double_t cbParameters[5];
  Double_t covarianceMatrix[5][5];
  
  cbParameters[0] = fFitTotal->GetParameter(3);
  cbParameters[1] = fFitTotal->GetParameter(4);
  cbParameters[2] = fFitTotal->GetParameter(5);
  cbParameters[3] = fFitTotal->GetParameter(6);
  cbParameters[4] = fFitTotal->GetParameter(7);
  
  for ( int iy = 0; iy < 5; ++iy )
  {
    for ( int ix = 0; ix < 5; ++ix )
    {
      covarianceMatrix[ix][iy] = (r->GetCovarianceMatrix())(ix+3,iy+3);
    }
  }
  
  double njpsi = fcb->Integral(m-n*s,m+n*s)/h.GetBinWidth(1);
  
  double nerr = fcb->IntegralError(m-n*s,m+n*s,&cbParameters[0],&covarianceMatrix[0][0])/h.GetBinWidth(1);
  
  Set("NofJpsi",njpsi,nerr);
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Result::FitJpsiPCBE(TH1& h)
{
  std::cout << "Fit with jpsi alone (pol2 + CB + expo)" << std::endl;
  
  const Double_t xmin(2.0);
  const Double_t xmax(5.0);
  
  fFitTotal = new TF1("fFitTotal",funcJpsiPCBE,xmin,xmax,10);
  fFitTotal->SetParNames("p0","p1","p2","N","alpha","n","mean","sigma","expocste","exposlope");
  
  fFitTotal->SetParLimits(3,0,h.GetMaximum()*2); // N

  const Double_t cbalpha(0.98);
  const Double_t cbn(5.2);
  
  fFitTotal->FixParameter(4,cbalpha);
  fFitTotal->FixParameter(5,cbn);
  
  fFitTotal->SetParLimits(6,2,4); // mean
  fFitTotal->SetParLimits(7,0.05,0.2); // sigma
  
  TF1* fpol2 = new TF1("pol2","pol2",xmin,xmax);
                       
  h.Fit(fpol2,"+","",2,2.8);
  
  fFitTotal->SetParameter(0,fpol2->GetParameter(0));
  fFitTotal->SetParameter(1,fpol2->GetParameter(1));
  fFitTotal->SetParameter(2,fpol2->GetParameter(2));

  TF1* fexpo = new TF1("expo","expo",xmin,xmax);
  
  h.Fit(fexpo,"+","",3.5,4.5);
  
  fFitTotal->SetParameter(8,fexpo->GetParameter(0));
  fFitTotal->SetParameter(9,fexpo->GetParameter(1));
    
  fFitTotal->SetParameter(3,h.GetMaximum()),
  fFitTotal->SetParameter(4,cbalpha);
  fFitTotal->SetParameter(5,cbn);
  fFitTotal->SetParameter(6,3.15);
  fFitTotal->SetParameter(7,0.1);
  
  h.Fit(fFitTotal,"+","",2.5,5);
    
  Set("MeanJpsi",fFitTotal->GetParameter(6),fFitTotal->GetParError(6));
  Set("SigmaJpsi",fFitTotal->GetParameter(7),fFitTotal->GetParError(7));
  
  double m = MeanJpsi();
  double s = SigmaJpsi();
  double n = 2.0;
  
  TF1* fcb = new TF1("fcb",funcCB,xmin,xmax,5);
  fcb->SetParameters(fFitTotal->GetParameter(3),
                     fFitTotal->GetParameter(4),
                     fFitTotal->GetParameter(5),
                     fFitTotal->GetParameter(6),
                     fFitTotal->GetParameter(7));
  
  fcb->SetLineColor(6);
  fcb->SetNpx(100);
  TLine* l1 = new TLine(m-n*s,0,m-n*s,fFitTotal->GetParameter(3));
  TLine* l2 = new TLine(m+n*s,0,m+n*s,fFitTotal->GetParameter(3));
  l1->SetLineColor(6);
  l2->SetLineColor(6);
  h.GetListOfFunctions()->Add(fcb);
  h.GetListOfFunctions()->Add(l1);
  h.GetListOfFunctions()->Add(l2);
  
  
  Set("NofJpsi",fFitTotal->Integral(m-n*s,m+n*s)/h.GetBinWidth(1),fFitTotal->IntegralError(m-n*s,m+n*s)/h.GetBinWidth(1));
  
  //  Set("NofJpsi",fFitTotal->Integral(0,10)/h.GetBinWidth(1),fFitTotal->IntegralError(0,10)/h.GetBinWidth(1));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Result::FitJpsiCBE(TH1& h)
{
  std::cout << "Fit with jpsi alone" << std::endl;
  
  const Double_t xmin(1.5);
  const Double_t xmax(8.0);
  
  fFitTotal = new TF1("fFitTotal",funcJpsiCBE,xmin,xmax,7);
  fFitTotal->SetParNames("N","alpha","n","mean","sigma","expocste","exposlope");
  
//  fFitTotal->SetParameters(h.GetMaximum(),1,5,3.0,0.07,1.5,3,1,0);

  fFitTotal->SetParameters(1,1.15,3.6,3.0,0.07,1,-1);

  fFitTotal->SetParLimits(0,0,h.GetMaximum()); // N
//  fFitTotal->SetParLimits(1,0.1,2); // alpha
  fFitTotal->FixParameter(1,0.98);
//  fFitTotal->SetParLimits(2,0.01,5); // n
  fFitTotal->FixParameter(2,5.2);
  fFitTotal->SetParLimits(3,2.8,3.5); // mean
  fFitTotal->SetParLimits(4,0.05,0.2); // sigma
  
  TF1* fexpo = new TF1("expo","expo",xmin,xmax);
  
  h.Fit(fexpo,"+","",2,3);
  
  fFitTotal->SetParameter(5,fexpo->GetParameter(0));
  fFitTotal->SetParameter(6,fexpo->GetParameter(1));
  
  h.Fit(fFitTotal,"+","",2,5);
  
  
  Set("MeanJpsi",fFitTotal->GetParameter(3),fFitTotal->GetParError(3));
  Set("SigmaJpsi",fFitTotal->GetParameter(4),fFitTotal->GetParError(4));
  
  double m = MeanJpsi();
  double s = SigmaJpsi();
  double n = 2.0;
  
  TF1* fcb = new TF1("fcb",funcCB,xmin,xmax,5);
  fcb->SetParameters(fFitTotal->GetParameter(0),
                     fFitTotal->GetParameter(1),
                     fFitTotal->GetParameter(2),
                     fFitTotal->GetParameter(3),
                     fFitTotal->GetParameter(4));

  fcb->SetLineColor(6);
  fcb->SetNpx(1000);
  TLine* l1 = new TLine(m-n*s,0,m-n*s,fFitTotal->GetParameter(0));
  TLine* l2 = new TLine(m+n*s,0,m+n*s,fFitTotal->GetParameter(0));
  l1->SetLineColor(6);
  l2->SetLineColor(6);
  h.GetListOfFunctions()->Add(fcb);
  h.GetListOfFunctions()->Add(l1);
  h.GetListOfFunctions()->Add(l2);
  
  
  Set("NofJpsi",fFitTotal->Integral(m-n*s,m+n*s)/h.GetBinWidth(1),fFitTotal->IntegralError(m-n*s,m+n*s)/h.GetBinWidth(1));
  
  //  Set("NofJpsi",fFitTotal->Integral(0,10)/h.GetBinWidth(1),fFitTotal->IntegralError(0,10)/h.GetBinWidth(1));
  
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Result::FitJpsiECBE(TH1& h)
{
  std::cout << "Fit with jpsi alone (expo + CB + expo)" << std::endl;
  
  const Double_t xmin(1.5);
  const Double_t xmax(8.0);
  
  fFitTotal = new TF1("fFitTotal",funcJpsiECBE,xmin,xmax,9);
  fFitTotal->SetParNames("e0","s0","N","alpha","n","mean","sigma","e1","s1");
  
  fFitTotal->SetParameters(1,-1,1,1.15,3.6,3.2,0.06,-1);

  fFitTotal->SetParLimits(0,0,h.GetMaximum()*2);
  
  fFitTotal->FixParameter(3,0.98); // alpha
  fFitTotal->FixParameter(4,5.2); // n
  fFitTotal->SetParLimits(5,2.8,3.5); // mean
  fFitTotal->SetParLimits(6,0.05,0.2); // sigma
  
  TF1* fexpo1 = new TF1("expo1","expo",xmin,xmax);
  TF1* fexpo2 = new TF1("expo2","expo",xmin,xmax);
  
  h.Fit(fexpo1,"","",1.5,3);
  
  fFitTotal->SetParameter(0,fexpo1->GetParameter(0));
  fFitTotal->SetParameter(1,fexpo1->GetParameter(1));

  h.Fit(fexpo2,"","",3.5,5.0);

  fFitTotal->SetParameter(7,fexpo2->GetParameter(0));
  fFitTotal->SetParameter(8,fexpo2->GetParameter(1));

  const char* fitOption = "SI+";
  
  TFitResultPtr r = h.Fit(fFitTotal,fitOption,"",2,5);
  
  Set("MeanJpsi",fFitTotal->GetParameter(5),fFitTotal->GetParError(5));
  Set("SigmaJpsi",fFitTotal->GetParameter(6),fFitTotal->GetParError(6));
  
  double m = MeanJpsi();
  double s = SigmaJpsi();
  double n = 3.0;
  
  TF1* fcb = new TF1("fcb",funcCB,xmin,xmax,5);
  fcb->SetParameters(fFitTotal->GetParameter(2),
                     fFitTotal->GetParameter(3),
                     fFitTotal->GetParameter(4),
                     fFitTotal->GetParameter(5),
                     fFitTotal->GetParameter(6));

  fcb->SetParError(0,fFitTotal->GetParError(2));
  fcb->SetParError(1,fFitTotal->GetParError(3));
  fcb->SetParError(2,fFitTotal->GetParError(4));
  fcb->SetParError(3,fFitTotal->GetParError(5));
  fcb->SetParError(4,fFitTotal->GetParError(6));
  
  fcb->SetLineColor(6);
  fcb->SetNpx(1000);
  TLine* l1 = new TLine(m-n*s,0,m-n*s,fFitTotal->GetParameter(2));
  TLine* l2 = new TLine(m+n*s,0,m+n*s,fFitTotal->GetParameter(2));
  l1->SetLineColor(6);
  l2->SetLineColor(6);
  h.GetListOfFunctions()->Add(fcb);
  h.GetListOfFunctions()->Add(l1);
  h.GetListOfFunctions()->Add(l2);
  
  
  Double_t cbParameters[5];
  Double_t covarianceMatrix[5][5];
  
  cbParameters[0] = fFitTotal->GetParameter(2);
  cbParameters[1] = fFitTotal->GetParameter(3);
  cbParameters[2] = fFitTotal->GetParameter(4);
  cbParameters[3] = fFitTotal->GetParameter(5);
  cbParameters[4] = fFitTotal->GetParameter(6);
  
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

//_____________________________________________________________________________
void AliAnalysisMuMu::Result::FitJpsi(TH1& h)
{
  std::cout << "Fit with jpsi alone" << std::endl;

  const Double_t xmin(1.5);
  const Double_t xmax(8.0);

  fFitTotal = new TF1("fFitTotal",funcCB2,xmin,xmax,7);
  fFitTotal->SetParNames("N","alpha","n","mean","sigma","alphaprime","nprime");
  fFitTotal->SetParameters(h.GetMaximum(),1,5,3.0,0.07,1.5,3);
  fFitTotal->SetParLimits(0,0,h.GetMaximum()*2); // N
  fFitTotal->SetParLimits(1,0,10); // alpha
  fFitTotal->SetParLimits(2,1,10); // n
  fFitTotal->SetParLimits(3,1,4); // mean
  fFitTotal->SetParLimits(4,0.01,1); // sigma
  fFitTotal->SetParLimits(5,0,10); // alpha
  fFitTotal->SetParLimits(6,1,10); // n
  
  h.Fit(fFitTotal,"+","",2,5);
  
  
  Set("MeanJpsi",fFitTotal->GetParameter(3),fFitTotal->GetParError(3));
  Set("SigmaJpsi",fFitTotal->GetParameter(4),fFitTotal->GetParError(4));

  double m = MeanJpsi();
  double s = SigmaJpsi();
  double n = 2.0;
  
  Set("NofJpsi",fFitTotal->Integral(m-n*s,m+n*s)/h.GetBinWidth(1),fFitTotal->IntegralError(m-n*s,m+n*s)/h.GetBinWidth(1));

//  Set("NofJpsi",fFitTotal->Integral(0,10)/h.GetBinWidth(1),fFitTotal->IntegralError(0,10)/h.GetBinWidth(1));

}

//_____________________________________________________________________________
void AliAnalysisMuMu::Result::FitUpsilon(TH1& h)
{
  std::cout << "Fit with upsilon alone" << std::endl;
  
  const Double_t xmin(6.0);
  const Double_t xmax(12.0);
  
  fFitTotal = new TF1("fFitTotal",funcCB2,xmin,xmax,7);
  fFitTotal->SetParNames("N","alpha","n","mean","sigma","alphaprime","nprime");
  fFitTotal->SetParameters(h.GetMaximum(),1,5,9.46,0.2,1.5,3);
  fFitTotal->SetParLimits(0,0,h.GetMaximum()*2); // N
  fFitTotal->SetParLimits(1,0,10); // alpha
  fFitTotal->SetParLimits(2,1,10); // n
  fFitTotal->SetParLimits(3,8,12); // mean
  fFitTotal->SetParLimits(4,0.01,1); // sigma
  fFitTotal->SetParLimits(5,0,10); // alpha
  fFitTotal->SetParLimits(6,1,10); // n
  
  h.Fit(fFitTotal,"+","",6,12);
  
  
  Set("MeanUpsilon",fFitTotal->GetParameter(3),fFitTotal->GetParError(3));
  Set("SigmaUpsilon",fFitTotal->GetParameter(4),fFitTotal->GetParError(4));
  
  double m = MeanUpsilon();
  double s = SigmaUpsilon();
  double n = 2.0;
  
  Set("NofUpsilon",fFitTotal->Integral(m-n*s,m+n*s)/h.GetBinWidth(1),fFitTotal->IntegralError(m-n*s,m+n*s)/h.GetBinWidth(1));  
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Result::Fit(Int_t nrebin)
{
  //  new TCanvas;
  
  static int n(0);
  
  fRebin = nrebin;
  
//  TH2* h2 = static_cast<TH2*>(fHC->Histo("MinvUSPt"));  
//  int binptmin = h2->GetXaxis()->FindBin(2.0);
//  TH1* hminv = static_cast<TH1*>(h2->ProjectionY("minv0",binptmin,h2->GetXaxis()->GetNbins()));
  
  TH1* hminv = fHC->Histo("MinvUSPt:py");
  
  if (!hminv || hminv->GetEntries()<100 ) return;
    
  fMinv = static_cast<TH1*>(hminv->Clone(Form("minv%d",n)));
  fMinv->Rebin(fRebin);
  fMinv->SetDirectory(0);
//  fMinv->SetStats(false);
  fMinv->SetTitle(TriggerClass());
                  
  if ( ( fFitType & kJpsi ) &&
      ( fFitType & kJpsiPrime ) )
  {
    FitJpsiJpsiPrimeCustom(*fMinv);
//    FitJpsiJpsiPrimeCB(*fMinv);
  }
  else if ( ( fFitType & kJpsi ) && ( fFitType & kPbPb2011 ) )
  {
    if ( fFitType & kMatchAny )
    {
      FitJpsiECBE(*fMinv); // for MATCH
    }
    else
    {
      //    FitJpsiPCBE(*fMinv);
      FitJpsiGCBE(*fMinv); // for MATCHLOW      
    }
    fMinv->GetXaxis()->SetRangeUser(2.0,5.0);
    fMinv->SetMarkerStyle(kFullDotMedium);
    fMinv->SetMarkerColor(kBlue);
    fMinv->SetLineColor(kBlue);
    fMinv->SetDrawOption("EP");
  }
  else
  {
    int bminjpsi = fMinv->GetXaxis()->FindBin(2);
    int bmaxjpsi = fMinv->GetXaxis()->FindBin(4);

    int bminupsilon = fMinv->GetXaxis()->FindBin(8);
    int bmaxupsilon = fMinv->GetXaxis()->FindBin(12);

    if ( fMinv->Integral(bminjpsi,bmaxjpsi) > 100.0 ) FitJpsi(*fMinv);
    if ( fMinv->Integral(bminupsilon,bmaxupsilon) > 100.0 ) FitUpsilon(*fMinv);
  }
  
  TH2* hpp = static_cast<TH2*>(fHC->Histo("MinvPPPt"));
  TH2* hmm = static_cast<TH2*>(fHC->Histo("MinvMMPt"));
  
  if ( hpp && hmm )
  {
    TH2* htmp = static_cast<TH2*>(hpp->Clone("htmp"));
    htmp->SetDirectory(0);
    htmp->Add(hmm);
    htmp->Scale(0.5);
  
    fMinvLS = htmp->ProjectionY();
    fMinvLS->SetDirectory(0);
    fMinvLS->Rebin(fRebin);
      
    delete htmp;
  }
  else
  {
    AliDebug(1,Form("hpp=%p hmm=%p",hpp,hmm)); // this might happen for small runs and
    // strict cuts on pairs
  }
  
//  Print();
}

//_____________________________________________________________________________
TMap* AliAnalysisMuMu::Result::Map() const
{
  if (!fMap)
  {
    fMap = new TMap;
    fMap->SetOwnerKeyValue(kTRUE,kTRUE);
  }
  return fMap;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Result::Set(const char* name, double value, double error)
{
  TObjArray* p = static_cast<TObjArray*>(Map()->GetValue(name));
  if (!p) 
  {
    TParameter<double>* v = new TParameter<double>(name,value);
    TParameter<double>* e = new TParameter<double>(name,error);
    p = new TObjArray;
    p->SetOwner(kTRUE);
    p->Add(v);
    p->Add(e);
    Map()->Add(new TObjString(name),p);
  }
  else
  {
    TParameter<double>* v = static_cast<TParameter<double>*>(p->At(0));
    TParameter<double>* e = static_cast<TParameter<double>*>(p->At(1));
    v->SetVal(value);
    e->SetVal(error);    
  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::Result::HasValue(const char* name) const
{
  return ( Map()->GetValue(name) != 0x0 );
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMu::Result::GetValue(const char* name) const
{
  TObjArray* p = static_cast<TObjArray*>(Map()->GetValue(name));
  if (p)
  {
    TParameter<double>* val = static_cast<TParameter<double>*>(p->At(0));
    return val->GetVal();
  }
  return 0.0;
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMu::Result::GetError(const char* name) const
{
  TObjArray* p = static_cast<TObjArray*>(Map()->GetValue(name));
  if (p)
  {
    TParameter<double>* val = static_cast<TParameter<double>*>(p->At(1));
    return val->GetVal();
  }
  return 0.0;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Result::Print(Option_t* /*opt*/) const
{
  std::cout << Form("%20s - %s %s %s - NRUNS %d - NTRIGGER %10d",
               GetName(),
               EventSelection(),
                    PairSelection(),
                    CentralitySelection(),
                    NofRuns(),
               NofTriggers()) << std::endl;
  
  if (NofJpsi())
  {
    std::cout << Form("\t\tNjpsi %7.2f +- %5.2f \n\t\tRatio (x10^4) %7.2f +- %5.2f CHI2/NDF %5.2f"
                 "\n\t\tPEAK %7.2f +- %5.2f Gev/c^2 SIGMA %7.2f +- %5.2f MeV/c^2",
                 NofJpsi(),ErrorOnNofJpsi(),
                 NofJpsi()*1E4/fNofTriggers,
                 ErrorOnNofJpsi()*1E4/fNofTriggers,
                 (fFitTotal ? fFitTotal->GetChisquare()/fFitTotal->GetNDF() : 0),
                 MeanJpsi(),ErrorOnMeanJpsi(),
                 1E3*SigmaJpsi(),1E3*ErrorOnSigmaJpsi()) << std::endl;

  }

  if (NofJpsiPrime())
  {
    std::cout << Form("\t\tNjpsiPrime %7.2f +- %5.2f \n\t\tRatio (x10^4) %7.2f +- %5.2f CHI2/NDF %5.2f"
                 "\n\t\tPEAK %7.2f +- %5.2f Gev/c^2 SIGMA %7.2f +- %5.2f MeV/c^2",
                 NofJpsiPrime(),ErrorOnNofJpsiPrime(),
                 NofJpsiPrime()*1E4/fNofTriggers,
                 ErrorOnNofJpsiPrime()*1E4/fNofTriggers,
                 (fFitTotal ? fFitTotal->GetChisquare()/fFitTotal->GetNDF() : 0),
                 MeanJpsiPrime(),ErrorOnMeanJpsiPrime(),
                 1E3*SigmaJpsiPrime(),1E3*ErrorOnSigmaJpsiPrime()) << std::endl;
    
  }
  
  if ( NofUpsilon() )
  {
    std::cout << Form("\t\tNupsilon %7.2f +- %5.2f \n\t\tRatio (x10^4) %7.2f +- %5.2f CHI2/NDF %5.2f"
                 "\n\t\tPEAK %7.2f +- %5.2f Gev/c^2 SIGMA %7.2f +- %5.2f MeV/c^2",
                 NofUpsilon(),ErrorOnNofUpsilon(),
                 NofUpsilon()*1E4/fNofTriggers,
                 ErrorOnNofUpsilon()*1E4/fNofTriggers,
                 (fFitTotal ? fFitTotal->GetChisquare()/fFitTotal->GetNDF() : 0),
                 MeanUpsilon(),ErrorOnMeanUpsilon(),
                 1E3*SigmaUpsilon(),1E3*ErrorOnSigmaUpsilon()) << std::endl;
  }
}


//_____________________________________________________________________________
TString FindTrigger(const AliHistogramCollection& hc,
                    const char* base,
                    const char* selection,
                    const char* paircut,
                    const char* centrality)
{
  /// find the trigger containing the MinvPt histograms
  
  std::vector<std::string> trigger2test;
  
  //  trigger2test.push_back(Form("%s5-B-NOPF-ALLNOTRD",base));
  //  trigger2test.push_back(Form("%s1-B-NOPF-ALLNOTRD",base));
  //  trigger2test.push_back(Form("%s1B-ABCE-NOPF-MUON",base));
  if ( TString(base).Contains("||") || TString(base).Contains("-") )
  {
    trigger2test.push_back(base);
  }
  else
  {
    trigger2test.push_back(Form("%s-B-NOPF-ALLNOTRD",base));
    trigger2test.push_back(Form("%s-B-NOPF-MUON",base));
    trigger2test.push_back(Form("%s-S-NOPF-ALLNOTRD",base));
    trigger2test.push_back(Form("%s-S-NOPF-MUON",base));
  }
  trigger2test.push_back("ANY");
  
  for ( std::vector<std::string>::size_type i = 0; i < trigger2test.size(); ++i ) 
  {
    std::string trigger = trigger2test[i];
    
    if ( hc.Histo(selection,trigger.c_str(),centrality,paircut,"MinvUSPt") )
    {
      return trigger.c_str();
    }
  }
  
//  AliWarningGeneral("FindTrigger",Form("DID NOT FIND TRIGGER base=%s selection=%s paircut=%s centrality=%s",
//                  base,selection,paircut,centrality));
//  for ( std::vector<std::string>::size_type i = 0; i < trigger2test.size(); ++i ) 
//  {
//    AliWarningGeneral("FindTrigger",Form("tested trigger = %s",trigger2test[i].c_str()));
//  }
  return "";
}

//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________

TString AliAnalysisMuMu::fgOCDBPath("raw://");

TString AliAnalysisMuMu::fgDefaultDimuonTriggers("CMUL7-S-NOPF-ALLNOTRD,CMUL7-S-NOPF-MUON,CMUL8-S-NOPF-MUON,CMUL7-B-NOPF-ALLNOTRD,CMUL7-B-NOPF-MUON,CMUU7-B-NOPF-ALLNOTRD,CMUU7-B-NOPF-MUON,CPBI1MUL-B-NOPF-MUON");

TString AliAnalysisMuMu::fgDefaultMuonTriggers("CMSL7-S-NOPF-MUON,CMSL7-S-NOPF-ALLNOTRD,CMSL8-S-NOPF-MUON,CMSL8-S-NOPF-ALLNOTRD,CMSL7-B-NOPF-MUON,CMUS1-B-NOPF-MUON,CMUS7-B-NOPF-MUON");

TString AliAnalysisMuMu::fgDefaultMinbiasTriggers("CINT7-B-NOPF-ALLNOTRD,CINT7-S-NOPF-ALLNOTRD,CINT8-B-NOPF-ALLNOTRD,CINT8-S-NOPF-ALLNOTRD,CINT1-B-NOPF-ALLNOTRD,CPBI2_B1-B-NOPF-ALLNOTRD");

TString AliAnalysisMuMu::fgDefaultEventSelectionList("ALL");

TString AliAnalysisMuMu::fgDefaultPairSelectionList("pMATCHLOWRABSDCABOTH");

//_____________________________________________________________________________
AliAnalysisMuMu::AliAnalysisMuMu(const char* filename) : TObject(),
fFilename(filename),
fHistogramCollection(0x0),
fCounterCollection(0x0),
fDimuonTriggers(fgDefaultDimuonTriggers),
fMuonTriggers(fgDefaultMuonTriggers),
fMinbiasTriggers(fgDefaultMinbiasTriggers),
fEventSelectionList(fgDefaultEventSelectionList),
fPairSelectionList(fgDefaultPairSelectionList)
{
  // ctor
  
  GetCollections(fFilename,fHistogramCollection,fCounterCollection);
}

//_____________________________________________________________________________
AliAnalysisMuMu::~AliAnalysisMuMu()
{
  // dtor
  gROOT->CloseFiles();
  delete fCounterCollection;
  delete fHistogramCollection;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::CentralityCheck(const char* filelist)
{
  // Check if we get correctly filled centrality
  
  TObjArray* files = ReadFileList(filelist);
  
  if (!files || files->IsEmpty() ) return;
  
  TIter next(files);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    AliHistogramCollection* hc(0x0);
    AliCounterCollection* cc(0x0);
    
    if (!GetCollections(str->String().Data(),hc,cc)) continue;
    
    int run = RunNumberFromFileName(str->String().Data());
    
    TH1* h = hc->Histo("/ALL/CPBI1MUL-B-NOPF-MUON/Centrality");
 
    float percent(0);
    
    if (h)
    {
      percent = 100*h->Integral(1,1)/h->Integral();
    }
    
    std::cout << Form("RUN %09d PERCENT %7.2f",run,percent) << std::endl;
    
    delete hc;
  }
  
  gROOT->CloseFiles();
  
}

//_____________________________________________________________________________
void AliAnalysisMuMu::BasicCounts(Bool_t detailTriggers,
                                  ULong64_t* totalNmb,
                                  ULong64_t* totalNmsl,
                                  ULong64_t* totalNmul)
{
  // Report of some basic numbers, like number of MB and MUON triggers, 
  // both before and after physics selection, and comparison with 
  // the total number of such triggers (taken from the OCDB scalers)
  // if requested.
  //
  // filename is assumed to be a root filecontaining a list containing
  //    an AliCounterCollection (or directly an AliCounterCollection)
  //
  // if detailTriggers is kTRUE, each kind of (MB,MUL,MSL) is counted separately
  //
  
  if (!fHistogramCollection || !fCounterCollection) return;
  
  TObjArray* runs = fCounterCollection->GetKeyWords("run").Tokenize(",");
  TIter nextRun(runs);

  TObjArray* triggers = fCounterCollection->GetKeyWords("trigger").Tokenize(",");
  TIter nextTrigger(triggers);

  TObjArray* events = fCounterCollection->GetKeyWords("event").Tokenize(",");

  Bool_t doPS = (events->FindObject("PSALL") != 0x0);
  
  TObjString* srun;
  TObjString* strigger;

  ULong64_t localNmb(0);
  ULong64_t localNmsl(0);
  ULong64_t localNmul(0);
  
  if ( totalNmb) *totalNmb = 0;
  if ( totalNmsl) *totalNmsl = 0;
  if ( totalNmul ) *totalNmul = 0;

  while ( ( srun = static_cast<TObjString*>(nextRun()) ) )
  {
    std::cout << Form("RUN %09d ",srun->String().Atoi());
    
    TString details;
    ULong64_t nmb(0);
    ULong64_t nmsl(0);
    ULong64_t nmul(0);
    
    nextTrigger.Reset();
    
    Int_t nofPS(0);
    
    while ( ( strigger = static_cast<TObjString*>(nextTrigger()) ) )
    {
      
      if ( !fgDefaultMinbiasTriggers.Contains(strigger->String().Data()) &&
           !fgDefaultMuonTriggers.Contains(strigger->String().Data()) &&
           !fgDefaultDimuonTriggers.Contains(strigger->String().Data()) ) continue;
          
      ULong64_t n = fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d",
                                                    strigger->String().Data(),"ALL",srun->String().Atoi()));

      details += TString::Format("\n%50s %10lld",strigger->String().Data(),n);
      

      ULong64_t nps = fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d",
                                                      strigger->String().Data(),"PSALL",srun->String().Atoi()));

      if ( doPS )
      {
        details += TString::Format(" PS %5.1f %%",nps*100.0/n);
      }

      if (nps)
      {
        ++nofPS;
      }
      
      if ( fMinbiasTriggers.Contains(strigger->String()) )
      {
        nmb += n;
        if ( totalNmb) (*totalNmb) += n;
        localNmb += n;
      }
      else if ( fMuonTriggers.Contains(strigger->String()) )
      {
        nmsl += n;
        if ( totalNmsl) (*totalNmsl) += n;
        localNmsl += n;
      }
      else if ( fDimuonTriggers.Contains(strigger->String()) )
      {
        nmul += n;
        if ( totalNmul ) (*totalNmul) += n;
        localNmul += n;
      }      
    }
    
    std::cout << Form("MB %10lld MSL %10lld MUL %10lld %s",
                 nmb,nmsl,nmul,(nofPS == 0 ? "(NO PS AVAIL)": ""));
    
    if ( detailTriggers )
    {
      std::cout << details.Data();
    }
    std::cout << std::endl;
  }

  if ( !totalNmul && !totalNmsl && !totalNmb )
  {
    std::cout << std::endl << Form("%13s MB %10lld MSL %10lld MUL %10lld ","TOTAL",
                                   localNmb,localNmsl,localNmul) << std::endl;
  }

  delete runs;
  delete triggers;
  delete events;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::BasicCountsEvolution(const char* filelist, Bool_t detailTriggers)
{
  // Report of some basic numbers, like number of MB and MUON triggers,
  // both before and after physics selection, and comparison with
  // the total number of such triggers (taken from the OCDB scalers)
  // if requested.
  //
  // if detailTriggers is kTRUE, each kind of (MB,MUL,MSL) is counted separately
  //
  // To change the list of (single muon, dimuon, MB) triggers, use
  // the SetDefault*TriggerList methods prior to call this one
  //
  
  TObjArray* files = ReadFileList(filelist);
  
  if (!files || files->IsEmpty() ) return;
  
  TIter next(files);
  TObjString* str;
  
  ULong64_t totalNmb(0);
  ULong64_t totalNmsl(0);
  ULong64_t totalNmul(0);

  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    AliAnalysisMuMu m(str->String().Data());
    
    ULong64_t nmb(0);
    ULong64_t nmsl(0);
    ULong64_t nmul(0);
    
    m.BasicCounts(detailTriggers,&nmb,&nmsl,&nmul);
    
    totalNmb += nmb;
    totalNmsl += nmsl;
    totalNmul += nmul;
  }
  
  std::cout << std::endl << Form("%13s MB %10lld MSL %10lld MUL %10lld ","TOTAL",
                                 totalNmb,totalNmsl,totalNmul) << std::endl;

}

//_____________________________________________________________________________
TObjArray* AliAnalysisMuMu::CompareJpsiPerCMUUWithBackground(const char* jpsiresults,
                                                                   const char* backgroundresults)
{
  TFile* fjpsi = FileOpen(jpsiresults);
  TFile* fbck = FileOpen(backgroundresults);
  
  if (!fjpsi || !fbck) return 0x0;
  
  TGraph* gjpsi = static_cast<TGraph*>(fjpsi->Get("jpsipercmuu"));
    
  std::vector<std::string> checks;

  checks.push_back("muminus-CMUU7-B-NOPF-ALLNOTRD");
  checks.push_back("muplus-CMUU7-B-NOPF-ALLNOTRD");
  checks.push_back("muminus-CMUSH7-B-NOPF-MUON");
  checks.push_back("muplus-CMUSH7-B-NOPF-MUON");
  
  if (!gjpsi) return 0x0;

  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);
  
  for ( std::vector<std::string>::size_type j = 0; j < checks.size(); ++j )
  {
    
    TGraph* gback = static_cast<TGraph*>(fbck->Get(checks[j].c_str()));
    
    if (!gback) continue;

    if ( gjpsi->GetN() != gback->GetN() )
    {
      AliErrorClass("graphs have different number of points !");
      continue;
    }
    
    TGraphErrors* g = new TGraphErrors(gjpsi->GetN());
    
    for ( int i = 0; i < gjpsi->GetN(); ++i ) 
    {
      double r1,r2,y1,y2;
      
      gjpsi->GetPoint(i,r1,y1);
      gback->GetPoint(i,r2,y2);
      
      if ( r1 != r2 ) 
      {
        AliWarningClass(Form("run[%d]=%d vs %d",i,(int)r1,(int)r2));
        continue;
      }
      
      g->SetPoint(i,y2,y1);
      //    g->SetPointError(i,gjpsi->GetErrorY(i),gback->GetErrorY(i));
    }
    
    g->SetMarkerStyle(25+j);
    g->SetMarkerSize(1.2);
    if (j==0)
    {
      g->Draw("ap");
    }
    else
    {
      g->Draw("p");
    }
    g->SetLineColor(j+1);
    g->SetMarkerColor(j+1);
    g->SetName(checks[j].c_str());
    a->AddLast(g);
  }
  
  return a;
}

//_____________________________________________________________________________
TGraph* AliAnalysisMuMu::CompareJpsiPerCMUUWithSimu(const char* realjpsiresults,
                                                             const char* simjpsiresults)
{
  TFile* freal = FileOpen(realjpsiresults);
  TFile* fsim = FileOpen(simjpsiresults);
  
  if (!freal || !fsim) return 0x0;
  
  TGraph* greal = static_cast<TGraph*>(freal->Get("jpsipercmuu"));
  TGraph* gsim = static_cast<TGraph*>(fsim->Get("jpsipercmuu"));
  
  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);
  
  if ( greal->GetN() != gsim->GetN() )
  {
    AliErrorClass("graphs have different number of points !");
    return 0x0;
  }
    
  TGraphErrors* g = new TGraphErrors(greal->GetN());
  TGraphErrors* gratio = new TGraphErrors(greal->GetN());
    
  for ( int i = 0; i < greal->GetN(); ++i ) 
  {
    double r1,r2,y1,y2;
    
    greal->GetPoint(i,r1,y1);
    gsim->GetPoint(i,r2,y2);
    
    if ( r1 != r2 ) 
    {
      AliWarningClass(Form("run[%d]=%d vs %d",i,(int)r1,(int)r2));
      continue;
    }
    
    double ratio(0.0);
    
    if ( TMath::Abs(y1)<1E-6 || TMath::Abs(y2)<1E-6)
    {
      g->SetPoint(i,0,0);
      g->SetPointError(i,0,0);
    }
    else
    {    
      g->SetPoint(i,y2,y1);
      g->SetPointError(i,greal->GetErrorY(i),gsim ->GetErrorY(i));
      ratio = y2/y1;
    }
    gratio->SetPoint(i,r1,ratio);
  }
    
  g->SetMarkerStyle(25);
  g->SetMarkerSize(1.2);

  new TCanvas;
  
  g->Draw("ap");

  g->SetLineColor(1);
  g->SetMarkerColor(1);
  g->SetName("jpsipercmuurealvssim");

  new TCanvas;
  
  greal->Draw("alp");
  gsim->SetLineColor(4);
  
  gsim->Draw("lp");

  new TCanvas;
  gratio->Draw("alp");
  
  return g;
}

//_____________________________________________________________________________
TObjArray* AliAnalysisMuMu::ComputeBackgroundEvolution(const char* filelist, 
                                                       const char* triggerList, 
                                                       const char* outputFile,
                                                       const char* outputMode)
{
  // triggerList is a list of complete trigger names, separated by space
  // of the triggers to consider : only the first one found in the list 
  // is used for each run.
  
  TObjArray* files = ReadFileList(filelist);
  
  if (!files || files->IsEmpty() ) return 0x0;
  
  TIter next(files);
  TObjString* str;
  
  const char* ps = "ALL";
  const char* centrality = "PP";
  const char* ts1 = "sMATCHLOWRABS";
  const char* ts2 = "sMATCHLOWRABSDCA";
  
  std::map<std::string, std::vector<float> > runs;
  std::map<std::string, std::vector<float> > errruns;
  std::map<std::string, std::vector<float> > yplus,erryplus;
  std::map<std::string, std::vector<float> > yminus,erryminus;
  
  TObjArray* triggers = TString(triggerList).Tokenize(",");
  
  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);

  Bool_t bothSigns(kFALSE);
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    AliInfoClass(str->String().Data());
    
    AliHistogramCollection* hc(0x0);
    AliCounterCollection* cc(0x0);
    
    if (!GetCollections(str->String().Data(),hc,cc)) continue;
    
    TIter nextHisto(hc->CreateIterator());
    TH1* h;
    int nplus(0), nminus(0);
    
    while ( ( h = static_cast<TH1*>(nextHisto()) ) )
    {
      if ( TString(h->GetName()).EndsWith("Plus") )
      {
        nplus++;
      }
      if ( TString(h->GetName()).EndsWith("Minus") )
      {
        nminus++;
      }
    }
    
    if  (nminus==nplus && nplus>0 ) 
    {
      bothSigns = kTRUE;
    }

    AliInfoClass(Form("Both signs = %d",bothSigns));
    
    TIter nextTrigger(triggers);
    TObjString* trigger;
    TH1* h1p(0x0);
    TH1* h1m(0x0);
    TH1* h2p(0x0);
    TH1* h2m(0x0);
    TString format;
    
    while ( ( trigger = static_cast<TObjString*>(nextTrigger()) ) )
    {
      if  (bothSigns)
      {
        format = "/%s/%s/%s/%s/PtEtaMuPlus:py";
      }
      else
      {
        format = "/%s/%s/%s/%s/PtEtaMu:py";
      }
      
      h1p = hc->Histo(Form(format.Data(),ps,trigger->String().Data(),centrality,ts1));
      
      if (!h1p)
      {
        continue;
      }
      
      AliInfoClass(Form("Will use trigger %s",trigger->String().Data()));
      
      h2p = hc->Histo(Form(format.Data(),ps,trigger->String().Data(),centrality,ts2));
      
      if ( bothSigns )
      {
        format.ReplaceAll("Plus","Minus");
        h1m = hc->Histo(Form(format.Data(),ps,trigger->String().Data(),centrality,ts1));
        h2m = hc->Histo(Form(format.Data(),ps,trigger->String().Data(),centrality,ts2));
      }
      else
      {
        h2m=h2p;
        h1m=h1p;
      }
      
      if (h1m && h2m && h1p && h2p)
      {
        runs[trigger->String().Data()].push_back(RunNumberFromFileName(str->String().Data()));
        errruns[trigger->String().Data()].push_back(0.5);
        
        double value = 100-h2m->Integral()*100/h1m->Integral();
        yminus[trigger->String().Data()].push_back(value);
        double e1 = 1.0/TMath::Sqrt(h1m->GetEntries());
        double e2 = 1.0/TMath::Sqrt(h2m->GetEntries());
        erryminus[trigger->String().Data()].push_back(TMath::Sqrt(e1*e1+e2*e2)*value);
        
        value = 100-h2p->Integral()*100/h1p->Integral();
        yplus[trigger->String().Data()].push_back(value);
        e1 = 1.0/TMath::Sqrt(h1p->GetEntries());
        e2 = 1.0/TMath::Sqrt(h2p->GetEntries());
        erryplus[trigger->String().Data()].push_back(TMath::Sqrt(e1*e1+e2*e2)*value);
      }
      else
      {
        std::cout << Form("Error : h1m %p h2m %p h1p %p h2p %p",h1m,h2m,h1p,h2p) << std::endl;
      }
    }
    
    delete hc;
    delete cc;
    TFile* f = static_cast<TFile*>(gROOT->GetListOfFiles()->FindObject(str->String().Data()));    
    delete f;
  }
  
  delete triggers;
  
  TFile* f = new TFile(outputFile,outputMode);
  
  std::map<std::string, std::vector<float> >::const_iterator it;
  
  for ( it = runs.begin(); it != runs.end(); ++it )
  {
    std::string triggerName = it->first;
    
  TGraphErrors* gp = new TGraphErrors(runs[triggerName].size(),&runs[triggerName][0],&yplus[triggerName][0],&errruns[triggerName][0],&erryplus[triggerName][0]);
  TGraphErrors* gm(0x0);
  
  if ( bothSigns ) 
  {
    gm = new TGraphErrors(runs[triggerName].size(),&runs[triggerName][0],&yminus[triggerName][0],&errruns[triggerName][0],&erryminus[triggerName][0]);
  }
  
  if ( bothSigns ) 
  {
    gp->Write(Form("muplus_%s",triggerName.c_str()),TObject::kOverwrite);
    gm->Write(Form("muminus_%s",triggerName.c_str()),TObject::kOverwrite);
  }
  else
  {
    gp->Write(Form("mu_%s",triggerName.c_str()),TObject::kOverwrite);
  }
  
  }
  
  delete f;
  
  return a;
}

//_____________________________________________________________________________
TMap*
AliAnalysisMuMu::ComputeJpsiEvolution(const char* filelist, const char* triggerList,
                                      const char* outputFile, const char* outputMode,
                                      Bool_t simulation)
{
  /// Compute some jpsi information for a list of files / trigger combinations
  
  TObjArray* files = ReadFileList(filelist);
  
  if (!files || files->IsEmpty() ) return 0x0;
  
  TMap results; // one TObjString->TObjArray per file
  results.SetOwnerKeyValue(kTRUE,kTRUE);
  
  TIter nextFile(files);
  TObjString* str;
  UInt_t fitType(0);
  
  while ( ( str = static_cast<TObjString*>(nextFile()) ) )
  {
    std::cout << str->String().Data() << std::endl;
    
    AliAnalysisMuMu m(str->String().Data());
    
    m.SetDimuonTriggerList(triggerList);
    m.SetEventSelectionList("ALL");
    
    TObjArray* array = m.Jpsi(simulation); // the array will contain results for all triggers in fDimuonTriggers variable
    
    if (!array)
    {
      AliWarningClass(Form("Got no jpsi for %s",str->String().Data()));
    }
    else
    {
      Result* r = static_cast<Result*>(array->First());
      if (!r) continue;
      fitType = r->FitType();
    }
    
    results.Add(new TObjString(str->String()), array);
  }
  
  if (!results.GetSize()) return 0x0;
  
  // compute the total over all files
  
  TMap* total = new TMap;
  total->SetOwnerKeyValue(kTRUE,kTRUE);
  
  nextFile.Reset();
  TObjArray* triggers = TString(triggerList).Tokenize(",");
  
  TIter nextTrigger(triggers);
  TObjString* trigger(0x0);
  
  while ( ( trigger = static_cast<TObjString*>(nextTrigger())))
  {
    nextFile.Reset();
    
    Int_t nruns(0);
    Int_t n(0);
    TList l;    
    Result* ref(0x0);
    AliHistogramCollection* hc(0x0);
    
    while ( ( str = static_cast<TObjString*>(nextFile()) ) )
    {
      TObjArray* a = static_cast<TObjArray*>(results.GetValue(str->String().Data()));
      
      Result* r(0x0);
      
      if (a)
      {
        r = static_cast<Result*>(a->FindObject(trigger->String().Data()));
        
        if (r)
        {
          if (!ref) ref = r;

          if ( !hc )
          {
            AliHistogramCollection* htmp = r->HC();
            if (!htmp)
            {
              continue;
            }
            hc = static_cast<AliHistogramCollection*>(htmp->Clone(Form("hc%d",0)));
          }
          else
          {
            l.Add(r->HC());
          }
          
          n += r->NofTriggers();
          ++nruns;
        }
      }
    }
    
    hc->Merge(&l);
    
    if (!ref) continue;
    
    Result* sum = new Result(ref->TriggerClass(),ref->EventSelection(),
                             ref->PairSelection(),ref->CentralitySelection(),
                             n,hc,1,fitType);

    sum->SetNofRuns(nruns);
    
    total->Add(new TObjString(trigger->String().Data()),sum);
    
  }

  StdoutToAliInfoClass(total->Print(););
  
  TFile* fout = new TFile(outputFile,outputMode);
  
  results.Write("rbr",TObject::kSingleKey|TObject::kOverwrite);
  
  total->Write("total",TObject::kSingleKey|TObject::kOverwrite);
  
  delete fout;
  
  AliInfoClass(Form("%d files analyzed",files->GetEntries()));
  
  return total;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::DecodeFileName(const char* filename,
                                             TString& period,
                                             int& esdpass,
                                             int& aodtrain,
                                             int& runnumber)
{
  esdpass=aodtrain=runnumber=-1;
  period="";
  
  TString sfile(gSystem->BaseName(filename));
  
  if (!sfile.BeginsWith("LHC") && !sfile.BeginsWith("SIM") ) 
  {
    std::cerr << Form("filename %s does not start with LHC or SIM",filename) << std::endl;
    return kFALSE;
  }
  
  int year;
  char p;
  Bool_t ok(kFALSE);
  
  if ( sfile.BeginsWith("LHC") ) 
  {
    if ( sfile.Contains("pass") && sfile.Contains("AODMUON") )
    {
      int pass;
      sscanf(sfile.Data(),"LHC%2d%c_pass%d_AODMUON%03d_%09d",&year,&p,&pass,&aodtrain,&runnumber); 
      period = TString::Format("LHC%2d%c",year,p);
      ok=kTRUE;
    }
    else if ( sfile.Contains("pass") && sfile.Contains("_muon_") && sfile.Contains("AOD000") )
    {
      // LHC11c_pass2_muon_AOD000_000152087.saf.root
      sscanf(sfile.Data(),"LHC%2d%c_pass%d_muon_AOD000_%09d",&year,&p,&esdpass,&runnumber); 
      period = TString::Format("LHC%2d%c",year,p);
      ok=kTRUE;      
    }
    else if ( sfile.Contains("_muon_calo_") && sfile.Contains("AODMUON000") )
    {
      //      LHC12h_muon_calo_AODMUON000_000190112.saf.root
      sscanf(sfile.Data(),"LHC%2d%c_muon_calo_AODMUON000_%09d",&year,&p,&runnumber);
      period = TString::Format("LHC%2d%c",year,p);
      ok=kTRUE;
    }
    else if ( sfile.Contains("_muon_calo_") && sfile.Contains("AOD000") )
    {
      //      LHC12h_muon_calo_AOD000_000190112.saf.root
      sscanf(sfile.Data(),"LHC%2d%c_muon_calo_AOD000_%09d",&year,&p,&runnumber);
      period = TString::Format("LHC%2d%c",year,p);
      ok=kTRUE;
    }
    else if ( sfile.Contains("AODMUON" ) )
    {    
      sscanf(sfile.Data(),"LHC%2d%c_AODMUON%03d_%09d",&year,&p,&aodtrain,&runnumber); 
      period = TString::Format("LHC%2d%c",year,p);
      ok=kTRUE;
    }
    else if ( sfile.Contains("AOD000") ) 
    {
      sscanf(sfile.Data(),"LHC%2d%c_muons_AOD000_%09d",&year,&p,&runnumber); 
      period = TString::Format("LHC%2d%c",year,p);
      ok=kTRUE;
    }
    else if ( sfile.Contains("ESD_OUTER000"))
    {
      sscanf(sfile.Data(),"LHC%2d%c_cpass1_ESD_OUTER000_%09d",&year,&p,&runnumber);
      ok=kTRUE;
    }
  }
  else if ( sfile.BeginsWith("SIM_JPSI3" ) )
  {
    sscanf(sfile.Data(),"SIM_JPSI3_%09d",&runnumber);
    ok = kTRUE;
  }
  else if ( sfile.BeginsWith("SIM_UPSILON" ) )
  {
    sscanf(sfile.Data(),"SIM_UPSILON_%09d",&runnumber);
    ok = kTRUE;
  }
  else if ( sfile.BeginsWith("SIM_JPSI" ) )
  {
    sscanf(sfile.Data(),"SIM_JPSI_LHC%2d%c_%09d",&year,&p,&runnumber);
    period = TString::Format("LHC%2d%c",year,p);
    ok = kTRUE;
  }
  
  if (!ok)
  {
    std::cerr << Form("Can not decode %s",filename) << std::endl;
    return kFALSE;
  }
  
  return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisMuMu::DrawFill(Int_t run1, Int_t run2, double ymin, double ymax, const char* label)
{
  AliDebugClass(1,Form("RUN1 %09d RUN2 %09d YMIN %e YMAX %e %s",
                       run1,run2,ymin,ymax,label));
  TBox* b = new TBox(run1*1.0,ymin,run2*1.0,ymax);
  b->SetFillColor(5);
  b->Draw();
  TText* text = new TText((run1+run2)/2.0,ymax*0.85,label);
  text->SetTextSize(0.025);
  text->SetTextFont(42);
  text->SetTextAlign(23);
  text->SetTextAngle(45);
  text->Draw();
}

//_____________________________________________________________________________
TString 
AliAnalysisMuMu::ExpandPathName(const char* file)
{
  // An expand method that lives alien URL as they are
  TString sfile;
  
  if ( !sfile.BeginsWith("alien://") )
  {
    return gSystem->ExpandPathName(file);
  }
  else
  {
    if (!gGrid) TGrid::Connect("alien://");
    if (!gGrid) return "";    
  }
  
  return file;
}

//_____________________________________________________________________________
TFile* 
AliAnalysisMuMu::FileOpen(const char* file)
{
  // Open a file after expansion of its name
  
  return TFile::Open(ExpandPathName(file).Data());
}

//_____________________________________________________________________________
ULong64_t AliAnalysisMuMu::GetTriggerScalerCount(const char* triggerList, Int_t runNumber)
{
  // Get the expected (from OCDB scalers) trigger count
  
  AliAnalysisTriggerScalers ts(runNumber,fgOCDBPath.Data());
  
  TObjArray* triggers = TString(triggerList).Tokenize(",");
  TObjString* trigger;
  TIter next(triggers);
  ULong64_t n(0);
  
  while ( ( trigger = static_cast<TObjString*>(next()) ) )
  {
    AliAnalysisTriggerScalerItem* item = ts.GetTriggerScaler(runNumber,"L2A",trigger->String().Data());
    if (item)
    {
      n += item->Value();
    }
    delete item;
  }
  delete triggers;
  
  return n;
}

//_____________________________________________________________________________
UInt_t AliAnalysisMuMu::GetSum(AliCounterCollection& cc, const char* triggerList, const char* eventSelection, Int_t runNumber)
{
  TObjArray* ktrigger = cc.GetKeyWords("trigger").Tokenize(",");
  TObjArray* kevent = cc.GetKeyWords("event").Tokenize(",");
  TObjArray* a = TString(triggerList).Tokenize(" ");
  TIter next(a);
  TObjString* str;
  
  UInt_t n(0);
  
  TString sEventSelection(eventSelection);
  sEventSelection.ToUpper();
  
  if ( kevent->FindObject(sEventSelection.Data()) ) 
  {
    while ( ( str = static_cast<TObjString*>(next()) ) )
    {
      if ( ktrigger->FindObject(str->String().Data()) )
      {
        if ( runNumber < 0 ) 
        {
          n +=  static_cast<UInt_t>(cc.GetSum(Form("trigger:%s/event:%s",str->String().Data(),eventSelection)));              
        }
        else
        {
          n +=  static_cast<UInt_t>(cc.GetSum(Form("trigger:%s/event:%s/run:%d",str->String().Data(),eventSelection,runNumber)));                        
        }
      }
    }
  }
  
  delete a;
  delete ktrigger;
  delete kevent;
  return n;
}

//_____________________________________________________________________________
TObjArray* 
AliAnalysisMuMu::Jpsi(Bool_t simulation)
{
  // Fit the J/psi (and psiprime) peaks for the triggers in fDimuonTriggers list
  
  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);
  
  if ( simulation )
  {
    const char* selection = "ALL";
    
    Add(a,GetResult(*fHistogramCollection,*fCounterCollection,"ANY",selection,"pMATCHLOWRABSBOTH","PP",Result::kJpsi));
  }
  else
  {
    TObjArray* triggerArray = fDimuonTriggers.Tokenize(",");
    TObjArray* eventTypeArray = fEventSelectionList.Tokenize(",");
    TObjArray* pairCutArray = fPairSelectionList.Tokenize(",");
    
    TIter nextTrigger(triggerArray);
    TIter nextEventType(eventTypeArray);
    TIter nextPairCut(pairCutArray);
    
    TObjString* trigger;
    TObjString* eventType;
    TObjString* pairCut;
    
    while ( ( trigger = static_cast<TObjString*>(nextTrigger())) )
    {
      AliDebug(1,Form("TRIGGER %s",trigger->String().Data()));
      
      nextEventType.Reset();
      
      while ( ( eventType = static_cast<TObjString*>(nextEventType())) )
      {
        AliDebug(1,Form("EVENTTYPE %s",eventType->String().Data()));
        
        nextPairCut.Reset();
        
        while ( ( pairCut = static_cast<TObjString*>(nextPairCut())) )
        {
          AliDebug(1,Form("PAIRCUT %s",pairCut->String().Data()));
          Add(a,GetResult(*fHistogramCollection,*fCounterCollection,
                          trigger->String().Data(),
                          eventType->String().Data(),
                          pairCut->String().Data(),
                          "PP",
                          Result::kJpsi | Result::kJpsiPrime,2));
        }
      }
    }
    
    delete triggerArray;
    delete eventTypeArray;
    delete pairCutArray;
  }
  
  if ( a->GetLast() < 0 )
  {
    delete a;
    a = 0;
  }
  return a;
  
}

//_____________________________________________________________________________
void
AliAnalysisMuMu::PlotJpsiEvolution(const char* resultFile, const char* triggerList, Bool_t fillBoundaries,
                                   const char* efficiencyFile)
{
  if ( efficiencyFile && strlen(efficiencyFile) > 0 )
  {
    std::ifstream in(gSystem->ExpandPathName(efficiencyFile));
    if (!in.bad())
    {
      char line[1024];
      int run;
      float eff, error;
      while ( in.getline(line,1023,'\n') )
      {
        sscanf(line,"%d %f ± %f",&run,&eff,&error);
        AliInfoClass(Form("%09d %8.6f ± %8.6f",run,eff,error));
      }
    }
    
    return;
  }
  
  TFile* f = TFile::Open(gSystem->ExpandPathName(resultFile));
  
  std::map<int, std::pair<int,int> > fills;
  
  TMap* m = static_cast<TMap*>(f->Get("rbr"));

  TIter next(m);
  TObjString* str;
  
  TObjArray files;
  files.SetOwner(kTRUE);
  
  while ( ( str = static_cast<TObjString*>(next())) )
  {
    files.Add(new TObjString(str->String()));
  }
  
  files.Sort();
  
  std::map<std::string, std::vector<float> > x_jpsirate;
  std::map<std::string, std::vector<float> > y_jpsirate;
  std::map<std::string, std::vector<float> > xerr_jpsirate;
  std::map<std::string, std::vector<float> > yerr_jpsirate;
  
  TIter nextTrigger(TString(triggerList).Tokenize(","));
  TObjString* trigger(0x0);
  
  int runMin(100000000);
  int runMax(0);

  TIter nextFile(&files);

  while ( ( trigger = static_cast<TObjString*>(nextTrigger())))
  {
    TString triggerClass(trigger->String());
    
    nextFile.Reset();
    
    
    while ( ( str = static_cast<TObjString*>(nextFile())) )
    {
      TObjArray* a = static_cast<TObjArray*>(m->GetValue(str->String().Data()));
      if (!a) continue;
      Result* r = static_cast<Result*>(a->FindObject(triggerClass.Data()));
      if (!r) continue;

      TString period;
      int aodtrain,esdpass,runnumber;

      if ( DecodeFileName(str->String().Data(),period,esdpass,aodtrain,runnumber) )
      {
        runMin = TMath::Min(runMin,runnumber);
        runMax = TMath::Max(runMax,runnumber);
        
        x_jpsirate[triggerClass.Data()].push_back(runnumber);
        xerr_jpsirate[triggerClass.Data()].push_back(0.5);

        if ( fillBoundaries )
        {
          AliAnalysisTriggerScalers ts(runnumber,fgOCDBPath.Data());
          int fill = ts.GetFillNumberFromRunNumber(runnumber);
          
          if (fills.count(fill))
          {
            std::pair<int,int>& p = fills[fill];
            p.first = TMath::Min(runnumber,p.first);
            p.second = TMath::Max(runnumber,p.second);
          }
          else
          {
            fills[fill] = std::make_pair<int,int>(runnumber,runnumber);
          }
        }
        
        Double_t y(0.0);
        Double_t yerr(0.0);
        
        if ( TMath::Finite(r->SigmaJpsi()) && r->NofTriggers() > 10 )
        {
          y = 100*r->NofJpsi()/r->NofTriggers();
          
          if ( r->NofJpsi() > 0 )
          {
            yerr = y * TMath::Sqrt( (r->ErrorOnNofJpsi()*r->ErrorOnNofJpsi())/(r->NofJpsi()*r->NofJpsi()) + 1.0/r->NofTriggers());
          }
        }
        
        y_jpsirate[triggerClass.Data()].push_back(y);
        yerr_jpsirate[triggerClass.Data()].push_back(yerr);
      }
    }
  }

  delete f;
  
  TCanvas* c = new TCanvas("cJpsiRateEvolution","cJpsiRateEvolution");
  
  c->Draw();
  
  Double_t ymin(0);
  Double_t ymax(2);
  
  TH2* h = new TH2F("h","h;RunNumber;J/#psi per CMUL (%)",100,runMin,runMax,100,ymin,ymax);
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  h->GetXaxis()->SetNoExponent();
  
  h->Draw();

  if (fillBoundaries)
  {
    std::map<int, std::pair<int,int> >::const_iterator it;
    
    for ( it = fills.begin(); it != fills.end(); ++it )
    {
      const std::pair<int,int>& p = it->second;
      TString fillnumber;
      fillnumber.Form("%d",it->first);
      DrawFill(p.first,p.second,ymin,ymax,fillnumber.Data());
    }
  }

  h->Draw("sameaxis");
  
  //c->RedrawAxis("g");

  nextTrigger.Reset();
  
  int i(0);
  int color[] = { 1,2,4,5,6 };
  int marker[] = { 20,23,25,21,22 };
  
  while ( ( trigger = static_cast<TObjString*>(nextTrigger())))
  {
    std::vector<float>& x = x_jpsirate[trigger->String().Data()];
    std::vector<float>& y = y_jpsirate[trigger->String().Data()];
    std::vector<float>& xerr = xerr_jpsirate[trigger->String().Data()];
    std::vector<float>& yerr = yerr_jpsirate[trigger->String().Data()];
    
    TGraphErrors* g = new TGraphErrors(x.size(),&x[0],&y[0],&xerr[0],&yerr[0]);
    
    g->SetLineColor(color[i]);
    g->SetMarkerColor(color[i]);
    g->SetMarkerStyle(marker[i]);
    g->GetXaxis()->SetNoExponent();
    g->Draw("LP");
//    g->Print();

    Double_t m2 = g->GetMean(2);
    
    TLine* line = new TLine(runMin,m2,runMax,m2);
    line->SetLineColor(color[i]);
    line->Draw();
    
    AliInfoClass(Form("TRIGGER %s MEAN %7.2f",trigger->String().Data(),m2));
    ++i;
  }
  
  
}

//_____________________________________________________________________________
AliAnalysisMuMu::Result*
AliAnalysisMuMu::GetResult(const AliHistogramCollection& hc,
                           AliCounterCollection& cc,
                           const char* base,
                           const char* selection,
                           const char* paircut,
                           const char* centrality,
                           UInt_t fitType,
                           Int_t nrebin)
{
  Result* r(0x0);

  TString trigger = FindTrigger(hc,base,selection,paircut,centrality);
  
  if ( trigger == "" )
  {
    return 0;
  }
  
  Int_t ntrigger = (Int_t)cc.GetSum(Form("trigger:%s/event:%s",trigger.Data(),selection));

//  new TCanvas;
  
  r = new Result(trigger.Data(),
                 selection,
                 paircut,
                 centrality,
                 ntrigger,
                 hc.Project(selection,trigger,centrality,paircut),
                 nrebin,
                 fitType);
  
  return r;
}

//_____________________________________________________________________________
Bool_t
AliAnalysisMuMu::GetCollections(const char* rootfile,
                                      AliHistogramCollection*& hc,
                                      AliCounterCollection*& cc)
{
  hc = 0x0;
  cc = 0x0;
  
  TFile* f = static_cast<TFile*>(gROOT->GetListOfFiles()->FindObject(rootfile));
  
  if (!f)
  {
    f = TFile::Open(rootfile);    
  }
  
  if ( !f || !f->IsOpen() ) 
  {
    return kFALSE;    
  }
  
  TList* list = static_cast<TList*>(f->Get("chist"));
  
  if (!list) return kFALSE;
  
  hc = static_cast<AliHistogramCollection*>(list->At(0));
  cc = static_cast<AliCounterCollection*>(list->At(1));
  
  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::PlotBackgroundEvolution(const char* gfile, const char* triggerList)
{
  // plot the graphs found in the file (and generated using the ComputeBackgroundEvolution() method)
  
  TFile* f = TFile::Open(ExpandPathName(gfile).Data());    
  
  if ( !f || !f->IsOpen() ) 
  {
    return;
  }
  
  SetColorScheme();
  

  TCanvas* c = new TCanvas("background-evolution","background-evolution");
  
  c->Draw();
  
  TLegend* l = new TLegend(0.4,0.6,0.97,0.97);
  l->SetFillColor(0);
  l->SetTextColor(AliAnalysisMuMu::kBlue);
  l->SetLineColor(AliAnalysisMuMu::kBlue);
  
  TObjArray* triggers = TString(triggerList).Tokenize(",");
  
  gStyle->SetOptTitle(0);
  
  TObjString* str(0x0);
  TIter next(triggers);
  Int_t i(0);
  Int_t run1(99999999);
  Int_t run2(0);
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TGraph* g = static_cast<TGraph*>(f->Get(Form("mu_%s",str->String().Data())));
    if (!g) continue;
    run1 = TMath::Min(run1,TMath::Nint(g->GetX()[0]));
    run2 = TMath::Max(run2,TMath::Nint(g->GetX()[g->GetN()-1]));
  }
  
  AliInfoClass(Form("run1 %d run2 %d",run1,run2));
  
  TH2* hframe = new TH2F("hframe","hframe",run2-run1+1,run1,run2,100,0,100);
  hframe->Draw();
  hframe->GetXaxis()->SetNoExponent();
  hframe->GetYaxis()->SetTitle("Background percentage");

  next.Reset();
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TGraph* g = static_cast<TGraph*>(f->Get(Form("mu_%s",str->String().Data())));
    if (!g)
    {
      AliErrorClass(Form("Graph mu_%s not found",str->String().Data()));
      continue;
    }
    
    Int_t color(i+1);
    
    if (i==0) color = AliAnalysisMuMu::kBlue;
    if (i==1) color = AliAnalysisMuMu::kOrange;
    
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(20+i);
    
    g->Draw("LP");
    
    TLegendEntry* le = l->AddEntry(g,str->String().Data(),"lp");
    le->SetTextColor(color);
    
    g->GetYaxis()->SetTitleColor(AliAnalysisMuMu::kBlue);
    g->GetXaxis()->SetTitleColor(AliAnalysisMuMu::kBlue);
//    g->GetXaxis()->SetNoExponent();
//    g->GetYaxis()->SetRangeUser(0,100);
//    g->GetYaxis()->SetTitle("Background percentage");
    g->Print();
    
    ++i;
  }
  
  l->Draw();
  delete triggers;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::SetColorScheme()
{
  new TColor(AliAnalysisMuMu::kBlue,4/255.0,44/255.0,87/255.0,"my blue");
  new TColor(AliAnalysisMuMu::kOrange,255/255.0,83/255.0,8/255.0,"my orange");
  new TColor(AliAnalysisMuMu::kGreen,152/255.0,202/255.0,52/255.0,"my green");
  
  gStyle->SetGridColor(AliAnalysisMuMu::kBlue);
  
  gStyle->SetFrameLineColor(AliAnalysisMuMu::kBlue);
  gStyle->SetAxisColor(AliAnalysisMuMu::kBlue,"xyz");
  gStyle->SetLabelColor(AliAnalysisMuMu::kBlue,"xyz");
  
  gStyle->SetTitleColor(AliAnalysisMuMu::kBlue);
  gStyle->SetTitleTextColor(AliAnalysisMuMu::kBlue);
  gStyle->SetLabelColor(AliAnalysisMuMu::kBlue);
  gStyle->SetStatTextColor(AliAnalysisMuMu::kBlue);
  
  gStyle->SetOptStat(0);
}

//_____________________________________________________________________________
void 
AliAnalysisMuMu::SinglePtPlot(const char* rootfile)
{
  AliHistogramCollection* histogramCollection(0x0);
  AliCounterCollection* counterCollection(0x0);
  
  if (!GetCollections(rootfile,histogramCollection,counterCollection))
  {
    return;
  }

  TCanvas* c1 = new TCanvas("singlept","singlept");
  
  gStyle->SetTextSize(1.0);
  
  c1->SetLogy();
  
  TLegend* l = new TLegend(0.12,0.12,0.4,0.3,"","NDC");
  l->SetFillStyle(0);
  l->SetLineWidth(0);
  l->SetLineColor(0);
  l->SetTextColor(kBlack);
  l->SetTextFont(42);
  l->SetTextSize(0.035);
  
  const char* cuts[] = { "ALL","ETA","ETARABS","ETARABSMATCH","ETARABSMATCHDCA" };
  const char* cutnames[] = {
    "all",
    "+ -4 < #eta < -2.5",
    "+ 171^{#circ} < #theta_{abs} < 178^{#circ}",
    "+ trigger matching",
    "+ PxDCA"
  };
  const int colors[] = { 1,2,3,4,6 };

  for ( int i = 0; i < 5; ++i )
  {
    TH1* minus = histogramCollection->Histo(Form("/PS/CPBI1MSH-B-NOPF-MUON/CENT80/s%s/PtEtaMuMinus:py",cuts[i]));
    TH1* plus = histogramCollection->Histo(Form("/PS/CPBI1MSH-B-NOPF-MUON/CENT80/s%s/PtEtaMuPlus:py",cuts[i]));
    if (!minus || !plus)
    {
      AliErrorClass(Form("Form cannot get histos for cut %s",cuts[i]));
      continue;
    }
    TH1* h = static_cast<TH1*>(minus->Clone(Form("h%d",i)));
    h->Add(plus);
    h->SetDirectory(0);
    h->SetLineColor(colors[i]);
    h->SetMinimum(1);
    h->SetMarkerSize(0);
    h->SetLineWidth(2);
    h->SetStats(0);
    h->SetBit(TH1::kNoTitle);
    if (i==0)
    {
      h->Draw("e");
    }
    else
    {
      h->Draw("esame");
    }
    h->GetYaxis()->SetTitle("dN/dp_{t} (counts/0.5 GeV/c)");
    h->GetXaxis()->SetTitle("p_{t} (GeV/c)");
    
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetTitleSize(0.04);

    h->GetXaxis()->SetRangeUser(0,30);
    l->AddEntry(h,cutnames[i]);
  }
  
  l->Draw();
  
  TList extralines;
  extralines.SetOwner(kTRUE);
  
  extralines.Add(new TObjString("PbPb #sqrt{s_{NN}}=2.76 TeV"));
  extralines.Add(new TObjString("MUON high p_{t} trigger events"));
  extralines.Add(new TObjString("w/ phys. sel."));
  extralines.Add(new TObjString("w/ reco. vertex"));
  extralines.Add(new TObjString("centrality 0-80 %"));
  
  ALICEseal(0.5,0.93,&extralines);
  
  delete histogramCollection;
  delete counterCollection;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::TriggerCountCoverage(const char* triggerList, Bool_t compact)
{
  // Give the fraction of triggers (in triggerList) relative 
  // to what is expected in the scalers
  
  TGrid::Connect("alien://"); // to insure the "Trying to connect to server... message does not pollute our output later on...
  
  AliLog::EType_t oldLevel = static_cast<AliLog::EType_t>(AliLog::GetGlobalLogLevel());
  
  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  
  if (!fHistogramCollection || !fCounterCollection) return;
  
  TObjArray* runs = fCounterCollection->GetKeyWords("run").Tokenize(",");
  TIter nextRun(runs);
  
  TObjArray* triggers = fCounterCollection->GetKeyWords("trigger").Tokenize(",");
  TIter nextTrigger(triggers);
  
  TObjString* srun;
  TObjString* strigger;
  
  TString striggerList(triggerList);
  
  ULong64_t total(0);
  ULong64_t totalExpected(0);
  
  while ( ( srun = static_cast<TObjString*>(nextRun()) ) )
  {
    std::cout << Form("RUN %09d ",srun->String().Atoi());
    
    if (!compact) std::cout << std::endl;
    
    nextTrigger.Reset();
    
    while ( ( strigger = static_cast<TObjString*>(nextTrigger()) ) )
    {
      if ( !striggerList.Contains(strigger->String().Data()) ) 
      {
        continue;
      }
      ULong64_t n = fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d",
                                                    strigger->String().Data(),"ALL",srun->String().Atoi()));
   
      ULong64_t expected = GetTriggerScalerCount(strigger->String().Data(),srun->String().Atoi());
    
      
      total += n;
      totalExpected += expected;
      
      std::cout << Form("%30s %9lld expected %9lld ",strigger->String().Data(),n,expected);
      
      if ( expected > 0 ) {
        std::cout << Form("fraction %5.1f %%",n*100.0/expected);
      }
      if (!compact)
      {
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;
  }
  
  std::cout << Form("TOTAL %lld expected %lld fraction %5.1f %%",
               total,totalExpected,totalExpected ? total*100.0/totalExpected : 0.0) << std::endl;
  
  AliLog::SetGlobalLogLevel(oldLevel);
  delete triggers;
  delete runs;
}

//_____________________________________________________________________________
void AnalyisResultLocation(const char* runlist, const char* basedir, const char* what)
{
    std::ifstream in(runlist);
  int run;
  while ( in >> run )
  {
    std::cout << Form("%s/%09d/%s",basedir,run,what) << std::endl;
  }
}

//_____________________________________________________________________________
void plot(const char* file="results.root", const char* pdfname="toto.pdf")
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  
  TFile* f = TFile::Open(file);
  
  TObjArray* a = static_cast<TObjArray*>(f->Get("results"));
  
  if (!a) return;
  
  TCanvas* c = new TCanvas("jpsiresults","jpsiresults");
  
  c->Draw();
  
  TLegend* l = new TLegend(0.5,0.7,0.9,0.9);
  l->SetFillColor(0);
  TH1* h0(0x0);
  
  for ( int i = 0; i <= a->GetLast(); ++i ) 
  {
    AliAnalysisMuMu::Result* r = static_cast<AliAnalysisMuMu::Result*>(a->At(i));
    r->Print();
    
    TH1* h = static_cast<TH1*>(r->Minv());
    
    h->SetStats(kFALSE);
    h->SetXTitle("M_{#mu^{+}#mu^{-}} (Gev/c^{2})");
    h->SetLineColor(i+1);
    
    h->SetMaximum(5E3);
    h->SetMarkerStyle(0);
    
    if ( i == 0 )
    {
      h0 = h;
      h->Draw("hist");
    }
    else
    {
      h->Draw("histsame");
    }
    
    TObjArray* n = TString(r->TriggerClass()).Tokenize("-");
    TObjString* nn = static_cast<TObjString*>(n->First());
    l->AddEntry(h,Form("%6d %s",r->NofTriggers(),nn->String().Data()));
    
    delete n;
  }
  
  h0->Draw("histsame");
  l->Draw();
  c->SetLogy();
  
  delete f;
  
  c->SaveAs(pdfname);
}

//_____________________________________________________________________________
TObjArray*
AliAnalysisMuMu::ReadFileList(const char* filelist)
{
  // 
  // read the filelist and try to order it by runnumber
  // 
  // filelist can either be a real filelist (i.e. a text file containing
  // root filenames) or a root file itself.
  //
  
  char line[1024];
  
  TObjArray* files = new TObjArray;
  files->SetOwner(kTRUE);
  
  TString sfilelist(ExpandPathName(filelist));
  
  if ( sfilelist.EndsWith(".root") )
  {
    files->Add(new TObjString(sfilelist.Data()));
    return files;
  }
  
  std::set<int> runnumbers;
  std::map<int,std::string> filemap;
  
    std::ifstream in(sfilelist.Data());
  
  TString period;
  int aodtrain,esdpass,runnumber;
  
  while ( in.getline(line,1022,'\n') )
  {
    DecodeFileName(line,period,esdpass,aodtrain,runnumber);
    
    AliDebugClass(1,Form("line %s => period %s esdpass %d aodtrain %d runnumber %09d",
                         line,period.Data(),esdpass,aodtrain,runnumber));
    
    filemap.insert(std::make_pair<int,std::string>(runnumber,line));
    runnumbers.insert(runnumber);      
  }
  
  in.close();
  
  std::set<int>::const_iterator it;
  
  for ( it = runnumbers.begin(); it != runnumbers.end(); ++it ) 
  {
    files->Add(new TObjString(filemap[*it].c_str()));
  }
  
  return files;
}

//_____________________________________________________________________________
Int_t 
AliAnalysisMuMu::RunNumberFromFileName(const char* filename)
{
  TString period;
  int esdpass,aodtrain,runnumber;
  Bool_t ok = DecodeFileName(filename,period,esdpass,aodtrain,runnumber);
  if ( ok ) return runnumber;
  return -1;
}
