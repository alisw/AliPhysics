#include "AliHBTFits.h"
//_________________________________________________
///////////////////////////////////////////////////////////////////////////////////
//
// class AliHBTFits
//
// Sets of methods for fittig correlation functions
//
//
//
// 
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TH1.h>
#include <TH3.h>
#include <TF1.h>
#include <TF3.h>
#include <TError.h>
#include <TDirectory.h>
#include <TMath.h>

ClassImp(AliHBTFits)

TF1*     AliHBTFits::fgF1 = 0x0;
TF1*     AliHBTFits::fgF2 = 0x0;

/******************************************************************************************/
AliHBTFits::~AliHBTFits()
{
 //dtor
 delete fgF1;
 delete fgF2;
}
/******************************************************************************************/

/******************************************************************************************/
/********    Qout-Qside-Qlong  Cyl Surf   *************************************************/
/******************************************************************************************/

void AliHBTFits::FitQOutQSideQLongCylSurf(const TString& hname,Option_t* fopt, Float_t xmax, Float_t ymax,Float_t zmax)
{
  //Fits 3D histo with QOutQSideQLongCylSurf
 TH3* h = dynamic_cast<TH3*>(gDirectory->Get(hname));
 if (h == 0x0)
  {
    ::Error("FitQOutCylSurf","Can not find histogram %s",hname.Data());
    return;
  }
  
 delete gROOT->GetListOfFunctions()->FindObject("BorisTomasikCylSurfQOutQSideQLong");

 delete  fgF1;
 delete  fgF2;
 // here x is phi
 //par 0: R
 //par 1: Phi
 //par 2: qOut
 //par 3: qSide
 fgF1 = new TF1("ff1","sin([0]*([2]*cos(x)+[3]*sin(x))/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF2 = new TF1("ff2","cos([0]*([2]*cos(x)+[3]*sin(x))/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF1->SetNpx(1000);
 fgF2->SetNpx(1000);
 
 if (xmax <= 0.0) xmax = h->GetXaxis()->GetXmax();
 if (ymax <= 0.0) ymax = h->GetYaxis()->GetXmax();
 if (zmax <= 0.0) zmax = h->GetZaxis()->GetXmax();
 
 TF3 *theory = new TF3("BorisTomasikCylSurfQOutQSideQLong",&AliHBTFits::QOutQSideQLongCylSurf,0.0,xmax,0.0,ymax,0.0,zmax,5);//par :aLambda, aR, aJ,aPhi,
 theory->SetNpx(1000);
 
 theory->SetParameter(0,.75);
 theory->SetParameter(1,6.0);
 theory->SetParameter(2,1.0);
 theory->SetParameter(3,2.0);
 theory->SetParameter(4,0.9);
 
 theory->SetParName(0,"\\lambda");
 theory->SetParName(1,"R");
 theory->SetParName(2,"J");
 theory->SetParName(3,"L");
 theory->SetParName(4,"\\Phi");
 
 h->Fit(theory,fopt,"E");

}

/******************************************************************************************/
Double_t AliHBTFits::QOutQSideQLongCylSurf(Double_t *x, Double_t *par)
{
  //QOutQSideQLongCylSurf Function
  Double_t qout = x[0];
  Double_t qside = x[1];
  Double_t qlong = x[2];
  
  Double_t aLambda = par[0];
  Double_t aR = par[1];
  Double_t aJ = par[2];
  Double_t aL = par[3];
  Double_t aPhi = par[4];
  
  static Double_t oldlam = aLambda;
  static Double_t oldPhi = aPhi;

  if (oldPhi != aPhi)
   {
    ::Info("QOutQSideQLongCylSurf","out=%f, side=%f, long=%f\n lambda=%f, R=%f, J=%f, L=%f Phi=%f",
                                 qout,qside,qlong,aLambda,aR,aJ,aL,aPhi);
    oldPhi = aPhi;
   }
  
  if (oldlam != aLambda)
   {
    ::Info("QOutQSideQLongCylSurf","out=%f, side=%f, long=%f\n lambda=%f, R=%f, J=%f, L=%f Phi=%f",
                                 qout,qside,qlong,aLambda,aR,aJ,aL,aPhi);
    oldlam = aLambda;
   }
  
  Double_t result=aLambda*TMath::Exp(-(qout*qout+qside*qside+qlong*qlong)*aJ*aJ/0.0389273);
  
  if ( (qlong != 0.0) && (aL != 0.0) )
   {
     Double_t qlL = qlong*aL/2.0;
     Double_t sinqlL = TMath::Sin(qlL);
     Double_t sin2qlL = sinqlL*sinqlL;
     Double_t longpart = sin2qlL/(qlL*qlL);
     result*= longpart;
   }
  
  if (aPhi > 0.0)
   {
    fgF1->SetParameter(0,aR);
    fgF2->SetParameter(0,aR);
  
    fgF1->SetParameter(1,aPhi);
    fgF2->SetParameter(1,aPhi);
  
    fgF1->SetParameter(2,qout);
    fgF2->SetParameter(2,qout);

    fgF1->SetParameter(3,qside);
    fgF2->SetParameter(3,qside);
    
    Double_t sp = fgF1->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-8);
    Double_t cp = fgF2->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-8);
    
//    ::Info("QSideCylSurf","sp=%f, cp=%f",cp,sp);
    
    result *= cp*cp + sp*sp;
//    ::Info("QSideCylSurf","2: result=%f",result);

    Double_t tmp = TMath::BesselI0(1./aPhi);
    tmp = tmp*tmp*TMath::TwoPi()*TMath::TwoPi();
    if (tmp == 0.0)
     {
       ::Error("QOutCylSurf","Div by 0");
       return 1.0;
     }
    result=result/tmp;
   }
   
  result+=1.; 

//  ::Info("QSideCylSurf","Final: result=%f",result);
  return result;
;  
}
/******************************************************************************************/
/********    Qout  Cyl Surf   *************************************************************/
/******************************************************************************************/

void AliHBTFits::FitQOutCylSurf(const TString& hname, Option_t* fopt, Float_t max)
{
 TH1D* h = dynamic_cast<TH1D*>(gDirectory->Get(hname));
 if (h == 0x0)
  {
    ::Error("FitQOutCylSurf","Can not find histogram %s",hname.Data());
    return;
  }
  
 delete gROOT->GetListOfFunctions()->FindObject("BorisTomasikCylSurfQOut");

 delete  fgF1;
 delete  fgF2;
 fgF1 = new TF1("ff1","sin([0]*[2]*cos(x)/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF2 = new TF1("ff2","cos([0]*[2]*cos(x)/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF1->SetNpx(1000);
 fgF2->SetNpx(1000);
 
 if (max <= 0.0) max = h->GetXaxis()->GetXmax();
 TF1 *theory = new TF1("BorisTomasikCylSurfQOut",&AliHBTFits::QOutCylSurf,0.0,max,3);//par : aR, aPhi
 theory->SetNpx(1000);
 theory->SetParameter(0,6.0);
 theory->SetParameter(1,1.0);
 theory->SetParameter(2,0.75);
 theory->SetParName(0,"R_{out}");
 theory->SetParName(1,"\\Phi");
 theory->SetParName(2,"\\lambda");
 
 h->Fit(theory,fopt,"E");
 
}
/******************************************************************************************/

Double_t AliHBTFits::QOutCylSurf(Double_t *x, Double_t *par)
{
  //Qout 
  Double_t qout = x[0];
  Double_t aR = par[0];
  Double_t aPhi = par[1];
  Double_t aLambda = par[2];

//  ::Info("QOutCylSurf","q=%f, lambda=%f, Rout=%f, Phi=%f",qout,aLambda,aR,aPhi);

  Double_t result=aLambda*TMath::Exp(-qout*qout/0.0389273);
  
//  ::Info("QOutCylSurf","1: result=%f",result);
  
  if (aPhi > 0.0)
   {
    fgF1->SetParameter(0,aR);
    fgF2->SetParameter(0,aR);
  
    fgF1->SetParameter(1,aPhi);
    fgF2->SetParameter(1,aPhi);
  
    fgF1->SetParameter(2,qout);
    fgF2->SetParameter(2,qout);
    
    Double_t sp = fgF1->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-8);
    Double_t cp = fgF2->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-8);
    
//    ::Info("QOutCylSurf","sp=%f, cp=%f",cp,sp);
    
    result = result * (cp*cp + sp*sp);
//    ::Info("QOutCylSurf","2: result=%f",result);

    Double_t tmp = TMath::BesselI0(1./aPhi);
    tmp = tmp*tmp*TMath::TwoPi()*TMath::TwoPi();
    if (tmp == 0.0)
     {
       ::Error("QOutCylSurf","Div by 0");
       return 1.0;
     }
    result=result/tmp;
   }
   
  result+=1.; 

//  ::Info("QOutCylSurf","Final: result=%f",result);
  return result;
}
/******************************************************************************************/
/********    Qside  Cyl Surf   ************************************************************/
/******************************************************************************************/

void AliHBTFits::FitQSideCylSurf(const TString& hname, Option_t* fopt, Float_t max)
{
//Fits QSide According to Boris Tomasik formula
 TH1D* h = dynamic_cast<TH1D*>(gDirectory->Get(hname));
 if (h == 0x0)
  {
    ::Error("FitQSideCylSurf","Can not find histogram %s",hname.Data());
    return;
  }
 delete gROOT->GetListOfFunctions()->FindObject("BorisTomasikCylSurfQSide");
 delete fgF1;
 delete fgF2;
 fgF1 = new TF1("ff1","sin([0]*[2]*sin(x)/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF2 = new TF1("ff2","cos([0]*[2]*sin(x)/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF1->SetNpx(1000);
 fgF2->SetNpx(1000);
 
 if (max <= 0.0) max = h->GetXaxis()->GetXmax();
 TF1 *theory = new TF1("BorisTomasikCylSurfQSide",&AliHBTFits::QSideCylSurf,0.0,max,3);//par : aR, aPhi
 theory->SetNpx(1000);
 theory->SetParameter(0,6.0);
 theory->SetParameter(1,1.0);
 theory->SetParameter(2,0.75);
 theory->SetParName(0,"R_{out}");
 theory->SetParName(1,"\\Phi");
 theory->SetParName(2,"\\lambda");
 
 h->Fit(theory,fopt,"E");
 
}
/******************************************************************************************/

Double_t AliHBTFits::QSideCylSurf(Double_t *x, Double_t *par)
{
  //Qout 
  Double_t qside = x[0];
  Double_t aR = par[0];
  Double_t aPhi = par[1];
  Double_t aLambda = par[2];

//  ::Info("QSideCylSurf","q=%f, lambda=%f, Rside=%f, Phi=%f",qside,aLambda,aR,aPhi);

  Double_t result=aLambda*TMath::Exp(-qside*qside/0.0389273);
  if (aPhi > 0.0)
   {
    fgF1->SetParameter(0,aR);
    fgF2->SetParameter(0,aR);
  
    fgF1->SetParameter(1,aPhi);
    fgF2->SetParameter(1,aPhi);
  
    fgF1->SetParameter(2,qside);
    fgF2->SetParameter(2,qside);
    
    Double_t sp = fgF1->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-7);
    Double_t cp = fgF2->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-7);

    result *= cp*cp + sp*sp;

    Double_t tmp = TMath::BesselI0(1./aPhi);
    tmp = tmp*tmp*TMath::TwoPi()*TMath::TwoPi();
    if (tmp == 0.0)
     {
       ::Error("QSideCylSurf","Div by 0");
       return 1.0;
     }
    result/=tmp;
   }
  return result + 1.;
}
/******************************************************************************************/




