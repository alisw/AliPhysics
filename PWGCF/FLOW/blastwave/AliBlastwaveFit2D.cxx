#include<stdio.h>

#include"AliBlastwaveFit2D.h"

#include"TMath.h"
#include"TH1.h"
#include"TGraphErrors.h"

ClassImp(AliBlastwaveFit2D);


TF2 *AliBlastwaveFit2D::fgFuncIntYield = NULL;
TF2 *AliBlastwaveFit2D::fgFuncIntV2 = NULL;

AliBlastwaveFit2D::AliBlastwaveFit2D(const char *name,Double_t mass) :
  AliBlastwaveFit(name,mass)
{
  if(!fgFuncIntYield){
     fgFuncIntYield = new TF2("fFuncInt2Yield",AliBlastwaveFit2D::FunctionIntYield,0,2*TMath::Pi(),0,1,7);
     fgFuncIntYield->SetNpx(100);
     fgFuncIntYield->SetNpy(100);
  }
  if(!fgFuncIntV2){
     fgFuncIntV2 = new TF2("fFuncInt2V2",AliBlastwaveFit2D::FunctionIntV2,0,2*TMath::Pi(),0,1,7);
     fgFuncIntV2->SetNpx(100);
     fgFuncIntV2->SetNpy(100);
  }

  Initialize();
}
//------------------------------------------------------------------------------
AliBlastwaveFit2D::AliBlastwaveFit2D() :
  AliBlastwaveFit()
{  
  if(!fgFuncIntYield){
     fgFuncIntYield = new TF2("fFuncInt2Yield",AliBlastwaveFit2D::FunctionIntYield,0,2*TMath::Pi(),0,1,7);
     fgFuncIntYield->SetNpx(100);
     fgFuncIntYield->SetNpy(100);
  }
  if(!fgFuncIntV2){
     fgFuncIntV2 = new TF2("fFuncInt2V2",AliBlastwaveFit2D::FunctionIntV2,0,2*TMath::Pi(),0,1,7);
     fgFuncIntV2->SetNpx(100);
     fgFuncIntV2->SetNpy(100);
  }
  Initialize();
}
//------------------------------------------------------------------------------
AliBlastwaveFit2D::~AliBlastwaveFit2D(){
}
//------------------------------------------------------------------------------
Double_t AliBlastwaveFit2D::FunctionIntYield(Double_t x[],Double_t par[]){
  // x[0] = phi
  // x[1] = r/R
  // par[0] = pt
  // par[1] = T
  // par[2] = s2
  // par[3] = rh0_0
  // par[4] = rho_a
  // par[5] = R_power
  // par[6] = mass 

  Double_t mt = TMath::Sqrt(par[0]*par[0]+par[6]*par[6]);
  Double_t rho = (par[3]+par[4]*TMath::Cos(2*x[0]))*TMath::Power(x[1],par[5]);
  Double_t alfat = (par[0]/par[1])*TMath::SinH(rho);
  Double_t betat = (mt/par[1])*TMath::CosH(rho);
  Double_t results;
  if(betat < 200) results = TMath::BesselI(0,alfat)*TMath::BesselK(1,betat)*(1+2*par[2]*TMath::Cos(2*x[0])) *x[1];
  else results = 0.5*TMath::Exp(alfat-betat)/sqrt(alfat*betat)*(1+2*par[2]*TMath::Cos(2*x[0])) *x[1];

  return results*par[0]*mt;
}
//------------------------------------------------------------------------------
Double_t AliBlastwaveFit2D::FunctionIntV2(Double_t x[],Double_t par[]){
  // x[0] = phi
  // x[1] = r/R
  // par[0] = pt
  // par[1] = T_fo
  // par[2] = s2
  // par[3] = rh0_0
  // par[4] = rho_a
  // par[5] = R_power
  // par[6] = mass

  Double_t mt = TMath::Sqrt(par[0]*par[0]+par[6]*par[6]);
  Double_t rho = (par[3]+par[4]*TMath::Cos(2*x[0]))*TMath::Power(x[1],par[5]);
  Double_t alfat = (par[0]/par[1])*TMath::SinH(rho);
  Double_t betat = (mt/par[1])*TMath::CosH(rho);
  Double_t results;
  if(betat < 200) results = TMath::Cos(2*x[0])*TMath::BesselI(2,alfat)*TMath::BesselK(1,betat)*(1+2*par[2]*TMath::Cos(2*x[0])) *x[1];
  else results = TMath::Cos(2*x[0])*0.5*TMath::Exp(alfat-betat)/sqrt(alfat*betat)*(1+2*par[2]*TMath::Cos(2*x[0])) *x[1];

  return results*par[0]*mt;
}
//------------------------------------------------------------------------------
Int_t AliBlastwaveFit2D::SetParameter(Int_t ipar,Double_t val){
  Bool_t status = kTRUE;

  if(ipar < 0 || ipar > 4) return status;

  if(fFunctionYield){
     fFunctionYield->SetParameter(ipar,val);
     status = kFALSE;
  }
  if(fFunctionV2){
     fFunctionV2->SetParameter(ipar,val);
     status = kFALSE;
  }
  return status;
}
//------------------------------------------------------------------------------
Int_t AliBlastwaveFit2D::SetNormalization(){
  Bool_t status = kTRUE;

  TH1 *h = (TH1 *) GetSpectrumObjCopy();
  if(fFunctionYield && h){
      fFunctionYield->SetParameter(6,1.0);
      Double_t intH=0;
      Double_t intF=0;
      for (Int_t ibin = 1; ibin <= h->GetNbinsX(); ibin++) {
	  Double_t pt = h->GetBinCenter(ibin);
	  if(pt < GetMaxPt() && pt > GetMinPt() && TMath::Abs(h->GetBinContent(ibin)) > h->GetBinError(ibin)+0.0001){
	      Double_t width = h->GetBinWidth(ibin);
	      intH += h->GetBinContent(ibin)*width;
	      intF += EvalYield(pt)*width;
//  	      for(Double_t x = -width*0.45;x < width*0.5;x+=width*0.1)
//  		  intF += EvalYield(pt+x)*width*0.1;
	  }
      }

      Double_t norm = 1;
      if(intF > 0) norm = intH/intF;
      fFunctionYield->SetParameter(6,norm);
      status = kFALSE;
  }
  return status;
}
//------------------------------------------------------------------------------
void AliBlastwaveFit2D::SetMass(Double_t mass) {
  fMass=mass;
  if(fFunctionV2) fFunctionV2->SetParameter(5,fMass);
  if(fFunctionYield) fFunctionYield->SetParameter(5,fMass);   
}
//------------------------------------------------------------------------------
void AliBlastwaveFit2D::Initialize(){
  char nameFuncYield[100];
  char nameFuncV2[100];
  snprintf(nameFuncYield,100,"%sFuncYield",GetName());
  snprintf(nameFuncV2,100,"%sFuncV2",GetName());

  if(fFunctionV2) delete fFunctionV2;
  if(fFunctionYield) delete fFunctionYield;
  fFunctionV2 = new TF1(nameFuncV2,V2,0,6,6);
  fFunctionV2->SetParameter(5,fMass);  
  fFunctionYield = new TF1(nameFuncYield,Pt,0,6,7);
  fFunctionYield->SetParameter(5,fMass);
  fFunctionV2->SetNpx(200);
  fFunctionYield->SetNpx(200);

  for(Int_t i =0;i < 6;i++) fFunctionV2->SetParName(i,fgParName[i]);
  for(Int_t i =0;i < 7;i++) fFunctionYield->SetParName(i,fgParName[i]);
}
//------------------------------------------------------------------------------
Double_t AliBlastwaveFit2D::V2(Double_t x[],Double_t par[]){
  // par[0] = T_fo
  // par[1] = s2
  // par[2] = rh0_mean
  // par[3] = rho_v2
  // par[4] = R_power
  // par[5] = mass
  Double_t maxRho = par[2]*0.5*(par[4]+2)/(1+2*par[1]*par[3]);

  fgFuncIntYield->SetParameter(0,x[0]);
  fgFuncIntYield->SetParameter(1,par[0]);
  fgFuncIntYield->SetParameter(2,par[1]);
  fgFuncIntYield->SetParameter(3,maxRho);
  fgFuncIntYield->SetParameter(4,par[3]*maxRho*2);
  fgFuncIntYield->SetParameter(5,par[4]);
  fgFuncIntYield->SetParameter(6,par[5]);

  fgFuncIntV2->SetParameter(0,x[0]);
  fgFuncIntV2->SetParameter(1,par[0]);
  fgFuncIntV2->SetParameter(2,par[1]);
  fgFuncIntV2->SetParameter(3,maxRho);
  fgFuncIntV2->SetParameter(4,par[3]*maxRho*2);
  fgFuncIntV2->SetParameter(5,par[4]);
  fgFuncIntV2->SetParameter(6,par[5]);

  Double_t res = fgFuncIntV2->Integral(0.,2*TMath::Pi(),0.,1.);
  Double_t den = fgFuncIntYield->Integral(0.,2*TMath::Pi(),0.,1.);
  
  if(den == 0) return 0.0;

  return res/den;
}
//------------------------------------------------------------------------------
Double_t AliBlastwaveFit2D::Pt(Double_t x[],Double_t par[]){
  // par[0] = T_fo
  // par[1] = s2
  // par[2] = rh0_mean
  // par[3] = rho_v2
  // par[4] = R_power
  // par[5] = mass

  Double_t maxRho = par[2]*0.5*(par[4]+2)/(1+2*par[1]*par[3]);

  fgFuncIntYield->SetParameter(0,x[0]);
  fgFuncIntYield->SetParameter(1,par[0]);
  fgFuncIntYield->SetParameter(2,par[1]);
  fgFuncIntYield->SetParameter(3,maxRho);
  fgFuncIntYield->SetParameter(4,par[3]*maxRho*2);
  fgFuncIntYield->SetParameter(5,par[4]);
  fgFuncIntYield->SetParameter(6,par[5]);

  fgFuncIntV2->SetParameter(0,x[0]);
  fgFuncIntV2->SetParameter(1,par[0]);
  fgFuncIntV2->SetParameter(2,par[1]);
  fgFuncIntV2->SetParameter(3,maxRho);
  fgFuncIntV2->SetParameter(4,par[3]*maxRho*2);
  fgFuncIntV2->SetParameter(5,par[4]);
  fgFuncIntV2->SetParameter(6,par[5]);

  Double_t res = par[6]*1E+5*fgFuncIntYield->Integral(0.,2*TMath::Pi(),0.,1.);

  return res;
}
//------------------------------------------------------------------------------
const Float_t AliBlastwaveFit2D::GetMeanBeta(){

  TF2 fbeta("fbeta","TMath::TanH(([0]+[1]*TMath::Cos(2*x))*TMath::Power(y,[2]))*y*(1+2*[3]*cos(2*x))",0,2*TMath::Pi(),0,1);
  fbeta.SetNpx(1000);
  fbeta.SetNpy(1000);

  fbeta.SetParameter(0,fgFuncIntYield->GetParameter(3));
  fbeta.SetParameter(1,fgFuncIntYield->GetParameter(4));
  fbeta.SetParameter(2,fgFuncIntYield->GetParameter(5));
  fbeta.SetParameter(3,fgFuncIntYield->GetParameter(2));

  return fbeta.Integral(0.,2*TMath::Pi(),0.,1.)/TMath::Pi();
}
//------------------------------------------------------------------------------
const Float_t AliBlastwaveFit2D::GetMeanBeta(Double_t par[]){
  // par[0] = T_fo
  // par[1] = s2
  // par[2] = rh0_mean
  // par[3] = rho_v2
  // par[4] = R_power
  // par[5] = mass

  Double_t maxRho = par[2]*0.5*(par[4]+2)/(1+2*par[1]*par[3]);

  TF2 fbeta("fbeta","TMath::TanH(([0]+[1]*TMath::Cos(2*x))*TMath::Power(y,[2]))*y*(1+2*[3]*cos(2*x))",0,2*TMath::Pi(),0,1);
  fbeta.SetNpx(1000);
  fbeta.SetNpy(1000);

  fbeta.SetParameter(0,maxRho);
  fbeta.SetParameter(1,2*par[3]*maxRho);
  fbeta.SetParameter(2,par[4]);
  fbeta.SetParameter(3,par[1]);

  return fbeta.Integral(0.,2*TMath::Pi(),0.,1.)/TMath::Pi();
}
//------------------------------------------------------------------------------
void AliBlastwaveFit2D::SwitchOffFlow(TMinuit *m) const{
    m->FixParameter(kParS2);
    m->FixParameter(kParRho2);

    printf("AliBlastwaveFit2D: No flow histos -> Flow parameters put to zero\n");
}
