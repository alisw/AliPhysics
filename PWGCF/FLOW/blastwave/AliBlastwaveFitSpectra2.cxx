#include<stdio.h>

#include"AliBlastwaveFitSpectra2.h"

#include"TMath.h"

ClassImp(AliBlastwaveFitSpectra2);


TF1 *AliBlastwaveFitSpectra2::fgFuncIntYield = NULL;

AliBlastwaveFitSpectra2::AliBlastwaveFitSpectra2(const char *name,Double_t mass) :
  AliBlastwaveFit(name,mass)
{
  if(!fgFuncIntYield){
     fgFuncIntYield = new TF1("fFuncInt2Yield",AliBlastwaveFitSpectra2::FunctionIntYield,0,1,5);
     fgFuncIntYield->SetNpx(100);
  }

  Initialize();
}
//------------------------------------------------------------------------------
AliBlastwaveFitSpectra2::AliBlastwaveFitSpectra2() :
  AliBlastwaveFit("BlastwaveFitSpectra2",0.0)
{  
  if(!fgFuncIntYield){
     fgFuncIntYield = new TF1("fFuncInt2Yield",AliBlastwaveFitSpectra2::FunctionIntYield,0,1,5);
     fgFuncIntYield->SetNpx(100);
  }

  Initialize();
}
//------------------------------------------------------------------------------
AliBlastwaveFitSpectra2::~AliBlastwaveFitSpectra2(){
}
//------------------------------------------------------------------------------
Double_t AliBlastwaveFitSpectra2::FunctionIntYield(Double_t x[],Double_t par[]){
  // x[0] = r/R
  // par[0] = pt
  // par[1] = T
  // par[2] = beta_av
  // par[3] = R_power
  // par[4] = mass 

  Double_t maxrho = par[2]*0.5*(par[3]+2);

  Double_t mt = TMath::Sqrt(par[0]*par[0]+par[4]*par[4]);
  Double_t rho = maxrho*TMath::Power(x[0],par[3]);
  if(rho > 0.9999999999999999) rho = 0.9999999999999999;
  rho = TMath::ATanH(rho);
  Double_t alfat = (par[0]/par[1])*TMath::SinH(rho);
  Double_t betat = (mt/par[1])*TMath::CosH(rho);
  Double_t results;
  if(betat < 200) results = TMath::BesselI(0,alfat)*TMath::BesselK(1,betat) *x[0];
  else results = 0.5*TMath::Exp(alfat-betat)/sqrt(alfat*betat) *x[0];

  return results*par[0]*mt;
}
//------------------------------------------------------------------------------
Int_t AliBlastwaveFitSpectra2::SetParameter(Int_t ipar,Double_t val){
  Bool_t status = kTRUE;

  if(ipar < 0 || ipar > 4) return status;

  if(fFunctionYield){
     fFunctionYield->SetParameter(ipar,val);
     status = kFALSE;
  }
  return status;
}
//------------------------------------------------------------------------------
Int_t AliBlastwaveFitSpectra2::SetNormalization(){
  Bool_t status = kTRUE;

  TH1 *h = (TH1 *) GetSpectrumObjCopy();
  if(fFunctionYield && h){
      fFunctionYield->SetParameter(4,1.0);
      Double_t intH=0;
      Double_t intF=0;
      for (Int_t ibin = 1; ibin <= h->GetNbinsX(); ibin++) {
	  Double_t pt = h->GetBinCenter(ibin);
	  if(pt < GetMaxPt() && pt > GetMinPt() && TMath::Abs(h->GetBinContent(ibin)) > h->GetBinError(ibin)+0.0001){
	      Double_t width = h->GetBinWidth(ibin);
	      intH += h->GetBinContent(ibin)*width;
	      for(Double_t x = -width*0.45;x < width*0.5;x+=width*0.1)
		  intF += EvalYield(pt+x)*width*0.1;
	  }
      }

      Double_t norm = 1;
      if(intF > 0) norm = intH/intF;
      fFunctionYield->SetParameter(4,norm);
      status = kFALSE;
  }

  return status;
}
//------------------------------------------------------------------------------
void AliBlastwaveFitSpectra2::SetMass(Double_t mass) {
  fMass=mass;
  if(fFunctionYield) fFunctionYield->SetParameter(3,fMass);   
}
//------------------------------------------------------------------------------
void AliBlastwaveFitSpectra2::Initialize(){
  char nameFuncYield[100];
  snprintf(nameFuncYield,100,"%sFuncYield",GetName());

  if(fFunctionYield) delete fFunctionYield;
  fFunctionYield = new TF1(nameFuncYield,Pt,0,6,7);
  fFunctionYield->SetParameter(3,fMass);
  fFunctionYield->SetNpx(200);

  for(Int_t i =0;i < 5;i++) fFunctionYield->SetParName(i,fgParName[i]);
}
//------------------------------------------------------------------------------
Double_t AliBlastwaveFitSpectra2::Pt(Double_t x[],Double_t par[]){
  // par[0] = T_fo
  // par[1] = beta_av
  // par[2] = R_power
  // par[3] = mass

  fgFuncIntYield->SetParameter(0,x[0]);
  fgFuncIntYield->SetParameter(1,par[0]);
  fgFuncIntYield->SetParameter(2,par[1]);
  fgFuncIntYield->SetParameter(3,par[2]);
  fgFuncIntYield->SetParameter(4,par[3]);


  Double_t res = par[4]*1E+5*fgFuncIntYield->Integral(0.,1.);

  return res;
}
//------------------------------------------------------------------------------
const Float_t AliBlastwaveFitSpectra2::GetMeanBeta(){
    return fgFuncIntYield->GetParameter(2);
}
//------------------------------------------------------------------------------
const Float_t AliBlastwaveFitSpectra2::GetMeanBeta(Double_t par[]){
    return par[1];
}
