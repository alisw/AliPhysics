#include<stdio.h>

#include"AliBlastwaveFitSpectra.h"

#include"TMath.h"

ClassImp(AliBlastwaveFitSpectra);


TF1 *AliBlastwaveFitSpectra::fgFuncIntYield = NULL;
const char *AliBlastwaveFitSpectra::fgParName[5] = {"T_{FO}","mean #rho","#gamma","mass","norm"};
Float_t AliBlastwaveFitSpectra::fgStartValues[3] = {0.1,1.2,1.1};
const Float_t AliBlastwaveFitSpectra::fgStepValues[3] = {0.001,0.001,0.001};
const Float_t AliBlastwaveFitSpectra::fgMinValues[3] = {0.0001,0.01,0.5};
const Float_t AliBlastwaveFitSpectra::fgMaxValues[3] = {1.0,10.0,5};

AliBlastwaveFitSpectra::AliBlastwaveFitSpectra(const char *name,Double_t mass) :
  AliBlastwaveFit(name,mass)
{
  if(!fgFuncIntYield){
     fgFuncIntYield = new TF1("fFuncInt2Yield",AliBlastwaveFitSpectra::FunctionIntYield,0,1,5);
     fgFuncIntYield->SetNpx(100);
  }

  Initialize();
}
//------------------------------------------------------------------------------
AliBlastwaveFitSpectra::AliBlastwaveFitSpectra() :
  AliBlastwaveFit("BlastwaveFitSpectra",0.0)
{  
  if(!fgFuncIntYield){
     fgFuncIntYield = new TF1("fFuncInt2Yield",AliBlastwaveFitSpectra::FunctionIntYield,0,1,5);
     fgFuncIntYield->SetNpx(100);
  }

  Initialize();
}
//------------------------------------------------------------------------------
AliBlastwaveFitSpectra::~AliBlastwaveFitSpectra(){
}
//------------------------------------------------------------------------------
Double_t AliBlastwaveFitSpectra::FunctionIntYield(Double_t x[],Double_t par[]){
  // x[0] = r/R
  // par[0] = pt
  // par[1] = T
  // par[2] = rho_av
  // par[3] = R_power
  // par[4] = mass 

  Double_t maxrho = par[2]*0.5*(par[3]+2);

  Double_t mt = TMath::Sqrt(par[0]*par[0]+par[4]*par[4]);
  Double_t rho = maxrho*TMath::Power(x[0],par[3]);
  Double_t alfat = (par[0]/par[1])*TMath::SinH(rho);
  Double_t betat = (mt/par[1])*TMath::CosH(rho);
  Double_t results;
  if(betat < 200) results = TMath::BesselI(0,alfat)*TMath::BesselK(1,betat) *x[0];
  else results = 0.5*TMath::Exp(alfat-betat)/sqrt(alfat*betat) *x[0];

  return results*par[0]*mt;
}
//------------------------------------------------------------------------------
Int_t AliBlastwaveFitSpectra::SetParameter(Int_t ipar,Double_t val){
  Bool_t status = kTRUE;

  if(ipar < 0 || ipar > 4) return status;

  if(fFunctionYield){
     fFunctionYield->SetParameter(ipar,val);
     status = kFALSE;
  }
  return status;
}
//------------------------------------------------------------------------------
Int_t AliBlastwaveFitSpectra::SetNormalization(){
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
void AliBlastwaveFitSpectra::SetMass(Double_t mass) {
  fMass=mass;
  if(fFunctionYield) fFunctionYield->SetParameter(3,fMass);   
}
//------------------------------------------------------------------------------
void AliBlastwaveFitSpectra::Initialize(){
  char nameFuncYield[100];
  snprintf(nameFuncYield,100,"%sFuncYield",GetName());

  if(fFunctionYield) delete fFunctionYield;
  fFunctionYield = new TF1(nameFuncYield,Pt,0,6,7);
  fFunctionYield->SetParameter(3,fMass);
  fFunctionYield->SetNpx(200);

  for(Int_t i =0;i < 5;i++) fFunctionYield->SetParName(i,fgParName[i]);
}
//------------------------------------------------------------------------------
Double_t AliBlastwaveFitSpectra::Pt(Double_t x[],Double_t par[]){
  // par[0] = T_fo
  // par[1] = rho_av
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
const Float_t AliBlastwaveFitSpectra::GetMeanBeta(){
    TF1 fbeta("fbeta","TMath::TanH([0]*TMath::Power(x,[1]))*x",0,1);
    fbeta.SetNpx(1000);

    Double_t maxrho = fgFuncIntYield->GetParameter(2)*0.5*(fgFuncIntYield->GetParameter(3)+2);

    fbeta.SetParameter(0,maxrho);
    fbeta.SetParameter(1,fgFuncIntYield->GetParameter(3));
    
    return fbeta.Integral(0.,1.)*2;
}
//------------------------------------------------------------------------------
const Float_t AliBlastwaveFitSpectra::GetMeanBeta(Double_t par[]){
  // par[0] = T_fo
  // par[1] = rho_av
  // par[2] = R_power

    TF1 fbeta("fbeta","TMath::TanH([0]*TMath::Power(x,[1]))*x",0,1);
    fbeta.SetNpx(1000);

    Double_t maxrho = par[1]*0.5*(par[2]+2);

    fbeta.SetParameter(0,maxrho);
    fbeta.SetParameter(1,par[2]);
    
    return fbeta.Integral(0.,1.)*2;
}
