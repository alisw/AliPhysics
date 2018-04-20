// -*- C++ -*-

#include <TF1.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TMinuitMinimizer.h>
#include <Math/Functor.h>

class AliVdMPileup {
public:
  AliVdMPileup()
    : fMinuit("minimize", 4)
    , fAnotC()
    , fCnotA() {}

  void DoFit(TGraphErrors* gAnotC,
             TGraphErrors* gCnotA,
             const Double_t par0[4])
  {
    fAnotC = gAnotC;
    fCnotA = gCnotA;

    fMinuit.UseStaticMinuit(kTRUE);
    ROOT::Math::Functor fcn(*this, 4);
    fMinuit.SetFunction(fcn);

    fMinuit.SetLimitedVariable( 0, "r_{A}",    par0[0], 0.01, 0, 1);
    fMinuit.SetLimitedVariable( 1, "r_{C}",    par0[1], 0.01, 0, 1);
    fMinuit.SetLimitedVariable( 2, "bkgd_{A}", par0[2], 5e-6, 5e-7, 1e-3);
    fMinuit.SetLimitedVariable( 3, "bkgd_{C}", par0[3], 5e-6, 5e-7, 1e-3);
    fMinuit.Minimize();
    fMinuit.PrintResults();
  }

  Double_t GetChi2() const { return fMinuit.MinValue(); }
  Double_t GetNDF() const { return fAnotC->GetN() + fCnotA->GetN() - fMinuit.NFree(); }

  Double_t operator()(const Double_t *par) {
    const Double_t eps = 1e-7;
    Double_t val = 0.0;
    const Double_t *x=NULL, *y=NULL, *ex=NULL, *ey=NULL;
    x  = fCnotA->GetX();
    y  = fCnotA->GetY();
    ex = fCnotA->GetEY();
    ey = fCnotA->GetEY();
    TF1 fcnCnotA("fcnCnotA", this, &AliVdMPileup::fcnCnotA,0,1,4);
    fcnCnotA.SetParameters(par);
    for (Int_t j=0, n=fCnotA->GetN(); j<n; ++j) {
      const Double_t mu = x[j];
      if (mu < 0) continue;
      val += TMath::Power(y[j]-fcnCnotA.Eval(mu), 2)/(ey[j]*ey[j] +
                                                      0*TMath::Power(ex[j]*fcnCnotA.Derivative(mu, nullptr, eps), 2));
    }
    x  = fAnotC->GetX();
    y  = fAnotC->GetY();
    ex = fAnotC->GetEX();
    ey = fAnotC->GetEY();
    TF1 fcnAnotC("fcnAnotC", this, &AliVdMPileup::fcnAnotC,0,1,4);
    fcnAnotC.SetParameters(par);
    for (Int_t j=0, n=fAnotC->GetN(); j<n; ++j) {
      const Double_t mu = x[j];
      if (mu < 0) continue;
      val += TMath::Power(y[j]-fcnAnotC.Eval(mu), 2)/(ey[j]*ey[j] +
                                                      0*TMath::Power(ex[j]*fcnAnotC.Derivative(mu, nullptr, eps), 2));
    }
    return val;
  }

  const Double_t* GetPar() const { return fMinuit.X(); }

  Double_t fcnAnotC(const Double_t *x, const Double_t *par) const {
    const Double_t mu = trova(x[0], par);
    return rAnotC(&mu,par)/rAC(&mu,par);
  }
  Double_t fcnCnotA(const Double_t *x, const Double_t *par) const {
    const Double_t mu = trova(x[0], par);
    return rCnotA(&mu,par)/rAC(&mu,par);
  }

  static Double_t trova(Double_t y, const Double_t *par) {
    Double_t x[2] = {0.0, 2.0}; // xmin,xmax
    const Double_t eps=1e-16;   // 1e-14
    for(Int_t i=0; i<100 && x[1]-x[0] > eps; ++i){
      const Double_t xmed = 0.5*(x[1]+x[0]);
      const Double_t val  = rAC(&xmed, par);
      x[val >= y] = xmed;
    }
    return 0.5*(x[1]+x[0]);
  }
private:
  static Double_t rAnotC(const Double_t* x, const Double_t* par) {
    const Double_t mu = x[0];
    return (1 - TMath::Exp(-par[0]*mu-par[2])) * TMath::Exp(-mu) * TMath::Exp(-par[1]*mu-par[3]);
  }
  static Double_t rCnotA(const Double_t* x, const Double_t* par) {
    const Double_t mu = x[0];
    return (1 - TMath::Exp(-par[1]*mu-par[3])) * TMath::Exp(-mu) * TMath::Exp(-par[0]*mu-par[2]);
  }
  static Double_t rAC(const Double_t* x, const Double_t* par) {
    const Double_t mu = x[0];
    return 1-TMath::Exp(-mu) + TMath::Exp(-mu)*(1 - TMath::Exp(-par[1]*mu-par[3]))*(1 - TMath::Exp(-par[0]*mu-par[2]));
  }

  TMinuitMinimizer fMinuit;
  TGraphErrors* fAnotC;
  TGraphErrors* fCnotA;
} ;
