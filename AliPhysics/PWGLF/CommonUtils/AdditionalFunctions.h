#ifndef ADDITIONAL_FUNCTIONS_H
#define ADDITIONAL_FUNCTIONS_H

#if !defined (__CINT__) || (defined(__MAKECINT__))
#include "TMath.h"
#include "TF1.h"
#include "Riostream.h"
#endif

double LevyTsallis_Func(const double *x, const double *p) {
  /* dN/dpt */

  double pt = x[0];
  double mass = p[0];
  double mt = TMath::Sqrt(pt * pt + mass * mass);
  double n = p[1];
  double C = p[2];
  double norm = p[3];

  double part1 = (n - 1.) * (n - 2.);
  double part2 = n * C * (n * C + mass * (n - 2.));
  double part3 = part1 / part2;
  double part4 = 1. + (mt - mass) / n / C;
  double part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 * LevyTsallis(const Char_t *name, double mass, double n = 5., double C = 0.1, double norm = 1.) {

  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e6);
  return fLevyTsallis;
}

double Hagedorn_Func(double *x, double *par) {
  Float_t pt =x[0];
  double mt=TMath::Sqrt(pt*pt+par[5]*par[5]);
  double f=-1.0*par[0]*mt-par[1]*mt*mt;
  f=TMath::Exp(f)+mt/par[2];
  f=TMath::Power(f,par[3]);
  f=pt*par[4]/f;
  return f;
}


TF1 * Hagedorn (const Char_t *name, double mass) {

  TF1 *fun = new TF1(name, Hagedorn_Func, 0, 50., 6);
  fun->SetParNames("a","b","m0","n","Scale","mass");
  fun->FixParameter(5,mass);
  return fun;

}
double PowerLawdNdptTimesPt(const double x, const double* p)
{
  return p[0]*x*TMath::Power((1.0 + x/p[1] ),(-1.0*p[2]));
}


double MTExpdNdptTimesPt(const double x, const double* p) {
  return p[0]*x*TMath::Exp(-1.0*TMath::Sqrt(x*x+p[2]*p[2])/p[1]);
}

double UA1_Func(const double * x, const double* p) {

  // "mass","p0star","pt0","n","T","norm"
  double mass   = p[0];
  double p0star = p[1];
  double pt0    = p[2];
  double n      = p[3];
  double temp   = p[4];
  double norm   = p[5];
  double xx = x[0];

  double parplaw[3]={0.0,0.0,0.0};
  parplaw[0]=norm;
  parplaw[1]=pt0;
  parplaw[2]=n;
  double parpexp[3]={0.0,0.0,0.0};
  parpexp[0]=norm;
  parpexp[1]=temp;
  parpexp[2]=mass;


  double normMT =MTExpdNdptTimesPt(p0star,parpexp) >0 ? PowerLawdNdptTimesPt(p0star,parplaw) / MTExpdNdptTimesPt(p0star, parpexp) *  parpexp[0]: 1;
  parpexp[0]=normMT;

  if (TMath:: Abs( MTExpdNdptTimesPt(p0star, parpexp)- PowerLawdNdptTimesPt(p0star,parplaw)) > 0.0001 ) {
    Printf("UA1 - Wrong norm") ;
    Printf(" p0* %f  NMT: %f  N: %f  PL: %f  MT: %f", p0star, normMT, norm, PowerLawdNdptTimesPt(p0star,parplaw),MTExpdNdptTimesPt(p0star, parpexp));
  }

  if (xx > p0star) return PowerLawdNdptTimesPt(xx,parplaw);
  else return MTExpdNdptTimesPt(xx, parpexp);

}

TF1 * UA1 (const Char_t *name, double mass) {
  TF1 *fun = new TF1(name, UA1_Func, 0, 50., 6);
  fun->FixParameter(0,mass);
  return fun;
}

double Bylinkin_Func(const double * x, const double* p) {
  double mass   = p[0];
  double a1 = p[1];
  double a2 = p[2];
  double t1  = p[3];
  double t2  = p[4];
  double n   = p[5];

  double pt=x[0];
  double mt=TMath::Sqrt(pt*pt+mass*mass);

  double f1=TMath::Exp(-1.0*(mt-mass)/t1);
  double f2=1.0+(pt*pt)/(t2*t2*n);
  f2=TMath::Power(f2,-1.0*n);
  return pt*a1*(f1+a2*f2);
}

TF1 * Bylinkin(const Char_t *name, double mass) {
  TF1 *fun = new TF1(name,Bylinkin_Func , 0, 50., 6);
  fun->SetParNames("mass","Norm_{global}","Norm_{n}","T_{exp}","T_{n}","n");
  fun->FixParameter(0,mass);
  return fun;
}

double Boltzmann_Func(const double *x, const double *p) {
  /* dN/dpt */
  double pt = x[0];
  double mass = p[0];
  double mt = TMath::Sqrt(pt * pt + mass * mass);
  double T = p[1];
  double norm = p[2];

  return pt * norm * mt * TMath::Exp(-mt / T);
}

#endif
