//--------------------------------------------------------------------------
//
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Module: EvtBtoXsgammaFermiUtil.cc
//
// Description:
//      Class to hold various fermi functions and their helper functions. The 
//      fermi functions are used in EvtBtoXsgammaKagan.
//
// Modification history:
//
//      Jane Tinslay       March 21, 2001       Module created
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "EvtGenModels/EvtBtoXsgammaFermiUtil.hh"
#include "EvtGenModels/EvtItgTwoCoeffFcn.hh"
#include "EvtGenModels/EvtBtoXsgammaRootFinder.hh"
#include "EvtGenModels/EvtItgFunction.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtReport.hh"

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <math.h>
using std::endl;

double EvtBtoXsgammaFermiUtil::FermiExpFunc(double y, const std::vector<double> &coeffs) {

  //coeffs: 1 = lambdabar, 2 = a, 3 = lam1, 4 = norm
  // report(Severity::Info,"EvtGen")<<coeffs[4]<<endl;
  return (pow(1. - (y/coeffs[1]),coeffs[2])*exp((-3.*pow(coeffs[1],2.)/coeffs[3])*y/coeffs[1]))/coeffs[4];

}

double EvtBtoXsgammaFermiUtil::FermiGaussFunc(double y, const std::vector<double> &coeffs) {

  //coeffs: 1 = lambdabar, 2 = a, 3 = c, 4 = norm
  return (pow(1. - (y/coeffs[1]),coeffs[2])*exp(-pow(coeffs[3],2.)*pow(1. - (y/coeffs[1]),2.)))/coeffs[4];

}

double EvtBtoXsgammaFermiUtil::FermiGaussFuncRoot(double lambdabar, double lam1, double mb, std::vector<double> &gammaCoeffs) {
 
  std::vector<double> coeffs1(3);
  std::vector<double> coeffs2(3);

  coeffs1[0]=0.2;
  coeffs1[1]=lambdabar;
  coeffs1[2]=0.0;
  
  coeffs2[0]=0.2;
  coeffs2[1]=lambdabar;
  coeffs2[2]=-lam1/3.;

  EvtItgTwoCoeffFcn *lhFunc = new EvtItgTwoCoeffFcn(&FermiGaussRootFcnA, -mb, lambdabar, coeffs1, gammaCoeffs);
  EvtItgTwoCoeffFcn *rhFunc = new EvtItgTwoCoeffFcn(&FermiGaussRootFcnB, -mb, lambdabar, coeffs2, gammaCoeffs);

  EvtBtoXsgammaRootFinder *rootFinder = new EvtBtoXsgammaRootFinder();
  
  double root = rootFinder->GetGaussIntegFcnRoot(lhFunc, rhFunc, 1.0e-4, 1.0e-4, 40, 40, -mb, lambdabar, 0.2, 0.4, 1.0e-6);

  delete rootFinder; rootFinder=0;
  delete lhFunc; lhFunc=0;
  delete rhFunc; rhFunc=0;
  return root;

}

double EvtBtoXsgammaFermiUtil::FermiGaussRootFcnA(double y, const std::vector<double> &coeffs1, const std::vector<double> &coeffs2) {

  
  //coeffs1: 0=ap, 1=lambdabar, coeffs2=gamma function coeffs
  double cp = Gamma((2.0 + coeffs1[0])/2., coeffs2)/Gamma((1.0 + coeffs1[0])/2., coeffs2);

  return (y*y)*pow((1. - (y/coeffs1[1])),coeffs1[0])*exp(-pow(cp,2)*pow((1.-(y/coeffs1[1])),2.));

}

double EvtBtoXsgammaFermiUtil::FermiGaussRootFcnB(double y, const std::vector<double> &coeffs1, const std::vector<double> &coeffs2) {

  //coeffs1: 0=ap, 1=lambdabar, coeffs2=gamma function coeffs
  double cp = Gamma((2.0 + coeffs1[0])/2., coeffs2)/Gamma((1.0 + coeffs1[0])/2., coeffs2);
  return pow((1. - (y/coeffs1[1])),coeffs1[0])*exp(-pow(cp,2)*pow((1.-(y/coeffs1[1])),2.));

}

double EvtBtoXsgammaFermiUtil::Gamma(double z, const std::vector<double> &coeffs) {

  //Lifted from Numerical Recipies in C
  double x, y, tmp, ser;
 
  int j;
  y = z;
  x = z;
  
  tmp = x + 5.5;
  tmp = tmp - (x+0.5)*log(tmp);
  ser=1.000000000190015;

  for (j=0;j<6;j++) {
    y = y +1.0;
    ser = ser + coeffs[j]/y;
  }

  return exp(-tmp+log(2.5066282746310005*ser/x));

}

double EvtBtoXsgammaFermiUtil::BesselK1(double x) {

  //Lifted from Numerical Recipies in C : Returns the modified Bessel
  //function K_1(x) for positive real x
  if (x<0.0) report(Severity::Info,"EvtGen") <<"x is negative !"<<endl;
  
  double y, ans;

  if (x <= 2.0) {
    y=x*x/4.0;
    ans = (log(x/2.0)*BesselI1(x))+(1.0/x)*(1.0+y*(0.15443144+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1+y*(-0.110404e-2+y*(-0.4686e-4)))))));
  }
  else {
    y=2.0/x;
    ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2+y*(0.325614e-2+y*(-0.68245e-3)))))));
  }
  return ans;

}

double EvtBtoXsgammaFermiUtil::BesselI1(double x) {
   //Lifted from Numerical Recipies in C : Returns the modified Bessel
  //function I_1(x) for any real x

  double ax, ans;
  double y;

  ax=fabs(x);
  if ( ax  < 3.75) {
    y=x/3.75;
    y*=y;
    ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  }
  else {
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1 -y*0.420059e-2));
    ans=0.398914228+y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans*=(exp(ax)/sqrt(ax));
  }
    return x < 0.0 ? -ans:ans;
}

double EvtBtoXsgammaFermiUtil::FermiRomanFuncRoot(double lambdabar, double lam1) {
 
  EvtItgFunction *lhFunc = new EvtItgFunction(&FermiRomanRootFcnA, -1.e-6, 1.e6);
  
  EvtBtoXsgammaRootFinder *rootFinder = new EvtBtoXsgammaRootFinder();
  double rhSide = 1.0 - (lam1/(3.0*lambdabar*lambdabar));

  double rho = rootFinder->GetRootSingleFunc(lhFunc, rhSide, 0.1, 0.4, 1.0e-6);
  //rho=0.250353;
  report(Severity::Info,"EvtGen")<<"rho/2 "<<rho/2.<<" bessel "<<BesselK1(rho/2.)<<endl;
  double pF = lambdabar*sqrt(EvtConst::pi)/(rho*exp(rho/2.)*BesselK1(rho/2.));
  report(Severity::Info,"EvtGen")<<"rho "<<rho<<" pf "<<pF<<endl;
  
  delete lhFunc; lhFunc=0;
  delete rootFinder; rootFinder=0;
  return rho;

}

double EvtBtoXsgammaFermiUtil::FermiRomanRootFcnA(double y) {

   return EvtConst::pi*(2. + y)*pow(y,-2.)*exp(-y)*pow(BesselK1(y/2.),-2.);

}
double EvtBtoXsgammaFermiUtil::FermiRomanFunc(double y, const std::vector<double> &coeffs) {
  if (y == (coeffs[1]-coeffs[2])) y=0.99999999*(coeffs[1]-coeffs[2]);

  //coeffs: 1 = mB, 2=mb, 3=rho, 4=lambdabar, 5=norm
  double pF = coeffs[4]*sqrt(EvtConst::pi)/(coeffs[3]*exp(coeffs[3]/2.)*BesselK1(coeffs[3]/2.));
  //  report(Severity::Info,"EvtGen")<<" pf "<<y<<" "<<pF<<" "<<coeffs[1]<<" "<<coeffs[2]<<" "<<coeffs[3]<<" "<<coeffs[4]<<" "<<coeffs[5]<<endl;
  //double pF=0.382533;

  //report(Severity::Info,"EvtGen")<<(coeffs[1]-coeffs[2])*(1./(sqrt(EvtConst::pi)*pF))<<endl;
  //report(Severity::Info,"EvtGen")<<(1.-y/(coeffs[1]-coeffs[2]))<<endl;
  //report(Severity::Info,"EvtGen")<<(coeffs[1]-coeffs[2])<<endl;
  //report(Severity::Info,"EvtGen")<<(coeffs[1]-coeffs[2])*(1.-y/(coeffs[1]-coeffs[2]))<<endl;

  //report(Severity::Info,"EvtGen")<<" "<<pF*coeffs[3]/((coeffs[1]-coeffs[2])*(1.-y/(coeffs[1]-coeffs[2])))<<endl;
  // report(Severity::Info,"EvtGen")<<" "<<((coeffs[1]-coeffs[2])/pF)*(1. -y/(coeffs[1]-coeffs[2]))<<endl;

  //report(Severity::Info,"EvtGen")<<"result "<<(coeffs[1]-coeffs[2])*(1./(sqrt(EvtConst::pi)*pF))*exp(-(1./4.)*pow(pF*(coeffs[3]/((coeffs[1]-coeffs[2])*(1.-y/(coeffs[1]-coeffs[2])))) - ((coeffs[1]-coeffs[2])/pF)*(1. -y/(coeffs[1]-coeffs[2])),2.))/coeffs[5];

  //report(Severity::Info,"EvtGen")<<"leaving"<<endl;
  return (coeffs[1]-coeffs[2])*(1./(sqrt(EvtConst::pi)*pF))*exp(-(1./4.)*pow(pF*(coeffs[3]/((coeffs[1]-coeffs[2])*(1.-y/(coeffs[1]-coeffs[2])))) - ((coeffs[1]-coeffs[2])/pF)*(1. -y/(coeffs[1]-coeffs[2])),2.))/coeffs[5];
 

}
