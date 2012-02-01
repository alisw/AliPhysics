/************************************************************************* 
* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. * 
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

//Author: Dariusz Miskowiec
//Date:   2010

//=============================================================================
// Coulomb correlation function
//=============================================================================

#include <TRandom2.h>
#include <TComplex.h>
#include <TFile.h>
#include <TH2.h>
#include "AliUnicorCoulomb.h"

ClassImp(AliUnicorCoulomb)

//=============================================================================
AliUnicorCoulomb::AliUnicorCoulomb(int zz, double mass, double R) : TGraph(1001) {

  // constructor
  // zz = +1 for pi+pi+, -1 for pi-pi-, etc. 
  // mass is the reduced mass: m1*m2/(m1+m2)

  SetMarkerStyle(20);
  int nstep = 10000;

  double rx = R;
  double ry = R;
  double rz = R;
  double r0 = 0;
  TRandom2 ran;

  if (zz>0) SetPoint(0,0,0);
  else SetPoint(0,0,1e9);
  for (int i=1; i<1000; i++) {
    double qinv = i/1000.0;
    double k = qinv/2.0;
    double sum = 0;
    for (int j=0; j<nstep; j++) {
      double x0 = ran.Gaus(0.0,rx);
      double y0 = ran.Gaus(0.0,ry);
      double z0 = ran.Gaus(0.0,rz);
      double t0 = ran.Gaus(0.0,r0);
      double x1 = ran.Gaus(0.0,rx);
      double y1 = ran.Gaus(0.0,ry);
      double z1 = ran.Gaus(0.0,rz);
      double t1 = ran.Gaus(0.0,r0);
      double dx = x1-x0;
      double dy = y1-y0;
      double dz = z1-z0;
      double dt = t1-t0;
      dz += dt*k/mass;
      sum += WaveFunction2(zz,mass,k,dx,dy,dz);
    }
    SetPoint(i,qinv,sum/nstep);
  }
  SetPoint(1000,1000,1);
}
//=============================================================================
double AliUnicorCoulomb::WaveFunction2(int zz, double mass, double k, double x, double y, double z) {

  // Square of normalized Coulomb wave function. 
  // Hypergeometrical function diverges for numerical reasons for large Qinv, z and x. 
  // Out of the convergence limit I will return 1.0

  double qinv12 = 1.2 * 2.0 * k;
  double p0 = -9.1569/qinv12;
  double p2 = 2.735e-2*qinv12;
  double wf2 = 1.0;
  int converg = (z>p0+p2*(x*x+y*y));
  if (converg) {
    TComplex co = WaveFunction(zz,mass,k,x,y,z);
    wf2 =  co.Rho2()*Gamow(zz,mass,k);
  }
  return wf2;
}
//=============================================================================
TComplex AliUnicorCoulomb::WaveFunction(int zz, double mass, double k, double x, double y, double z) {

  // Coulomb wave function in parabolic coordinates
  // Merzbacher, Quantum Mechanics, p.248
  // Landau-Lifshitz, Quantum Mechanics (Non-Relativistic Theory), p.518

  double r = sqrt(x*x+y*y+z*z);
  //  double ksi = r+z;
  double eta = r-z;
  TComplex jjj(0,1);
  TComplex a = -zz*mass/k/137.036*jjj;
  TComplex b(1.0,0.0);
  TComplex c = jjj*k*eta/0.197327;
  TComplex wf = TComplex::Exp(jjj*k*z/0.197327)*F1(a,b,c);
  //printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n",a.Re(),a.Im(),b.Re(),b.Im(),c.Re(),c.Im(),F1(a,b,c).Rho());
  return wf;
}
//=============================================================================
double AliUnicorCoulomb::Gamow(int zz, double m, double k) {

  // Squared Coulomb wave function in infinity = Gamow
  // zz - product of charges
  // m - reduced mass
  // k = Qinv/2 for equal masses

  double n = -zz*m/k/137.036;
  double a = 2*TMath::Pi()*n;
  double ga = a/(1.0-exp(-a));
  return ga;
}
//=============================================================================
TComplex AliUnicorCoulomb::F1(TComplex alpha, TComplex gamma, TComplex z) {

  // confluent hypergeometric function 1F1

  double prec = 0.0000001;
  TComplex term(1.0,0.0);
  TComplex sum(1.0,0.0);
  for (int n=1; n<1000; n++) {
    double u = (double) n;
    TComplex v = (alpha+u-1.0)/(gamma+u-1.0);
    TComplex w = z/u;
    term *= v;
    term *= w;
    sum += term;
    //    printf("%10d %10.3f %10.3f %10.3f %10.3f\n", n, v.Rho(), w.Rho(), term.Rho(), sum.Rho());
    if (TComplex::Abs(term/sum) < prec) return sum;
  }
  printf("F1 Maximum number of iterations reached\n");
  return sum;
}
//=============================================================================
 void AliUnicorCoulomb::Makehist(int zz, double m, const char *outfil) {

  // Make a two-dim histogram coulomb(R,Q) and store it on outfil.
  // Later to be used via TH2::Interpolate(). 
  // zz and m are the product of charges and the reduced mass. 

  TH2D *hi = new TH2D("co","coulomb",11,-0.5,10.5,999,0.0005,0.9995);
  hi->SetXTitle("R (fm)");
  hi->SetYTitle("Q (GeV/c)");
  AliUnicorCoulomb *co = 0;
  for (int i=1; i<=hi->GetNbinsX(); i++) {
    double R = hi->GetXaxis()->GetBinCenter(i);
    printf("R = %.1f fm\n",R); 
    co = new AliUnicorCoulomb(zz, m, R);
    for (int j=1; j<=hi->GetNbinsY(); j++) {
      double Q = hi->GetYaxis()->GetBinCenter(j);
      hi->SetBinContent(i,j,co->Eval(Q));
    }
    delete co;
  }
  TFile::Open(outfil,"recreate");
  hi->Write();
  gFile->Close();
}
//=============================================================================
