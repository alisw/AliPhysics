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

/*
$Log$
*/

#include "AliBoost.h"
 
ClassImp(AliBoost) // Class implementation to enable ROOT I/O
 
AliBoost::AliBoost()
{
// Creation of a Lorentz boost object and initialisation of parameters
 fGamma=1;
 fBeta2=0;
 Double_t a[3]={0,0,0};
 fBeta.SetVector(a,"sph");
}
///////////////////////////////////////////////////////////////////////////
AliBoost::~AliBoost()
{
// Default destructor
}
///////////////////////////////////////////////////////////////////////////
void AliBoost::SetBeta(Ali3Vector b)
{
// Setting of boost parameters on basis of beta 3-vector
 fBeta2=b.Dot(b);
 fBeta=b;

 if (fBeta2 > 1.)
 {
  cout << " *AliBoost::SetBeta* beta > 1." << endl;
 }
 Double_t test=1.-fBeta2;
 fGamma=0;
 if (test > 0.) fGamma=sqrt(1./test);
}
///////////////////////////////////////////////////////////////////////////
void AliBoost::SetGamma(Double_t g,Ali3Vector v)
{
// Setting of boost parameters on basis of gamma and direction 3-vector
 if (g >= 1.)
 {
  fGamma=g;
  fBeta2=1.-(1./(fGamma*fGamma));
  fBeta=v*sqrt(fBeta2);
 }
 else
 {
  cout << " *AliBoost::SetGamma* Invalid input gamma = " << g << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliBoost::Set4Momentum(Ali4Vector& p)
{
// Setting of boost parameters on basis of momentum 4-vector data
 Double_t E=p.GetScalar();
 if (E <= 0.)
 {
  cout << " *AliBoost::Set4Momentum* Unphysical situation." << endl;
  p.Info();
 }
 else
 {
  Ali3Vector b=p.Get3Vector();
  b=b/E;
  SetBeta(b);
 }
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector AliBoost::GetBetaVector()
{
// Provide the the beta 3-vector
 return fBeta;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliBoost::GetBeta()
{
// Provide the norm of the beta 3-vector
 return sqrt(fBeta2);
}
///////////////////////////////////////////////////////////////////////////
Double_t AliBoost::GetGamma()
{
// Provide the gamma factor
 return fGamma;
}
///////////////////////////////////////////////////////////////////////////
void AliBoost::Info(TString f)
{
// Printing of the boost parameter info in coordinate frame f

 cout << " *AliBoost::Info* beta = " << sqrt(fBeta2) << " gamma = " << fGamma << endl
      << "  Beta";
 fBeta.Info(f);
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector AliBoost::Boost(Ali4Vector& v)
{
// Perform the Lorentz boost on the 4-vector v
 if (fBeta2 > 1.e-20)
 {
  Double_t E=v.GetScalar();
  Ali3Vector p=v.Get3Vector();
  Double_t pdotbeta=p.Dot(fBeta);

  Double_t Eprim;
  Eprim=fGamma*(E-pdotbeta);

  Ali3Vector term1,term2,pprim;
  term1=fBeta*((fGamma-1.)*pdotbeta/fBeta2);
  term2=fBeta*(fGamma*E);
  pprim=p+term1-term2;

  Ali4Vector w;
  w.SetVector(Eprim,pprim);

  return w;
 }
 else
 {
  return v;
 }
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector AliBoost::Inverse(Ali4Vector& vprim)
{
// Perform the inverse Lorentz boost on the 4-vector vprim
 if (fBeta2 > 1.e-20)
 {
  Double_t Eprim=vprim.GetScalar();
  Ali3Vector pprim=vprim.Get3Vector();
  Double_t pprimdotbeta=pprim.Dot(fBeta);

  Double_t E;
  E=fGamma*(Eprim+pprimdotbeta);

  Ali3Vector term1,term2,p;
  term1=fBeta*((fGamma-1.)*pprimdotbeta/fBeta2);
  term2=fBeta*(fGamma*Eprim);
  p=pprim+term1+term2;

  Ali4Vector w;
  w.SetVector(E,p);

  return w;
 }
 else
 {
  return vprim;
 }
}
///////////////////////////////////////////////////////////////////////////
