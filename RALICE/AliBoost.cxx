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
Revision 1.2  1999/09/29 09:24:28  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////
// Class AliBoost
// Perform various Lorentz transformations.
//
// Example :
// =========
//
// Float_t a[3]={0.1,0.2,0.3};
// Float_t ea[3]={0.01,0.02,0.03};
// Ali3Vector beta;
// beta.SetVector(a,"car");
// beta.SetErrors(ea,"car");
//
// AliBoost b1;
// b1.SetBeta(beta);
// b1.Info();
//
// Float_t b[4]={14,1,2,3};
// Float_t eb[4]={1.4,0.1,0.2,0.3};
// Ali4Vector p;
// p.SetVector(b,"car");
// p.SetErrors(eb,"car");
// Ali4Vector pprim=b1.Boost(p);
// p.Info();
// pprim.Info();
//
// p=b1.Inverse(pprim);
// pprim.Info();
// p.Info();
//
// Float_t c[4]={5,0,0,4};
// Float_t ec[4]={0.5,0,0,0.4};
// Ali4Vector q;
// q.SetVector(c,"car");
// q.SetErrors(ec,"car");
//
// AliBoost b2;
// b2.Set4Momentum(q);
// b2.Info("sph");
//
//--- Author: Nick van Eijndhoven 14-may-1996 UU-SAP Utrecht
//- Modified: NvE 24-oct-1999 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliBoost.h"
 
ClassImp(AliBoost) // Class implementation to enable ROOT I/O
 
AliBoost::AliBoost()
{
// Creation of a Lorentz boost object and initialisation of parameters.
// Beta is set to (0,0,0) and consequently Gamma=1. 
// All errors are initialised to 0. 
 Double_t a[3]={0,0,0};
 fBeta.SetVector(a,"sph");
 fGamma=1;
 fDgamma=0;
 fDresult=0;
}
///////////////////////////////////////////////////////////////////////////
AliBoost::~AliBoost()
{
// Default destructor.
}
///////////////////////////////////////////////////////////////////////////
void AliBoost::SetBeta(Ali3Vector b)
{
// Setting of boost parameters on basis of beta 3-vector.
// The errors on the beta 3-vector are taken from the input 3-vector.
// The gamma value and its error are calculated accordingly.
 fBeta=b;
 Double_t beta2=fBeta.Dot(fBeta);
 Double_t dbeta2=fBeta.GetResultError();

 if (beta2 > 1.)
 {
  cout << " *AliBoost::SetBeta* beta > 1." << endl;
 }
 fGamma=0;
 fDgamma=0;
 Double_t temp=1.-beta2;
 if (temp > 0.)
 {
  fGamma=sqrt(1./temp);
  fDgamma=fabs(dbeta2/(2.*pow(temp,1.5)));
 }
}
///////////////////////////////////////////////////////////////////////////
void AliBoost::Set4Momentum(Ali4Vector& p)
{
// Setting of boost parameters on basis of momentum 4-vector data.
// The errors of the input 4-vector are used to calculate the
// errors on the beta 3-vector and the gamma factor.
 Double_t E=p.GetScalar();
 Double_t dE=p.GetResultError();
 if (E <= 0.)
 {
  cout << " *AliBoost::Set4Momentum* Unphysical situation." << endl;
  p.Info();
 }
 else
 {
  Ali3Vector b=p.Get3Vector();
  Double_t vb[3],eb[3];
  b.GetVector(vb,"car");
  b.GetErrors(eb,"car");
  b=b/E;
  for (Int_t i=0; i<3; i++)
  {
   eb[i]=sqrt(pow(eb[i]/E,2)+pow(vb[i]*dE/(E*E),2));
  }
  b.SetErrors(eb,"car");
  SetBeta(b);
 }
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector AliBoost::GetBetaVector()
{
// Provide the beta 3-vector.
 return fBeta;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliBoost::GetBeta()
{
// Provide the norm of the beta 3-vector.
// The error on the value can be obtained via GetResultError().
 Double_t norm=fBeta.GetNorm();
 fDresult=fBeta.GetResultError();
 return norm;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliBoost::GetGamma()
{
// Provide the gamma factor.
// The error on the value can be obtained via GetResultError().
 fDresult=fDgamma;
 return fGamma;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliBoost::GetResultError()
{
// Provide the error on the result of an operation yielding a scalar.
// E.g. GetBeta() or GetGamma()
 return fDresult;
}
///////////////////////////////////////////////////////////////////////////
void AliBoost::Info(TString f)
{
// Printing of the boost parameter info in coordinate frame f.
 Double_t beta=fBeta.GetNorm();
 Double_t dbeta=fBeta.GetResultError();
 cout << " *AliBoost::Info* beta : " << beta << " error : " << dbeta
      << " gamma : " << fGamma << " error : " << fDgamma << endl;
 cout << "  Beta"; 
 fBeta.Info(f);
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector AliBoost::Boost(Ali4Vector& v)
{
// Perform the Lorentz boost on the 4-vector v.
// Error propagation is performed automatically.
// Note : As an approximation Beta and p.Dot(Beta) are considered as
//        independent quantities.

 Double_t beta=fBeta.GetNorm();
 Double_t dbeta=fBeta.GetResultError();

 Double_t beta2=pow(beta,2);

 if (beta > 1.e-10)
 {
  Double_t E=v.GetScalar();
  Double_t dE=v.GetResultError();

  Ali3Vector p=v.Get3Vector();

  Double_t pdotbeta=p.Dot(fBeta);
  Double_t dpdotbeta=p.GetResultError();

  // Determine the new vector components
  Double_t Eprim=fGamma*(E-pdotbeta);

  Double_t z=((fGamma-1.)*pdotbeta/beta2)-fGamma*E;
  Ali3Vector add=fBeta*z;

  // Determine errors on the new vector components
  Double_t dEprim=sqrt(pow((E-pdotbeta)*fDgamma,2)+pow(fGamma*dE,2)
                      +pow(fGamma*dpdotbeta,2));
  Double_t dz=sqrt( pow(((fGamma-1.)/beta2)*dpdotbeta,2) + pow(fGamma*dE,2)
                   +pow((
    ((2./beta)-(4.*pow(beta,3)-6.*pow(beta,5))/(2.*pow((pow(beta,4)-pow(beta,6)),1.5)))*pdotbeta
    +beta*E/pow(fGamma,3))*dbeta,2) );

  Double_t vb[3],eb[3];
  fBeta.GetVector(vb,"car");
  fBeta.GetErrors(eb,"car");
  for (Int_t i=0; i<3; i++)
  {
   eb[i]=sqrt(pow(z*eb[i],2)+pow(vb[i]*dz,2));
  }
  add.SetErrors(eb,"car");

  // Calculate the new 3-vector
  Ali3Vector pprim=p+add;

  // Set the components and errors of the new 4-vector 
  Ali4Vector w;
  w.SetVector(Eprim,pprim);
  w.SetScalarError(dEprim);

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
// Perform the inverse Lorentz boost on the 4-vector vprim.
// Error propagation is performed automatically.
// Note : As an approximation Beta and pprim.Dot(Beta) are considered as
//        independent quantities.

 Double_t beta=fBeta.GetNorm();
 Double_t dbeta=fBeta.GetResultError();

 Double_t beta2=pow(beta,2);

 if (beta > 1.e-10)
 {
  Double_t Eprim=vprim.GetScalar();
  Double_t dEprim=vprim.GetResultError();

  Ali3Vector pprim=vprim.Get3Vector();

  Double_t pprimdotbeta=pprim.Dot(fBeta);
  Double_t dpprimdotbeta=pprim.GetResultError();

  // Determine the new vector components
  Double_t E=fGamma*(Eprim+pprimdotbeta);

  Double_t z=((fGamma-1.)*pprimdotbeta/beta2)+fGamma*Eprim;
  Ali3Vector add=fBeta*z;

  // Determine errors on the prime-vector components
  Double_t dE=sqrt(pow((Eprim+pprimdotbeta)*fDgamma,2)+pow(fGamma*dEprim,2)
                      +pow(fGamma*dpprimdotbeta,2));
  Double_t dz=sqrt( pow(((fGamma-1.)/beta2)*dpprimdotbeta,2) + pow(fGamma*dEprim,2)
                   +pow((
    ((2./beta)-(4.*pow(beta,3)-6.*pow(beta,5))/(2.*pow((pow(beta,4)-pow(beta,6)),1.5)))*pprimdotbeta
    -beta*Eprim/pow(fGamma,3))*dbeta,2) );

  Double_t vb[3],eb[3];
  fBeta.GetVector(vb,"car");
  fBeta.GetErrors(eb,"car");
  for (Int_t i=0; i<3; i++)
  {
   eb[i]=sqrt(pow(z*eb[i],2)+pow(vb[i]*dz,2));
  }
  add.SetErrors(eb,"car");

  // Calculate the new 3-vector
  Ali3Vector p=pprim+add;

  // Set the components and errors of the new 4-vector 
  Ali4Vector w;
  w.SetVector(E,p);
  w.SetScalarError(dE);

  return w;
 }
 else
 {
  return vprim;
 }
}
///////////////////////////////////////////////////////////////////////////
