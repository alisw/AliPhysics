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
// Class Ali3Vector
// Handling of 3-vectors in various reference frames.
//
// This class is meant to serve as a base class for ALICE objects
// that have 3-dimensional vector characteristics.
// Error propagation is performed automatically. 
//
// Note :
// ------
// Vectors (v), Errors (e) and reference frames (f) are specified via
// SetVector(Float_t* v,TString f)
// SetErrors(Float_t* e,TString f)
// under the following conventions :
//
// f="car" ==> v in Cartesian coordinates   (x,y,z)
// f="sph" ==> v in Spherical coordinates   (r,theta,phi)
// f="cyl" ==> v in Cylindrical coordinates (rho,phi,z)
//
// All angles are in radians.
//
// Example :
// ---------
//
// Ali3Vector a;
// Float_t v[3]={-1,25,7};
// Float_t e[3]={0.03,0.5,0.21};
// a.SetVector(v,"car");
// a.SetErrors(e,"car");
// a.Info();
//
// Float_t vec[3];
// Float_t err[3];
// a.GetVector(vec,"sph");
// a.GetErrors(vec,"sph");
//
// Ali3Vector b;
// Float_t v2[3]={6,-18,33};
// Float_t e2[3]={0.19,0.45,0.93};
// b.SetVector(v2,"car");
// b.SetErrors(e2,"car");
//
// Float_t dotpro=a.Dot(b);
// Float_t doterror=a.GetResultError();
//
// Ali3Vector c=a.Cross(b);
// c.Info("sph");
// c.GetVector(vec,"cyl");
// c.GetErrors(err,"cyl");
//
// Float_t norm=c.GetNorm();
// Float_t normerror=c.GetResultError();
//
// c=a+b;
// c=a-b;
// c=a*5;
//
//--- Author: Nick van Eijndhoven 30-mar-1999 UU-SAP Utrecht
//- Modified: NvE 25-oct-1999 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "Ali3Vector.h"
 
ClassImp(Ali3Vector) // Class implementation to enable ROOT I/O
 
Ali3Vector::Ali3Vector()
{
// Creation of an Ali3Vector object and initialisation of parameters
// All attributes initialised to 0
 fV=0;
 fTheta=0;
 fPhi=0;
 fDx=0;
 fDy=0;
 fDz=0;
 fDresult=0;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector::~Ali3Vector()
{
// Destructor to delete dynamically allocated memory
}
///////////////////////////////////////////////////////////////////////////
void Ali3Vector::SetVector(Double_t* v,TString f)
{
// Store vector according to reference frame f
// All errors will be reset to 0
 fDx=0;
 fDy=0;
 fDz=0;
 fDresult=0;

 Double_t pi=acos(-1.);

 Int_t frame=0;
 if (f == "car") frame=1;
 if (f == "sph") frame=2;
 if (f == "cyl") frame=3;

 Double_t x,y,z,rho,phi;

 switch (frame)
 {
  case 1: // Cartesian coordinates
   x=v[0];
   y=v[1];
   z=v[2];
   fV=sqrt(x*x+y*y+z*z);
   fTheta=0;
   if (fV && fabs(z/fV)<=1.)
   {
    fTheta=acos(z/fV);
   }
   else
   {
    if (z<0.) fTheta=pi;
   }
   if (fTheta<0.) fTheta+=2.*pi;
   fPhi=0;
   if (x || y) fPhi=atan2(y,x);
   if (fPhi<0.) fPhi+=2.*pi;
   break;

  case 2: // Spherical coordinates
   fV=v[0];
   fTheta=v[1];
   fPhi=v[2];
   break;

  case 3: // Cylindrical coordinates
   rho=v[0];
   phi=v[1];
   z=v[2];
   fV=sqrt(rho*rho+z*z);
   fPhi=phi;
   if (fPhi<0.) fPhi+=2.*pi;
   fTheta=0;
   if (fV && fabs(z/fV)<=1.)
   {
    fTheta=acos(z/fV);
   }
   else
   {
    if (z<0.) fTheta=pi;
   }
   if (fTheta<0.) fTheta+=2.*pi;
   break;

  default: // Unsupported reference frame
   cout << "*Ali3Vector::SetVector* Unsupported frame : " << f << endl
        << " Possible frames are 'car', 'sph' and 'cyl'." << endl; 
   fV=0;
   fTheta=0;
   fPhi=0;
   break;
 }
}
///////////////////////////////////////////////////////////////////////////
void Ali3Vector::GetVector(Double_t* v,TString f)
{
// Provide vector according to reference frame f
 Int_t frame=0;
 if (f == "car") frame=1;
 if (f == "sph") frame=2;
 if (f == "cyl") frame=3;

 switch (frame)
 {
  case 1: // Cartesian coordinates
   v[0]=fV*sin(fTheta)*cos(fPhi);
   v[1]=fV*sin(fTheta)*sin(fPhi);
   v[2]=fV*cos(fTheta);
   break;

  case 2: // Spherical coordinates
   v[0]=fV;
   v[1]=fTheta;
   v[2]=fPhi;
   break;

  case 3: // Cylindrical coordinates
   v[0]=fV*sin(fTheta);
   v[1]=fPhi;
   v[2]=fV*cos(fTheta);
   break;

  default: // Unsupported reference frame
   cout << "*Ali3Vector::GetVector* Unsupported frame : " << f << endl
        << " Possible frames are 'car', 'sph' and 'cyl'." << endl; 
   for (Int_t i=0; i<3; i++)
   {
    v[i]=0;
   }
   break;
 }
}
///////////////////////////////////////////////////////////////////////////
void Ali3Vector::SetVector(Float_t* v,TString f)
{
// Store vector according to reference frame f
// All errors will be reset to 0
 Double_t vec[3];
 for (Int_t i=0; i<3; i++)
 {
  vec[i]=v[i];
 }
 SetVector(vec,f);
}
///////////////////////////////////////////////////////////////////////////
void Ali3Vector::GetVector(Float_t* v,TString f)
{
// Provide vector according to reference frame f
 Double_t vec[3];
 GetVector(vec,f);
 for (Int_t i=0; i<3; i++)
 {
  v[i]=vec[i];
 }
}
///////////////////////////////////////////////////////////////////////////
void Ali3Vector::SetErrors(Double_t* e,TString f)
{
// Store errors according to reference frame f
// The error on scalar results is reset to 0
 fDresult=0;

 Int_t frame=0;
 if (f == "car") frame=1;
 if (f == "sph") frame=2;
 if (f == "cyl") frame=3;

 Double_t dx2,dy2,dz2,rho;

 switch (frame)
 {
  case 1: // Cartesian coordinates
   fDx=fabs(e[0]);
   fDy=fabs(e[1]);
   fDz=fabs(e[2]);
   break;

  case 2: // Spherical coordinates
   dx2=pow((cos(fPhi)*sin(fTheta)*e[0]),2)+pow((fV*cos(fTheta)*cos(fPhi)*e[1]),2)
       +pow((fV*sin(fTheta)*sin(fPhi)*e[2]),2);
   dy2=pow((sin(fPhi)*sin(fTheta)*e[0]),2)+pow((fV*cos(fTheta)*sin(fPhi)*e[1]),2)
       +pow((fV*sin(fTheta)*cos(fPhi)*e[2]),2);
   dz2=pow((cos(fTheta)*e[0]),2)+pow((fV*sin(fTheta)*e[1]),2);
   fDx=sqrt(dx2);
   fDy=sqrt(dy2);
   fDz=sqrt(dz2);
   break;

  case 3: // Cylindrical coordinates
   rho=fV*sin(fTheta);
   dx2=pow((cos(fPhi)*e[0]),2)+pow((rho*sin(fPhi)*e[1]),2);
   dy2=pow((sin(fPhi)*e[0]),2)+pow((rho*cos(fPhi)*e[1]),2);
   fDx=sqrt(dx2);
   fDy=sqrt(dy2);
   fDz=fabs(e[2]);
   break;

  default: // Unsupported reference frame
   cout << "*Ali3Vector::SetErrors* Unsupported frame : " << f << endl
        << " Possible frames are 'car', 'sph' and 'cyl'." << endl; 
   fDx=0;
   fDy=0;
   fDz=0;
   break;
 }
}
///////////////////////////////////////////////////////////////////////////
void Ali3Vector::GetErrors(Double_t* e,TString f)
{
// Provide errors according to reference frame f
 Int_t frame=0;
 if (f == "car") frame=1;
 if (f == "sph") frame=2;
 if (f == "cyl") frame=3;

 Double_t dr2,dtheta2,dphi2,rho,drho2;
 Double_t v[3]; 

 switch (frame)
 {
  case 1: // Cartesian coordinates
   e[0]=fDx;
   e[1]=fDy;
   e[2]=fDz;
   break;

  case 2: // Spherical coordinates
   GetVector(v,"car");
   if (fV) 
   {
    dr2=(pow((v[0]*fDx),2)+pow((v[1]*fDy),2)+pow((v[2]*fDz),2))/(fV*fV);
   }
   else
   {
    dr2=0;
   }
   if (v[2]-fV)
   {
    dtheta2=(v[2]*v[2]/(pow(fV,4)-pow(v[2],2)*pow(fV,2)))*dr2
            +pow(fDz,2)/(pow(fV,2)-pow(v[2],2));
   }
   else
   {
//    dr2=fDz*fDz;
    dtheta2=0;
   }
   if (v[0] || v[1])
   {
    dphi2=(pow((v[1]*fDx),2)+pow((v[0]*fDy),2))/(pow(v[0],2)+pow(v[1],2));
   }
   else
   {
    dphi2=0;
   }
   e[0]=sqrt(dr2);
   e[1]=sqrt(dtheta2);
   e[2]=sqrt(dphi2);
   break;

  case 3: // Cylindrical coordinates
   GetVector(v,"car");
   rho=fV*sin(fTheta);
   if (rho) 
   {
    drho2=(pow((v[0]*fDx),2)+pow((v[1]*fDy),2))/(rho*rho);
   }
   else
   {
    drho2=0;
   }
   if (v[0] || v[1])
   {
    dphi2=(pow((v[1]*fDx),2)+pow((v[0]*fDy),2))/(pow(v[0],2)+pow(v[1],2));
   }
   else
   {
    dphi2=0;
   }
   e[0]=sqrt(drho2);
   e[1]=sqrt(dphi2);
   e[2]=fDz;
   break;

  default: // Unsupported reference frame
   cout << "*Ali3Vector::GetErrors* Unsupported frame : " << f << endl
        << " Possible frames are 'car', 'sph' and 'cyl'." << endl; 
   for (Int_t i=0; i<3; i++)
   {
    e[i]=0;
   }
   break;
 }
}
///////////////////////////////////////////////////////////////////////////
void Ali3Vector::SetErrors(Float_t* e,TString f)
{
// Store errors according to reference frame f
// The error on scalar results is reset to 0
 Double_t vec[3];
 for (Int_t i=0; i<3; i++)
 {
  vec[i]=e[i];
 }
 SetErrors(vec,f);
}
///////////////////////////////////////////////////////////////////////////
void Ali3Vector::GetErrors(Float_t* e,TString f)
{
// Provide errors according to reference frame f
 Double_t vec[3];
 GetErrors(vec,f);
 for (Int_t i=0; i<3; i++)
 {
  e[i]=vec[i];
 }
}
///////////////////////////////////////////////////////////////////////////
void Ali3Vector::Info(TString f)
{
// Print vector components according to reference frame f
 if (f=="car" || f=="sph" || f=="cyl")
 {
  Double_t vec[3],err[3];
  GetVector(vec,f);
  GetErrors(err,f);
  cout << " Vector in " << f << " coordinates : "
       << vec[0] << " " << vec[1] << " " << vec[2] << endl; 
  cout << "   Err. in " << f << " coordinates : "
       << err[0] << " " << err[1] << " " << err[2] << endl; 
 }
 else
 {
  cout << " *Ali3Vector::Info* Unsupported frame : " << f << endl
       << "  Possible frames are 'car', 'sph' and 'cyl'." << endl; 
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t Ali3Vector::GetNorm()
{
// Provide the norm of the current vector
// The error on the scalar result (norm) is updated accordingly
 Double_t e[3];
 GetErrors(e,"sph");
 fDresult=e[0]; 
 return fV;
}
///////////////////////////////////////////////////////////////////////////
Double_t Ali3Vector::GetPseudoRapidity()
{
// Provide the pseudo-rapidity w.r.t. the z-axis.
// In other words : eta=-log(tan(theta/2))
// The error on the scalar result (pseudo-rap.) is updated accordingly
 Double_t v[3];
 GetVector(v,"sph");
 Double_t thetahalf=v[1]/2.;
 Double_t arg=tan(thetahalf);
 Double_t eta=0;
 if (arg>0) eta=-log(arg);
 Double_t e[3];
 GetErrors(e,"sph");
 Double_t prod=cos(thetahalf)*sin(thetahalf);
 fDresult=0;
 if (prod) fDresult=fabs(e[1]/2.*prod);
 return eta;
}
///////////////////////////////////////////////////////////////////////////
Double_t Ali3Vector::Dot(Ali3Vector& q)
{
// Provide the dot product of the current vector with vector q
// The error on the scalar result (dotproduct) is updated accordingly

 Double_t dotpro=0;

 if ((this) == &q) // Check for special case v.Dot(v)
 {
  Double_t norm=GetNorm();
  Double_t dnorm=GetResultError();
  dotpro=pow(norm,2);
  fDresult=2.*norm*dnorm;
 }
 else
 {
  Double_t a[3],b[3];
  Double_t ea[3],eb[3];
  Double_t d2=0;

  GetVector(a,"car");
  GetErrors(ea,"car");
  q.GetVector(b,"car");
  q.GetErrors(eb,"car");
  for (Int_t i=0; i<3; i++)
  {
   dotpro+=a[i]*b[i];
   d2+=pow(b[i]*ea[i],2)+pow(a[i]*eb[i],2);
  }
  fDresult=sqrt(d2);
 }

 return dotpro;
}
///////////////////////////////////////////////////////////////////////////
Double_t Ali3Vector::GetResultError()
{
// Provide the error on the result of an operation yielding a scalar
// E.g. GetNorm() or Dot()
 return fDresult;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector Ali3Vector::Cross(Ali3Vector& q)
{
// Provide the cross product of the current vector with vector q
// Error propagation is performed automatically
 Double_t a[3],b[3],c[3];
 Double_t ea[3],eb[3],ec[3],d2;

 GetVector(a,"car");
 GetErrors(ea,"car");
 q.GetVector(b,"car");
 q.GetErrors(eb,"car");

 c[0]=a[1]*b[2]-a[2]*b[1];
 c[1]=a[2]*b[0]-a[0]*b[2];
 c[2]=a[0]*b[1]-a[1]*b[0];

 d2=pow(b[2]*ea[1],2)+pow(a[1]*eb[2],2)
   +pow(b[1]*ea[2],2)+pow(a[2]*eb[1],2);
 ec[0]=sqrt(d2);

 d2=pow(b[0]*ea[2],2)+pow(a[2]*eb[0],2)
   +pow(b[2]*ea[0],2)+pow(a[0]*eb[2],2);
 ec[1]=sqrt(d2);

 d2=pow(b[1]*ea[0],2)+pow(a[0]*eb[1],2)
   +pow(b[0]*ea[1],2)+pow(a[1]*eb[0],2);
 ec[2]=sqrt(d2);

 Ali3Vector v;
 v.SetVector(c,"car");
 v.SetErrors(ec,"car");
  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector Ali3Vector::operator+(Ali3Vector& q)
{
// Add vector q to the current vector
// Error propagation is performed automatically
 Double_t a[3],b[3],ea[3],eb[3];

 GetVector(a,"car");
 GetErrors(ea,"car");
 q.GetVector(b,"car");
 q.GetErrors(eb,"car");

 for (Int_t i=0; i<3; i++)
 {
  a[i]+=b[i];
  ea[i]=sqrt(pow(ea[i],2)+pow(eb[i],2));
 }

 Ali3Vector v;
 v.SetVector(a,"car");
 v.SetErrors(ea,"car");
  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector Ali3Vector::operator-(Ali3Vector& q)
{
// Subtract vector q from the current vector
// Error propagation is performed automatically
 Double_t a[3],b[3],ea[3],eb[3];

 GetVector(a,"car");
 GetErrors(ea,"car");
 q.GetVector(b,"car");
 q.GetErrors(eb,"car");

 for (Int_t i=0; i<3; i++)
 {
  a[i]-=b[i];
  ea[i]=sqrt(pow(ea[i],2)+pow(eb[i],2));
 }

 Ali3Vector v;
 v.SetVector(a,"car");
 v.SetErrors(ea,"car");
  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector Ali3Vector::operator*(Double_t s)
{
// Multiply the current vector with a scalar s.
// Error propagation is performed automatically.
 Double_t a[3],ea[3];

 GetVector(a,"car");
 GetErrors(ea,"car");

 for (Int_t i=0; i<3; i++)
 {
  a[i]*=s;
  ea[i]*=s;
 }

 Ali3Vector v;
 v.SetVector(a,"car");
 v.SetErrors(ea,"car");
  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector Ali3Vector::operator/(Double_t s)
{
// Divide the current vector by a scalar s
// Error propagation is performed automatically

 if (fabs(s)<1.e-20) // Protect against division by 0
 {
  cout << " *Ali3Vector::/* Division by 0 detected. No action taken." << endl;
  return *this;
 }
 else
 {
  Double_t a[3],ea[3];

  GetVector(a,"car");
  GetErrors(ea,"car");

  for (Int_t i=0; i<3; i++)
  {
   a[i]/=s;
   ea[i]/=s;
  }

  Ali3Vector v;
  v.SetVector(a,"car");
  v.SetErrors(ea,"car");
  
  return v;
 }
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector& Ali3Vector::operator+=(Ali3Vector& q)
{
// Add vector q to the current vector
// Error propagation is performed automatically
 Double_t a[3],b[3],ea[3],eb[3];

 GetVector(a,"car");
 GetErrors(ea,"car");
 q.GetVector(b,"car");
 q.GetErrors(eb,"car");

 for (Int_t i=0; i<3; i++)
 {
  a[i]+=b[i];
  ea[i]=sqrt(pow(ea[i],2)+pow(eb[i],2));
 }

 SetVector(a,"car");
 SetErrors(ea,"car");
  
 return *this;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector& Ali3Vector::operator-=(Ali3Vector& q)
{
// Subtract vector q from the current vector
// Error propagation is performed automatically
 Double_t a[3],b[3],ea[3],eb[3];

 GetVector(a,"car");
 GetErrors(ea,"car");
 q.GetVector(b,"car");
 q.GetErrors(eb,"car");

 for (Int_t i=0; i<3; i++)
 {
  a[i]-=b[i];
  ea[i]=sqrt(pow(ea[i],2)+pow(eb[i],2));
 }

 SetVector(a,"car");
 SetErrors(ea,"car");
  
 return *this;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector& Ali3Vector::operator*=(Double_t s)
{
// Multiply the current vector with a scalar s
// Error propagation is performed automatically
 Double_t a[3],ea[3];

 GetVector(a,"car");
 GetErrors(ea,"car");

 for (Int_t i=0; i<3; i++)
 {
  a[i]*=s;
  ea[i]*=s;
 }

 SetVector(a,"car");
 SetErrors(ea,"car");
  
 return *this;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector& Ali3Vector::operator/=(Double_t s)
{
// Divide the current vector by a scalar s
// Error propagation is performed automatically

 if (fabs(s)<1.e-20) // Protect against division by 0
 {
  cout << " *Ali3Vector::/=* Division by 0 detected. No action taken." << endl;
  return *this;
 }
 else
 {
  Double_t a[3],ea[3];

  GetVector(a,"car");
  GetErrors(ea,"car");

  for (Int_t i=0; i<3; i++)
  {
   a[i]/=s;
   ea[i]/=s;
  }

  SetVector(a,"car");
  SetErrors(ea,"car");
  
  return *this;
 }
}
///////////////////////////////////////////////////////////////////////////
