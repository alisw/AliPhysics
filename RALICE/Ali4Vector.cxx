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
// Class Ali4Vector
// Handling of Lorentz 4-vectors in various reference frames.
//
// This class is meant to serve as a base class for ALICE objects
// that have Lorentz 4-vector characteristics.
// Error propagation is performed automatically.
//
// All 4-vectors are treated in the contravariant form and the convention
// for the metric and the 4-vector components is according to the one
// used in the book "Classical Electrodynamics" by J.D. Jackson.
//
// A 4-vector is said to have a scalar part and a 3-vector part,
// which is indicated by the notation
//
//    x^i = (x^0,x^1,x^2,x^3)
// 
// The scalar part   = x^0
// The 3-vector part = (x^1,x^2,x^3)
//
// In view of accuracy and the fact that e.g. particle identity (mass)
// is preserved in many physics processes, the Lorentz invariant 
// (x^i*x_i) is internally saved together with the scalar part.
//
// This allows the following two modes of functionality :
//
// Scalar mode    : The scalar part and the 3-vector part are considered
//                  as basic quantities and the invariant with its error
//                  is derived from these.  
// Invariant mode : The invariant and the 3-vector part are considered
//                  as basic quantities and the scalar with its error
//                  is derived from these.
//
// The philosophy followed here is the following :
// ===============================================
//
// 1) Invokation of SetVector() sets the scalar and 3-vector parts 
//    and the invariant is calculated from these.
//    Automatically the scalar mode is selected and invokation of
//    SetErrors() will calculate the error on the invariant.
//
// 2) In case the scalar part is modified via SetScalar(), scalar mode is
//    automatically selected and the Lorentz invariant (x^i*x_i) and its
//    error are updated accordingly.
//    The 3-vector part is NOT modified.
//    This situation arises when one e.g. precisely determines the time
//    or energy (x^0).     
//
// 3) In case the Lorentz invariant (x^i*x_i) is modified via SetInvariant(),
//    invariant mode is selected automatically and the scalar part and its
//    error are updated accordingly.
//    The 3-vector part is NOT modified.
//    This situation arises when one e.g. precisely determines the mass.     
//
// 4) In case the vector part is modified via Set3Vector(), then the 
//    current mode determines whether the scalar or the invariant is updated. 
//    Scalar mode    : The Lorentz invariant (x^i*x_i) and its error are updated;
//                     the scalar part and its error are NOT modified. 
//                     This situation arises when one e.g. improves the 3-position
//                     vector for a particle with a very precise timing.     
//    Invariant mode : The scalar part and its error are updated;
//                     the Lorentz invariant (x^i*x_i) and its error are NOT modified.
//                     This situation arises when one e.g. improves the 3-momentum
//                     vector for a particle with known mass.     
//
// The dotproduct is defined such that p.Dot(p) yields the Lorentz invariant
// scalar of the 4-vector p (i.e. m**2 in case p is a 4-momentum).   
//
// Note :
// ------
// Vectors (v), Errors (e) and reference frames (f) are specified via
// SetVector(Float_t* v,TString f)
// SetErrors(Float_t* e,TString f)
// under the following conventions :
//
// f="car" ==> 3-vector part of v in Cartesian coordinates   (x,y,z)
// f="sph" ==> 3-vector part of v in Spherical coordinates   (r,theta,phi)
// f="cyl" ==> 3-vector part of v in Cylindrical coordinates (rho,phi,z)
//
// All angles are in radians.
//
// Example :
// ---------
//
// Ali4Vector a;
//
// Float_t v[4]={25,-1,3,7};
// a.SetVector(v,"car");
//
// Float_t vec[4];
// a.GetVector(vec,"sph");
//
// Ali4Vector b;
// Float_t v2[4]={33,6,-18,2};
// b.SetVector(v2,"car");
//
// Float_t dotpro=a.Dot(b);
//
// Float_t x0=16;
// Ali3Vector x;
// Float_t vec2[3]={1,2,3};
// x.SetVector(vec2,"car");
//
// Ali4Vector c;
// c.SetVector(x0,x);
// c.GetVector(vec,"car");
// c.Info("cyl");
// c=a+b;
// c=a-b;
// c=a*5;
//
//--- Author: Nick van Eijndhoven 01-apr-1999 UU-SAP Utrecht
//- Modified: NvE 15-oct-1999 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "Ali4Vector.h"
 
ClassImp(Ali4Vector) // Class implementation to enable ROOT I/O
 
Ali4Vector::Ali4Vector()
{
// Creation of a contravariant 4-vector and initialisation of parameters.
// All values are initialised to 0.
// Scalar mode is initially selected.
 fScalar=1;
 fV2=0;
 fDv2=0;
 fV0=0;
 fDv0=0;
 fDresult=0;
 Double_t a[3]={0,0,0};
 fV.SetVector(a,"sph");
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector::~Ali4Vector()
{
// Destructor to delete dynamically allocated memory
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::SetVector(Double_t v0,Ali3Vector v)
{
// Store contravariant vector.
// The error on the scalar part is initialised to 0.
// The errors on the vector part are taken from the input Ali3Vector.
// Scalar mode is automatically selected. 
// The error on scalar result operations is reset to 0.
 fScalar=1;
 fV0=v0;
 fV=v;
 fV2=pow(v0,2)-fV.Dot(fV);
 SetScalarError(0);
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::SetVector(Double_t* v,TString f)
{
// Store vector according to reference frame f.
// All errors are initialised to 0.
// Scalar mode is automatically selected. 
// The error on scalar result operations is reset to 0.
 fScalar=1;
 Double_t a[3];
 for (Int_t i=0; i<3; i++)
 {
  a[i]=v[i+1];
 }
 fV0=v[0];
 fV.SetVector(a,f);
 fV2=pow(fV0,2)-fV.Dot(fV);
 fDv2=0;
 fDv0=0;
 fDresult=0;
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::GetVector(Double_t* v,TString f)
{
// Provide 4-vector components according to reference frame f
// and according to the current mode.
// Scalar mode    : The scalar part is directly returned via v[0].
// Invariant mode : The scalar part is re-calculated via the value
//                  of the Lorentz invariant and then returned via v[0].
 if (fScalar)
 {
  v[0]=fV0;
 }
 else
 {
  v[0]=sqrt(fV.Dot(fV)+fV2);
 } 
 Double_t a[3];
 fV.GetVector(a,f);
 for (Int_t i=0; i<3; i++)
 {
  v[i+1]=a[i];
 }
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::SetVector(Float_t* v,TString f)
{
// Store vector according to reference frame f.
// All errors are initialised to 0.
// Scalar mode is automatically selected. 
// The error on scalar result operations is reset to 0.
 Double_t vec[4];
 for (Int_t i=0; i<4; i++)
 {
  vec[i]=v[i];
 }
 SetVector(vec,f);
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::GetVector(Float_t* v,TString f)
{
// Provide 4-vector components according to reference frame f
// and according to the current mode.
// Scalar mode    : The scalar part is directly returned via v[0].
// Invariant mode : The scalar part is re-calculated via the value
//                  of the Lorentz invariant and then returned via v[0].
 Double_t vec[4];
 GetVector(vec,f);
 for (Int_t i=0; i<4; i++)
 {
  v[i]=vec[i];
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t Ali4Vector::GetScalar()
{
// Provide the scalar part.
// The error on the scalar value is available via GetResultError()
// after invokation of GetScalar().
 if (fScalar)
 {
  fDresult=fDv0;
  return fV0;
 }
 else
 {
  Double_t dot=fV.Dot(fV);
  Double_t ddot=fV.GetResultError();
  Double_t v02=dot+fV2;
  Double_t dv02=sqrt(pow(ddot,2)+pow(fDv2,2));
  Double_t v0=sqrt(fabs(v02));
  Double_t dv0=0;
  if (v0) dv0=dv02/(2.*v0);
  fDresult=dv0;
  return v0;
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t Ali4Vector::GetResultError()
{
// Provide the error on the result of an operation yielding a scalar
// E.g. GetScalar(), GetInvariant() or Dot()
 return fDresult;
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::SetScalar(Double_t v0,Double_t dv0)
{
// Modify the scalar part (v0) and its error (dv0).
// The default value for dv0 is 0.
// The vector part is not modified. 
// Scalar mode is automatically selected
// ==> Lorentz invariant and its error are updated. 
// The error on scalar result operations is reset to 0.
 fScalar=1;
 fV0=v0;
 fV2=pow(v0,2)-fV.Dot(fV);
 SetScalarError(dv0);
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::SetScalarError(Double_t dv0)
{
// Set the error on the scalar part.
// If in scalar mode, update error on the invariant accordingly.
// The error on scalar result operations is reset to 0.
 fDv0=dv0;
 if (fScalar)
 {
  Double_t norm=fV.GetNorm();
  Double_t dnorm=fV.GetResultError();
  fDv2=sqrt(pow(2.*fV0*fDv0,2)+pow(2.*norm*dnorm,2));
 } 
 fDresult=0;
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::Set3Vector(Ali3Vector v)
{
// Set the 3-vector part, the errors are taken from the input Ali3Vector
// Scalar mode    : The scalar part and its error are not modified,
//                  the Lorentz invariantand its error are re-calculated.
// Invariant mode : The Lorentz invariant and its error are not modified,
//                  the scalar part and its error are re-calculated.
// The error on scalar result operations is reset to 0.
 fV=v;
 if (fScalar)
 {
  SetScalar(fV0,fDv0);
 }
 else
 {
  SetInvariant(fV2,fDv2);
 }
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::Set3Vector(Double_t* v,TString f)
{
// Set the 3-vector part according to reference frame f
// The errors on the vector part are initialised to 0
// Scalar mode    : The scalar part and its error are not modified,
//                  the Lorentz invariantand its error are re-calculated.
// Invariant mode : The Lorentz invariant and its error are not modified,
//                  the scalar part and its error are re-calculated.
// The error on scalar result operations is reset to 0.
 Double_t a[3];
 for (Int_t i=0; i<3; i++)
 {
  a[i]=v[i];
 }
 fV.SetVector(a,f);

 if (fScalar)
 {
  SetScalar(fV0,fDv0);
 }
 else
 {
  SetInvariant(fV2,fDv2);
 }
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::Set3Vector(Float_t* v,TString f)
{
// Set the 3-vector part according to reference frame f
// The errors on the vector part are initialised to 0
// The Lorentz invariant is not modified
// The error on scalar result operations is reset to 0.
 Double_t vec[3];
 for (Int_t i=0; i<3; i++)
 {
  vec[i]=v[i];
 }
 Set3Vector(vec,f);
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::SetInvariant(Double_t v2,Double_t dv2)
{
// Modify the Lorentz invariant (v2) quantity v^i*v_i and its error (dv2).
// The default value for the error dv2 is 0.
// The vector part is not modified.
// Invariant mode is automatically selected
// ==> the scalar part and its error are updated.
// The error on scalar result operations is reset to 0.
//
 fScalar=0;
 fV2=v2;
 fDv2=dv2;
 fV0=GetScalar();
 fDv0=GetResultError();
 fDresult=0;
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::SetInvariantError(Double_t dv2)
{
// Set the error on the Lorentz invariant.
// If in invariant mode, update error on the scalar part accordingly.
// The error on scalar result operations is reset to 0.
 fDv2=dv2;
 if (!fScalar)
 {
  fDv0=GetResultError();
 }
 fDresult=0; 
}
///////////////////////////////////////////////////////////////////////////
Double_t Ali4Vector::GetInvariant()
{
// Provide the Lorentz invariant v^i*v_i.
// The error on the Lorentz invariant is available via GetResultError()
// after invokation of GetInvariant().
 if (!fScalar)
 {
  fDresult=fDv2;
  return fV2;
 }
 else
 {
  Double_t inv=Dot(*this);
  return inv;
 }
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector Ali4Vector::Get3Vector()
{
// Provide the 3-vector part
 return fV;
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::SetErrors(Double_t* e,TString f)
{
// Store errors for vector v^i according to reference frame f
// If in scalar mode, update error on the invariant accordingly.
// The error on scalar result operations is reset to 0.
 Double_t a[3];
 for (Int_t i=0; i<3; i++)
 {
  a[i]=e[i+1];
 }
 SetScalarError(e[0]);
 fV.SetErrors(a,f);
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::SetErrors(Float_t* e,TString f)
{
// Store errors for vector v^i according to reference frame f
// If in scalar mode, update error on the invariant accordingly.
// The error on scalar result operations is reset to 0.
 Double_t a[4];
 for (Int_t i=0; i<4; i++)
 {
  a[i]=e[i];
 }
 SetErrors(a,f);
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::GetErrors(Double_t* e,TString f)
{
// Provide errors for vector v^i according to reference frame f
// and according to the current mode.
// Scalar mode    : The error on the scalar part is directly returned via e[0].
// Invariant mode : The error on the scalar part is re-calculated via the error
//                  value on the Lorentz invariant and then returned via e[0].
 Double_t a[3];
 fV.GetErrors(a,f);

 e[0]=GetResultError();
 for (Int_t i=0; i<3; i++)
 {
  e[i+1]=a[i];
 }
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::GetErrors(Float_t* e,TString f)
{
// Provide errors for vector v^i according to reference frame f
// and according to the current mode.
// Scalar mode    : The error on the scalar part is directly returned via e[0].
// Invariant mode : The error on the scalar part is re-calculated via the error
//                  value on the Lorentz invariant and then returned via e[0].
 Double_t a[4];
 GetErrors(a,f);
 for (Int_t i=0; i<4; i++)
 {
  e[i]=a[i];
 }
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::Info(TString f)
{
// Print contravariant vector components and errors according to
// reference frame f and according to the current mode.
// Scalar mode    : The scalar part and its error are directly returned.
// Invariant mode : The scalar part and its error are re-calculated via the
//                  value (and error) of the Lorentz invariant.

 if (f=="car" || f=="sph" || f=="cyl")
 {
  Double_t vec[4],err[4];
  GetVector(vec,f);
  GetErrors(err,f);
  Double_t inv=GetInvariant();
  Double_t dinv=GetResultError();
  cout << " Contravariant vector in " << f << " coordinates : "
       << vec[0] << " " << vec[1] << " " << vec[2] << " " << vec[3] << endl; 
  cout << " ------------- Errors in " << f << " coordinates : "
       << err[0] << " " << err[1] << " " << err[2] << " " << err[3] << endl; 
  cout << " --- Lorentz invariant (v^i*v_i) : " << inv << " error : " << dinv << endl;
 }
 else
 {
  cout << " *Ali4Vector::Info* Unsupported frame : " << f << endl
       << "  Possible frames are 'car', 'sph' and 'cyl'." << endl; 
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t Ali4Vector::Dot(Ali4Vector& q)
{
// Provide the dot product of the current vector with vector q
 Double_t dotpro=0;
 Double_t a0=GetScalar();
 Double_t da0=GetResultError();
 if ((this) == &q) // Check for special case v.Dot(v)
 {
  Double_t norm=fV.GetNorm();
  Double_t dnorm=fV.GetResultError();
  dotpro=pow(a0,2)-pow(norm,2);
  fDresult=sqrt(pow(2.*a0*da0,2)+pow(2.*norm*dnorm,2));
 }
 else
 {
  Double_t b0=q.GetScalar();
  Double_t db0=q.GetResultError();
  Ali3Vector b=q.Get3Vector();

  Double_t dot=fV.Dot(b);
  Double_t ddot=fV.GetResultError();

  dotpro=a0*b0-dot;

  fDresult=sqrt(pow(b0*da0,2)+pow(a0*db0,2)+pow(ddot,2));
 }

 return dotpro;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector Ali4Vector::operator+(Ali4Vector& q)
{
// Add 4-vector q to the current 4-vector
// Error propagation is performed automatically
 Double_t a0=GetScalar();
 Double_t da0=GetResultError();
 Ali3Vector a=Get3Vector();
 Double_t b0=q.GetScalar();
 Double_t db0=q.GetResultError();
 Ali3Vector b=q.Get3Vector();

 Double_t c0=a0+b0;
 Ali3Vector c=a+b;
 Double_t dc0=sqrt(pow(da0,2)+pow(db0,2));

 Ali4Vector v;
 v.SetVector(c0,c);
 v.SetScalarError(dc0);  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector Ali4Vector::operator-(Ali4Vector& q)
{
// Subtract 4-vector q from the current 4-vector
// Error propagation is performed automatically
 Double_t a0=GetScalar();
 Double_t da0=GetResultError();
 Ali3Vector a=Get3Vector();
 Double_t b0=q.GetScalar();
 Double_t db0=q.GetResultError();
 Ali3Vector b=q.Get3Vector();

 Double_t c0=a0-b0;
 Ali3Vector c=a-b;
 Double_t dc0=sqrt(pow(da0,2)+pow(db0,2));

 Ali4Vector v;
 v.SetVector(c0,c);
 v.SetScalarError(dc0);  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector Ali4Vector::operator*(Double_t s)
{
// Multiply the current 4-vector with a scalar s
// Error propagation is performed automatically
 Double_t a0=GetScalar();
 Double_t da0=GetResultError();
 Ali3Vector a=Get3Vector();

 a0*=s;
 a*=s;
 da0*=s;

 Ali4Vector v;
 v.SetVector(a0,a);
 v.SetScalarError(da0);  
  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector Ali4Vector::operator/(Double_t s)
{
// Divide the current vector by a scalar s
// Error propagation is performed automatically

 if (fabs(s)<1.e-20) // Protect against division by 0
 {
  cout << " *Ali4Vector::/* Division by 0 detected. No action taken." << endl;
  return *this;
 }
 else
 {
  Double_t a0=GetScalar();
  Double_t da0=GetResultError();
  Ali3Vector a=Get3Vector();

  a0/=s;
  a/=s;
  da0/=s;

  Ali4Vector v;
  v.SetVector(a0,a);
  v.SetScalarError(da0);  
  
  return v;
 }
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector& Ali4Vector::operator+=(Ali4Vector& q)
{
// Add 4-vector q to the current 4-vector
// Error propagation is performed automatically
 Double_t a0=GetScalar();
 Double_t da0=GetResultError();
 Ali3Vector a=Get3Vector();
 Double_t b0=q.GetScalar();
 Double_t db0=q.GetResultError();
 Ali3Vector b=q.Get3Vector();

 Double_t c0=a0+b0;
 Ali3Vector c=a+b;
 Double_t dc0=sqrt(pow(da0,2)+pow(db0,2));

 SetVector(c0,c);
 SetScalarError(dc0);  
  
 return *this;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector& Ali4Vector::operator-=(Ali4Vector& q)
{
// Subtract 4-vector q from the current 4-vector
// Error propagation is performed automatically
 Double_t a0=GetScalar();
 Double_t da0=GetResultError();
 Ali3Vector a=Get3Vector();
 Double_t b0=q.GetScalar();
 Double_t db0=q.GetResultError();
 Ali3Vector b=q.Get3Vector();

 Double_t c0=a0-b0;
 Ali3Vector c=a-b;
 Double_t dc0=sqrt(pow(da0,2)+pow(db0,2));

 SetVector(c0,c);
 SetScalarError(dc0);  

 return *this;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector& Ali4Vector::operator*=(Double_t s)
{
// Multiply the current 4-vector with a scalar s
// Error propagation is performed automatically
 Double_t a0=GetScalar();
 Double_t da0=GetResultError();
 Ali3Vector a=Get3Vector();

 a0*=s;
 a*=s;
 da0*=s;

 SetVector(a0,a);
 SetScalarError(da0);  

 return *this;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector& Ali4Vector::operator/=(Double_t s)
{
// Divide the current vector by a scalar s
// Error propagation is performed automatically

 if (fabs(s)<1.e-20) // Protect against division by 0
 {
  cout << " *Ali4Vector::/* Division by 0 detected. No action taken." << endl;
  return *this;
 }
 else
 {
  Double_t a0=GetScalar();
  Double_t da0=GetResultError();
  Ali3Vector a=Get3Vector();

  a0/=s;
  a/=s;
  da0/=s;

  SetVector(a0,a);
  SetScalarError(da0);  
  
  return *this;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t Ali4Vector::GetScalarFlag()
{
// Provide the value of the fScalar flag (for internal use only).
 return fScalar;
}
///////////////////////////////////////////////////////////////////////////
