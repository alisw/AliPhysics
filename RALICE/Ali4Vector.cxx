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

#include "Ali4Vector.h"
 
ClassImp(Ali4Vector) // Class implementation to enable ROOT I/O
 
Ali4Vector::Ali4Vector()
{
// Creation of a contravariant 4-vector and initialisation of parameters
 fV0=0;
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
// Store contravariant vector
 fV0=v0;
 fV=v;
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::SetVector(Double_t* v,TString f)
{
// Store vector according to reference frame f
 fV0=v[0];
 Double_t a[3];
 for (Int_t i=0; i<3; i++)
 {
  a[i]=v[i+1];
 }
 fV.SetVector(a,f);
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::GetVector(Double_t* v,TString f)
{
// Provide vector according to reference frame f
 v[0]=fV0;
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
// Store vector according to reference frame f
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
// Provide vector according to reference frame f
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
// Provide the scalar part
 return fV0;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector Ali4Vector::Get3Vector()
{
// Provide the 3-vector part
 return fV;
}
///////////////////////////////////////////////////////////////////////////
void Ali4Vector::Info(TString f)
{
// Print contravariant vector components according to reference frame f
 if (f=="car" || f=="sph" || f=="cyl")
 {
  Double_t vec[4];
  GetVector(vec,f);
  cout << " Contravariant vector in " << f << " coordinates : "
       << vec[0] << " " << vec[1] << " " << vec[2] << " " << vec[3] << endl; 
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
 Double_t a[4],b[4];
 Double_t dotpro;

 GetVector(a,"car");
 q.GetVector(b,"car");

 dotpro=a[0]*b[0];
 for (Int_t i=1; i<4; i++)
 {
  dotpro-=a[i]*b[i];
 }
 
 return dotpro;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector Ali4Vector::operator+(Ali4Vector& q)
{
// Add vector q to the current vector
 Double_t a[4],b[4];

 GetVector(a,"car");
 q.GetVector(b,"car");

 for (Int_t i=0; i<4; i++)
 {
  a[i]+=b[i];
 }

 Ali4Vector v;
 v.SetVector(a,"car");
  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector Ali4Vector::operator-(Ali4Vector& q)
{
// Subtract vector q from the current vector
 Double_t a[4],b[4];

 GetVector(a,"car");
 q.GetVector(b,"car");

 for (Int_t i=0; i<4; i++)
 {
  a[i]-=b[i];
 }

 Ali4Vector v;
 v.SetVector(a,"car");
  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector Ali4Vector::operator*(Double_t s)
{
// Multiply the current vector with a scalar s
 Double_t a[4];

 GetVector(a,"car");

 for (Int_t i=0; i<4; i++)
 {
  a[i]*=s;
 }

 Ali4Vector v;
 v.SetVector(a,"car");
  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector Ali4Vector::operator/(Double_t s)
{
// Divide the current vector by a scalar s

 if (fabs(s)<1.e-20) // Protect against division by 0
 {
  cout << " *Ali4Vector::/* Division by 0 detected. No action taken." << endl;
  return *this;
 }
 else
 {
  Double_t a[4];

  GetVector(a,"car");

  for (Int_t i=0; i<4; i++)
  {
   a[i]/=s;
  }

  Ali4Vector v;
  v.SetVector(a,"car");
  
  return v;
 }
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector& Ali4Vector::operator+=(Ali4Vector& q)
{
// Add vector q to the current vector
 Double_t a[4],b[4];

 GetVector(a,"car");
 q.GetVector(b,"car");

 for (Int_t i=0; i<4; i++)
 {
  a[i]+=b[i];
 }

 SetVector(a,"car");
  
 return *this;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector& Ali4Vector::operator-=(Ali4Vector& q)
{
// Subtract vector q from the current vector
 Double_t a[4],b[4];

 GetVector(a,"car");
 q.GetVector(b,"car");

 for (Int_t i=0; i<4; i++)
 {
  a[i]-=b[i];
 }

 SetVector(a,"car");
  
 return *this;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector& Ali4Vector::operator*=(Double_t s)
{
// Multiply the current vector with a scalar s
 Double_t a[4];

 GetVector(a,"car");

 for (Int_t i=0; i<4; i++)
 {
  a[i]*=s;
 }

 SetVector(a,"car");
  
 return *this;
}
///////////////////////////////////////////////////////////////////////////
Ali4Vector& Ali4Vector::operator/=(Double_t s)
{
// Divide the current vector by a scalar s

 if (fabs(s)<1.e-20) // Protect against division by 0
 {
  cout << " *Ali4Vector::/=* Division by 0 detected. No action taken." << endl;
  return *this;
 }
 else
 {
  Double_t a[4];

  GetVector(a,"car");

  for (Int_t i=0; i<4; i++)
  {
   a[i]/=s;
  }

  SetVector(a,"car");
  
  return *this;
 }
}
///////////////////////////////////////////////////////////////////////////
