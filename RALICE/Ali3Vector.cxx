#include "Ali3Vector.h"
 
ClassImp(Ali3Vector) // Class implementation to enable ROOT I/O
 
Ali3Vector::Ali3Vector()
{
// Creation of an Ali3Vector object and initialisation of parameters
 fV=0;
 fTheta=0;
 fPhi=0;
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
void Ali3Vector::Info(TString f)
{
// Print vector components according to reference frame f
 if (f=="car" || f=="sph" || f=="cyl")
 {
  Double_t vec[3];
  GetVector(vec,f);
  cout << " Vector in " << f << " coordinates : "
       << vec[0] << " " << vec[1] << " " << vec[2] << endl; 
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
 return fV;
}
///////////////////////////////////////////////////////////////////////////
Double_t Ali3Vector::Dot(Ali3Vector& q)
{
// Provide the dot product of the current vector with vector q
 Double_t a[3],b[3];
 Double_t dotpro=0;

 GetVector(a,"car");
 q.GetVector(b,"car");
 for (Int_t i=0; i<3; i++)
 {
  dotpro+=a[i]*b[i];
 }
 
 return dotpro;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector Ali3Vector::Cross(Ali3Vector& q)
{
// Provide the cross product of the current vector with vector q
 Double_t a[3],b[3],c[3];

 GetVector(a,"car");
 q.GetVector(b,"car");

 c[0]=a[1]*b[2]-a[2]*b[1];
 c[1]=a[2]*b[0]-a[0]*b[2];
 c[2]=a[0]*b[1]-a[1]*b[0];

 Ali3Vector v;
 v.SetVector(c,"car");
  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector Ali3Vector::operator+(Ali3Vector& q)
{
// Add vector q to the current vector
 Double_t a[3],b[3];

 GetVector(a,"car");
 q.GetVector(b,"car");

 for (Int_t i=0; i<3; i++)
 {
  a[i]+=b[i];
 }

 Ali3Vector v;
 v.SetVector(a,"car");
  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector Ali3Vector::operator-(Ali3Vector& q)
{
// Subtract vector q from the current vector
 Double_t a[3],b[3];

 GetVector(a,"car");
 q.GetVector(b,"car");

 for (Int_t i=0; i<3; i++)
 {
  a[i]-=b[i];
 }

 Ali3Vector v;
 v.SetVector(a,"car");
  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector Ali3Vector::operator*(Double_t s)
{
// Multiply the current vector with a scalar s
 Double_t a[3];

 GetVector(a,"car");

 for (Int_t i=0; i<3; i++)
 {
  a[i]*=s;
 }

 Ali3Vector v;
 v.SetVector(a,"car");
  
 return v;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector Ali3Vector::operator/(Double_t s)
{
// Divide the current vector by a scalar s

 if (fabs(s)<1.e-20) // Protect against division by 0
 {
  cout << " *Ali3Vector::/* Division by 0 detected. No action taken." << endl;
  return *this;
 }
 else
 {
  Double_t a[3];

  GetVector(a,"car");

  for (Int_t i=0; i<3; i++)
  {
   a[i]/=s;
  }

  Ali3Vector v;
  v.SetVector(a,"car");
  
  return v;
 }
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector& Ali3Vector::operator+=(Ali3Vector& q)
{
// Add vector q to the current vector
 Double_t a[3],b[3];

 GetVector(a,"car");
 q.GetVector(b,"car");

 for (Int_t i=0; i<3; i++)
 {
  a[i]+=b[i];
 }

 SetVector(a,"car");
  
 return *this;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector& Ali3Vector::operator-=(Ali3Vector& q)
{
// Subtract vector q from the current vector
 Double_t a[3],b[3];

 GetVector(a,"car");
 q.GetVector(b,"car");

 for (Int_t i=0; i<3; i++)
 {
  a[i]-=b[i];
 }

 SetVector(a,"car");
  
 return *this;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector& Ali3Vector::operator*=(Double_t s)
{
// Multiply the current vector with a scalar s
 Double_t a[3];

 GetVector(a,"car");

 for (Int_t i=0; i<3; i++)
 {
  a[i]*=s;
 }

 SetVector(a,"car");
  
 return *this;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector& Ali3Vector::operator/=(Double_t s)
{
// Divide the current vector by a scalar s

 if (fabs(s)<1.e-20) // Protect against division by 0
 {
  cout << " *Ali3Vector::/=* Division by 0 detected. No action taken." << endl;
  return *this;
 }
 else
 {
  Double_t a[3];

  GetVector(a,"car");

  for (Int_t i=0; i<3; i++)
  {
   a[i]/=s;
  }

  SetVector(a,"car");
  
  return *this;
 }
}
///////////////////////////////////////////////////////////////////////////
