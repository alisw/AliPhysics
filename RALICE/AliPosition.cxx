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

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class AliPosition
// Handling of positions in various reference frames.
//
// This class is meant to serve as a base class for ALICE objects
// that have a unique position in 3-dimensional space.
//
// Note :
// ------
// Positions (r), errors (e), reference frames (f) and angular units (u)
// are specified via
//
//    SetPosition(Float_t* r,TString f,TString u)
//    SetPositionErrors(Float_t* e,TString f,TString u)
//
// under the following conventions :
//
// f="car" ==> r in Cartesian coordinates   (x,y,z)
// f="sph" ==> r in Spherical coordinates   (r,theta,phi)
// f="cyl" ==> r in Cylindrical coordinates (rho,phi,z)
//
// u="rad" ==> angles in radians
// u="deg" ==> angles in degrees
//
// The "f" and "u" facilities only serve as a convenient user interface.
// Internally the actual storage of the various components is performed
// in a unique way. This allows setting/retrieval of vector components in a
// user selected frame/unit convention at any time. 
//
// The metric unit scale for the coordinates can be defined by the user
// via the SetUnitScale() memberfunction.
// This enables standardised expressions using numerical values of
// physical constants by means of the GetUnitScale() memberfunction.
// By default the unit scale is set to meter, corresponding to invokation
// of SetUnitScale(1).
// The user can specify a certain required metric unit scale in retreiving
// position components and/or distances.
// Please refer to the corresponding member functions for further details.
//   
//
// Example :
// ---------
//
// AliPosition q;
// Float_t pos[3]={-1,25,7};
// Float_t err[3]={0.08,1.85,0.5};
// q.SetPosition(pos,"car");
// q.SetPositionErrors(pos,"car");
// Float_t loc[3],dloc[3];
// q.GetPosition(loc,"sph","deg");
// q.GetPositionErrors(dloc,"sph","deg");
//
//--- Author: Nick van Eijndhoven 06-feb-1999 UU-SAP Utrecht
//- Modified: NvE $Date$ UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliPosition.h"
#include "Riostream.h"
 
ClassImp(AliPosition) // Class implementation to enable ROOT I/O
 
AliPosition::AliPosition()
{
// Creation of an AliPosition object and initialisation of parameters.
// The unit scale for position coordinates is initialised to cm.
 fScale=1;
 fTstamp=0;
}
///////////////////////////////////////////////////////////////////////////
AliPosition::~AliPosition()
{
// Destructor to delete dynamically allocated memory
 if (fTstamp)
 {
  delete fTstamp;
  fTstamp=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliPosition::AliPosition(const AliPosition& p) : Ali3Vector(p)
{
// Copy constructor
 fScale=p.fScale;
 fTstamp=0;
 if (p.fTstamp) fTstamp=new AliTimestamp(*(p.fTstamp));
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::SetPosition(Double_t* r,TString f,TString u)
{
// Store position according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".

 SetVector(r,f,u);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::GetPosition(Double_t* r,TString f,TString u,Float_t scale) const
{
// Provide position according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".
//
// By default the coordinates will be provided in the metric unit scale as
// stored in the AliPosition object.
// However, the user can select a different metric unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the position coordinates in cm.

 Ali3Vector v=(Ali3Vector)(*this);
 if (scale>0) v*=fScale/scale;
 v.GetVector(r,f,u);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::SetPosition(Float_t* r,TString f,TString u)
{
// Store position according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".

 SetVector(r,f,u);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::GetPosition(Float_t* r,TString f,TString u,Float_t scale) const
{
// Provide position according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".
//
// By default the coordinates will be provided in the metric unit scale as
// stored in the AliPosition object.
// However, the user can select a different metric unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the position coordinates in cm.

 Ali3Vector v=(Ali3Vector)(*this);
 if (scale>0) v*=fScale/scale;
 v.GetVector(r,f,u);
}
///////////////////////////////////////////////////////////////////////////
AliPosition& AliPosition::GetPosition()
{
// Provide position
 return (*this);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::SetPosition(Ali3Vector& r)
{
// Set position
 Double_t a[3];
 r.GetVector(a,"sph");
 SetVector(a,"sph");
 r.GetErrors(a,"car");
 SetErrors(a,"car");
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::SetPositionErrors(Double_t* r,TString f,TString u)
{
// Store position errors according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".

 SetErrors(r,f,u);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::GetPositionErrors(Double_t* r,TString f,TString u,Float_t scale) const
{
// Provide position errors according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".
//
// By default the coordinate errors will be provided in the metric unit scale as
// stored in the AliPosition object.
// However, the user can select a different metric unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the position coordinate errors in cm.

 Ali3Vector v=(Ali3Vector)(*this);
 if (scale>0) v*=fScale/scale;
 v.GetErrors(r,f,u);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::SetPositionErrors(Float_t* r,TString f,TString u)
{
// Store position errors according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".

 SetErrors(r,f,u);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::GetPositionErrors(Float_t* r,TString f,TString u,Float_t scale) const
{
// Provide position errors according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".
//
// By default the coordinate errors will be provided in the metric unit scale as
// stored in the AliPosition object.
// However, the user can select a different metric unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the position coordinate errors in cm.

 Ali3Vector v=(Ali3Vector)(*this);
 if (scale>0) v*=fScale/scale;
 v.GetErrors(r,f,u);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::ResetPosition()
{
// Reset the position and corresponding errors to 0.
 Double_t r[3]={0,0,0};
 SetVector(r,"sph");
 SetErrors(r,"car");
}
///////////////////////////////////////////////////////////////////////////
Double_t AliPosition::GetDistance(AliPosition& p,Float_t scale)
{
// Provide distance of the current AliPosition to position p.
// The error on the result can be obtained as usual by invoking
// GetResultError() afterwards. 
//
// By default the distance will be provided in the metric unit scale of
// the current AliPosition.
// This implies that the results of r1.GetDistance(r2) and r2.GetDistance(r1)
// may be numerically different in case r1 and r2 have different metric units.
// However, the user can specify a required metric unit scale by specification
// of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the distance in cm.
// As such it is possible to obtain a correctly computed distance even in case
// the position coordinates have a different unit scale.
// However, it is recommended to work always with one single unit scale.
//
 Ali3Vector d=(Ali3Vector)p;
 Float_t pscale=p.GetUnitScale();
 if ((pscale/fScale > 1.1) || (fScale/pscale > 1.1)) d=d*(pscale/fScale);
 Ali3Vector q=(Ali3Vector)(*this);
 d=d-q;
 Double_t dist=d.GetNorm();
 fDresult=d.GetResultError();

 if (scale>0)
 {
  dist*=fScale/scale;
  fDresult*=fScale/scale;
 }
 return dist;
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::SetUnitScale(Float_t s)
{
// Set the unit scale for the position coordinates.
// The scale is normalised w.r.t. the meter, so setting the unit scale
// to 0.01 means that all position coordinates are in cm.
// By default the unit scale is set to meter in the AliPosition constructor.
// It is recommended to use one single unit scale throughout a complete
// analysis and/or simulation project.
//
// Note : This memberfunction does not modify the numerical values of
//        the position coordinates.
//        It only specifies their numerical meaning.
// 
 if (s>0.)
 {
  fScale=s;
 }
 else
 {
  cout << " *AliPosition::SetUnitScale* Invalid argument s = " << s << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliPosition::GetUnitScale() const
{
// Provide the unit scale for the position coordinates.
// The scale is normalised w.r.t. the meter, so a unit scale of 0.01
// means that all position coordinates are in cm.
 return fScale;
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::SetTimestamp(AliTimestamp& t)
{
// Store the timestamp for this position.
 if (fTstamp) delete fTstamp;
 fTstamp=new AliTimestamp(t);
}
///////////////////////////////////////////////////////////////////////////
AliTimestamp* AliPosition::GetTimestamp()
{
// Provide the timestamp of this position.
 return fTstamp;
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::RemoveTimestamp()
{
// Remove the timestamp from this postion.
 if (fTstamp)
 {
  delete fTstamp;
  fTstamp=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::Data(TString f,TString u) const
{
// Provide all position/time information within the coordinate frame f.
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 Ali3Vector::Data(f,u);
 cout << "   Metric unit : " << fScale << " meter" << endl;
 if (fTstamp) fTstamp->Date(1);
} 
///////////////////////////////////////////////////////////////////////////
