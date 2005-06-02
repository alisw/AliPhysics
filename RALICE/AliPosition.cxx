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
// Positions (r), errors (e) and reference frames (f) are specified via
//
//    SetPosition(Float_t* r,TString f)
//    SetPositionErrors(Float_t* e,TString f)
//
// under the following conventions :
//
// f="car" ==> r in Cartesian coordinates   (x,y,z)
// f="sph" ==> r in Spherical coordinates   (r,theta,phi)
// f="cyl" ==> r in Cylindrical coordinates (rho,phi,z)
//
// All angles are in radians.
//
// The unit scale for the coordinates can be defined by the user
// via the SetUnitScale() memberfunction.
// This enables standardised expressions using numerical values of
// physical constants by means of the GetUnitScale() memberfunction.
// By default the unit scale is set to cm, corresponding to invokation
// of SetUnitScale(0.01).
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
// q.GetPosition(loc,"sph");
// q.GetPositionErrors(dloc,"sph");
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
 fScale=0.01;
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
void AliPosition::SetPosition(Double_t* r,TString f)
{
// Store position according to reference frame f
 SetVector(r,f);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::GetPosition(Double_t* r,TString f) const
{
// Provide position according to reference frame f
 GetVector(r,f);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::SetPosition(Float_t* r,TString f)
{
// Store position according to reference frame f
 SetVector(r,f);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::GetPosition(Float_t* r,TString f) const
{
// Provide position according to reference frame f
 GetVector(r,f);
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
void AliPosition::SetPositionErrors(Double_t* r,TString f)
{
// Store position errors according to reference frame f
 SetErrors(r,f);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::GetPositionErrors(Double_t* r,TString f) const
{
// Provide position errors according to reference frame f
 GetErrors(r,f);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::SetPositionErrors(Float_t* r,TString f)
{
// Store position errors according to reference frame f
 SetErrors(r,f);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::GetPositionErrors(Float_t* r,TString f) const
{
// Provide position errors according to reference frame f
 GetErrors(r,f);
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
Double_t AliPosition::GetDistance(AliPosition& p)
{
// Provide distance of the current AliPosition to position p.
// The error on the result can be obtained as usual by invoking
// GetResultError() afterwards. 
//
// In the case of two positions with different unit scales, the distance
// will be provided in the unit scale of the current AliPosition.
// This implies that in such cases the results of r.GetDistance(q) and
// q.GetDistance(r) will be numerically different.
// As such it is possible to obtain a correctly computed distance between
// positions which have different unit scales.
// However, it is recommended to work always with one single unit scale.
//
 Ali3Vector d=(Ali3Vector)p;
 Float_t pscale=p.GetUnitScale();
 if ((pscale/fScale > 1.1) || (fScale/pscale > 1.1)) d=d*(pscale/fScale);
 Ali3Vector q=(Ali3Vector)(*this);
 d=d-q;
 Double_t dist=d.GetNorm();
 fDresult=d.GetResultError();
 return dist;
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::SetUnitScale(Float_t s)
{
// Set the unit scale for the position coordinates.
// The scale is normalised w.r.t. the meter, so setting the unit scale
// to 0.01 means that all position coordinates are in cm.
// By default the unit scale is set to cm in the AliPosition constructor.
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
void AliPosition::Data(TString f) const
{
// Provide all position/time information within the coordinate frame f.
 Ali3Vector::Data(f);
 if (fTstamp) fTstamp->Date(1);
} 
///////////////////////////////////////////////////////////////////////////
