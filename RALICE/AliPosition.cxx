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
// Creation of an AliPosition object and initialisation of parameters
}
///////////////////////////////////////////////////////////////////////////
AliPosition::~AliPosition()
{
// Destructor to delete dynamically allocated memory
}
///////////////////////////////////////////////////////////////////////////
AliPosition::AliPosition(const AliPosition& p) : Ali3Vector(p)
{
// Copy constructor
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::SetPosition(Double_t* r,TString f)
{
// Store position according to reference frame f
 SetVector(r,f);
}
///////////////////////////////////////////////////////////////////////////
void AliPosition::GetPosition(Double_t* r,TString f)
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
void AliPosition::GetPosition(Float_t* r,TString f)
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
void AliPosition::GetPositionErrors(Double_t* r,TString f)
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
void AliPosition::GetPositionErrors(Float_t* r,TString f)
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
// Provide distance to position p.
// The error on the result can be obtained as usual by invoking
// GetResultError() afterwards. 
 Ali3Vector d=(Ali3Vector)((*this)-p);
 Double_t dist=d.GetNorm();
 fDresult=d.GetResultError();
 return dist;
}
///////////////////////////////////////////////////////////////////////////
