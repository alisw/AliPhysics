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

///////////////////////////////////////////////////////////////////////////
// Class AliPositionObj
// Handling of positions in various reference frames.
//
// This class is meant to provide an AliPosition object which is derived
// from TObject such that it can be stored in e.g. TObjArray etc...
// and that it can be written out using the ROOT I/O machinery.
//
// Example :
// =========
//
// Float_t a[3]={1,2,3};
// Float_t ea[3]={0.01,0.02,0.03};
// Float_t b[3]={4,5,6};
// Float_t eb[3]={0.04,0.05,0.06};
//
// AliPosition r1,r2;
//
// r1.SetPosition(a,"car");
// r1.SetPositionErrors(ea,"car");
// r2.SetPosition(b,"car");
// r2.SetPositionErrors(eb,"car");
//
// Ali3Vector sum=r1+r2;
// Ali3Vector rel=r1-r2;
//
// AliPositionObj rr1(r1);
// AliPositionObj rr2;
// rr2.Load(r2);
// AliPositionObj ssum(r1+r2);
//
// rr1.Info();
// rr2.Info();
// ssum.Info();
//
//--- Author: Nick van Eijndhoven 18-oct-1999 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliPositionObj.h"
 
ClassImp(AliPositionObj) // Class implementation to enable ROOT I/O
 
AliPositionObj::AliPositionObj()
{
// Creation of an AliPositionObj object and initialisation of parameters.
// All attributes initialised to 0.
}
///////////////////////////////////////////////////////////////////////////
AliPositionObj::AliPositionObj(Ali3Vector& q)
{
// Creation of an AliPositionObj object and initialisation of parameters.
// All attributes are initialised to the values of the input Ali3Vector.
 Load(q);
}
///////////////////////////////////////////////////////////////////////////
AliPositionObj::~AliPositionObj()
{
// Destructor to delete dynamically allocated memory.
}
///////////////////////////////////////////////////////////////////////////
void AliPositionObj::Load(Ali3Vector& q)
{
// Load all attributes of the input Ali3Vector into this AliPositionObj object.
 Double_t temp=q.GetResultError();
 Double_t a[3];
 q.GetVector(a,"sph");
 SetPosition(a,"sph");
 q.GetErrors(a,"car");
 SetPositionErrors(a,"car");
 fDresult=temp;
}
///////////////////////////////////////////////////////////////////////////
