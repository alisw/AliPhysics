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
// Class Ali4VectorObj
// Handling of Lorentz 4-vectors in various reference frames.
//
// This class is meant to provide an Ali4Vector object which is derived
// from TObject such that it can be stored in e.g. TObjArray etc...
// and that it can be written out using the ROOT I/O machinery.
//
// Example :
// =========
//
// Float_t a[4]={5,1,2,3};
// Float_t ea[4]={0.05,0.01,0.02,0.03};
// Float_t b[4]={10,4,5,6};
// Float_t eb[4]={0.1,0.04,0.05,0.06};
//
// Ali4Vector v,w;
//
// v.SetVector(a,"car");
// v.SetErrors(ea,"car");
// w.SetVector(b,"car");
// w.SetErrors(eb,"car");
//
// Ali4Vector add=v+w;
//
// Ali4Vector sub=v-w;
//
// Ali4VectorObj vec1(add);
//
// Ali4VectorObj vec2;
// vec2.Load(sub);
//
// vec1.Info();
// vec2.Info();
//
//--- Author: Nick van Eijndhoven 18-oct-1999 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "Ali4VectorObj.h"
 
ClassImp(Ali4VectorObj) // Class implementation to enable ROOT I/O
 
Ali4VectorObj::Ali4VectorObj()
{
// Creation of an Ali4VectorObj object and initialisation of parameters.
// All attributes initialised to 0.
}
///////////////////////////////////////////////////////////////////////////
Ali4VectorObj::Ali4VectorObj(Ali4Vector& q)
{
// Creation of an Ali3VectorObj object and initialisation of parameters.
// All attributes are initialised to the values of the input Ali3Vector.
 Load(q);
}
///////////////////////////////////////////////////////////////////////////
Ali4VectorObj::~Ali4VectorObj()
{
// Destructor to delete dynamically allocated memory.
}
///////////////////////////////////////////////////////////////////////////
void Ali4VectorObj::Load(Ali4Vector& q)
{
// Load all attributes of the input Ali4Vector into this Ali4VectorObj object.
 Int_t temp1=q.GetScalarFlag();
 Double_t temp2=q.GetResultError();
 Double_t a[4];
 q.GetVector(a,"sph");
 SetVector(a,"sph");
 q.GetErrors(a,"car");
 SetErrors(a,"car");
 fScalar=temp1;
 fDresult=temp2;
}
///////////////////////////////////////////////////////////////////////////
