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
// Class Ali3VectorObj
// Handling of 3-vectors in various reference frames.
//
// This class is meant to provide an Ali3Vector object which is derived
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
// Ali3Vector v,w;
//
// v.SetVector(a,"car");
// v.SetErrors(ea,"car");
// w.SetVector(b,"car");
// w.SetErrors(eb,"car");
//
// Ali3Vector cross=v.Cross(w);
//
// Ali3Vector add=v+w;
//
// Ali3VectorObj vec1(cross);
//
// Ali3VectorObj vec2;
// vec2.Load(add);
//
// vec1.Info();
// vec2.Info();
//
//--- Author: Nick van Eijndhoven 18-oct-1999 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "Ali3VectorObj.h"
 
ClassImp(Ali3VectorObj) // Class implementation to enable ROOT I/O
 
Ali3VectorObj::Ali3VectorObj()
{
// Creation of an Ali3VectorObj object and initialisation of parameters.
// All attributes initialised to 0.
}
///////////////////////////////////////////////////////////////////////////
Ali3VectorObj::Ali3VectorObj(Ali3Vector& q)
{
// Creation of an Ali3VectorObj object and initialisation of parameters.
// All attributes are initialised to the values of the input Ali3Vector.
 Load(q);
}
///////////////////////////////////////////////////////////////////////////
Ali3VectorObj::~Ali3VectorObj()
{
// Destructor to delete dynamically allocated memory.
}
///////////////////////////////////////////////////////////////////////////
void Ali3VectorObj::Load(Ali3Vector& q)
{
// Load all attributes of the input Ali3Vector into this Ali3VectorObj object.
 Double_t temp=q.GetResultError();
 Double_t a[3];
 q.GetVector(a,"sph");
 SetVector(a,"sph");
 q.GetErrors(a,"car");
 SetErrors(a,"car");
 fDresult=temp;
}
///////////////////////////////////////////////////////////////////////////
