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

#include "AliPosition.h"
 
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
}
///////////////////////////////////////////////////////////////////////////
