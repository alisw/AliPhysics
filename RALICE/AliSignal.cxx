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

#include "AliSignal.h"
 
ClassImp(AliSignal) // Class implementation to enable ROOT I/O
 
AliSignal::AliSignal()
{
// Creation of an AliSignal object and initialisation of parameters
 Reset();
}
///////////////////////////////////////////////////////////////////////////
AliSignal::~AliSignal()
{
// Destructor to delete dynamically allocated memory
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::Reset()
{
// Reset all values
 Float_t r[3]={0,0,0};
 SetPosition(r,"sph");
 fSignal=0;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetSignal(Float_t sig)
{
// Store signal value
 fSignal=sig;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSignal::GetSignal()
{
// Provide signal value
 return fSignal;
}
///////////////////////////////////////////////////////////////////////////
