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
Revision 1.3  1999/11/03 14:23:17  fca
New version of RALICE introduced

Revision 1.2  1999/09/29 09:24:28  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////
// Class AliCalmodule
// Description of a module in a calorimeter system.
// A matrix geometry is assumed, such that a module
// is identified by (row,col) and contains a certain signal.
// Note : row and col start counting at 1.
//
//--- Author: Nick van Eijndhoven 13-jun-1997 UU-SAP Utrecht
//- Modified: NvE 18-jan-2000 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliCalmodule.h"
 
ClassImp(AliCalmodule) // Class implementation to enable ROOT I/O
 
AliCalmodule::AliCalmodule()
{
// Default constructor, all module data is set to 0
 fRow=0;
 fCol=0;
 fSigc=0;
 fDead=0;
 fGain=1;
}
///////////////////////////////////////////////////////////////////////////
AliCalmodule::~AliCalmodule()
{
// Default destructor
}
///////////////////////////////////////////////////////////////////////////
AliCalmodule::AliCalmodule(Int_t row,Int_t col,Float_t sig)
{
// Module constructor with initialisation of module data
 fRow=row;
 fCol=col;
 AliSignal::SetSignal(sig);
 fSigc=sig;
 fDead=0;
 fGain=1;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::SetRow(Int_t i)
{
// Set the row number for this module
 fRow=i;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::SetColumn(Int_t i)
{
// Set the column number for this module
 fCol=i;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::SetSignal(Int_t row,Int_t col,Float_t sig)
{
// Set or change the data of the module
 fRow=row;
 fCol=col;
 AliSignal::SetSignal(sig);
 fSigc=sig;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::AddSignal(Int_t row,Int_t col,Float_t sig)
{
// Add or change the data of the module
 fRow=row;
 fCol=col;
 AliSignal::AddSignal(sig);
 fSigc+=sig;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::SetClusteredSignal(Float_t sig)
{
// Set or change the signal of the module after clustering
 fSigc=sig;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::SetDead()
{
// Indicate the module as dead
 fDead=1;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::SetAlive()
{
// Indicate the module as dead
 fDead=0;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::SetGain(Float_t gain)
{
// Set the gain value of the readout system
 fGain=gain;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalmodule::GetRow()
{
// Provide the row number of the module
 return fRow;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalmodule::GetColumn()
{
// Provide the column number of the module
 return fCol;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalmodule::GetClusteredSignal()
{
// Provide the signal of the module after clustering
 if (!fDead)
 {
  return fSigc;
 }
 else
 {
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalmodule::GetDeadValue()
{
// Provide the value of the dead indicator (1=dead 0=alive)
 return fDead;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalmodule::GetGain()
{
// Provide the gain value of the readout system
 return fGain;
}
///////////////////////////////////////////////////////////////////////////
