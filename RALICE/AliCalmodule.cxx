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

#include "AliCalmodule.h"
 
ClassImp(AliCalmodule) // Class implementation to enable ROOT I/O
 
AliCalmodule::AliCalmodule()
{
// Default constructor, all module data is set to 0
 fRow=0;
 fCol=0;
 fSigc=0;
 fEdge=0;
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
 fSignal=sig;
 fSigc=sig;
 fEdge=0;
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
 fSignal=sig;
 fSigc=sig;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::AddSignal(Int_t row,Int_t col,Float_t sig)
{
// Add or change the data of the module
 fRow=row;
 fCol=col;
 fSignal+=sig;
 fSigc+=sig;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::SetClusteredSignal(Float_t sig)
{
// Set or change the signal of the module after clustering
 fSigc=sig;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::SetEdgeOn()
{
// Indicate the module as edge module
 fEdge=1;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::SetEdgeOff()
{
// Indicate the module as non-edge module
 fEdge=0;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::EdgeUp()
{
// Increase the edge value by 1
// This simplifies treatment of edge modules around temp. dead modules
 fEdge+=1;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::EdgeDown()
{
// Decrease the edge value by 1
// This simplifies treatment of edge modules around temp. dead modules
 if (fEdge > 0) fEdge-=1;
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
Float_t AliCalmodule::GetSignal()
{
// Provide the signal value of the module
 if (!fDead)
 {
  return fSignal;
 }
 else
 {
  return 0;
 }
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
Int_t AliCalmodule::GetEdgeValue()
{
// Provide the value of the edge indicator (1=edge 0=no-edge)
 return fEdge;
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
