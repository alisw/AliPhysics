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
// Class AliCalmodule
// Description of a module in a calorimeter system.
// A matrix geometry is assumed, such that a module
// is identified by (row,col) and contains a certain signal.
// Note : row and col start counting at 1.
//
//--- Author: Nick van Eijndhoven 13-jun-1997 UU-SAP Utrecht
//- Modified: NvE $Date$ UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliCalmodule.h"
#include "Riostream.h"
 
ClassImp(AliCalmodule) // Class implementation to enable ROOT I/O
 
AliCalmodule::AliCalmodule() : AliSignal()
{
// Default constructor, all module data is set to 0
 fRow=0;
 fCol=0;
 fSigc=0;
}
///////////////////////////////////////////////////////////////////////////
AliCalmodule::~AliCalmodule()
{
// Default destructor
}
///////////////////////////////////////////////////////////////////////////
AliCalmodule::AliCalmodule(const AliCalmodule& m) : AliSignal(m)
{
// Copy constructor
 fRow=m.fRow;
 fCol=m.fCol;
 fSigc=m.fSigc;
}
///////////////////////////////////////////////////////////////////////////
AliCalmodule::AliCalmodule(Int_t row,Int_t col,Double_t sig) : AliSignal()
{
// Module constructor with initialisation of module data
 fRow=row;
 fCol=col;
 AliSignal::SetSignal(sig);
 fSigc=sig;
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
void AliCalmodule::SetSignal(Double_t sig,Int_t j)
{
// Set or change the data of the module.
// This is an extension of AliSignal::SetSignal in view of the clustered signal.
 AliSignal::SetSignal(sig,j);
 if (j==1) fSigc=sig;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::AddSignal(Double_t sig,Int_t j)
{
// Add or change the data of the module
// This is an extension of AliSignal::AddSignal in view of the clustered signal.
 AliSignal::AddSignal(sig,j);
 if (j==1) fSigc+=sig;
}
///////////////////////////////////////////////////////////////////////////
void AliCalmodule::SetClusteredSignal(Double_t sig)
{
// Set or change the signal of the module after clustering
 fSigc=sig;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalmodule::GetRow() const
{
// Provide the row number of the module
 return fRow;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalmodule::GetColumn() const
{
// Provide the column number of the module
 return fCol;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalmodule::GetClusteredSignal() const
{
// Provide the signal of the module after clustering.
 Int_t dead=GetDeadValue();
 if (!dead)
 {
  return fSigc;
 }
 else
 {
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
TObject* AliCalmodule::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers like AliCalorimeter when adding objects in case the
// container owns the objects. This feature allows e.g. AliCalorimeter
// to store either AliCalmodule objects or objects derived from AliCalmodule
// via tha AddSignal memberfunction, provided these derived classes also have
// a proper Clone memberfunction. 

 AliCalmodule* m=new AliCalmodule(*this);
 if (name)
 {
  if (strlen(name)) m->SetName(name);
 }
 return m;
}
///////////////////////////////////////////////////////////////////////////
