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
//
// Class AliMUONVGeometryBuilder
// -----------------------------
// Abstract base class for geometry construction per chamber(s).
// Author: Ivana Hrivnacova, IPN Orsay
// 23/01/2004

#include <TObjArray.h>

#include "AliMUONVGeometryBuilder.h"
#include "AliMUONChamber.h"

ClassImp(AliMUONVGeometryBuilder)

//______________________________________________________________________________
AliMUONVGeometryBuilder::AliMUONVGeometryBuilder(
                                AliMUONChamber* ch1, AliMUONChamber* ch2,
                                AliMUONChamber* ch3, AliMUONChamber* ch4,
                                AliMUONChamber* ch5, AliMUONChamber* ch6)
 : TObject(),
   fChambers(0)
 {
// Standard constructor

  // Create the chambers array
  fChambers = new TObjArray();
  
  if (ch1) fChambers->Add(ch1);
  if (ch2) fChambers->Add(ch2);
  if (ch3) fChambers->Add(ch3);
  if (ch4) fChambers->Add(ch4);
  if (ch5) fChambers->Add(ch5);
  if (ch6) fChambers->Add(ch6);

}


//______________________________________________________________________________
AliMUONVGeometryBuilder::AliMUONVGeometryBuilder()
 : TObject(),
   fChambers(0)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONVGeometryBuilder::AliMUONVGeometryBuilder(const AliMUONVGeometryBuilder& rhs)
  : TObject(rhs)
{
// Protected copy constructor

  Fatal("Copy constructor", 
        "Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONVGeometryBuilder::~AliMUONVGeometryBuilder() {
//
  if (fChambers) {
    fChambers->Clear(); // Sets pointers to 0 sinche it is not the owner
    delete fChambers;
  }
}

//______________________________________________________________________________
AliMUONVGeometryBuilder& 
AliMUONVGeometryBuilder::operator = (const AliMUONVGeometryBuilder& rhs) 
{
// Protected assignement operator

  // check assignement to self
  if (this == &rhs) return *this;

  Fatal("operator=", 
        "Assignment operator is not implemented.");
    
  return *this;  
}

//
// public methods
//

//______________________________________________________________________________
AliMUONChamber*  AliMUONVGeometryBuilder::GetChamber(Int_t chamberId) const
{
// Returns the chamber specified by chamberId
// ---

  for (Int_t i=0; i<fChambers->GetEntriesFast(); i++) {
    AliMUONChamber* chamber = (AliMUONChamber*)fChambers->At(i);
    if ( chamber->GetId() == chamberId) return chamber;
  }   
  
  return 0;
}  
