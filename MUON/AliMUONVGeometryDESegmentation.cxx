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

//
// Class AliVMUONGeometryDESegmentation
// ----------------------------------
// Extension for AliSegmentation interface,
// added functions:
//  Bool_t  HasPad(Float_t x, Float_t y, Float_t z);
//  Bool_t  HasPad(Int_t ix, Int_t iy);
//
// Author:Ivana Hrivnacova, IPN Orsay

#include "AliLog.h"

#include "AliMUONVGeometryDESegmentation.h"

ClassImp(AliMUONVGeometryDESegmentation)


//______________________________________________________________________________
AliMUONVGeometryDESegmentation::AliMUONVGeometryDESegmentation() 
: AliSegmentation()
{
// Normal/default constructor
}

//______________________________________________________________________________
AliMUONVGeometryDESegmentation::AliMUONVGeometryDESegmentation(
                                  const AliMUONVGeometryDESegmentation& rhs) 
  : AliSegmentation(rhs)
{
// Copy constructor
  AliFatal("Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONVGeometryDESegmentation::~AliMUONVGeometryDESegmentation() {
// Destructor
} 

//
// operators
//

//______________________________________________________________________________
AliMUONVGeometryDESegmentation& 
AliMUONVGeometryDESegmentation::operator=(const AliMUONVGeometryDESegmentation& rhs)
{
// Copy operator 

  // check assignement to self
  if (this == &rhs) return *this;

  AliFatal("Assignment operator is not implemented.");
    
  return *this;  
}
