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
// Class AliMUONVGeometryDEIndexing
// --------------------------------
// The abstract singleton base class for definition of
// the conversion between the detection element Ids and 
// the indexing in a simple array.
//
// Author: Ivana Hrivnacova, IPN Orsay

#include "AliMUONVGeometryDEIndexing.h"

ClassImp(AliMUONVGeometryDEIndexing)

//______________________________________________________________________________
AliMUONVGeometryDEIndexing::AliMUONVGeometryDEIndexing()
 : TObject()
{ 
// Standard/default constructor
}

//______________________________________________________________________________
AliMUONVGeometryDEIndexing::~AliMUONVGeometryDEIndexing() {
//
}

//______________________________________________________________________________
Int_t AliMUONVGeometryDEIndexing::GetModuleId(Int_t detElemId)
{
// Get module Id from detection element Id
// ---

  return detElemId/100 - 1;
}  

