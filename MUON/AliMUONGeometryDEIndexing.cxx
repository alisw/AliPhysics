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
// Class AliMUONGeometryDEIndexing
// -------------------------------
// The class that provides conversion between the detection element Id
// and the index in the array.
// Used in storing DE transformations and segmentations.
// (See more in the header file.) 
//
// Author: Ivana Hrivnacova, IPN Orsay

#include <Riostream.h>
#include <TGeoMatrix.h>
#include <TObjString.h>

#include "AliLog.h"

#include "AliMUONGeometryDEIndexing.h"
#include "AliMUONConstants.h"

ClassImp(AliMUONGeometryDEIndexing)

//______________________________________________________________________________
AliMUONGeometryDEIndexing::AliMUONGeometryDEIndexing(
                              Int_t moduleId, Int_t nofDetElements)
 : AliMUONVGeometryDEIndexing(),
   fModuleId(moduleId),
   fNofDetElements(nofDetElements)
{ 
// Standard constructor
}

//______________________________________________________________________________
AliMUONGeometryDEIndexing::AliMUONGeometryDEIndexing()
 : AliMUONVGeometryDEIndexing(),
   fModuleId(0),
   fNofDetElements(0)
{ 
// Ddefault constructor
}

//______________________________________________________________________________
AliMUONGeometryDEIndexing::~AliMUONGeometryDEIndexing() {
//
}

//
// private methods
//

//______________________________________________________________________________
Int_t AliMUONGeometryDEIndexing::GetFirstDetElemId() const
{
// Get first detection element Id for chamber specified by moduleId
// ---

  return (fModuleId+1)*100;
}  

//
// public methods
//

//______________________________________________________________________________
Int_t AliMUONGeometryDEIndexing::GetDetElementIndex(Int_t detElemId) const
{
// Returns the index of detector element specified by detElemId
// ---

  if ( fNofDetElements == 0 ) {
    AliFatal("The number of detection elements has not been set.");
    return 0;
  }  

  return detElemId - GetFirstDetElemId();
}  

//______________________________________________________________________________
Int_t AliMUONGeometryDEIndexing::GetDetElementId(Int_t detElemIndex) const
{
// Returns the ID of detector element specified by index
// ---

  if ( fNofDetElements == 0 ) {
    AliFatal("The number of detection elements has not been set.");
    return 0;
  }  

  return GetFirstDetElemId() + detElemIndex;
}  
