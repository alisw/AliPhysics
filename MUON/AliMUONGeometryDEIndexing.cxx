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

Int_t const AliMUONGeometryDEIndexing::fgkSeparator = 100; 

//
// static methods
//

//______________________________________________________________________________
Int_t AliMUONGeometryDEIndexing::GetModuleId(Int_t detElemId)
{
// Get module Id from detection element Id
// ---

  return detElemId/fgkSeparator - 1;
}  

//______________________________________________________________________________
Int_t AliMUONGeometryDEIndexing::GetDEIndex(Int_t detElemId)
{
/// Returns the index of detector element specified by detElemId

  return detElemId - detElemId/fgkSeparator*fgkSeparator;
 }  

//______________________________________________________________________________
Int_t AliMUONGeometryDEIndexing::GetDEId(Int_t moduleId, Int_t detElemIndex)
{
/// Returns the ID of detector element specified by index

  return ( moduleId + 1 ) * fgkSeparator + detElemIndex;
}  
