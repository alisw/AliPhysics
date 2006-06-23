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
// Class AliMUONGeometryStore
// -------------------------------------
// The class contains the array of the detection elements,
// which are sorted using the AliMUONVGeometryDEIndexing class.
// The class provides fast access to detection element via DetElemId.
//
// Author: Ivana Hrivnacova, IPN Orsay

#include <Riostream.h>
#include <TGeoMatrix.h>

#include "AliLog.h"

#include "AliMUONGeometryStore.h"

const Int_t AliMUONGeometryStore::fgkInitSize = 100;
const Int_t AliMUONGeometryStore::fgkCoefficient = 100;

/// \cond CLASSIMP
ClassImp(AliMUONGeometryStore)
/// \endcond

//
// static methods
//

//______________________________________________________________________________
Int_t AliMUONGeometryStore::GetModuleId(Int_t detElemId)
{
/// Get module Id from detection element Id

  return detElemId/fgkCoefficient - 1;
}  

//
// Constructor/destructor
//

//______________________________________________________________________________
AliMUONGeometryStore::AliMUONGeometryStore(Bool_t isOwner)
 : TObject(),
   fObjects(fgkInitSize)
{ 
/// Standard constructor
  
  fObjects.SetOwner(isOwner);
  for (Int_t i=0; i<fgkInitSize; i++) fObjects[i] = 0;
}

//______________________________________________________________________________
AliMUONGeometryStore::AliMUONGeometryStore()
 : TObject(),
   fObjects()
{
/// Default constructor
}

//______________________________________________________________________________
AliMUONGeometryStore::~AliMUONGeometryStore() 
{
/// Destructor

}

//
// private methods
//

//______________________________________________________________________________
Int_t AliMUONGeometryStore::GetDEIndex(Int_t detElemId) const
{
/// Return the index of detector element specified by detElemId

  return detElemId - detElemId/fgkCoefficient*fgkCoefficient;
 }  

//
// public methods
//

//______________________________________________________________________________
void AliMUONGeometryStore::Add(Int_t objectId, TObject* object)
{
/// Add detection element in the array
/// if detection element with the same Id is not yet present. 
 
  // Expand array if the init size has been reached
  Int_t index = GetDEIndex(objectId);
  while ( index >= fObjects.GetSize() ) {
    Int_t size = fObjects.GetSize();
    fObjects.Expand(size + fgkInitSize);
    for (Int_t i=size; i<fObjects.GetSize(); i++) fObjects[i] = 0;
  }  

  // Add to the map  
  if ( !Get(objectId, false) )
    fObjects.AddAt(object, index);
  else 
    AliWarning(Form("The detection element %d is already present", objectId));  
}		      

//______________________________________________________________________________
TObject* 
AliMUONGeometryStore::Get(Int_t objectId, Bool_t warn) const
{
/// Return the object for the specified detector element Id

  Int_t index = GetDEIndex(objectId);
  
  if ( index >= 0 && index < fObjects.GetEntriesFast() )
    return (TObject*) fObjects.At(index);
  else {
    if (warn) AliWarning(Form("Index %d out of limits", index));
    return 0;  
  }  
}  

//______________________________________________________________________________
Int_t  AliMUONGeometryStore::GetNofEntries() const
{
/// Return number of entries.
/// \todo Add check if the array is already filled

  return fObjects.GetEntriesFast();
}  
 

//______________________________________________________________________________
TObject* 
AliMUONGeometryStore::GetEntry(Int_t index) const
{
/// Return entry at specified index.
/// \todo Add check if the array is already filled
  
  return (TObject*) fObjects.At(index);
}  
