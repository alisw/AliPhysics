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
// Class AliMUONGeometryModule
// -----------------------------
// Class for definition of the detector module parameters
// (the transformations of detection elements, mapping between
//  sensitive volumes and detection elements).
//
// Author: Ivana Hrivnacova, IPN Orsay

#include <TVirtualMC.h>
#include <TGeoMatrix.h>
#include <TObjArray.h>
#include <TArrayI.h>
#include <Riostream.h>

#include "AliLog.h"	

#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryEnvelope.h"
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONGeometryDetElement.h"	
#include "AliMUONGeometryStore.h"	
#include "AliMUONGeometrySVMap.h"	
#include "AliMUONGeometryDEIndexing.h"

ClassImp(AliMUONGeometryModule)

//______________________________________________________________________________
AliMUONGeometryModule::AliMUONGeometryModule(Int_t moduleId)
 : TObject(),
   fModuleId(moduleId),
   fMotherVolume("ALIC"),
   fNofSVs(0),
   fSVVolumeIds(0),
   fTransformation(0),
   fEnvelopes(0),
   fDEIndexing(0),
   fDetElements(0),
   fSVMap(0)
{
// Standard constructor

  // Chamber transformation
  fTransformation = new TGeoCombiTrans("");

  // Arrays of volumes Ids
  fSVVolumeIds = new TArrayI(20);

  // Sensitive volumes map
  fSVMap = new AliMUONGeometrySVMap(100);

  // Get indexing
  fDEIndexing = new AliMUONGeometryDEIndexing(fModuleId, 0);

  // Det elements transformation stores
  fDetElements = new AliMUONGeometryStore(fDEIndexing);
    
  // Envelope store
  fEnvelopes = new AliMUONGeometryEnvelopeStore(fDetElements);  
}


//______________________________________________________________________________
AliMUONGeometryModule::AliMUONGeometryModule()
 : TObject(),
   fModuleId(0),
   fMotherVolume(),
   fNofSVs(0),
   fSVVolumeIds(0),
   fTransformation(0),
   fEnvelopes(0),
   fDEIndexing(0),
   fDetElements(0),
   fSVMap(0)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONGeometryModule::AliMUONGeometryModule(const AliMUONGeometryModule& rhs)
  : TObject(rhs)
{
  AliFatal("Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONGeometryModule::~AliMUONGeometryModule() {
//

  delete fTransformation;
  delete fSVVolumeIds;
  delete fEnvelopes;
  delete fDEIndexing;
  delete fDetElements;
  delete fSVMap;
}

//______________________________________________________________________________
AliMUONGeometryModule& 
AliMUONGeometryModule::operator = (const AliMUONGeometryModule& rhs) 
{
  // check assignement to self
  if (this == &rhs) return *this;

  AliFatal("Assignment operator is not implemented.");
    
  return *this;  
}

//
// private methods
//

//______________________________________________________________________________
Int_t AliMUONGeometryModule::GetSVIndex(Int_t svVolId) const
{
// Returns the index of the volume specified by volId
// if it is present in the list of sensitive volumes 
// (or -1 if not present).
 
  for (Int_t i=0; i<fNofSVs; i++) {
      if (fSVVolumeIds->At(i) == svVolId) return i;
  }
  return -1;
}

//
// public methods
//

//______________________________________________________________________________
void  AliMUONGeometryModule::Global2Local(Int_t detElemId,
                                  Float_t xg, Float_t yg, Float_t zg, 
                                  Float_t& xl, Float_t& yl, Float_t& zl) const
{
// Transforms point from the global reference frame (ALIC)
// to the local reference frame of the detection element specified
// by detElemId.
// ---

  // Get detection element
  AliMUONGeometryDetElement* detElement = GetDetElement(detElemId);
  if (!detElement) return;
   
  // Transform point
  detElement->Global2Local(xg, yg, zg, xl, yl, zl);
}
				  
//______________________________________________________________________________
void  AliMUONGeometryModule::Global2Local(Int_t detElemId,
                                  Double_t xg, Double_t yg, Double_t zg, 
                                  Double_t& xl, Double_t& yl, Double_t& zl) const
{
// Transforms point from the global reference frame (ALIC)
// to the local reference frame of the detection element specified
// by detElemId.
// ---

   // Get detection element
   AliMUONGeometryDetElement* detElement = GetDetElement(detElemId);
   if (!detElement) return;
   
   // Transform point
   detElement->Global2Local(xg, yg, zg, xl, yl, zl);
}
				  
//______________________________________________________________________________
void  AliMUONGeometryModule::Local2Global(Int_t detElemId,
                 Float_t xl, Float_t yl, Float_t zl, 
                 Float_t& xg, Float_t& yg, Float_t& zg) const
{
// Transforms point from the local reference frame of the detection element 
// specified by detElemId to the global reference frame (ALIC).
// ---

  // Get detection element
  AliMUONGeometryDetElement* detElement = GetDetElement(detElemId);
  if (!detElement) return;
   
   // Transform point
  detElement->Local2Global(xl, yl, zl, xg, yg, zg);  
}

//______________________________________________________________________________
void  AliMUONGeometryModule::Local2Global(Int_t detElemId,
                 Double_t xl, Double_t yl, Double_t zl, 
                 Double_t& xg, Double_t& yg, Double_t& zg) const
{
// Transforms point from the local reference frame of the detection element 
// specified by detElemId to the global reference frame (ALIC).
// ---

   // Get detection element
   AliMUONGeometryDetElement* detElement = GetDetElement(detElemId);
   if (!detElement) return;
   
   // Transform point
   detElement->Local2Global(xl, yl, zl, xg, yg, zg); 
}

//______________________________________________________________________________
void  AliMUONGeometryModule::SetTranslation(const TGeoTranslation& translation)
{
// Sets the module position wrt world.
// ---

  fTransformation
    ->SetTranslation(const_cast<Double_t*>(translation.GetTranslation()));
}  

//______________________________________________________________________________
void  AliMUONGeometryModule::SetRotation(const TGeoRotation& rotation)
{
// Sets the module rotation wrt ALIC.
// ---

  TGeoRotation* rot = new TGeoRotation();
  rot->SetMatrix(const_cast<Double_t*>(rotation.GetRotationMatrix()));

  fTransformation->SetRotation(rot);
}  

//______________________________________________________________________________
void  AliMUONGeometryModule::SetSensitiveVolume(Int_t svVolId)
{
// Adds the volume specified by volId to the list of sensitive
// volumes
// ---
  
  // Resize TArrayI if needed
  if (fSVVolumeIds->GetSize() == fNofSVs) fSVVolumeIds->Set(2*fNofSVs);

  fSVVolumeIds->AddAt(svVolId, fNofSVs++);
}      

//______________________________________________________________________________
void  AliMUONGeometryModule::SetSensitiveVolume(const TString& volName)
{
// Adds the volume specified by volName to the list of sensitive
// volumes
// ---

  SetSensitiveVolume(gMC->VolId(volName));
}      

//______________________________________________________________________________
void  AliMUONGeometryModule::SetAlign(Bool_t align)
{
// Sets alignement option to enevelope store.
// ---
  
  fEnvelopes->SetAlign(align);
}  

//______________________________________________________________________________
AliMUONGeometryDetElement* 
AliMUONGeometryModule::FindBySensitiveVolume(const TString& sensVolume) const
{
// Finds TGeoCombiTrans for the detector element Id specified by aligned volume 
// ---

  Int_t detElemId = fSVMap->GetDetElemId(sensVolume);

  if (!detElemId) return 0; 
        // The specified sensitive volume is not in the map   
  
  return (AliMUONGeometryDetElement*)fDetElements->Get(detElemId);
}  

//______________________________________________________________________________
Bool_t AliMUONGeometryModule::IsSensitiveVolume(Int_t volId) const
{
// Checks if the volume specified by volId is present in the list
// of sensitive volumes.
// ---

  for (Int_t i=0; i<fNofSVs; i++) {
      if (fSVVolumeIds->At(i) == volId) return kTRUE;
  }
  return kFALSE;
}

//______________________________________________________________________________
Bool_t AliMUONGeometryModule::IsSensitiveVolume(const TString& volName) const
{
// Checks if the volume specified by volName  is present in the list
// of sensitive volumes.
// ---

  return IsSensitiveVolume(gMC->VolId(volName));
}

//______________________________________________________________________________
AliMUONGeometryDetElement*
AliMUONGeometryModule::GetDetElement(Int_t detElemId) const
{
// Returns thethe detection element specified by detElemId.
// Gives error if detection element is not defined.
// ---

   // Get detection element
   AliMUONGeometryDetElement* detElement
     = (AliMUONGeometryDetElement*) fDetElements->Get(detElemId);

   if (!detElement) {
     AliError(Form("Detection element %d not found", detElemId));
     return 0;
   }  

   return detElement;
}

/*				  
//______________________________________________________________________________
Int_t  AliMUONGeometryModule::GetNofDetElements() const
{
// Returns the number of detection elements
// ---

  return fDEIndexing->GetNofDetElements();
}   

//______________________________________________________________________________
Int_t  AliMUONGeometryModule::GetDetElemId(Int_t i) const
{
// Returns the i-th detection element id
// ---

  return fDEIndexing->GetDetElementId(i);
}   
*/
