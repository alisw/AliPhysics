// $Id$
//
// Class AliMUONChamberGeometry
// -----------------------------
// Class for definititon of the MUON chamber positions in ALIC
// Author: Ivana Hrivnacova, IPN Orsay
// 23/01/2004

#include <TVirtualMC.h>
#include <TGeoMatrix.h>
#include <TObjArray.h>
#include <TArrayI.h>
#include <Riostream.h>

#include "AliMUONChamberGeometry.h"
#include "AliMUONGeometryEnvelope.h"
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONGeometryTransformStore.h"	
#include "AliMUONConstants.h"
#include "AliLog.h"

ClassImp(AliMUONChamberGeometry)

//______________________________________________________________________________
AliMUONChamberGeometry::AliMUONChamberGeometry(Int_t chamberId)
 : TObject(),
   fChamberId(chamberId),
   fMotherVolume("ALIC"),
   fNofSVs(0),
   fSVVolumeIds(0),
   fTransformation(0),
   fDETransforms(0),
   fEnvelopes(0),
   fSVMap(0)
{
// Standard constructor

  // Chamber transformation
  fTransformation = new TGeoCombiTrans("");

  // Arrays of volumes Ids
  fSVVolumeIds = new TArrayI(20);

  // Sensitive volumes map
  fSVMap = new AliMUONGeometrySVMap(100);

  // Det elements transformation store
  fDETransforms = new AliMUONGeometryTransformStore(
                         AliMUONConstants::GetFirstDetElemId(chamberId), 
			 AliMUONConstants::NofDetElements(chamberId),
			 fSVMap);  
  // Envelope store
  fEnvelopes = new AliMUONGeometryEnvelopeStore(fDETransforms);  
}


//______________________________________________________________________________
AliMUONChamberGeometry::AliMUONChamberGeometry()
 : TObject(),
   fChamberId(0),
   fMotherVolume(),
   fNofSVs(0),
   fSVVolumeIds(0),
   fTransformation(0),
   fDETransforms(0),
   fEnvelopes(0),
   fSVMap(0)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONChamberGeometry::AliMUONChamberGeometry(const AliMUONChamberGeometry& rhs)
  : TObject(rhs)
{
  AliFatal("Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONChamberGeometry::~AliMUONChamberGeometry() {
//

  delete fTransformation;
  delete fSVVolumeIds;
  delete fEnvelopes;
  delete fDETransforms;
  delete fSVMap;
}

//______________________________________________________________________________
AliMUONChamberGeometry& 
AliMUONChamberGeometry::operator = (const AliMUONChamberGeometry& rhs) 
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
Int_t AliMUONChamberGeometry::GetSVIndex(Int_t svVolId) const
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
void  AliMUONChamberGeometry::SetTranslation(const TGeoTranslation& translation)
{
// Sets the chamber position wrt ALIC.
// ---

  fTransformation
    ->SetTranslation(const_cast<Double_t*>(translation.GetTranslation()));
}  

//______________________________________________________________________________
void  AliMUONChamberGeometry::SetRotation(const TGeoRotation& rotation)
{
// Sets the chamber rotation wrt ALIC.
// ---

  TGeoRotation* rot = new TGeoRotation();
  rot->SetMatrix(const_cast<Double_t*>(rotation.GetRotationMatrix()));

  fTransformation->SetRotation(rot);
}  

//______________________________________________________________________________
void  AliMUONChamberGeometry::SetSensitiveVolume(Int_t svVolId)
{
// Adds the volume specified by volId to the list of sensitive
// volumes
// ---
  
  // Resize TArrayI if needed
  if (fSVVolumeIds->GetSize() == fNofSVs) fSVVolumeIds->Set(2*fNofSVs);

  fSVVolumeIds->AddAt(svVolId, fNofSVs++);
}      

//______________________________________________________________________________
void  AliMUONChamberGeometry::SetSensitiveVolume(const TString& volName)
{
// Adds the volume specified by volName to the list of sensitive
// volumes
// ---

  SetSensitiveVolume(gMC->VolId(volName));
}      

//______________________________________________________________________________
void  AliMUONChamberGeometry::SetAlign(Bool_t align)
{
// Sets alignement option to enevelope store.
// ---
  
  fEnvelopes->SetAlign(align);
}  

//______________________________________________________________________________
Bool_t AliMUONChamberGeometry::IsSensitiveVolume(Int_t volId) const
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
Bool_t AliMUONChamberGeometry::IsSensitiveVolume(const TString& volName) const
{
// Checks if the volume specified by volName  is present in the list
// of sensitive volumes.
// ---

  return IsSensitiveVolume(gMC->VolId(volName));
}
