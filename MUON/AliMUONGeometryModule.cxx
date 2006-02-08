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

#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryEnvelope.h"
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONGeometryDetElement.h"	
#include "AliMUONGeometryStore.h"	
#include "AliMUONStringIntMap.h"	

#include "AliLog.h"	

#include <TVirtualMC.h>
#include <TGeoMatrix.h>
#include <TObjArray.h>
#include <TArrayI.h>
#include <Riostream.h>

ClassImp(AliMUONGeometryModule)

//______________________________________________________________________________
AliMUONGeometryModule::AliMUONGeometryModule(Int_t moduleId)
 : TObject(),
   fIsVirtual(true),
   fNofSVs(0),
   fSVVolumeIds(0),
   fEnvelopes(0),
   fSVMap(0),
   fTransformer(0)
{
/// Standard constructor

  // Arrays of volumes Ids
  fSVVolumeIds = new TArrayI(20);

  // Sensitive volumes map
  fSVMap = new AliMUONStringIntMap();

  // Geometry parametrisation
  fTransformer = new AliMUONGeometryModuleTransformer(moduleId);
    
  // Envelope store
  fEnvelopes = new AliMUONGeometryEnvelopeStore(
                             fTransformer->GetDetElementStore());  
}


//______________________________________________________________________________
AliMUONGeometryModule::AliMUONGeometryModule()
 : TObject(),
   fIsVirtual(true),
   fNofSVs(0),
   fSVVolumeIds(0),
   fEnvelopes(0),
   fSVMap(0),
   fTransformer(0)
{
/// Default constructor
}


//______________________________________________________________________________
AliMUONGeometryModule::AliMUONGeometryModule(const AliMUONGeometryModule& rhs)
  : TObject(rhs)
{
/// Protected copy constructor

  AliFatal("Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONGeometryModule::~AliMUONGeometryModule() 
{
/// Destructor

  delete fSVVolumeIds;
  delete fEnvelopes;
  delete fSVMap;
  delete fTransformer;
}

//______________________________________________________________________________
AliMUONGeometryModule& 
AliMUONGeometryModule::operator = (const AliMUONGeometryModule& rhs) 
{
/// Protected assignement operator

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
/// Return the index of the volume specified by volId
/// if it is present in the list of sensitive volumes 
/// (or -1 if not present).
 
  for (Int_t i=0; i<fNofSVs; i++) {
      if (fSVVolumeIds->At(i) == svVolId) return i;
  }
  return -1;
}

//
// public methods
//

//______________________________________________________________________________
void  AliMUONGeometryModule::SetTransformation(const TGeoCombiTrans& transform)
{
/// Set the module position wrt world.

  fTransformer->SetTransformation(transform);
}  

//______________________________________________________________________________
void AliMUONGeometryModule::SetVolumePath(const TString& volumePath)
{ 
/// Set the volume path to transformer

  fTransformer->SetVolumePath(volumePath);
}

//______________________________________________________________________________
void  AliMUONGeometryModule::SetSensitiveVolume(Int_t svVolId)
{
/// Add the volume specified by volId to the list of sensitive
/// volumes
  
  // Resize TArrayI if needed
  if (fSVVolumeIds->GetSize() == fNofSVs) fSVVolumeIds->Set(2*fNofSVs);

  fSVVolumeIds->AddAt(svVolId, fNofSVs++);
}      

//______________________________________________________________________________
void  AliMUONGeometryModule::SetSensitiveVolume(const TString& volName)
{
/// Add the volume specified by volName to the list of sensitive
/// volumes

  SetSensitiveVolume(gMC->VolId(volName));
}      

//______________________________________________________________________________
void  AliMUONGeometryModule::SetAlign(Bool_t align)
{
/// Set alignement option to envelope store.
  
  fEnvelopes->SetAlign(align);
}  

//______________________________________________________________________________
AliMUONGeometryDetElement* 
AliMUONGeometryModule::FindBySensitiveVolume(const TString& sensVolume) const
{
/// Find TGeoCombiTrans for the detector element Id specified by aligned volume 

  Int_t detElemId = fSVMap->Get(sensVolume);

  if (!detElemId) return 0; 
        // The specified sensitive volume is not in the map   
  
  return fTransformer->GetDetElement(detElemId);
}  

//______________________________________________________________________________
Bool_t AliMUONGeometryModule::IsSensitiveVolume(Int_t volId) const
{
/// Check if the volume specified by volId is present in the list
/// of sensitive volumes.

  for (Int_t i=0; i<fNofSVs; i++) {
      if (fSVVolumeIds->At(i) == volId) return kTRUE;
  }
  return kFALSE;
}

//______________________________________________________________________________
Bool_t AliMUONGeometryModule::IsSensitiveVolume(const TString& volName) const
{
/// Check if the volume specified by volName  is present in the list
/// of sensitive volumes.

  return IsSensitiveVolume(gMC->VolId(volName));
}
