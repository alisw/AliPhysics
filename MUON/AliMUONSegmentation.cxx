/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *      SigmaEffect_thetadegrees                                                                  *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
//
// ----------------------------
// Class AliMUONSegmentation
// ----------------------------
// Manager class for geometry construction via geometry builders.
// Author: Ivana Hrivnacova, IPN Orsay

#include <iostream>

#include <TObjArray.h>

#include "AliMUONSegmentation.h"
#include "AliMUONVGeometryDESegmentation.h"
#include "AliMUONGeometrySegmentation.h"
\
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"

#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONSegmentation)
/// \endcond
 
//______________________________________________________________________________
AliMUONSegmentation::AliMUONSegmentation(Int_t nofModules)
  : TObject(),
    fMpSegmentations(0),
    fDESegmentations(0)
{
/// Standard constructor

  // Create array for mapping segmentations
  fMpSegmentations = new TObjArray();
  fMpSegmentations->SetOwner(kTRUE);

  // Create array for DE segmentations
  fDESegmentations = new TObjArray();
  fDESegmentations->SetOwner(kTRUE);

  // Create array for modules segmentations
  for (Int_t cathod = 0; cathod < 2; cathod++) {
    fModuleSegmentations[cathod] = new TObjArray(nofModules);
    fModuleSegmentations[cathod]->SetOwner(true);
    
    for (Int_t i=0; i<nofModules; i++)
      fModuleSegmentations[cathod]->AddAt(0, i);
  }  

  AliDebug(1, Form("ctor this = %p", this) ); 
}

//______________________________________________________________________________
AliMUONSegmentation::AliMUONSegmentation() 
  : TObject(),
    fMpSegmentations(0),
    fDESegmentations(0)
{
/// Default constructor

  fModuleSegmentations[0] = 0;
  fModuleSegmentations[1] = 0;

  AliDebug(1, Form("default (empty) ctor this = %p", this));
} 

//______________________________________________________________________________
AliMUONSegmentation::~AliMUONSegmentation()
{
/// Destructor

  AliDebug(1, Form("dtor this = %p", this));

  delete fMpSegmentations;
  delete fDESegmentations;
  delete fModuleSegmentations[0];
  delete fModuleSegmentations[1];
}

//
// private functions
//

//_____________________________________________________________________________
AliMUONGeometrySegmentation* 
AliMUONSegmentation::GetModuleSegmentation(
                        Int_t moduleId, Int_t cathod, Bool_t warn) const
{
/// Return the geometry module segmentation specified by moduleId

  if (cathod < 0 || cathod >= 2) {
    if (warn) {
      AliWarningStream() 
        << "Cathod: " << cathod << " outside limits" << std::endl;
    }			 
    return 0;  
  }  

  if (moduleId < 0 || moduleId >= fModuleSegmentations[cathod]->GetEntriesFast()) {
    if (warn) {
      AliWarningStream() 
        << "Index: " << moduleId << " outside limits" << std::endl;
    }			 
    return 0;  
  }  

  return (AliMUONGeometrySegmentation*) 
            fModuleSegmentations[cathod]->At(moduleId);
}    

//
// public functions
//

//_____________________________________________________________________________
void AliMUONSegmentation::AddMpSegmentation(AliMpVSegmentation* segmentation)
{
/// Add the mapping segmentation to the array if not present

  Bool_t isPresent = false;
  for (Int_t i=0; i<fMpSegmentations->GetEntries(); i++) 
    if ( (AliMpVSegmentation*)fMpSegmentations->At(i) == segmentation ) {
      isPresent = true;
      break;
    }  

  if (!isPresent) fMpSegmentations->Add(segmentation);
}

//_____________________________________________________________________________
void AliMUONSegmentation::AddDESegmentation(
                                AliMUONVGeometryDESegmentation* segmentation)
{
/// Add the DE segmentation to the array

  fDESegmentations->Add(segmentation);
  
  // Deregister the mapping segmentation contained in DE segmentation
  // from fMpSegmentations, if present
  const AliMpVSegmentation* kmpSeg = segmentation->GetMpSegmentation();
  
  for (Int_t i=0; i<fMpSegmentations->GetEntries(); i++) 
    if ( (const AliMpVSegmentation*)fMpSegmentations->At(i) == kmpSeg ) {
      fMpSegmentations->RemoveAt(i);
      break;
    }  
}

//_____________________________________________________________________________
void AliMUONSegmentation::AddModuleSegmentation(Int_t moduleId, Int_t cathod,
                             AliMUONGeometrySegmentation* segmentation)
{
/// Add the module segmentation to the array

  if (cathod < 0 || cathod >= 2) {
      AliErrorStream() 
        << "Cathod: " << cathod << " outside limits" << std::endl;
      return;	
  }			 

  if (moduleId < 0 || moduleId >= fModuleSegmentations[cathod]->GetEntriesFast()) {
    AliErrorStream() 
        << "Module Id: " << moduleId << " outside limits" << std::endl;
    return;  
  }  

  fModuleSegmentations[cathod]->AddAt(segmentation, moduleId);
}

//_____________________________________________________________________________
void  AliMUONSegmentation::Init()
{
/// Initialize all segmentations
 
  for (Int_t cathod = 0; cathod < 2; cathod++) {
    for (Int_t i = 0; i < fModuleSegmentations[cathod]->GetEntriesFast(); i++) {
    
      AliMUONGeometrySegmentation* moduleSegmentation
        = (AliMUONGeometrySegmentation*)fModuleSegmentations[cathod]->At(i);
    
      if (moduleSegmentation) moduleSegmentation->Init(i); 
    }  
  }  
}
			    
//_____________________________________________________________________________
AliMUONGeometrySegmentation* 
AliMUONSegmentation::GetModuleSegmentationByDEId(
                        Int_t detElemId, Int_t cathod, Bool_t warn) const
{
/// Return the geometry module specified by detElemId/cathod

  // Get moduleId 
  Int_t moduleId = AliMpDEManager::GetGeomModuleId(detElemId);

  return GetModuleSegmentation(moduleId, cathod, warn);
}    

//_____________________________________________________________________________
const AliMUONVGeometryDESegmentation* 
AliMUONSegmentation::GetDESegmentation(
                        Int_t detElemId, Int_t cathod, Bool_t warn) const
{
/// Return the DE segmentation specified by detElemId/cathod

  // Get geometry segmentation 
  AliMUONGeometrySegmentation* moduleSegmentation
    = GetModuleSegmentationByDEId(detElemId, cathod, warn);
    
  if ( !moduleSegmentation ) return 0; 
  
  return moduleSegmentation->GetDESegmentation(detElemId, warn);
}    

//_____________________________________________________________________________
const AliMpVSegmentation*
AliMUONSegmentation::GetMpSegmentation(
                     Int_t detElemId, Int_t cathod, Bool_t warn) const
{		     
/// Return the mapping segmentation specified by detElemId/cathod


  // Get DE segmentation 
  const AliMUONVGeometryDESegmentation* kdeSegmentation
    = GetDESegmentation(detElemId, cathod, warn);
    
  if ( !kdeSegmentation ) return 0; 
  
  return kdeSegmentation->GetMpSegmentation();
}    

//_____________________________________________________________________________
Bool_t 
AliMUONSegmentation::HasDE(Int_t detElemId, Int_t cathod) const
{
/// Return true if segmentation for detElemId and cathod is defined.

  const AliMUONVGeometryDESegmentation* kdeSegmentation
    = GetDESegmentation(detElemId, cathod, false);
    
  return ( kdeSegmentation != 0 ); 
  
}

//_____________________________________________________________________________
TString 
AliMUONSegmentation::GetDEName(Int_t detElemId, Int_t cathod) const
{
/// Get detection element name 

  AliMUONGeometrySegmentation* moduleSegmentation
    = GetModuleSegmentationByDEId(detElemId, cathod, true);
    
  return  moduleSegmentation->GetDEName(detElemId); 
}


