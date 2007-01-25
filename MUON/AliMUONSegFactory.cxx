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

////////////////////////////////////////////////////////////
//  Factory for muon chambers, segmentations and response //
////////////////////////////////////////////////////////////

/* $Id$ */

// --------------------------
// Class AliMUONSegFactory
// --------------------------
// New factory for building segmentations at all levels
// Authors: Ivana Hrivnacova, IPN Orsay

#include "AliMUONSegFactory.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONSegmentation.h"
#include "AliMUONGeometrySegmentation.h"
#include "AliMUONSt12QuadrantSegmentation.h"
#include "AliMUONSt345SlatSegmentation.h"
#include "AliMUONTriggerSegmentation.h"

#include "AliMpDEManager.h"
#include "AliMpDEIterator.h"
#include "AliMpSegmentation.h"
#include "AliMpDetElement.h"
#include "AliMpCathodType.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TMap.h>

/// \cond CLASSIMP
ClassImp(AliMUONSegFactory)
/// \endcond

//______________________________________________________________________________
AliMUONSegFactory::AliMUONSegFactory(const AliMUONGeometryTransformer* geometry)
    : TObject(),
      fDESegmentations(),
      fSegmentation(0),
      fkTransformer(geometry)
{  
/// Standard constructor

}

//______________________________________________________________________________
AliMUONSegFactory::AliMUONSegFactory(const TString& volPathsFileName,
                                     const TString& transformsFileName)
    : TObject(),
      fDESegmentations(),
      fSegmentation(0),
      fkTransformer(0)
{  
/// Standard constructor

  // Transformer
  AliMUONGeometryTransformer* transformer = new AliMUONGeometryTransformer(true);
  transformer->ReadGeometryData(volPathsFileName, transformsFileName);
  fkTransformer = transformer;  
}

//______________________________________________________________________________
  AliMUONSegFactory::AliMUONSegFactory()
    : TObject(),      
      fDESegmentations(),
      fSegmentation(0),
      fkTransformer(0)
{
/// Default constructor
}

//______________________________________________________________________________

AliMUONSegFactory::~AliMUONSegFactory()
{
/// Destructor

  //delete fSegmentation;
        // The segmentation is supposed to be deleted in the client code
}

//
// Private methods
//

//______________________________________________________________________________
Bool_t AliMUONSegFactory::IsGeometryDefined(Int_t ichamberId)
{
/// Return true, if det elements for the chamber with the given ichamber Id
/// are defined in geometry (the geometry builder for this chamber was activated)

  AliMpDEIterator it;
  it.First(ichamberId);
  Int_t firstDE = it.CurrentDEId(); 
  
  if ( ! fkTransformer ||
       ! fkTransformer->GetModuleTransformerByDEId(firstDE, false) )
       
    return kFALSE;
  
  return kTRUE;
}  

//__________________________________________________________________________
AliMUONSegmentation* AliMUONSegFactory::Segmentation()
{ 
/// Return the segmentation container, create it if it does not yet exist

  if ( ! fSegmentation ) 
    fSegmentation = new AliMUONSegmentation(AliMUONConstants::NGeomModules());

  return fSegmentation; 
}

//
// public methods
//

//______________________________________________________________________________
AliMUONVGeometryDESegmentation*  
AliMUONSegFactory::CreateDESegmentation(Int_t detElemId, AliMp::CathodType cath)
{ 
/// Create DE segmentation, operating in local DE reference frame

  // Check detElemId & cath  
  if ( ! AliMpDEManager::IsValidDetElemId(detElemId, true) ) return 0;
  
  // Check if transformer is defined
  if ( ! fkTransformer) {
    AliErrorStream() << "Geometry transformer not defined" << endl;
    return 0;
  }  

  // Only return it, if DE segmentation for this detElemId and cath
  // was already defined
  //
  const AliMUONVGeometryDESegmentation* kdeSegmentation
    = Segmentation()->GetDESegmentation(detElemId, cath, false);
  if ( kdeSegmentation ) 
    return const_cast<AliMUONVGeometryDESegmentation*>(kdeSegmentation);
  
  // Get module, create it if it does not exist 
  //
  AliMUONGeometrySegmentation* moduleSegmentation
    = Segmentation()->GetModuleSegmentationByDEId(detElemId, cath, false);
  if (! moduleSegmentation) {
    Int_t moduleId = AliMpDEManager::GetGeomModuleId(detElemId);	       
    moduleSegmentation 
      = new AliMUONGeometrySegmentation(
               fkTransformer->GetModuleTransformer(moduleId));
    Segmentation()->AddModuleSegmentation(moduleId, cath, moduleSegmentation);
  }        
   
  // Get DE segmentation for this DE type, create it if it does not exist 
  // 
  AliMUONVGeometryDESegmentation* deSegmentation = 0;
  AliMpDetElement* detElement = AliMpDEManager::GetDetElement(detElemId);
  TString deSegName = detElement->GetSegName(cath);
  TObject* objSegmentation = fDESegmentations.Get(deSegName);
  if ( objSegmentation ) 
    deSegmentation = (AliMUONVGeometryDESegmentation*)objSegmentation;  
    
  if ( !deSegmentation ) {

    

    // Get/Create mapping segmentation via mapping manager
    const AliMpVSegmentation* kmpSegmentation 
      = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, cath);
    AliMpVSegmentation* mpSegmentation 
      = const_cast<AliMpVSegmentation*>(kmpSegmentation);

    AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
    AliMp::PlaneType planeType = AliMpDEManager::GetPlaneType(detElemId, cath);
    
    switch (stationType) {

      case AliMp::kStation1:
      case AliMp::kStation2:
        deSegmentation = new AliMUONSt12QuadrantSegmentation(
	                         mpSegmentation, stationType, planeType);
        //cout << "   new AliMUONSt12QuadrantSegmentation "
	//     << StationTypeName(stationType) << "  "  
	//     << PlaneTypeName(planeType) << "  "
	//     << deName << endl;
	      				  
        break;
        
      case AliMp::kStation345:  	          
        deSegmentation = new AliMUONSt345SlatSegmentation(
	                         mpSegmentation, detElemId, planeType); 
        //cout << "   new AliMUONSt345SlatSegmentationV2 "			  
	//     << StationTypeName(stationType) << "  "  
	//     << PlaneTypeName(planeType) << "  "
	//     << deName << endl;				  
        break;
    
      case AliMp::kStationTrigger:  	          
        deSegmentation = new AliMUONTriggerSegmentation(
	                         mpSegmentation, detElemId, planeType); 
        //cout << "   new AliMUONTriggerSegmentation "			  
	//     << StationTypeName(stationType) << "  "  
	//     << PlaneTypeName(planeType) << "  "			  
	//     << deName << endl;				  
        break;
    }
    
    // Map new DE segmentation
    fDESegmentations.Add(deSegName, deSegmentation);
    Segmentation()->AddDESegmentation(deSegmentation);
  }
  
  // Add  DE segmentation to module
  //
  moduleSegmentation->Add(detElemId, deSegName, deSegmentation);  
  
  return deSegmentation;
}        
  

//______________________________________________________________________________
void
AliMUONSegFactory::CreateModuleSegmentations(Int_t chamberId, AliMp::CathodType cath)
{ 
/// Create module segmentation(s), operating in the global reference frame
/// Detection elements are defined via DE names map.

  // Check module Id 
  if ( ! AliMpDEManager::IsValidChamberId(chamberId, true) ) return;

  AliMpDEIterator it;
  for ( it.First(chamberId); ! it.IsDone(); it.Next() )
    CreateDESegmentation(it.CurrentDEId(), cath);
}  
    
//______________________________________________________________________________
AliMUONSegmentation*  
AliMUONSegFactory::CreateSegmentation()
{
/// Create segmentations on all levels and return their container.

  for (Int_t chamberId = 0; chamberId<AliMUONConstants::NCh(); chamberId++)
    for (Int_t cath = AliMp::kCath0; cath <= AliMp::kCath1; cath++) {
      if ( IsGeometryDefined(chamberId) )
        CreateModuleSegmentations( chamberId, AliMp::GetCathodType(cath));
  }

  return Segmentation();
}

