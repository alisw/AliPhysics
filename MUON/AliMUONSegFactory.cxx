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
#include "AliMUONSt345SlatSegmentationV2.h"
#include "AliMUONTriggerSegmentation.h"
#include "AliMUONTriggerSegmentationV2.h"
#include "AliMUONTriggerConstants.h"

#include "AliMpDEManager.h"
#include "AliMpDEIterator.h"

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
      fMpSegFactory(),
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
      fMpSegFactory(),
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
      fMpSegFactory(),
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
  Int_t firstDE = it.CurrentDE(); 
  
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
AliMpVSegmentation* 
AliMUONSegFactory::CreateMpSegmentation(Int_t detElemId, Int_t cath)
{
/// Create mapping segmentation for given detElemId and cath
/// using mapping manager

  AliMpVSegmentation* mpSegmentation 
    = fMpSegFactory.CreateMpSegmentation(detElemId, cath);

  Segmentation()->AddMpSegmentation(mpSegmentation);
  
  return mpSegmentation;
} 
    
//______________________________________________________________________________
AliMUONVGeometryDESegmentation*  
AliMUONSegFactory::CreateDESegmentation(Int_t detElemId, Int_t cath)
{ 
/// Create DE segmentation, operating in local DE reference frame

  // Check detElemId & cath  
  if ( ! AliMpDEManager::IsValid(detElemId, cath, true) ) return 0;
  
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
  TString deName = AliMpDEManager::GetDEName(detElemId, cath);
  TObject* objSegmentation = fDESegmentations.Get(deName);
  if ( objSegmentation ) 
    deSegmentation = (AliMUONVGeometryDESegmentation*)objSegmentation;  
    
  if ( !deSegmentation ) {

    // Get/Create mapping segmentation via mapping manager
    AliMpVSegmentation* mpSegmentation 
      = CreateMpSegmentation(detElemId, cath);
 
    AliMpStationType stationType = AliMpDEManager::GetStationType(detElemId);
    AliMpPlaneType planeType = AliMpDEManager::GetPlaneType(detElemId, cath);
    
    switch (stationType) {

      case kStation1:
      case kStation2:
        deSegmentation = new AliMUONSt12QuadrantSegmentation(
	                         mpSegmentation, stationType, planeType);
        //cout << "   new AliMUONSt12QuadrantSegmentation "
	//     << StationTypeName(stationType) << "  "  
	//     << PlaneTypeName(planeType) << "  "
	//     << deName << endl;
	      				  
        break;
        
      case kStation345:  	          
        deSegmentation = new AliMUONSt345SlatSegmentationV2(
	                         mpSegmentation, detElemId, planeType); 
        //cout << "   new AliMUONSt345SlatSegmentationV2 "			  
	//     << StationTypeName(stationType) << "  "  
	//     << PlaneTypeName(planeType) << "  "
	//     << deName << endl;				  
        break;
    
      case kStationTrigger:  	          
        deSegmentation = new AliMUONTriggerSegmentationV2(
	                         mpSegmentation, detElemId, planeType); 
        //cout << "   new AliMUONTriggerSegmentationV2 "			  
	//     << StationTypeName(stationType) << "  "  
	//     << PlaneTypeName(planeType) << "  "			  
	//     << deName << endl;				  
        break;
    }
    
    // Map new DE segmentation
    fDESegmentations.Add(deName, deSegmentation);
    Segmentation()->AddDESegmentation(deSegmentation);
  }
  
  // Add  DE segmentation to module
  //
  moduleSegmentation->Add(detElemId, deName, deSegmentation);  
  
  return deSegmentation;
}        
  

//______________________________________________________________________________
void
AliMUONSegFactory::CreateModuleSegmentations(Int_t chamberId, Int_t cath)
{ 
/// Create module segmentation(s), operating in the global reference frame
/// Detection elements are defined via DE names map.

  // Check cathod & module Id 
  if ( ! AliMpDEManager::IsValidCathod(cath, true) || 
       ! AliMpDEManager::IsValidChamberId(chamberId, true) ) return;

  AliMpDEIterator it;
  for ( it.First(chamberId); ! it.IsDone(); it.Next() )
    CreateDESegmentation(it.CurrentDE(), cath);
}  
    
//______________________________________________________________________________
AliMUONSegmentation*  
AliMUONSegFactory::CreateSegmentation(const TString& option)
{
/// Create segmentations on all levels and return their container.

  // Check options
  if ( option != "default"   && 
       option != "FactoryV2" && 
       option != "FactoryV3" &&
       option != "FactoryV4" &&
       option != "new") {

    AliErrorStream() << "Option " << option << " not defined." << endl;
    return 0;
  }         
 
  if ( option == "FactoryV2" || option == "FactoryV3" ) { 

    AliErrorStream() 
      << "Option " << option << " not supported anymore." << endl;
    return 0;
  }         

  for (Int_t chamberId = 0; chamberId<AliMUONConstants::NCh(); chamberId++)
    for (Int_t cath = 0; cath < 2; cath++) {
      if ( IsGeometryDefined(chamberId) )
        CreateModuleSegmentations( chamberId, cath);
  }

  return Segmentation();
}

