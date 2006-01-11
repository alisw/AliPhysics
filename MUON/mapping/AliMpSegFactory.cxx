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

#include "AliMpSegFactory.h"
#include "AliMpDEManager.h"
#include "AliMpSector.h"
#include "AliMpSectorReader.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpSt345Reader.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpTrigger.h"
#include "AliMpTriggerReader.h"
#include "AliMpTriggerSegmentation.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TMap.h>

ClassImp(AliMpSegFactory)

//______________________________________________________________________________
AliMpSegFactory::AliMpSegFactory()
    : TObject(),
      fMpSegmentations()
{  
/// Standard constructor
}

//______________________________________________________________________________
AliMpSegFactory::AliMpSegFactory(const AliMpSegFactory& rhs)
 : TObject(rhs)
{
/// Protected copy constructor

  AliFatal("Not implemented.");
}

//______________________________________________________________________________

AliMpSegFactory::~AliMpSegFactory()
{
/// Destructor

  // The segmentations is supposed to be deleted in the client code
}

//______________________________________________________________________________
AliMpSegFactory&  AliMpSegFactory::operator=(const AliMpSegFactory& rhs)
{
  // Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}    
          
//
// public methods
//

//______________________________________________________________________________
AliMpVSegmentation* 
AliMpSegFactory::CreateMpSegmentation(Int_t detElemId, Int_t cath)
{
/// Create mapping segmentation for given detElemId and cath
/// or return it if it was already built

  // Check detElemId & cath  
  if ( ! AliMpDEManager::IsValid(detElemId, cath, true) ) return 0;

  // If segmentation is already built, just return it
  //
  TString deName = AliMpDEManager::GetDEName(detElemId, cath);
  TObject* object = fMpSegmentations.Get(deName);
  if ( object ) return (AliMpVSegmentation*)object;

  // Read mapping data and create segmentation
  //
  AliMpStationType stationType = AliMpDEManager::GetStationType(detElemId);
  AliMpPlaneType planeType = AliMpDEManager::GetPlaneType(detElemId, cath);
  TString deTypeName = AliMpDEManager::GetDETypeName(detElemId, cath);

  AliMpVSegmentation* mpSegmentation = 0;

  if ( stationType == kStation1 || stationType == kStation2 ) {
    AliMpSectorReader reader(stationType, planeType);
    AliMpSector* sector = reader.BuildSector();
    mpSegmentation = new AliMpSectorSegmentation(sector);
  }
  else if ( stationType == kStation345 ) { 
    AliMpSlat* slat = AliMpSt345Reader::ReadSlat(deTypeName, planeType);
    mpSegmentation =  new AliMpSlatSegmentation(slat);
  }
  else if ( stationType == kStationTrigger ) {
    AliMpTrigger* trigger = AliMpTriggerReader::ReadSlat(deTypeName, planeType);
    mpSegmentation = new AliMpTriggerSegmentation(trigger);
  }
  else   
    AliErrorStream() << "Unknown station type" << endl;

  fMpSegmentations.Add(deName, mpSegmentation); 
  return mpSegmentation;
} 
    
//______________________________________________________________________________
void AliMpSegFactory::Delete(const Option_t* /*opt*/)
{
/// Delete all segmentations created with this manager

  fMpSegmentations.Clear();
}

