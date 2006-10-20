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
// $MpId: AliMpSegFactory.cxx,v 1.7 2006/05/24 13:58:34 ivana Exp $
// Category: management

// -----------------------
// Class AliMpSegFactory
// -----------------------
// The factory for building mapping segmentations
// Authors: Ivana Hrivnacova, IPN Orsay
//          Laurent Aphecetche, SUBATECH

#include "AliMpSegFactory.h"

#include "AliMpDEManager.h"
#include "AliMpExMap.h"
#include "AliMpSector.h"
#include "AliMpSectorReader.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpSt345Reader.h"
#include "AliMpTrigger.h"
#include "AliMpTriggerReader.h"
#include "AliMpTriggerSegmentation.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TMap.h>
#include <TObjString.h>
#include <TSystem.h>

/// \cond CLASSIMP
ClassImp(AliMpSegFactory)
/// \endcond

Int_t AliMpSegFactory::fgNumberOfInstances(0);

//______________________________________________________________________________
AliMpSegFactory::AliMpSegFactory()
: TObject(),
  fMpSegmentations(),
  fMpMap(new AliMpExMap(true))
{  
    /// Standard constructor
    AliDebug(1,"");
    fMpMap->SetOwner(true);
    ++fgNumberOfInstances;
}

//______________________________________________________________________________

AliMpSegFactory::~AliMpSegFactory()
{
/// Destructor

  // The segmentations is supposed to be deleted in the client code
  AliDebug(1,"");
  --fgNumberOfInstances;
}

//
// private methods
//

//_____________________________________________________________________________
AliMpExMap*
AliMpSegFactory::FillMpMap(Int_t detElemId)
{
/// Fill the map of electronic cards IDs to segmentations for
/// given detElemId

  AliDebugStream(2) << "detElemId=" << detElemId << endl;;
  
  AliMpExMap* mde = new AliMpExMap(true);
  mde->SetOwner(kFALSE);
  fMpMap->Add(detElemId,mde);
  
  AliMpVSegmentation* seg[2];
  TArrayI ecn[2];
  
  // Do it in 2 steps to be able to set the AliMpExMap size once for all,
  // to avoid annoying warning message in case of dynamical resizing.
  // (not critical).
  for ( Int_t cathode = 0; cathode < 2; ++cathode )
  {
    seg[cathode] = CreateMpSegmentation(detElemId,cathode);
    seg[cathode]->GetAllElectronicCardIDs(ecn[cathode]);
  }
  
  mde->SetSize(ecn[0].GetSize()+ecn[1].GetSize());
  
  for ( Int_t cathode = 0; cathode < 2; ++ cathode )
  {
    for ( Int_t i = 0; i < ecn[cathode].GetSize(); ++i )
    {
      mde->Add(ecn[cathode][i],const_cast<AliMpVSegmentation*>(seg[cathode]));
    }
  }
  
  return mde;
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

  AliDebugStream(3)
    << "Creating segmentation for detElemId=" << detElemId 
    << " cath=" << cath << endl;
  
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

//_____________________________________________________________________________
AliMpVSegmentation* 
AliMpSegFactory::CreateMpSegmentationByElectronics(Int_t detElemId,
                                                   Int_t ecId)
{
/// Create mapping segmentation for given detElemId and electronic card Id
/// (motif position Id) or return it if it was already built

  AliMpExMap* m = static_cast<AliMpExMap*>(fMpMap->GetValue(detElemId));
  
  if (!m)
  {
    m = FillMpMap(detElemId);
  }
  
  return static_cast<AliMpVSegmentation*>(m->GetValue(ecId));
}
    
//______________________________________________________________________________
void AliMpSegFactory::DeleteSegmentations()
{
/// Delete all segmentations created with this manager

  AliDebug(1,"deleting mpSegmentations");
  fMpSegmentations.Clear();
  AliDebug(1,"deleting mpMap");
  delete fMpMap;
  fMpMap = 0; 
  AliDebug(1,"done");
}

