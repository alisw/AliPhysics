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
// $MpId: AliMpSegmentation.cxx,v 1.7 2006/05/24 13:58:34 ivana Exp $
// Category: management

// -----------------------
// Class AliMpSegmentation
// -----------------------
// Singleton container class for mapping segmentations
// Authors: Ivana Hrivnacova, IPN Orsay
//          Laurent Aphecetche, SUBATECH

#include "AliMpSegmentation.h"

#include "AliMpDetElement.h"
#include "AliMpDEManager.h"
#include "AliMpDEIterator.h"
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
#include "AliMpCathodType.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TMap.h>
#include <TObjString.h>
#include <TSystem.h>

/// \cond CLASSIMP
ClassImp(AliMpSegmentation)
/// \endcond

AliMpSegmentation* AliMpSegmentation::fgInstance = 0;

//
// static methods
//

//______________________________________________________________________________
AliMpSegmentation* AliMpSegmentation::Instance()
{
/// Create the sementation if it does not yet exist
/// and return its instance

  if ( ! fgInstance )
    fgInstance = new AliMpSegmentation();
    
  return fgInstance;
}    

//
// ctors, dtor
//

//______________________________________________________________________________
AliMpSegmentation::AliMpSegmentation()
: TObject(),
  fMpSegmentations(true),
  fElCardsMap(true),
  fSlatMotifMap()
{  
/// Standard constructor

  AliDebug(1,"");
  fElCardsMap.SetOwner(true);

  // Create mapping segmentations for all detection elements
  for ( Int_t cath = AliMp::kCath0; cath <= AliMp::kCath1; cath ++ ) { 
    AliMpDEIterator it;
    for ( it.First(); ! it.IsDone(); it.Next() ) {
      CreateMpSegmentation(it.CurrentDEId(), AliMp::GetCathodType(cath));
    } 
  }  

  // Fill el cards map for all detection elements
  // of tracking chambers
  AliMpDEIterator it;
  for ( it.First(); ! it.IsDone(); it.Next() ) { 
    if ( AliMpDEManager::GetStationType(it.CurrentDEId()) != AliMp::kStationTrigger ) {
      FillElCardsMap(it.CurrentDEId());
    }
  }  
}

//______________________________________________________________________________
AliMpSegmentation::AliMpSegmentation(TRootIOCtor* /*ioCtor*/)
: TObject(),
  fMpSegmentations(),
  fElCardsMap(),
  fSlatMotifMap()
{  
/// Constructor for IO

  AliDebug(1,"");

  fgInstance = this;
}

//______________________________________________________________________________
AliMpSegmentation::~AliMpSegmentation()
{
/// Destructor

  AliDebug(1,"");

  // Segmentations are deleted with fMpSegmentations 
  // El cards arrays are deleted with fElCardsMap
  
  fgInstance = 0;
}

//
// private methods
//

//______________________________________________________________________________
AliMpVSegmentation* 
AliMpSegmentation::CreateMpSegmentation(Int_t detElemId, AliMp::CathodType cath)
{
/// Create mapping segmentation for given detElemId and cath
/// or return it if it was already built

  // Check detElemId & cath  
  if ( ! AliMpDEManager::IsValidDetElemId(detElemId, true) ) return 0;

  // If segmentation is already built, just return it
  //
  AliMpDetElement* detElement = AliMpDEManager::GetDetElement(detElemId);
  TString deSegName = detElement->GetSegName(cath);
  TObject* object = fMpSegmentations.Get(deSegName);
  if ( object ) return (AliMpVSegmentation*)object;

  AliDebugStream(3)
    << "Creating segmentation for detElemId=" << detElemId 
    << " cath=" << cath << endl;
  
  // Read mapping data and create segmentation
  //
  AliMp::StationType stationType = detElement->GetStationType();
  AliMp::PlaneType planeType = detElement->GetPlaneType(cath);
  TString deTypeName = detElement->GetSegType();

  AliMpVSegmentation* mpSegmentation = 0;

  if ( stationType == AliMp::kStation1 || stationType == AliMp::kStation2 ) {
    AliMpSectorReader reader(stationType, planeType);
    AliMpSector* sector = reader.BuildSector();
    mpSegmentation = new AliMpSectorSegmentation(sector, true);
  }
  else if ( stationType == AliMp::kStation345 ) { 
    AliMpSt345Reader reader(fSlatMotifMap);
    AliMpSlat* slat = reader.ReadSlat(deTypeName, planeType);
    mpSegmentation =  new AliMpSlatSegmentation(slat, true);
  }
  else if ( stationType == AliMp::kStationTrigger ) {
    AliMpTriggerReader reader(fSlatMotifMap);
    AliMpTrigger* trigger = reader.ReadSlat(deTypeName, planeType);
    mpSegmentation = new AliMpTriggerSegmentation(trigger, true);
  }
  else   
    AliErrorStream() << "Unknown station type" << endl;

  fMpSegmentations.Add(deSegName, mpSegmentation); 
  
  StdoutToAliDebug(3, fSlatMotifMap.Print(););
  
  return mpSegmentation;
} 

//_____________________________________________________________________________
AliMpExMap*
AliMpSegmentation::FillElCardsMap(Int_t detElemId)
{
/// Fill the map of electronic cards IDs to segmentations for
/// given detElemId

  AliDebugStream(2) << "detElemId=" << detElemId << endl;;
  
  AliMpExMap* mde = new AliMpExMap(true);
  mde->SetOwner(kFALSE);
  fElCardsMap.Add(detElemId,mde);

  const AliMpVSegmentation* seg[2];
  TArrayI ecn[2];
  
  // Do it in 2 steps to be able to set the AliMpExMap size once for all,
  // to avoid annoying warning message in case of dynamical resizing.
  // (not critical).
  for ( Int_t cathode = AliMp::kCath0; cathode <= AliMp::kCath1; ++cathode )
  {
    seg[cathode] = GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));
    seg[cathode]->GetAllElectronicCardIDs(ecn[cathode]);
  }
  
  mde->SetSize(ecn[0].GetSize()+ecn[1].GetSize());
  
  for ( Int_t cathode = AliMp::kCath0; cathode <= AliMp::kCath1; ++ cathode )
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
const AliMpVSegmentation* 
AliMpSegmentation::GetMpSegmentation(
                      Int_t detElemId, AliMp::CathodType cath, Bool_t warn) const
{
/// Return mapping segmentation for given detElemId and cath

  // Check detElemId & cath  
  if ( ! AliMpDEManager::IsValidDetElemId(detElemId, false) ) {
    
    if ( warn ) {
      AliWarningStream() 
        << "Invalid detElemId " << detElemId << endl;
    }	
    return 0;
  }  

  // If segmentation is already built, just return it
  //
  AliMpDetElement* detElement = AliMpDEManager::GetDetElement(detElemId);
  TString deSegName = detElement->GetSegName(cath);
  TObject* object = fMpSegmentations.Get(deSegName);
  if ( ! object ) {
    // Should never happen
    AliErrorStream() 
      << "Segmentation for detElemId/cathod " 
	<< detElemId << ", " << cath << " not defined" << endl;
    return 0;
  }  
  
  return static_cast<AliMpVSegmentation*>(object);
} 

//_____________________________________________________________________________
const AliMpVSegmentation* 
AliMpSegmentation::GetMpSegmentationByElectronics(
                      Int_t detElemId, Int_t ecId, Bool_t warn) const
{
/// Return mapping segmentation for given detElemId and electronic card Id
/// (motif position Id)

  AliMpExMap* m = static_cast<AliMpExMap*>(fElCardsMap.GetValue(detElemId));
  
  if (!m) {
    // Should never happen
    AliErrorStream() 
      << "Cannot find the el cards map for detElemId " << detElemId << endl;
    return 0;
  }  

  TObject* object = m->GetValue(ecId);
  if ( ! object ) {
    if ( warn ) {
      AliErrorStream() 
        << "Segmentation for electronic card " 
	<< ecId << " not found" << endl;
    }	
    return 0;
  }  
   
  return static_cast<AliMpVSegmentation*>(object);
}
