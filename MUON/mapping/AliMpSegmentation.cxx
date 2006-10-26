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
  fElCardsMap(true)
{  
/// Standard constructor

  AliDebug(1,"");
  fElCardsMap.SetOwner(true);

  // Create mapping segmentations for all detection elements
  for ( Int_t cath = 0; cath < 2; cath ++ ) { 
    AliMpDEIterator it;
    for ( it.First(); ! it.IsDone(); it.Next() ) {
      //if ( AliMpDEManager::GetChamberId(it.CurrentDE()) >= 4 ) break;
      CreateMpSegmentation(it.CurrentDE(), cath);
    } 
  }  

  // Fill el cards map for all detection elements
  // of tracking chambers
  AliMpDEIterator it;
  for ( it.First(); ! it.IsDone(); it.Next() ) { 
    //if ( AliMpDEManager::GetChamberId(it.CurrentDE()) >= 4 ) break;
    if ( AliMpDEManager::GetChamberId(it.CurrentDE()) >= 10 ) break;
    FillElCardsMap(it.CurrentDE());
  }  
}

//______________________________________________________________________________
AliMpSegmentation::AliMpSegmentation(TRootIOCtor* /*ioCtor*/)
: TObject(),
  fMpSegmentations(),
  fElCardsMap()
{  
/// Constructor for IO

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
AliMpSegmentation::CreateMpSegmentation(Int_t detElemId, Int_t cath)
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
    mpSegmentation = new AliMpSectorSegmentation(sector, true);
  }
  else if ( stationType == kStation345 ) { 
    AliMpSlat* slat = AliMpSt345Reader::ReadSlat(deTypeName, planeType);
    mpSegmentation =  new AliMpSlatSegmentation(slat, true);
  }
  else if ( stationType == kStationTrigger ) {
    AliMpTrigger* trigger = AliMpTriggerReader::ReadSlat(deTypeName, planeType);
    mpSegmentation = new AliMpTriggerSegmentation(trigger, true);
  }
  else   
    AliErrorStream() << "Unknown station type" << endl;

  fMpSegmentations.Add(deName, mpSegmentation); 
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
  for ( Int_t cathode = 0; cathode < 2; ++cathode )
  {
    seg[cathode] = GetMpSegmentation(detElemId,cathode);
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
const AliMpVSegmentation* 
AliMpSegmentation::GetMpSegmentation(
                      Int_t detElemId, Int_t cath, Bool_t warn) const
{
/// Return mapping segmentation for given detElemId and cath

  // Check detElemId & cath  
  if ( ! AliMpDEManager::IsValid(detElemId, cath, false) ) {
    
    if ( warn ) {
      AliWarningStream() 
        << "Invalid detElemId/cathod (" 
	<< detElemId << ", " << cath << ")" << endl;
    }	
    return 0;
  }  

  TString deName = AliMpDEManager::GetDEName(detElemId, cath);
  TObject* object = fMpSegmentations.Get(deName);
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
