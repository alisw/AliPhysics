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

//-----------------------------------------------------------------------------
// Class AliMpSegmentation
// -----------------------
// Singleton container class for mapping segmentations
// Authors: Ivana Hrivnacova, IPN Orsay
//          Laurent Aphecetche, SUBATECH
//-----------------------------------------------------------------------------

#include "AliMpSegmentation.h"

#include "AliMpDataStreams.h"
#include "AliMpDetElement.h"
#include "AliMpDEStore.h"
#include "AliMpDEManager.h"
#include "AliMpDEIterator.h"
#include "AliMpExMap.h"
#include "AliMpFastSegmentation.h"
//#include "AliMpFastSegmentationV2.h"
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
#include "AliMpSlatMotifMap.h"


#include "AliLog.h"

#include <Riostream.h>
#include <TMap.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TClass.h>

#include <cassert>

/// \cond CLASSIMP
ClassImp(AliMpSegmentation)
/// \endcond

AliMpSegmentation* AliMpSegmentation::fgInstance = 0;

//
// static methods
//

//______________________________________________________________________________
AliMpSegmentation* AliMpSegmentation::Instance(Bool_t warn)
{
/// Return its instance

  if ( ! fgInstance && warn ) {
    AliWarningClass("Segmentation has not been loaded");
  }  
    
  return fgInstance;
}    

//______________________________________________________________________________
AliMpSegmentation* AliMpSegmentation::ReadData(const AliMpDataStreams& dataStreams,
                                               Bool_t warn)
{
/// Load the sementation from ASCII data files
/// and return its instance

  if ( fgInstance ) {
    if ( warn )
      AliWarningClass("Segmentation has been already loaded");
    return fgInstance;
  }  
  
  if ( dataStreams.GetReadFromFiles() )
    AliInfoClass("Reading segmentation from ASCII files.");

  fgInstance = new AliMpSegmentation(dataStreams);
  return fgInstance;
}    

//
// ctors, dtor
//

//______________________________________________________________________________
AliMpSegmentation::AliMpSegmentation(const AliMpDataStreams& dataStreams)
: TObject(),
  fkDataStreams(dataStreams),
  fDetElements(0),
  fMpSegmentations(true),
  fElCardsMap(),
  fSlatMotifMap(new AliMpSlatMotifMap)
{  
/// Standard constructor - segmentation is loaded from ASCII data files

  AliDebug(1,"");
  
  fElCardsMap.SetOwner(kTRUE);
  
  // Load DE data
  if ( ! AliMpDEStore::Instance(false) )  
    AliMpDEStore::ReadData(dataStreams);
  fDetElements = AliMpDEStore::Instance();  

  // Create mapping segmentations for all detection elements
  AliMpDEIterator it;
  for ( it.First(); ! it.IsDone(); it.Next() ) 
  {
    Int_t n(0);
    
    for ( Int_t cath = AliMp::kCath0; cath <= AliMp::kCath1; ++cath ) 
    { 
      if ( CreateMpSegmentation(it.CurrentDEId(), AliMp::GetCathodType(cath)) ) ++n;
    }
     
    if ( n == 2 &&  // should always be the case except when we play with the CreateMpSegmentation method...
        AliMpDEManager::GetStationType(it.CurrentDEId()) != AliMp::kStationTrigger ) // only for tracker
    {
        // Fill el cards map for all detection elements
        // of tracking chambers
        FillElCardsMap(it.CurrentDEId());
    }      
  } 
}

//______________________________________________________________________________
AliMpSegmentation::AliMpSegmentation(TRootIOCtor* ioCtor)
: TObject(),
  fkDataStreams(ioCtor),
  fDetElements(0),
  fMpSegmentations(),
  fElCardsMap(ioCtor),
  fSlatMotifMap(0)
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

  delete fDetElements;

  // Segmentations are deleted with fMpSegmentations 
  // El cards arrays are deleted with fElCardsMap
  
  delete fSlatMotifMap;
  
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

  if ( stationType == AliMp::kStation12 ) {
    AliMq::Station12Type station12Type = detElement->GetStation12Type();
    AliMpSectorReader reader(fkDataStreams, station12Type, planeType);
    AliMpSector* sector = reader.BuildSector();
    mpSegmentation = new AliMpFastSegmentation(new AliMpSectorSegmentation(sector, true));
  }
  else if ( stationType == AliMp::kStation345 ) { 
    AliMpSt345Reader reader(fkDataStreams,fSlatMotifMap);
    AliMpSlat* slat = reader.ReadSlat(deTypeName, planeType);
    mpSegmentation = new AliMpFastSegmentation(new AliMpSlatSegmentation(slat, true));
  }
  else if ( stationType == AliMp::kStationTrigger ) {
    AliMpTriggerReader reader(fkDataStreams,fSlatMotifMap);
    AliMpTrigger* trigger = reader.ReadSlat(deTypeName, planeType);
    mpSegmentation = new AliMpTriggerSegmentation(trigger, true);
  }
  else   
    AliErrorStream() << "Unknown station type" << endl;

  if ( mpSegmentation ) fMpSegmentations.Add(deSegName, mpSegmentation); 
  
//  StdoutToAliDebug(3, fSlatMotifMap.Print(););
  
  return mpSegmentation;
} 

//_____________________________________________________________________________
AliMpExMap*
AliMpSegmentation::FillElCardsMap(Int_t detElemId)
{
/// Fill the map of electronic cards IDs to segmentations for
/// given detElemId

  AliDebugStream(2) << "detElemId=" << detElemId << endl;;
  
  AliMpExMap* mde = new AliMpExMap;
  mde->SetOwner(kFALSE);
  fElCardsMap.Add(detElemId,mde);

  const AliMpVSegmentation* seg[2];
  TArrayI ecn[2];
  
  for ( Int_t cathode = AliMp::kCath0; cathode <= AliMp::kCath1; ++cathode )
  {
    seg[cathode] = GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));
    seg[cathode]->GetAllElectronicCardIDs(ecn[cathode]);
    for ( Int_t i = 0; i < ecn[cathode].GetSize(); ++i )
    {
      mde->Add(ecn[cathode][i],const_cast<AliMpVSegmentation*>(seg[cathode]));
    }
  }
  
  assert( mde->GetSize() > 0 );
  
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

//_____________________________________________________________________________
const AliMpSector*  
AliMpSegmentation::GetSector(const AliMpVSegmentation* kSegmentation, 
                             Bool_t warn) const
{
/// Return sector for given mapping segmentation.
/// If segmentation is not of sector type, zero is returned 
/// and an Error is issued if warn is set true (default). 

  if ( ! kSegmentation ) return 0;

  if ( kSegmentation->StationType() != AliMp::kStation12 ) {
    if ( warn ) {
      AliErrorStream() 
        << "Segmentation is not of sector type" << endl;
     }   
     return 0;
  }
    
  // If fast segmentation
  const AliMpFastSegmentation* fseg 
    = dynamic_cast<const AliMpFastSegmentation*>(kSegmentation);
  if ( fseg ) {   
    return 
      static_cast<const AliMpSectorSegmentation*>(fseg->GetHelper())->GetSector();
  }
  
  // If fast segmentation V2
/*  
  const AliMpFastSegmentationV2* fsegV2 
    = dynamic_cast<const AliMpFastSegmentationV2*>(kSegmentation);
  if ( fsegV2 ) {   
    return 
      static_cast<const AliMpSectorSegmentation*>(fsegV2->GetHelper())->GetSector();
  }
*/  
  // If sector segmentation
  const AliMpSectorSegmentation* sseg 
    = dynamic_cast<const AliMpSectorSegmentation*>(kSegmentation);
  if ( sseg ) {   
    return sseg->GetSector();
  }
  
  // Should not get to this line
  AliErrorStream() << "Segemntation type not identified." << endl;
  return 0;         
}
                             
//_____________________________________________________________________________
const AliMpSector*  
AliMpSegmentation::GetSector(Int_t detElemId, AliMp::CathodType cath, 
                             Bool_t warn) const
{
/// Return sector for given detElemId and cath.
/// If segmentation is not of sector type, zero is returned 
/// and an Error is issued if warn is set true (default). 

  return GetSector(GetMpSegmentation(detElemId, cath, warn), warn);
}    
      
//_____________________________________________________________________________
const AliMpSector*  
AliMpSegmentation::GetSectorByElectronics(Int_t detElemId, Int_t elCardID, 
                             Bool_t warn) const
{
/// Return sector for given detElemId and elCardID.
/// If segmentation is not of sector type, zero is returned 
/// and an Error is issued if warn is set true (default). 

  return GetSector(GetMpSegmentationByElectronics(detElemId, elCardID, warn), warn);
}    
      
//_____________________________________________________________________________
const AliMpSlat*    
AliMpSegmentation::GetSlat(const AliMpVSegmentation* kSegmentation, 
                           Bool_t warn) const
{                           
/// Return slat for given mapping segmentation.
/// If segmentation is not of slat type, zero is returned 
/// and an Error is issued if warn is set true (default). 

  if ( ! kSegmentation ) return 0;
 
  if ( kSegmentation->StationType() != AliMp::kStation345 ) {
     if ( warn ) {
       AliErrorStream() 
         << "Segmentation is not of slat type" << endl;
     }    
     return 0;
  }
  
  // If fast segmentation
  const AliMpFastSegmentation* fseg 
    = dynamic_cast<const AliMpFastSegmentation*>(kSegmentation);
  if ( fseg ) {   
    return 
      static_cast<const AliMpSlatSegmentation*>(fseg->GetHelper())->Slat();
  }
  
  // If fast segmentation V2
/*
  const AliMpFastSegmentationV2* fsegV2 
    = dynamic_cast<const AliMpFastSegmentationV2*>(kSegmentation);
  if ( fsegV2 ) {   
    return 
      static_cast<const AliMpSlatSegmentation*>(fsegV2->GetHelper())->Slat();
  }
*/  
  // If slat segmentation
  const AliMpSlatSegmentation* sseg 
    = dynamic_cast<const AliMpSlatSegmentation*>(kSegmentation);
    
  if ( sseg ) {   
    return sseg->Slat();
  }
  
  // Should not get to this line
  AliErrorStream() << "Segemntation type not identified." << endl;
  return 0;         
}    
                           
//_____________________________________________________________________________
const AliMpSlat*  
AliMpSegmentation::GetSlat(Int_t detElemId, AliMp::CathodType cath, 
                           Bool_t warn) const
{
/// Return slat for given detElemId and cath.
/// If segmentation is not of slat type, zero is returned 
/// and an Error is issued if warn is set true (default). 

  return GetSlat(GetMpSegmentation(detElemId, cath, warn), warn);
}    

//_____________________________________________________________________________
const AliMpSlat*  
AliMpSegmentation::GetSlatByElectronics(Int_t detElemId, Int_t elCardID, 
                           Bool_t warn) const
{
/// Return slat for given detElemId and elCardID.
/// If segmentation is not of slat type, zero is returned 
/// and an Error is issued if warn is set true (default). 

  return GetSlat(GetMpSegmentationByElectronics(detElemId, elCardID, warn), warn);
}    

//_____________________________________________________________________________
const AliMpTrigger* 
AliMpSegmentation::GetTrigger(const AliMpVSegmentation* kSegmentation, 
                              Bool_t warn) const
{                               
/// Return trigger for given mapping segmentation.
/// If segmentation is not of trigger type, zero is returned 
/// and an Error is issued if warn is set true (default). 

  if ( ! kSegmentation ) return 0;
 
  if ( kSegmentation->StationType() != AliMp::kStationTrigger ) {
     if ( warn ) {
       AliErrorStream() 
         << "Segmentation is not of trigger type" << endl;
     }    
     return 0;
  }
  
  // If slat segmentation
  const AliMpTriggerSegmentation* tseg 
    = dynamic_cast<const AliMpTriggerSegmentation*>(kSegmentation);
  if ( tseg ) {   
    return tseg->Slat();
  }
  
  // If fast segmentation
  const AliMpFastSegmentation* fseg 
    = dynamic_cast<const AliMpFastSegmentation*>(kSegmentation);
    
  if ( fseg ) {   
    return 
      static_cast<const AliMpTriggerSegmentation*>(fseg->GetHelper())->Slat();
  }
  
  // If fast segmentation V2
/*  
  const AliMpFastSegmentationV2* fsegV2 
    = dynamic_cast<const AliMpFastSegmentationV2*>(kSegmentation);
    
  if ( fsegV2 ) {   
    return 
      static_cast<const AliMpTriggerSegmentation*>(fsegV2->GetHelper())->Slat();
  }
*/  
  
  // Should not get to this line
  AliErrorStream() << "Segemntation type not identified." << endl;
  return 0;         
}    

//_____________________________________________________________________________
const AliMpTrigger*  
AliMpSegmentation::GetTrigger(Int_t detElemId, AliMp::CathodType cath, 
                              Bool_t warn) const
{
/// Return trigger for given detElemId and cath.
/// If segmentation is not of trigger type, zero is returned 
/// and an Error is issued if warn is set true (default). 

  return GetTrigger(GetMpSegmentation(detElemId, cath, warn), warn);
}    

//_____________________________________________________________________________
const AliMpTrigger*  
AliMpSegmentation::GetTriggerByElectronics(Int_t detElemId, Int_t elCardID, 
                              Bool_t warn) const
{
/// Return trigger for given detElemId and elCardID.
/// If segmentation is not of trigger type, zero is returned 
/// and an Error is issued if warn is set true (default). 

  return GetTrigger(GetMpSegmentationByElectronics(detElemId, elCardID, warn), warn);
}    
