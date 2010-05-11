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
// $MpId: AliMpDetElement.cxx,v 1.4 2006/05/24 13:58:34 ivana Exp $
// Category: management

//-----------------------------------------------------------------------------
// Class AliMpDetElement
// --------------------
// The class defines the electronics properties of detection element
// Authors: Ivana Hrivnacova, IPN Orsay
//          Laurent Aphecetche, Christian Finck, SUBATECH Nantes
//-----------------------------------------------------------------------------

#include "AliMpDetElement.h"

#include "AliMpArrayI.h"
#include "AliMpConstants.h"
#include "AliMpDEManager.h"
#include "AliMpDCSNamer.h"
#include "AliMpHVUID.h"
#include "AliMpManuUID.h"
#include "AliMpPadUID.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"

#include "AliCodeTimer.h"
#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpDetElement)
/// \endcond

const char  AliMpDetElement::fgkNameSeparator = '_'; 

//______________________________________________________________________________
AliMpDetElement::AliMpDetElement(Int_t id, const TString& name,
                    const TString& segType, AliMp::PlaneType planeType)
  : TObject(),
    fId(id),
    fDdlId(-1),
    fName(name),
    fSegType(segType),
    fPlaneType(planeType),
    fBusPatchIds(false),
    fManuList(),
    fTrackerChannels(),
    fHVmanus(),
    fNofChannels(0)
{
/// Standard constructor
}

//______________________________________________________________________________
AliMpDetElement::AliMpDetElement(TRootIOCtor* ioCtor)
  : TObject(),
    fId(0),
    fDdlId(-1),
    fName(),
    fSegType(),
    fPlaneType(),
    fBusPatchIds(),
    fManuList(),
    fTrackerChannels(),
    fHVmanus(ioCtor),
    fNofChannels()
{
/// Root IO constructor
}

//______________________________________________________________________________
AliMpDetElement::~AliMpDetElement()
{
/// Destructor
}


//
// public methods
//

//______________________________________________________________________________
Bool_t AliMpDetElement::AddBusPatch(Int_t busPatchId)
{
/// Add bus patch Id if a bus patch with the same Id is not yet present;
/// return false if bus patch was not added

  if ( HasBusPatchId(busPatchId) ) {
    AliWarningStream() 
      << "Bus patch Id = " << busPatchId << " already present."
      << endl;
    return false;
  } 

  fBusPatchIds.Add(busPatchId); 
  return true;
}  
 
//______________________________________________________________________________
TString AliMpDetElement::GetSegName(AliMp::CathodType cathType) const
{
/// Return the segmentation name for the given catod type

  return fSegType + fgkNameSeparator + PlaneTypeName(GetPlaneType(cathType));
}	  

//______________________________________________________________________________
AliMp::PlaneType  AliMpDetElement::GetPlaneType(AliMp::CathodType cath) const 
{
/// Return plane type                                                      \n

  if ( cath == AliMp::kCath0 ) return fPlaneType;
  else                         return AliMp::OtherPlaneType(fPlaneType); 
}    

//______________________________________________________________________________
AliMp::CathodType AliMpDetElement::GetCathodType(AliMp::PlaneType planeType) const
{
/// Return cathod type for given planeType

  if ( fPlaneType == planeType ) return AliMp::kCath0;
  else                           return AliMp::kCath1;
}

//______________________________________________________________________________
AliMp::CathodType AliMpDetElement::GetCathodTypeFromManuId(Int_t manuId) const
{
/// Return cathod type for given manuId

  AliMp::PlaneType planeType = AliMp::kBendingPlane;
  if ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) ) 
  {
    planeType = AliMp::kNonBendingPlane;
  }
  return GetCathodType(planeType);
}

//______________________________________________________________________________
AliMp::StationType AliMpDetElement::GetStationType() const
{
/// Return station type                                                      \n
/// Failure causes Fatal error - as AliMp::StationType has no possibility
/// to return undefined value

  Int_t chamberId = AliMpDEManager::GetChamberId(fId, false);
  if ( ! AliMpDEManager::IsValidChamberId(chamberId, true) ) {
    AliFatal("Cannot return AliMp::StationType value.");
    return AliMp::kStation12;
  }  
  
  if ( chamberId >=  0 && chamberId <=  3 )  return AliMp::kStation12;
  if ( chamberId >=  4 && chamberId <=  9 )  return AliMp::kStation345;
  if ( chamberId >= 10 && chamberId <= 13 )  return AliMp::kStationTrigger;

  // Should never get to this line
  AliFatal("Cannot return AliMp::StationType value.");
  return AliMp::kStation12;
}

//______________________________________________________________________________
AliMq::Station12Type AliMpDetElement::GetStation12Type() const
{
/// Return station12 type                                                      \n
/// Failure causes Fatal error - as AliMp::Station12Type has no possibility
/// to return undefined value

  Int_t chamberId = AliMpDEManager::GetChamberId(fId, false);
  if ( ! AliMpDEManager::IsValidChamberId(chamberId, true) ) {
    AliFatal("Cannot return AliMp::StationType value.");
    return AliMq::kNotSt12;
  }  
  
  if ( chamberId ==  0 || chamberId ==  1 )  return AliMq::kStation1;
  if ( chamberId ==  2 || chamberId ==  3 )  return AliMq::kStation2;
  if ( chamberId >=  4 || chamberId <= 13 )  return AliMq::kNotSt12;

  // Should never get to this line
  AliFatal("Cannot return AliMp::StationType value.");
  return AliMq::kNotSt12;
}

//______________________________________________________________________________
Int_t AliMpDetElement::GetNofBusPatches() const
{
/// Return the number of bus patches in this detection element

  return fBusPatchIds.GetSize();
}  

//______________________________________________________________________________
Int_t AliMpDetElement::GetBusPatchId(Int_t index) const
{
/// Return the index-th bus patch

  if ( index < 0 || index >= GetNofBusPatches() ) {
    AliErrorStream()
      << "In DE = " << fId << ": Index " << index << " outside limits." << endl;
    return 0;
  }     

  return  fBusPatchIds.GetValue(index);
}   


//______________________________________________________________________________
Bool_t  AliMpDetElement::HasBusPatchId(Int_t busPatchId) const
{  
/// Return true if the bus patch Id is present

  return fBusPatchIds.HasValue(busPatchId);; 
}

//______________________________________________________________________________
Int_t 
AliMpDetElement::NofChannelsInManu(Int_t manuId) const
{
  /// Return the number of channels in a given manu
  
  Long_t uid = AliMpManuUID::BuildUniqueID(fId,manuId);
  
  return (Int_t)(fManuList.GetValue(uid));
}

//______________________________________________________________________________
Bool_t 
AliMpDetElement::IsExistingChannel(Int_t manuId, Int_t manuChannel) const
{
  /// Whether or not the channel is a valid one (does not tell if it is
  /// connected or not
  
  if ( NofChannelsInManu(manuId) > 0 && 
       manuChannel >= 0 && 
       manuChannel < AliMpConstants::ManuNofChannels() ) 
  {
    return kTRUE;
  }
  else
  {
    return kFALSE;
  }
}

//______________________________________________________________________________
Bool_t 
AliMpDetElement::IsConnectedChannel(Int_t manuId, Int_t manuChannel) const
{
  /// Whether or not the channel is a *connected* one (i.e. it is valid plus
  /// it corresponds to a real pad)
  
  return ( fTrackerChannels.GetValue(AliMpPadUID::BuildUniqueID(fId,manuId,manuChannel)) > 0 );
}

//______________________________________________________________________________
void
AliMpDetElement::AddManu(Int_t manuId)
{
  /// Fills the fManuList and fTrackerChannels
  AliMp::StationType stationType = AliMpDEManager::GetStationType(fId);
  
  if ( stationType == AliMp::kStationTrigger ) return;
    
  AliCodeTimerAuto("",0)
  
  AliDebug(1,Form("DE %4d Manu %4d",fId,manuId));

  AliCodeTimerStart(Form("%s",AliMp::StationTypeName(stationType).Data()));
  
  if ( fHVmanus.GetSize() == 0 ) 
  {
    fHVmanus.SetOwner(kTRUE); // to be 100% explicit
    
    // get the size, to avoid resizing when adding later on
    Int_t nmanus(0);
    
    AliMp::CathodType cathodes[] = { AliMp::kCath0, AliMp::kCath1 };
    
    for ( Int_t i = 0; i < 2; ++i ) 
    {
      const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentation(fId,cathodes[i]);
      
      TArrayI manus;
      
      seg->GetAllElectronicCardIDs(manus);

      nmanus += manus.GetSize();
    }
    
    fHVmanus.SetSize(nmanus);
  }
  
  const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(fId,manuId);
  
  Int_t n(0);
  
  for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i ) 
  {
    if ( seg->PadByLocation(manuId,i,kFALSE).IsValid() )
    {
      ++n;
      fTrackerChannels.Add((Long_t)AliMpPadUID::BuildUniqueID(fId,manuId,i),
                           (Long_t)1);
    }
  }
  
  fManuList.Add(AliMpManuUID::BuildUniqueID(fId,manuId),(Long_t)n);

  fNofChannels += n;
  
  AliMpDCSNamer hvNamer("TRACKER");
  
  Int_t index = hvNamer.ManuId2Index(fId,manuId);
                            
  UInt_t hvuid = AliMpHVUID::BuildUniqueID(fId,index);
  
  AliMpArrayI* hv = static_cast<AliMpArrayI*>(fHVmanus.GetValue(hvuid));
  
  if (!hv)
  {
    Bool_t sort(kFALSE);
    hv = new AliMpArrayI(sort);
    fHVmanus.Add(hvuid,hv);
  }
  
  hv->Add(manuId,kFALSE);        
  
  AliCodeTimerStop(Form("%s",AliMp::StationTypeName(stationType).Data()));
}
  
//______________________________________________________________________________
const AliMpArrayI* 
AliMpDetElement::ManusForHV(Int_t hvIndex) const
{
  /// Return the list of manus sharing a hv channel
  return static_cast<AliMpArrayI*>(fHVmanus.GetValue(AliMpHVUID::BuildUniqueID(fId,hvIndex)));
}

//______________________________________________________________________________
Int_t 
AliMpDetElement::NofManus() const
{
  /// Return the number of manus in this detection element
  return fManuList.GetSize();
}

