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
/// \class AliMUONPadStatusMapMaker
/// 
/// Convert a pad status container into a pad status *map* container
/// 
/// A pad status is one 32-bits word describing whether this pad pedestal, gains
/// hv is correct or not.
///
/// A pad status *map* is one 32-bits (of which 24 only are used)
/// word describing whether this pad neighbours are ok or not
/// (whether a pad is ok or not is determined by applying a given
/// bitmask to the pad status word). Each bit in this word is related to one
/// neighbour, assuming the pad itself is at bit 0
///
/// ----------------
/// |  3 |  5 |  8 |
/// ----------------
/// |  2 |  0 |  7 |
/// ----------------
/// |  1 |  4 |  6 |
/// ----------------
///
/// Note that for instance in NonBending plane of slats, at the boundaries
/// between two pad densities, the pictures is a bit different, e.g.
/// (bits in () are always zero)
///
/// so some care must be taken when designing a mask to be tested ;-) if you
/// want to go farther than immediate neighbours...
///
/// If a pad is at a physical boundary, is will for sure have some bits at 1
/// (i.e. a non-existing neighbour is considered = bad).
///
// author Laurent Aphecetche

#include "AliMUONPadStatusMapMaker.h"

#include "AliLog.h"
#include "AliMUONCalibParam1I.h"
#include "AliMUONObjectPair.h"
#include "AliMUONV2DStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDataIterator.h"
#include "AliMpArea.h"
#include "AliMpConstants.h"
#include "AliMpDEManager.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpVPadIterator.h"
#include "AliMpVSegmentation.h"
#include "Riostream.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TStopwatch.h"

#include <map>
#include <utility>

ClassImp(AliMUONPadStatusMapMaker)

Int_t AliMUONPadStatusMapMaker::fgkSelfDead = 1;

namespace
{
  Bool_t IsZero(Double_t x)
  {
    return TMath::Abs(x) < AliMpConstants::LengthTolerance();
  }
}

//_____________________________________________________________________________
AliMUONPadStatusMapMaker::AliMUONPadStatusMapMaker() 
: TObject(),
fStatus(0x0),
fMask(0),
fSegmentation(0x0),
fTimerComputeStatusMap(0x0)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONPadStatusMapMaker::~AliMUONPadStatusMapMaker()
{
  /// dtor
  delete fTimerComputeStatusMap;
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMapMaker::ComputeStatusMap(const TObjArray& neighbours,
                                           Int_t detElemId) const
{
  /// Given a list of neighbours of one pad (which includes the pad itself)
  /// compute the status map (aka deadmap) for that pad.
  
  fTimerComputeStatusMap->Start(kFALSE);
  
  Int_t statusMap(0);
  
 //Compute the statusmap related to the status of neighbouring
 //pads. An invalid pad means "outside of edges".
  Int_t i(0);
  TIter next(&neighbours);
  AliMpPad* p;
  
  while ( ( p = static_cast<AliMpPad*>(next()) ) )
  {
    Int_t status = 0;
    if ( !p->IsValid() )
    {
      status = -1;
    }
    else
    {
      status = GetPadStatus(detElemId,*p);
    }
    if ( ( fMask==0 && status !=0 ) || ( (status & fMask) != 0 ) )
    {
      statusMap |= (1<<i);
    }
    ++i;
  }
  
  fTimerComputeStatusMap->Stop();
  return statusMap;
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMapMaker::GetPadStatus(Int_t detElemId,
                                      const AliMpPad& pad) const
{
  /// Get the pad status
  Int_t manuId = pad.GetLocation().GetFirst();
  Int_t manuChannel = pad.GetLocation().GetSecond();
  AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fStatus->Get(detElemId,manuId));
  return param->ValueAsInt(manuChannel);
}

//_____________________________________________________________________________
Bool_t
AliMUONPadStatusMapMaker::IsValid(const AliMpPad& pad, 
                                  const TVector2& shift) const
{
  /// Whether pad.Position()+shift is within the detector
  TVector2 testPos = pad.Position() - pad.Dimensions() + shift;
  AliMpPad p = fSegmentation->PadByPosition(testPos,kFALSE);
  return p.IsValid();
}

//_____________________________________________________________________________
AliMUONV2DStore*
AliMUONPadStatusMapMaker::MakePadStatusMap(const AliMUONV2DStore& status,
                                           Int_t mask)
{
  /// Given the status store for all pads, compute a status map store
  /// for all pads. 
  /// @param mask is the status mask to be tested to tell if a pad is ok or not
  
  fStatus = &status;
  fMask = mask;
  
  AliMpExMap chamberTimers(kTRUE);
  fTimerComputeStatusMap = new TStopwatch;
  fTimerComputeStatusMap->Start(kTRUE);
  fTimerComputeStatusMap->Stop();
  
  TStopwatch timer;
  
  timer.Start(kTRUE);
  
  AliMUONV2DStore* statusMap = status.CloneEmpty();
  
  AliMUONVDataIterator* it = status.Iterator();
  AliMUONObjectPair* pair;
  
  while ( ( pair = static_cast<AliMUONObjectPair*>(it->Next()) ) )
  {
    AliMpIntPair* ip = static_cast<AliMpIntPair*>(pair->First());

    Int_t detElemId = ip->GetFirst();
    
    Int_t manuId = ip->GetSecond();
    Int_t chamber = AliMpDEManager::GetChamberId(detElemId);
    
    TStopwatch* chTimer = static_cast<TStopwatch*>(chamberTimers.GetValue(chamber));
    if (!chTimer)
    {
      chTimer = new TStopwatch;
      chTimer->Start(kTRUE);
      chTimer->Stop();
      chamberTimers.Add(chamber,chTimer);
    }
    
    chTimer->Start(kFALSE);
    
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
    fSegmentation = seg;
    
    AliMUONVCalibParam* statusEntry = static_cast<AliMUONVCalibParam*>(pair->Second());
    
    AliMUONVCalibParam* statusMapEntry = static_cast<AliMUONVCalibParam*>
      (statusMap->Get(detElemId,manuId));

    if (!statusMapEntry)
    {
      statusMapEntry = new AliMUONCalibParam1I(64,0);
      statusMap->Set(detElemId,manuId,statusMapEntry,false);
    }
    
    for ( Int_t manuChannel = 0; manuChannel < statusEntry->Size(); ++manuChannel ) 
    {
      // Loop over channels and for each channel loop on its immediate neighbours
      // to produce a statusMap word for this channel.
      
      AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId,manuChannel),kFALSE);
      
      Int_t statusMapValue(0);
      
      if ( pad.IsValid() )
      {
        TObjArray neighbours;
        neighbours.SetOwner(kTRUE);
        Int_t n = fSegmentation->GetNeighbours(pad,neighbours,true,true);
        statusMapValue = ComputeStatusMap(neighbours,detElemId);      
      }
      else
      {
        statusMapValue = fgkSelfDead;
      }
      statusMapEntry->SetValueAsInt(manuChannel,0,statusMapValue);
    }
    chTimer->Stop();
  }
  
  delete it;
  
  TExMapIter cit = chamberTimers.GetIterator();
  
  Long_t key, value;
  
  while ( cit.Next(key,value) ) 
  {
    TStopwatch* t = reinterpret_cast<TStopwatch*>(value);
    cout << Form("Chamber %2ld CPU time/manu %5.0f ms ",key,t->CpuTime()*1e3/t->Counter());
    t->Print();
  }

  cout << "ComputeStatusMap timer : ";
  fTimerComputeStatusMap->Print();
  cout<< endl;
  
  cout << "MakePadStatusMap total timer : ";
  timer.Print();
  cout << endl;
  
  return statusMap;
}

