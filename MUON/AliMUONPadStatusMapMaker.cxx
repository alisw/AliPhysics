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

//-----------------------------------------------------------------------------
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
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUONPadStatusMapMaker.h"

#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONVStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMpConstants.h"
#include "AliMpIntPair.h"
#include "AliMpManuList.h"
#include <Riostream.h>
#include <TList.h>
#include <TStopwatch.h>

/// \cond CLASSIMP
ClassImp(AliMUONPadStatusMapMaker)
/// \endcond

Int_t AliMUONPadStatusMapMaker::fgkSelfDead = 1;

//_____________________________________________________________________________
AliMUONPadStatusMapMaker::AliMUONPadStatusMapMaker(const AliMUONCalibrationData& calibData) 
: TObject(),
fStatus(0x0),
fMask(0),
fCalibrationData(calibData)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONPadStatusMapMaker::~AliMUONPadStatusMapMaker()
{
  /// dtor
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMapMaker::ComputeStatusMap(const AliMUONVCalibParam& neighbours,
                                          Int_t manuChannel,
                                          Int_t detElemId) const
{
  /// Given a list of neighbours of one pad (which includes the pad itself)
  /// compute the status map (aka deadmap) for that pad.
  
  Int_t statusMap(0);

  //Compute the statusmap related to the status of neighbouring
  //pads. An invalid pad means "outside of edges".

  Int_t n = neighbours.Dimension();
  for ( Int_t i = 0; i < n; ++i )
  {
    Int_t x = neighbours.ValueAsInt(manuChannel,i);
    Int_t m,c;
    neighbours.UnpackValue(x,m,c);
    if ( c < 0 ) continue;
    Int_t status = 0;
    if ( !m )
    {
      status = -1;
    }
    else
    {
      status = GetPadStatus(detElemId,m,c);
    }
    if ( ( fMask==0 && status !=0 ) || ( (status & fMask) != 0 ) )
    {
      statusMap |= (1<<i);
    }
  }
  return statusMap;
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMapMaker::GetPadStatus(Int_t detElemId, 
                                       Int_t manuId, Int_t manuChannel) const
                                      
{
  /// Get the pad status
  AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fStatus->FindObject(detElemId,manuId));
  return param->ValueAsInt(manuChannel);
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONPadStatusMapMaker::MakeEmptyPadStatusMap()
{
  /// Make an empty (but complete) statusMap
    
    AliMUONVStore* store = new AliMUON2DMap(true);
    
    TList* list = AliMpManuList::ManuList();
    
    AliMpIntPair* pair;
    
    TIter next(list);
    
    while ( ( pair = static_cast<AliMpIntPair*>(next()) ) ) 
    {
      Int_t detElemId = pair->GetFirst();
      Int_t manuId = pair->GetSecond();
      AliMUONVCalibParam* param = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),
                                                          detElemId,manuId,
                                                          0);
      store->Add(param);
    }
    
    delete list;
    
    return store;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONPadStatusMapMaker::MakePadStatusMap(const AliMUONVStore& status,
                                           Int_t mask)
{
  /// Given the status store for all pads, compute a status map store
  /// for all pads. 
  /// \param status
  /// \param mask is the status mask to be tested to tell if a pad is ok or not
  
  fStatus = &status;
  fMask = mask;
  
  TStopwatch timer;  
  timer.Start(kTRUE);
  
  AliMUONVStore* neighbourStore = fCalibrationData.Neighbours();
  
  AliMUONVStore* statusMap = status.Create();
  
  TIter next(status.CreateIterator());
  AliMUONVCalibParam* statusEntry;
  
  while ( ( statusEntry = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    Int_t detElemId = statusEntry->ID0();
    Int_t manuId = statusEntry->ID1();
        
    AliMUONVCalibParam* statusMapEntry = static_cast<AliMUONVCalibParam*>
      (statusMap->FindObject(detElemId,manuId));

    if (!statusMapEntry)
    {
      statusMapEntry = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),
                                               detElemId,manuId,0);
      statusMap->Add(statusMapEntry);
    }
    
    AliMUONVCalibParam* neighbours = static_cast<AliMUONVCalibParam*>
      (neighbourStore->FindObject(detElemId,manuId));
    
    if (!neighbours)
    {
      AliFatal(Form("Could not find neighbours for DE %d manuId %d",
                    detElemId,manuId));
      continue;
    }
    
    for ( Int_t manuChannel = 0; manuChannel < statusEntry->Size(); ++manuChannel ) 
    {
      // Loop over channels and for each channel loop on its immediate neighbours
      // to produce a statusMap word for this channel.
      
      Int_t statusMapValue(0);

      Int_t x = neighbours->ValueAsInt(manuChannel,0);
      
      if ( x > 0 )
      { 
        // channel is a valid one (i.e. (manuId,manuChannel) is an existing pad)
        statusMapValue = ComputeStatusMap(*neighbours,manuChannel,detElemId);
      }
      else
      {
        statusMapValue = fgkSelfDead;
      }
      
      statusMapEntry->SetValueAsInt(manuChannel,0,statusMapValue);
    }
  }
  timer.Stop();
  
  StdoutToAliInfo(
                  cout << "MakePadStatusMap total timer : ";
                  timer.Print();
                  cout << endl;
                  );

  return statusMap;
}

