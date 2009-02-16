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
/// Convert a pad statuses into pad status maps.
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
#include "AliMUONPadStatusMaker.h"
#include "AliMUONVStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMpConstants.h"
#include <Riostream.h>
#include <TList.h>
#include "AliCodeTimer.h"

/// \cond CLASSIMP
ClassImp(AliMUONPadStatusMapMaker)
/// \endcond

Int_t AliMUONPadStatusMapMaker::fgkSelfDead = 1;

//_____________________________________________________________________________
AliMUONPadStatusMapMaker::AliMUONPadStatusMapMaker(const AliMUONPadStatusMaker& padStatusMaker,
                                                   Int_t mask,
                                                   Bool_t deferredInitialization) 
: TObject(),
fkStatusMaker(padStatusMaker),
fMask(mask),
fStatusMap(new AliMUON2DMap(true))
{
  /// ctor
  if (!deferredInitialization)
  {
    AliCodeTimerAuto("Computing complete status map at once");
    AliMUONVStore* neighboursStore = padStatusMaker.NeighboursStore();
    AliMUONVCalibParam* param;
    TIter next(neighboursStore->CreateIterator());
    while ( ( param = static_cast<AliMUONVCalibParam*>(next()) ) )
    {
      Int_t detElemId = param->ID0();
      Int_t manuId = param->ID1();
      ComputeStatusMap(detElemId,manuId);
    }
  }
}

//_____________________________________________________________________________
AliMUONPadStatusMapMaker::~AliMUONPadStatusMapMaker()
{
  /// dtor
  delete fStatusMap;
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONPadStatusMapMaker::ComputeStatusMap(Int_t detElemId, Int_t manuId) const
{
  /// Compute the status map for a given manu, and add it to our internal
  /// fStatusMap internal storage
  
  AliCodeTimerAuto("(Int_t,Int_t)")
    
  AliMUONVCalibParam* param = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),
                                                      detElemId,manuId,-1);    
                                    
  Bool_t ok = fStatusMap->Add(param);
  if (!ok)
  {
    AliFatal(Form("Could not add manu %d of de %d",manuId,detElemId));
  }
                                  
  AliMUONVCalibParam* neighbours = fkStatusMaker.Neighbours(detElemId,manuId);
  
  AliMUONVCalibParam* statusParam = fkStatusMaker.PadStatus(detElemId,manuId);
  
  Int_t n = neighbours->Dimension();
  
  for ( Int_t manuChannel = 0; manuChannel < param->Size(); ++manuChannel )
  {
    Int_t statusMap(0);
    
    Int_t x = neighbours->ValueAsIntFast(manuChannel,0);
    if ( x < 0 ) 
    {
      // channel is not a valid one (i.e. (manuId,manuChannel) is not an existing pad)
      statusMap = -1;//fgkSelfDead;
      continue;
    }
        
    for ( Int_t i = 0; i < n; ++i )
    {
      // Compute the statusmap related to the status of neighbouring
      // pads. An invalid pad means "outside of edges".
            
      Int_t y = neighbours->ValueAsIntFast(manuChannel,i);      
      Int_t m,c;
      neighbours->UnpackValue(y,m,c);
      if ( c < 0 ) continue;
      Int_t status = 0;
      if ( !m )
      {
        status = -1;
      }
      else
      {
        status = statusParam->ValueAsIntFast(c); //fkStatusMaker.PadStatus(detElemId,m,c);
      }
      if ( ( fMask==0 && status !=0 ) || ( (status & fMask) != 0 ) )
      {
        statusMap |= (1<<i);
      }
    }    
    param->SetValueAsIntFast(manuChannel,0,statusMap);
  }
  return param;
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMapMaker::StatusMap(Int_t detElemId, Int_t manuId, 
                                    Int_t manuChannel) const
                                      
{
  /// Get the pad status map
  
  AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fStatusMap->FindObject(detElemId,manuId));
  if (!param)
  {
    // not yet computed, so do it now
    param = ComputeStatusMap(detElemId,manuId);
  }
  return param->ValueAsInt(manuChannel);
}
