/* $Id$ */

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

/**
 * >> Flat structure representing an ESDTrack <<
 *
 * To be used in the online and offline calibration schema.
 *
 * Class provides interface methods for 
 *   - Filling from AliESDtrack and AliExternalTrackParam, as well 
 *     as clusters from ESD friends (if requested)
 *   - HLT Filling to be added
 * 
 *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli
 *
 **************************************************************************/


#include "AliFlatESDFriendTrack.h"
#include "AliExternalTrackParam.h"
#include "AliESDfriendTrack.h"
#include "AliTPCseed.h"
#include "Riostream.h"

// _______________________________________________________________________________________________________
 AliFlatESDFriendTrack::AliFlatESDFriendTrack()
:
 AliVfriendTrack()
 ,fContentSize(0),
 fTPCOutPointer(-1),
 fITSOutPointer(-1),
 fTRDInPointer(-1),
 fTPCseedPointer(-1),
 fBitFlags(0)
{
  // Default constructor
  fContent[0]=0;
}

#pragma GCC diagnostic ignored "-Weffc++" 
AliFlatESDFriendTrack::AliFlatESDFriendTrack( AliVConstructorReinitialisationFlag f ) 
  :
  AliVfriendTrack( f )
{
  // constructor for reinitialisation of vtable
  if( fTPCseedPointer >= 0 ){
    AliFlatTPCseed *fp = reinterpret_cast< AliFlatTPCseed* >( fContent + fTPCseedPointer );
    fp->Reinitialize();
  }
}
#pragma GCC diagnostic warning "-Weffc++" 

void AliFlatESDFriendTrack::Reset()
{
  // reset
  fContentSize = 0;
  fTPCOutPointer = -1;
  fITSOutPointer = -1;
  fTRDInPointer = -1;
  fTPCseedPointer = -1;
  fBitFlags = 0; 
}

Int_t AliFlatESDFriendTrack::SetFromESDfriendTrack( const AliESDfriendTrack* track, size_t allocatedMemory )
{
  if( allocatedMemory < EstimateSize() ) return -1;
  Reset();
  if( !track ) return 0;
  SetSkipBit(track->TestSkipBit() );
  AddTrackParamTPCOut( track->GetTPCOut() );
  AddTrackParamITSOut( track->GetITSOut() );
  AddTrackParamTRDIn( track->GetTRDIn() );
  const AliTPCseed* seedP = NULL;
  {
    TObject* calibObject = NULL;
    for (Int_t idx = 0; (calibObject = track->GetCalibObject(idx)); ++idx) {
      if ((seedP = dynamic_cast<const AliTPCseed*>(calibObject))) {
	break;
      }
    }
  }
  AddTPCseed( seedP );
  return 0;
}
