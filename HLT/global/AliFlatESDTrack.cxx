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

#include "Rtypes.h"
#include "AliFlatESDTrack.h"
#include "AliFlatExternalTrackParam.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "Riostream.h"




// _______________________________________________________________________________________________________


// _______________________________________________________________________________________________________
Int_t AliFlatESDTrack::SetFromESDTrack(const AliESDtrack* track)
{
  // Fill external track parameters 
  fTrackParamMask = 0;
  fNTPCClusters = 0;
  fNITSClusters = 0;
  fContentSize = 0;
  
  if( !track ) return 0;

  Int_t iResult = SetExternalTrackParam( track,
					 track->GetInnerParam(),
					 track->GetTPCInnerParam(),
					 track->GetOuterParam(),
					 track->GetConstrainedParam(), NULL );
  fNITSClusters = track->GetTPCNcls();

  return iResult;
}

// _______________________________________________________________________________________________________
Int_t AliFlatESDTrack::SetExternalTrackParam( 
					     const AliExternalTrackParam* refittedParam,
					     const AliExternalTrackParam* innerParam,
					     const AliExternalTrackParam* innerTPC,
					     const AliExternalTrackParam* outerParam,
					     const AliExternalTrackParam* constrainedParam,
					     const AliExternalTrackParam* outerITS
					      ){
  // Fill external track parameters 

  fTrackParamMask = 0;
  fNTPCClusters = 0;
  fContentSize = 0;

  Int_t iResult = 0;

  Byte_t flag = 0x1;
  iResult = FillExternalTrackParam(refittedParam, flag);

  flag = 0x2;
  iResult = FillExternalTrackParam(innerParam, flag);
  
  flag = 0x4;
  iResult = FillExternalTrackParam(innerTPC, flag);
  
  flag = 0x8;
  iResult = FillExternalTrackParam(outerParam, flag);

  flag = 0x10;
  iResult = FillExternalTrackParam(constrainedParam, flag);

  flag = 0x20;
  iResult = FillExternalTrackParam(outerITS, flag);

  return iResult;
}

// _______________________________________________________________________________________________________
Int_t AliFlatESDTrack::FillExternalTrackParam(const AliExternalTrackParam* param, UShort_t flag) {
  // Fill external track parameters

  if (!param) return -1;

  //Printf("  DEBUG: CONTENT %d >> %p + 0x%07llx = %p", flag, fContent, fContentSize, fContent + fContentSize);

  AliFlatExternalTrackParam * current = reinterpret_cast<AliFlatExternalTrackParam*> (fContent + fContentSize);
  current->SetExternalTrackParam( param );    
  fTrackParamMask |= flag;
  fContentSize += sizeof(AliFlatExternalTrackParam);

  return 0;
}


// _______________________________________________________________________________________________________
Bool_t AliFlatESDTrack::GetXYZ(Double_t *p) const {
  const AliFlatExternalTrackParam *f = GetFlatTrackParam();
  if (!f) { return kFALSE; }
  p[0]=f->GetX();
  p[1]=f->GetY();
  p[2]=f->GetZ();
  return kTRUE;
}
