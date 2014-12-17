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
#include "AliTPCseed.h"
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
					 track->GetConstrainedParam() );
  fNTPCClusters = track->GetTPCNcls();
  fNITSClusters = track->GetITSNcls();

  track->GetImpactParameters(fImp, fImp+2 );
  fImp[5] = track->GetConstrainedChi2();

  track->GetImpactParametersTPC(fImpTPC, fImpTPC+2 );
  fImpTPC[5] = track->GetConstrainedChi2TPC();

  return iResult;
}

void  AliFlatESDTrack::GetESDTrack( AliESDtrack* esdTrack ) const
{
  // get esd track out of flat track

  if( !esdTrack ) return;

  AliTPCseed p;
  p.SetNumberOfClusters( GetNumberOfTPCClusters() );

  if( GetTrackParamOp( p )>=0 ){
    esdTrack->UpdateTrackParams( &p, AliESDtrack::kTPCout );
  }
  if( GetTrackParamTPCInner( p )>=0 ){
    esdTrack->UpdateTrackParams( &p, AliESDtrack::kTPCin );
  }
  if( GetTrackParamRefitted( p )>=0 ){
    esdTrack->UpdateTrackParams( &p, AliESDtrack::kTPCrefit );
  }
  if( GetTrackParamIp( p )>=0 ){
    p.SetNumberOfClusters( GetNumberOfITSClusters() );
   esdTrack->UpdateTrackParams( &p, AliESDtrack::kITSin );
  }

 if( GetTrackParamCp( p )>=0 ){
   esdTrack->SetImpactParameters( fImp, fImp+2, fImp[5], &p);
 } else {
   esdTrack->SetImpactParameters( fImp, fImp+2, fImp[5], NULL);
 }
 esdTrack->SetImpactParametersTPC( fImpTPC, fImpTPC+2, fImpTPC[5] );

}



// _______________________________________________________________________________________________________
Int_t AliFlatESDTrack::SetExternalTrackParam( 
					     const AliExternalTrackParam* refittedParam,
					     const AliExternalTrackParam* innerParam,
					     const AliExternalTrackParam* innerTPC,
					     const AliExternalTrackParam* outerParam,
					     const AliExternalTrackParam* constrainedParam
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
  //return the global track position
  const AliFlatExternalTrackParam *f = GetFlatTrackParam();
  if (!f) { return kFALSE; }

  p[0]=f->GetX();
  p[1]=f->GetY();
  p[2]=f->GetZ();
  return Local2GlobalPosition(p,f->GetAlpha());
}
