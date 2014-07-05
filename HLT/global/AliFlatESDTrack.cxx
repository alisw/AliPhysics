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

#include "TObject.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDfriendTrack.h"
#include "AliExternalTrackParam.h"
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"

#include "AliFlatESDEvent.h"
#include "AliFlatESDTrack.h"
#include "Riostream.h"

// _______________________________________________________________________________________________________
AliFlatESDTrack::AliFlatESDTrack() :
  // Default constructor
  fTrackParamMask(0),
  fNTPCClusters(0),
  fNITSClusters(0),

  fSize(0),
  fContent() {

}

// _______________________________________________________________________________________________________
AliFlatESDTrack::AliFlatESDTrack(AliFlatESDSpecialConstructorFlag f)
{
  //special contructor
  //use to restore the vtable pointer
  
  if(f == AliFlatESDReinitialize){
		AliFlatExternalTrackParam* trackParam = GetTrackParamRefitted();
		if (trackParam) { new (trackParam) AliFlatExternalTrackParam(f); }
		trackParam = GetTrackParamIp();
		if (trackParam) { new (trackParam) AliFlatExternalTrackParam(f); }
		trackParam = GetTrackParamTPCInner();
		if (trackParam) { new (trackParam) AliFlatExternalTrackParam(f); }
		trackParam = GetTrackParamOp();
		if (trackParam) { new (trackParam) AliFlatExternalTrackParam(f); }
		trackParam = GetTrackParamCp();
		if (trackParam) { new (trackParam) AliFlatExternalTrackParam(f); }
		trackParam = GetTrackParamITSOut();
		if (trackParam) { new (trackParam) AliFlatExternalTrackParam(f); }
		
		AliFlatTPCCluster* clusterTPC = GetTPCClusters();
		for (Int_t i=0; i<fNTPCClusters; i++){
		  new (clusterTPC) AliFlatTPCCluster(f);
		  clusterTPC++;
		 }
	}
	else AliFlatESDTrack();
}

// _______________________________________________________________________________________________________
AliFlatESDTrack::AliFlatESDTrack(const AliESDtrack* track, AliESDfriendTrack* friendTrack) :
  // Constructor
  fTrackParamMask(0),
  fNTPCClusters(0),
  fNITSClusters(0),
  fSize(0),
  fContent() {
  
  Fill(track, friendTrack);
}

// _______________________________________________________________________________________________________
AliFlatESDTrack::~AliFlatESDTrack() {
  // Destructor
  
}

// _______________________________________________________________________________________________________
ULong64_t AliFlatESDTrack::EstimateSize(Bool_t useESDFriends, Int_t nTPCClusters ) {
  // Estimate upper limit of the object size
  // -> Added objects have to be added here as well

  ULong64_t size = sizeof(AliFlatESDTrack) + (6*sizeof(AliFlatExternalTrackParam));

  if (useESDFriends){
    size += nTPCClusters*sizeof(AliFlatTPCCluster);
  }
  return size;
}


// _______________________________________________________________________________________________________
Int_t AliFlatESDTrack::Fill(const AliESDtrack* track, AliESDfriendTrack* friendTrack){
  // Fill external track parameters and friendTrack

  fTrackParamMask = 0;
  fNTPCClusters = 0;
  fNITSClusters = 0;
  fSize = 0;
  
  if( !track ) return 0;

  const AliExternalTrackParam *itsOut = 0;
  if (friendTrack) itsOut = friendTrack->GetITSOut();

  Int_t iResult = FillExternalTrackParam( track,
					  track->GetInnerParam(),
					  track->GetTPCInnerParam(),
					  track->GetOuterParam(),
					  track->GetConstrainedParam(),
					  itsOut );
  fNITSClusters = track->GetNcls(0);

  // -- Fill clusters from friend track
  // -------------------------------------------------------
  if (friendTrack) {
    //    Printf("DEBUG: Now filling clusters information for the current track");    

    // -- Get seed object
    TObject* calibObject = NULL;
    AliTPCseed* seed = NULL;
    /*
    for (Int_t idx = 0; (calibObject = friendTrack->GetCalibObject(idx)); ++idx) {      
      cout<<"Calibration object:"<<endl;
      std::cout<<calibObject->GetName()<<std::endl;
      calibObject->Print();
      cout<<"----------"<<endl;
    }
    friendTrack->Print();
    cout<<"ITS track: "<<(void*)friendTrack->GetITStrack()<<endl;
    cout<<"TRD track: "<<(void*)friendTrack->GetTRDtrack()<<endl;
    cout<<"ITS OUT track: "<<(void*)friendTrack->GetITSOut()<<endl;
    cout<<"ITS indices: ";
    if( friendTrack->GetITSindices() ){ 
      for( int i=0; i<friendTrack->GetMaxITScluster(); i++ ) cout<<friendTrack->GetITSindices()[i]<<" "; 
    }
    cout<<endl;
    */
    for (Int_t idx = 0; (calibObject = friendTrack->GetCalibObject(idx)); ++idx) {
      if ((seed = dynamic_cast<AliTPCseed*>(calibObject))) break;
    }

    // -- Fill cluster
    if (seed) {
      for (Int_t idxRow = 0; idxRow < 160; idxRow++){
	AliTPCclusterMI* currentCl = seed->GetClusterPointer(idxRow);
	if (currentCl) {
	  AliFlatTPCCluster* tmpCl = GetNextTPCClusterPointer();
	  new(tmpCl) AliFlatTPCCluster;
	  tmpCl->SetX(currentCl->GetX());
	  tmpCl->SetY(currentCl->GetY());
	  tmpCl->SetZ(currentCl->GetZ());	  
	 // tmpCl->SetPadRow(idxRow); // TO BE CHECKED IF THIS NEEDED or currentCl->GetRow();
	  tmpCl->SetPadRow(currentCl->GetRow());
	  tmpCl->SetSigmaY2(currentCl->GetSigmaY2());
	  tmpCl->SetSigmaZ2(currentCl->GetSigmaZ2());
	  tmpCl->SetCharge(currentCl->GetQ());
	  tmpCl->SetQMax(currentCl->GetMax());
	  StoreLastTPCCluster();
	}
	//	else
	//	  Printf("DEBUG: No cluster for row %d", idxRow);
      }
    }

    /*
    AliTPCseed* seed = NULL;
    for (Int_t idx = 0; (calibObject = friendTrack->GetCalibObject(idx)); ++idx) {
      if ((seed = dynamic_cast<AliTPCseed*>(calibObject))) break;
    }

    // -- Fill cluster
    if (seed) {
      for (Int_t idxRow = 0; idxRow < 160; idxRow++){
	AliTPCclusterMI* currentCl = seed->GetTPCClusterPointer(idxRow);
	if (currentCl) {
	  AliFlatTPCCluster &tmpCl = *GetNexTPCClusterPointer();
	  tmpCl.fX = currentCl->GetX();
	  tmpCl.fY = currentCl->GetY();
	  tmpCl.fZ = currentCl->GetZ();	  
	  tmpCl.fPadRow  = idxRow; // TO BE CHECKED IF THIS NEEDED or currentCl->GetRow();
	  tmpCl.fSigmaY2 = currentCl->GetSigmaY2();
	  tmpCl.fSigmaZ2 = currentCl->GetSigmaZ2();
	  tmpCl.fCharge  = currentCl->GetQ();
	  tmpCl.fQMax    = currentCl->GetMax();
	  StoreLastTPCCluster();
	}
	//	else
	//	  Printf("DEBUG: No cluster for row %d", idxRow);
      }
    }
    */

    //    else
    //      Printf("DEBUG: No seed object");

    //    Printf("DEBUG: Number of clusters for track = %d", fNTPCClusters);

    // -- Sorting clusters according to user defined function (increasing pad row numbering)
  //  std::sort(GetTPCClusters(), GetTPCClusters()+fNTPCClusters, AliFlatTPCCluster::SortClusters);
  }

  return iResult;
}

// _______________________________________________________________________________________________________
Int_t AliFlatESDTrack::FillExternalTrackParam( 
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
  fSize = 0;

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

  if (!param) 
    return -1;

  //Printf("  DEBUG: CONTENT %d >> %p + 0x%07llx = %p", flag, fContent, fSize, fContent + fSize);

  AliFlatExternalTrackParam * current = reinterpret_cast<AliFlatExternalTrackParam*> (fContent + fSize);
  new (current) AliFlatExternalTrackParam;
  current->SetAlpha(param->GetAlpha());
  current->SetX(param->GetX());
  current->SetY(param->GetY());
  current->SetZ(param->GetZ());
  current->SetSnp(param->GetSnp());
  current->SetTgl(param->GetTgl());
  current->SetSigned1Pt(param->GetSigned1Pt());
  
  const Double_t *cov = param->GetCovariance();
  for (Int_t idx = 0; idx <15; ++idx)
    current->fC[idx] = cov[idx];
    
  fTrackParamMask |= flag;
  fSize += sizeof(AliFlatExternalTrackParam);

  return 0;
}

// _______________________________________________________________________________________________________
UInt_t AliFlatESDTrack::CountBits(Byte_t field, UInt_t mask) {
  // Count bits in field
  UInt_t count = 0; 
  UInt_t reg = 0x0; 
  
  reg |= field;   
  reg &= mask;
  
  for (count = 0; reg; count++)
    reg &= reg - 1; 

  return count;
}
