// @(#) $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCEsdWriterComponent.cxx
    @author Matthias Richter
    @date   
    @brief  Writer component to store tracks of the HLT TPC conformal
            mapping tracker in the AliESD format

                                                                          */
#include "AliHLTTPCEsdWriterComponent.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "TTree.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCDefinitions.h"

/** global instance for component registration */
AliHLTTPCEsdWriterComponent gTPCEsdWriter;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCEsdWriterComponent)

AliHLTTPCEsdWriterComponent::AliHLTTPCEsdWriterComponent()
  :
  fTree(NULL),
  fESD(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCEsdWriterComponent::~AliHLTTPCEsdWriterComponent()
{
  // see header file for class documentation
}

int AliHLTTPCEsdWriterComponent::InitWriter()
{
  // see header file for class documentation
  int iResult=0;
  fESD = new AliESDEvent;
  if (fESD) {
    fESD->CreateStdContent();
    fTree = new TTree("esdTree", "Tree with HLT ESD objects");
    if (fTree) {
      fESD->WriteToTree(fTree);
    }
  }
  if (fTree==NULL) {
    iResult=-ENOMEM;
  }
  return iResult;
}

int AliHLTTPCEsdWriterComponent::CloseWriter()
{
  // see header file for class documentation
  int iResult=0;
  if (fTree) {
    WriteObject(kAliHLTVoidEventID, fTree);
    TTree* pTree=fTree;
    fTree=NULL;
    delete pTree;
  } else {
    HLTWarning("not initialized");
  }
  iResult=AliHLTRootFileWriterComponent::CloseWriter();
  return iResult;
}

int AliHLTTPCEsdWriterComponent::DumpEvent( const AliHLTComponentEventData& evtData,
					    const AliHLTComponentBlockData* blocks, 
					    AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;
  TTree* pTree=fTree;
  if (pTree) {
    if (fESD) {
      AliESDEvent* pESD=fESD;

      const AliHLTComponentBlockData* iter = NULL;
      AliHLTTPCTrackletData* inPtr=NULL;
      int bIsTrackSegs=0;
 
      for (int ndx=0; ndx<(int)evtData.fBlockCnt && iResult>=0; ndx++) {
	iter = blocks+ndx;
	if ( (bIsTrackSegs=(iter->fDataType == AliHLTTPCDefinitions::fgkTrackSegmentsDataType))==1 ||
	     iter->fDataType == AliHLTTPCDefinitions::fgkTracksDataType ) {
	  Int_t minslice=AliHLTTPCDefinitions::GetMinSliceNr(iter->fSpecification);
	  Int_t maxslice=AliHLTTPCDefinitions::GetMaxSliceNr(iter->fSpecification);
	  if (bIsTrackSegs==0) {
	    // slice parameter and data specification ignored, tracks already in global coordinates
	    minslice=-1;
	    maxslice=-1;
	  }
	  //HLTDebug("dataspec %#x minslice %d", iter->fSpecification, minslice);
	  if (minslice >=-1 && minslice<36) {
	    if (minslice!=maxslice) {
	      HLTWarning("data from multiple sectors in one block: "
			 "possible missmatch in treatment of local coordinate system");
	    }
	    AliHLTTPCTrackArray tracks;
	    inPtr=(AliHLTTPCTrackletData*)iter->fPtr;
	    HLTDebug("reading block %d (slice %d): %d tracklets", ndx, minslice, inPtr->fTrackletCnt);
	    tracks.FillTracks(inPtr->fTrackletCnt, inPtr->fTracklets, minslice, 0/*don't rotate*/);
	    if ((iResult=Tracks2ESD(&tracks, pESD))>=0) {
	    }
	  } else {
	    HLTError("invalid sector number");
	    iResult=-EBADF;
	  }
	}
      }
      if (iResult>=0) {
	pTree->Fill();
      }

      fESD->Reset();
    } else {
      iResult=-ENOMEM;
    }
  }
  return iResult;
}

int AliHLTTPCEsdWriterComponent::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=AliHLTRootFileWriterComponent::ScanArgument(argc, argv);
  return iResult;
}

int AliHLTTPCEsdWriterComponent::Tracks2ESD(AliHLTTPCTrackArray* pTracks, AliESDEvent* pESD)
{
  // see header file for class documentation
  int iResult=0;
  if (pTracks && pESD) {
    HLTDebug("converting %d tracks from track array", pTracks->GetNTracks());
    for (int i=0; i<pTracks->GetNTracks() && iResult>=0; i++) {
      AliHLTTPCTrack* pTrack=(*pTracks)[i];
      if (pTrack) {
	//HLTDebug("convert track %d", i);
	//pTrack->Print();
	int iLocal=pTrack->Convert2AliKalmanTrack();
	if (iLocal>=0) {
	AliESDtrack iotrack;
	iotrack.UpdateTrackParams(pTrack,AliESDtrack::kTPCin);
	iotrack.SetTPCPoints(pTrack->GetPoints());
	pESD->AddTrack(&iotrack);
	} else {
	  HLTError("conversion to AliKalmanTrack failed for track %d of %d", i, pTracks->GetNTracks());
	}
      } else {
	HLTError("internal missmatch in array");
	iResult=-EFAULT;
      }
    }
    
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}
