// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jacek Otwinowski <Jacek.Otwinowski@gsi.de>            *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTGlobalTrackMergerComponent.cxx
    @author Jacek Otwinowski
    @date   
    @brief  HLT global track merger component.
*/

using namespace std;
#include <climits>
#include <cassert>
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCGlobalMerger.h"
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCVertexData.h"
#include "AliHLTTPCTrack.h"

//#include "AliHLTTPCSpacePointData.h"
//#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTRDDefinitions.h"
#include <cstdlib>
#include <cerrno>

#include "AliESDEvent.h"
#include "AliTracker.h"
#include "AliTRDtrackV1.h"
#include "AliHLTGlobalTrackMerger.h"

#include "AliHLTGlobalTrackMergerComponent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalTrackMergerComponent);

//_____________________________________________________________________________
AliHLTGlobalTrackMergerComponent::AliHLTGlobalTrackMergerComponent() : AliHLTProcessor(), 
  fGlobalTrackMerger(0), 
  fESD(0)
{
}

//_____________________________________________________________________________
AliHLTGlobalTrackMergerComponent::~AliHLTGlobalTrackMergerComponent()
{
  // see header file for class documentation
}

//_____________________________________________________________________________
const char* AliHLTGlobalTrackMergerComponent::GetComponentID()
{
  // see header file for class documentation
  return "GlobalTrackMerger";
}

//_____________________________________________________________________________
void AliHLTGlobalTrackMergerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( AliHLTTPCDefinitions::fgkTracksDataType );
  list.push_back( AliHLTTRDDefinitions::fgkTRDSATracksDataType );
}

//_____________________________________________________________________________
AliHLTComponentDataType AliHLTGlobalTrackMergerComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeESDObject|kAliHLTDataOriginTPC;
}

//_____________________________________________________________________________
void AliHLTGlobalTrackMergerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 20000;
  inputMultiplier = 1.0;
}

//_____________________________________________________________________________
AliHLTComponent* AliHLTGlobalTrackMergerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTGlobalTrackMergerComponent;
}

//_____________________________________________________________________________
void AliHLTGlobalTrackMergerComponent::SetMergerParameters(Double_t maxy,Double_t maxz,Double_t maxkappa,Double_t maxpsi,Double_t maxtgl)
{
  // see header file for class documentation
  fGlobalTrackMerger->SetParameter( maxy, maxz, maxkappa, maxpsi, maxtgl );
}

//_____________________________________________________________________________
int AliHLTGlobalTrackMergerComponent::DoInit( int /*argc*/, const char** /*argv*/ )
{
  // see header file for class documentation
  int iResult = 0;
  
  // Init merger
  fGlobalTrackMerger = new AliHLTGlobalTrackMerger();
  
  // output of the component
  fESD = new AliESDEvent();
  if (fESD) {
     fESD->CreateStdContent();
  }

  if (!fGlobalTrackMerger || !fESD ) {
     HLTError("failed creating internal objects");
     iResult=-ENOMEM;
     return iResult;
  }

  SetMergerParameters();

  return iResult;
}

//_____________________________________________________________________________
int AliHLTGlobalTrackMergerComponent::DoDeinit()
{
  // see header file for class documentation
  if(fGlobalTrackMerger) delete fGlobalTrackMerger; fGlobalTrackMerger =0;
  if(fESD) delete fESD; fESD = 0;
  return 0;
}

//_____________________________________________________________________________
int AliHLTGlobalTrackMergerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, 
					      AliHLTUInt32_t& /*size*/, AliHLTComponentBlockDataList& /*outputBlocks*/ )
{
  //
  // global track merger function
  // takes TRD and TPC tracks and merges them
  //
  HLTInfo("DoEvent processing data");

  // see header file for class documentation
  int iResult=0;

  if(!fGlobalTrackMerger || !fESD) {
    HLTError("component not initialized");
    iResult=-ENOMEM;
    return iResult;
  }

  if(!blocks) {
    HLTError("no blocks");
    iResult=-EINVAL;
    return iResult;
  }

  const AliHLTComponentBlockData* iter=0;
  AliHLTTPCTrackletData* inPtr=0;
  Bool_t bIsTRDTrackDataBlock=kFALSE;
  Bool_t bIsTPCTrackDataBlock=kFALSE;
  TClonesArray *aTRDTracks=0;
  Int_t minSlice = INT_MAX, maxSlice = 0;
  Int_t slice;

  unsigned long ndx;
  Int_t nTRDDataBlocks = 0;
  Int_t nTPCDataBlocks = 0;
  for ( ndx = 0; ndx < evtData.fBlockCnt && iResult>=0; ndx++ )
  {
      iter = blocks+ndx;
      bIsTRDTrackDataBlock=kFALSE;
      bIsTPCTrackDataBlock=kFALSE;

      // check if TPC or TRD tracks
      if(!(bIsTRDTrackDataBlock=(iter->fDataType==AliHLTTRDDefinitions::fgkTRDSATracksDataType)) &&
         !(bIsTPCTrackDataBlock=(iter->fDataType==AliHLTTPCDefinitions::fgkTracksDataType))) {
	continue;
      }

      // collect TRD tracks from all SM
      // one TClonesArray of tracks per SM
      if(bIsTRDTrackDataBlock) 
      {
        nTRDDataBlocks++;
	if(nTRDDataBlocks>1) continue;
        for (TObject *pObj = (TObject *)GetFirstInputObject(AliHLTTRDDefinitions::fgkTRDSATracksDataType,"TClonesArray",0);
	  pObj !=0 && iResult>=0;
	  pObj = (TObject *)GetNextInputObject(0)) {
          aTRDTracks = dynamic_cast<TClonesArray*>(pObj);
          if (!aTRDTracks) continue;

          HLTInfo("reading block %d, trdTracks %d", ndx, aTRDTracks->GetEntriesFast());

	  // load TRD tracks
	  if (fGlobalTrackMerger->LoadTracks(aTRDTracks,fESD) == kFALSE) {
             HLTError("Cannot load TRD tracks");
             iResult=-ENOMEM;
             return iResult;
	  }
	}
	aTRDTracks->Delete();
      } 
      
      // collect TPC tracks from whole TPC
      if (bIsTPCTrackDataBlock) 
      {
        nTPCDataBlocks++;
	if(nTPCDataBlocks>1) continue;

        slice=AliHLTTPCDefinitions::GetMinSliceNr(iter->fSpecification);
	if(slice<minSlice) minSlice = slice;
	if(slice>maxSlice) maxSlice = slice;
	
        //minslice=AliHLTTPCDefinitions::GetMinSliceNr(iter->fSpecification);
        //maxslice=AliHLTTPCDefinitions::GetMaxSliceNr(iter->fSpecification);

        AliHLTTPCTrackArray tracks;
        inPtr=(AliHLTTPCTrackletData*)iter->fPtr;

        HLTInfo("reading block %d (slice %d): %d tracklets", ndx, slice, inPtr->fTrackletCnt);

        // read TPC track segments from memory
        if((iResult=tracks.FillTracksChecked(inPtr->fTracklets, inPtr->fTrackletCnt, iter->fSize, -1/*global track*/, 0/*don't rotate*/))>=0) 
	{
          // load TPC tracks
          if (fGlobalTrackMerger->LoadTracks(&tracks,fESD) == kFALSE) {
             HLTError("Cannot load TPC tracks");
             iResult=-ENOMEM;
             return iResult;
	  }
        }
      }
   }

   // set magnetic field 
   fESD->SetMagneticField(AliTracker::GetBz());

   // merge tracks
   Bool_t isMerged = fGlobalTrackMerger->Merge(fESD);
   if(!isMerged) {
     HLTInfo("No merged tracks");
   }

   // try to propagate all tracks to DCA to primary vertex
   fGlobalTrackMerger->PropagateTracksToDCA(fESD);

   // calculate specification
   // AliHLTUInt32_t iSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( minSlice, maxSlice, 0, 5 );
   // HLTInfo("minSlice %d, maxSlice %d", minSlice, maxSlice);

   // send output data
   //PushBack(fESD, kAliHLTDataTypeESDObject|kAliHLTDataOriginTPC, iSpecification);
   PushBack(fESD, kAliHLTDataTypeESDObject|kAliHLTDataOriginTPC);

   // clean ESD event content
   fESD->Reset();
  
return iResult;
}

	
