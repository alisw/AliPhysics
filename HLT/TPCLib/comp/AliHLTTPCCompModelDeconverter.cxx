// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
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

/** @file   AliHLTTPCCompModelDeconverter.cxx
    @author Timm Steinbeck
    @date   
    @brief  A copy processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCCompModelDeconverter.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCModelTrack.h"
#include "AliHLTTPCCompDataCompressorHelper.h"
#include <cerrno>

AliHLTTPCCompModelDeconverter::AliHLTTPCCompModelDeconverter():
    fInputTrackArray("AliHLTTPCModelTrack"),
    fTrackClusterModelData(0),
    fTrackClusterModelDataSize(0),
    fRemainingClustersModelData(0),
    fRemainingClustersModelDataSize(0)
    {
      // see header file for class documentation
      Init();
    }

AliHLTTPCCompModelDeconverter::~AliHLTTPCCompModelDeconverter()
    {
      // see header file for class documentation
    }

int AliHLTTPCCompModelDeconverter::Init()
    {
      // see header file for class documentation
      fInputTrackArray.Reset();
      fTrackClusterModelData = 0;
      fTrackClusterModelDataSize = 0;
      fRemainingClustersModelData = 0;
      fRemainingClustersModelDataSize = 0;
      return 0;
    }


int AliHLTTPCCompModelDeconverter::SetTrackClusterModelInputData( AliHLTUInt8_t* data, UInt_t size )
    {
      // see header file for class documentation
      fTrackClusterModelData = data;
      fTrackClusterModelDataSize = size;
      
      AliHLTUInt8_t* inputPtr = fTrackClusterModelData;
      AliHLTUInt8_t* inputEndPtr = fTrackClusterModelData+fTrackClusterModelDataSize;
      
      AliHLTUInt32_t version = *(AliHLTUInt32_t*)inputPtr;
      inputPtr += sizeof(AliHLTUInt32_t);
      if ( version != 0 )
	{
	  HLTError( "Unsupported version %hu. Only version 0 supported currently.", version );
	  return EIO;
	}
      
      while ( inputPtr+sizeof(AliHLTTPCTrackModel)+AliHLTTPCTransform::GetNRows()*sizeof(AliHLTTPCClusterModel) <= inputEndPtr )
	{
	  AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)fInputTrackArray.NextTrack();
	  if ( !track )
	    {
	      HLTError( "Error obtaining model track from track array." );
	      return EIO;
	    }
	  Int_t slice = ((AliHLTTPCClusterModel*)( inputPtr+sizeof(AliHLTTPCTrackModel) ))[ AliHLTTPCTransform::GetNRows()-1 ].fSlice;
	  track->Init(slice,-1);
	  AliHLTTPCTrackModel *model = track->GetModel();
	  AliHLTTPCClusterModel *clusters = track->GetClusters();
	  memcpy( model, inputPtr, sizeof(AliHLTTPCTrackModel) );
	  memcpy( clusters, inputPtr+sizeof(AliHLTTPCTrackModel), AliHLTTPCTransform::GetNRows()*sizeof(AliHLTTPCClusterModel) );
	  track->FillTrack();
	  
	  // validation test
	  //HLTInfo("track->GetNClusters() %d  GetNPresentClusters() %d", track->GetNClusters(), track->GetNPresentClusters());
	  
	  inputPtr += sizeof(AliHLTTPCTrackModel)+AliHLTTPCTransform::GetNRows()*sizeof(AliHLTTPCClusterModel);
	}
      
      if ( inputPtr!=inputEndPtr )
	{
	  HLTError( "Data format inconsistency on reading track model data." );
	  return EIO;
	}
      
      return 0;
    }

int AliHLTTPCCompModelDeconverter::SetRemainingClustersModelInputData( AliHLTUInt8_t* data, UInt_t size )
    {
      // see header file for class documentation
      fRemainingClustersModelData = data;
      fRemainingClustersModelDataSize = size;
      
      AliHLTUInt32_t version = *(AliHLTUInt32_t*)data;
      if ( version != 0 )
	{
	  HLTError( "Unsupported version %hu. Only version 0 supported currently.", version );
	  return EIO;
	}
      return 0;
    }

int AliHLTTPCCompModelDeconverter::DeconvertTracks( AliHLTUInt8_t* data, UInt_t& size )
    {
      // see header file for class documentation
      AliHLTTPCTrackletData* outPtr = (AliHLTTPCTrackletData*)data;
      UInt_t blockSize = fInputTrackArray.WriteTracks( outPtr->fTrackletCnt, outPtr->fTracklets );
      if ( blockSize >= size )
	{
	  HLTError( "Output size is bigger than allowed (%u instead of %u)",
		    (unsigned)blockSize, (unsigned)size );
	  size = 0;
	  return ENOBUFS;
	}
      size = blockSize;
      
      
      UInt_t clusterCounts[36][6];
      memset( clusterCounts, 0, sizeof(UInt_t)*36*6 );
      
      AliHLTTPCTrackSegmentData* tracklet = outPtr->fTracklets;
      
      // validation test
      //HLTInfo("fTrackletCnt = %d", outPtr->fTrackletCnt);
      //HLTInfo("tracklet Points = %d", tracklet->fNPoints);
      
      for ( UInt_t ii = 0; ii<outPtr->fTrackletCnt; ii++ )
	{
	  // validation test
	  //HLTInfo("Tracklet points %d", tracklet->fNPoints);
	  for ( UInt_t jj=0; jj<tracklet->fNPoints; jj++ )
	    {
	      //validation test
	      //HLTInfo("Hello World %d !",jj);
	      UInt_t slice, partition;
	      slice = AliHLTTPCSpacePointData::GetSlice(tracklet->fPointIDs[jj]);
	      partition = AliHLTTPCSpacePointData::GetPatch(tracklet->fPointIDs[jj]);
	      tracklet->fPointIDs[jj] = AliHLTTPCSpacePointData::GetID(slice,partition,clusterCounts[slice][partition]);
	      clusterCounts[slice][partition]++;
	    }
	  tracklet = (AliHLTTPCTrackSegmentData*) ( ((AliHLTUInt8_t*)tracklet)+sizeof(AliHLTTPCTrackSegmentData)+tracklet->fNPoints*sizeof(UInt_t) );
	}
      
      return 0;
    }

int AliHLTTPCCompModelDeconverter::DeconvertClusters( UInt_t slice, UInt_t patch, AliHLTUInt8_t* data, UInt_t& size )
    {
      // see header file for class documentation
      if ( size<sizeof(AliHLTTPCClusterData) )
	return ENOBUFS;
      Int_t charge;
      Float_t pad,time,sigmaY2,sigmaZ2;
      AliHLTTPCClusterData* clusterData = (AliHLTTPCClusterData*)data;
      clusterData->fSpacePointCnt = 0;
      AliHLTTPCSpacePointData* clusters = clusterData->fSpacePoints;
      unsigned long outSize = sizeof(AliHLTTPCClusterData);
      for(Int_t i=0; i<fInputTrackArray.GetNTracks(); i++)
	{
	  AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)fInputTrackArray.GetCheckedTrack(i);
	  if(!track)
	    continue;
	  for(Int_t padrow=AliHLTTPCTransform::GetFirstRow(patch); padrow <= AliHLTTPCTransform::GetLastRow(patch); padrow++)
	    {
	      if(!track->IsPresent(padrow))
		continue;
	      UInt_t thisSlice = track->GetClusterModel(padrow)->fSlice;
	      if ( thisSlice != slice )
		continue;
	      if ( clusterData->fSpacePointCnt >= (1<<22) )
		{
		  HLTError( "Too many clusters for slice %d patch %d", slice, patch );
		  break;
		}
	      if ( size<=outSize+sizeof(AliHLTTPCSpacePointData) )
		{
		  HLTError( "Not enough output space (%u bytes)", (unsigned)size );
		  return ENOBUFS;
		}
	      track->GetPad(padrow,pad);
	      track->GetTime(padrow,time);
	      track->GetClusterCharge(padrow,charge);
	      //NEW: Get parameters to create parSigma correctly
	      //track->GetCrossingAngleLUT(padrow);
	      //track->CalculateClusterWidths(padrow,kTRUE); // calculates parSigmas (with parametrisation) in raw coordinates
	      //HLTInfo("dangle %f", track->GetCrossingAngleLUT(padrow));
	      //HLTInfo("dparsigma %f",track->GetParSigmaY2(padrow));
	      track->GetSigmaY2(padrow,sigmaY2);
	      //AliHLTTPCClusterModel* test1 = track->GetClusterModel(padrow);
	      //HLTInfo("DsigmaY deconv. : %f",test1->fDSigmaY);
	      //HLTInfo("sigmaY2 deconv.: %f",sigmaY2);
	      track->GetSigmaZ2(padrow,sigmaZ2);
	      Float_t xyz[3];
	      AliHLTTPCTransform::RawHLT2Local( xyz, slice, padrow, pad, time );
	      clusters[clusterData->fSpacePointCnt].fX = xyz[0];
	      clusters[clusterData->fSpacePointCnt].fY = xyz[1];
	      clusters[clusterData->fSpacePointCnt].fZ = xyz[2];
	      clusters[clusterData->fSpacePointCnt].fPadRow = padrow;
	      clusters[clusterData->fSpacePointCnt].fSigmaY2 = sigmaY2*pow(AliHLTTPCTransform::GetPadPitchWidth(patch),2);
	      clusters[clusterData->fSpacePointCnt].fSigmaZ2 = sigmaZ2*pow(AliHLTTPCTransform::GetZWidth(),2);;
	      clusters[clusterData->fSpacePointCnt].fCharge = charge;
	      clusters[clusterData->fSpacePointCnt].SetUsed(kTRUE);
	      clusters[clusterData->fSpacePointCnt].SetTrackNumber(i);
	      clusters[clusterData->fSpacePointCnt].SetID(slice,patch,clusterData->fSpacePointCnt);

	      clusterData->fSpacePointCnt++;
	      outSize += sizeof(AliHLTTPCSpacePointData);
	    }
	}
      
    if ( fRemainingClustersModelDataSize )
      {
	AliHLTUInt8_t* inputPtr = fRemainingClustersModelData+sizeof(AliHLTUInt32_t);
	AliHLTUInt8_t* inputEndPtr = inputPtr+fRemainingClustersModelDataSize;
	for ( UInt_t thisSlice=0; thisSlice<36; thisSlice++ )
	  {
	    for ( UInt_t thisPatch=0; thisPatch<6; thisPatch++ )
	      {
		AliHLTUInt8_t rowCount = *inputPtr;
		inputPtr++;
		if ( !rowCount )
		  {
		    if ( thisSlice==slice && thisPatch==patch )
		      break;
		    continue;
		  }
		for ( UInt_t jj=0; jj < rowCount; jj++ )
		    {
		      AliHLTTPCRemainingRow *thisRow = (AliHLTTPCRemainingRow*)inputPtr;
		      if ( inputPtr+sizeof(AliHLTTPCRemainingRow)>inputEndPtr )
			{
			  HLTError( "Corrupt input data, cannot read row data for row %u of slice %u, partition %u", (unsigned)jj, (unsigned)thisSlice, (unsigned)thisPatch );
			  return EIO;
			}
		      AliHLTTPCRemainingCluster *cl = thisRow->fClusters;
		      if ( inputPtr+sizeof(AliHLTTPCRemainingRow)+thisRow->fNClusters*sizeof(AliHLTTPCRemainingCluster)>inputEndPtr )
			{
			  HLTError( "Corrupt input data, unable to read clusters for row %u, slice %u, partition %u", (unsigned)jj, (unsigned)thisSlice, (unsigned)thisPatch );
			  return EIO;
			}
		      Int_t padrow = thisRow->fPadRow;
		      if ( slice==thisSlice && patch==thisPatch )
			{
			  for ( UInt_t ii=0; ii<thisRow->fNClusters; ii++ )
			    {
			      Float_t xyz[3];
			      AliHLTTPCTransform::RawHLT2Local( xyz, slice, padrow, cl[ii].fPad, cl[ii].fTime );
			      clusters[clusterData->fSpacePointCnt].fX = xyz[0];
			      clusters[clusterData->fSpacePointCnt].fY = xyz[1];
			      clusters[clusterData->fSpacePointCnt].fZ = xyz[2];
			      clusters[clusterData->fSpacePointCnt].fPadRow = padrow;
			      clusters[clusterData->fSpacePointCnt].fSigmaY2 = cl[ii].fSigmaY2*pow(AliHLTTPCTransform::GetPadPitchWidth(patch),2);
			      clusters[clusterData->fSpacePointCnt].fSigmaZ2 = cl[ii].fSigmaZ2*pow(AliHLTTPCTransform::GetZWidth(),2);;
			      clusters[clusterData->fSpacePointCnt].fCharge = cl[ii].fCharge;
			      clusters[clusterData->fSpacePointCnt].SetUsed(kFALSE);
			      clusters[clusterData->fSpacePointCnt].SetTrackNumber(-1);
			      clusters[clusterData->fSpacePointCnt].SetID(slice,patch,clusterData->fSpacePointCnt);
			      clusterData->fSpacePointCnt++;
			      outSize += sizeof(AliHLTTPCSpacePointData);
			    }
			}
		      inputPtr += sizeof(AliHLTTPCRemainingRow)+thisRow->fNClusters*sizeof(AliHLTTPCRemainingCluster);
		    }
		
	      }
	    if ( thisSlice==slice )
	      break;
	  }
      }
    
    size = outSize;
    return 0;
    }
