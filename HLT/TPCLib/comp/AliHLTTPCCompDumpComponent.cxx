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

/** @file   AliHLTTPCCompDumpComponent.cxx
    @author Timm Steinbeck
    @date   10-08-2006
    @brief  A copy processing component for the HLT
            that writes the results of the Vestbo compression
            components to humanly readable files 
*/

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCCompDumpComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCModels.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCCompDataCompressorHelper.h"
#include <stdlib.h>
#include <errno.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCCompDumpComponent)
    
AliHLTTPCCompDumpComponent::AliHLTTPCCompDumpComponent()
  :
  fBitDataCurrentWord(0),
  fBitDataCurrentPosInWord(0),
  fBitDataCurrentInput(NULL),
  fBitDataCurrentInputStart(NULL),
  fBitDataCurrentInputEnd(NULL)
    {
      // see header file for class documentation
    }

AliHLTTPCCompDumpComponent::~AliHLTTPCCompDumpComponent()
    {
      // see header file for class documentation
    }

const char* AliHLTTPCCompDumpComponent::GetComponentID() 
    {
      // see header file for class documentation
    return "TPCCompDump"; // The ID of this component
    }

void AliHLTTPCCompDumpComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list) 
    {
      // see header file for class documentation
      list.clear(); // We do not have any requirements for our input data type(s).
      list.push_back( AliHLTTPCDefinitions::fgkClustersDataType );
      list.push_back( AliHLTTPCDefinitions::fgkTrackSegmentsDataType );
      list.push_back( AliHLTTPCDefinitions::fgkTracksDataType );
      list.push_back( AliHLTTPCDefinitions::fgkClusterTracksModelDataType );
      list.push_back( AliHLTTPCDefinitions::fgkRemainingClustersModelDataType );
    }

AliHLTComponent_DataType AliHLTTPCCompDumpComponent::GetOutputDataType()
    {
      // see header file for class documentation
      return AliHLTTPCDefinitions::fgkClusterTracksModelDataType;
    }

void AliHLTTPCCompDumpComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) 
    {
      // see header file for class documentation
      constBase = 0;
      inputMultiplier = 1.;
      //#warning Adapt input Multiplier to something more realistic
    }


// Spawn function, return new instance of this class
AliHLTComponent* AliHLTTPCCompDumpComponent::Spawn()
    {
      // see header file for class documentation
      return new AliHLTTPCCompDumpComponent;
    }

void AliHLTTPCCompDumpComponent::InitBitDataInput( AliHLTUInt8_t* input, UInt_t inputSize )
    {
      // see header file for class documentation
      fBitDataCurrentWord = 0;
      fBitDataCurrentPosInWord = 7;
      fBitDataCurrentInput = fBitDataCurrentInputStart = input;
      fBitDataCurrentInputEnd = input+inputSize;
      fBitDataCurrentWord = *fBitDataCurrentInput;
    }

 bool AliHLTTPCCompDumpComponent::InputBit( AliHLTUInt8_t & value )
    {
      // see header file for class documentation
      if ( fBitDataCurrentInput>=fBitDataCurrentInputEnd )
	return false;
      value = (fBitDataCurrentWord >> fBitDataCurrentPosInWord) & 1;
      if ( fBitDataCurrentPosInWord )
	fBitDataCurrentPosInWord--;
      else
	{
	  fBitDataCurrentInput++;
	  if ( fBitDataCurrentInput<fBitDataCurrentInputEnd )
	    {
	      fBitDataCurrentWord = *fBitDataCurrentInput;
	      fBitDataCurrentPosInWord = 7;
	    }
	}
      return true;
    }

bool AliHLTTPCCompDumpComponent::InputBits( AliHLTUInt8_t & value, UInt_t const & bitCount )
    {
      // see header file for class documentation
      if ( bitCount>8 )
	{
	  HLTFatal( "Internal error: Attempt to write more than 32 bits (%u)", (unsigned)bitCount );
	  return false;
	}
      AliHLTUInt64_t temp;
      if ( !InputBits( temp, bitCount ) )
	return false;
      value = (AliHLTUInt8_t)( temp & (AliHLTUInt64_t)0xFFFFFFFFULL );
      return true;
    }

bool AliHLTTPCCompDumpComponent::InputBits( AliHLTUInt16_t & value, UInt_t const & bitCount )
   {
     // see header file for class documentation
     if ( bitCount>16 )
       {
	 HLTFatal( "Internal error: Attempt to write more than 32 bits (%u)", (unsigned)bitCount );
	 return false;
       }
     AliHLTUInt64_t temp;
     if ( !InputBits( temp, bitCount ) )
       return false;
     value = (AliHLTUInt16_t)( temp & (AliHLTUInt64_t)0xFFFFFFFFULL );
     return true;
   }

bool AliHLTTPCCompDumpComponent::InputBits( AliHLTUInt32_t & value, UInt_t const & bitCount )
   {
     // see header file for class documentation
     if ( bitCount>32 ) 	 
       { 	 
	 HLTFatal( "Internal error: Attempt to write more than 32 bits (%u)", (unsigned)bitCount ); 	 
	 return false; 	 
       } 	 
     AliHLTUInt64_t temp; 	 
     if ( !InputBits( temp, bitCount ) ) 	 
       return false; 	 
     value = (AliHLTUInt32_t)( temp & (AliHLTUInt64_t)0xFFFFFFFFULL ); 	 
     return true; 	 
   }

bool AliHLTTPCCompDumpComponent::InputBits( Int_t & value, UInt_t const & bitCount )
   {
     // see header file for class documentation
     if ( bitCount>32 )
       {
	 HLTFatal( "Internal error: Attempt to write more than 32 bits (%u)", (unsigned)bitCount );
	 return false;
       }
     AliHLTUInt64_t temp;
     if ( !InputBits( temp, bitCount ) )
       return false;
     value = (Int_t)( temp & (AliHLTUInt64_t)0xFFFFFFFFULL );
     return true;
   }

bool AliHLTTPCCompDumpComponent::InputBits( AliHLTUInt64_t & value, UInt_t const & bitCount )
   {
     // see header file for class documentation
     if ( bitCount>64 )
       {
	 HLTFatal( "Internal error: Attempt to write more than 64 bits (%u)", (unsigned)bitCount );
	 return false;
       }
     UInt_t bitsToRead=bitCount;
     UInt_t curBitCount;
     value = 0;
     while ( bitsToRead>0 )
       {
	 if ( fBitDataCurrentInput>=fBitDataCurrentInputEnd )
	   return false;
	 if ( bitsToRead >= fBitDataCurrentPosInWord+1 )
	   curBitCount = fBitDataCurrentPosInWord+1;
	 else
	   curBitCount = bitsToRead;
	 value = (value << curBitCount) | ( (fBitDataCurrentWord >> (fBitDataCurrentPosInWord-curBitCount+1)) & ((1 << curBitCount)-1) );
	 if ( fBitDataCurrentPosInWord < curBitCount )
	   {
	     fBitDataCurrentInput++;
	     if ( fBitDataCurrentInput<fBitDataCurrentInputEnd )
	       {
		 fBitDataCurrentWord = *fBitDataCurrentInput;
		 fBitDataCurrentPosInWord = 7;
	       }
	   }
	 else
	   fBitDataCurrentPosInWord -= curBitCount;
	 bitsToRead -= curBitCount;
       }
     return true;
   }

void AliHLTTPCCompDumpComponent::Pad8Bits()
   {
     // see header file for class documentation
     if ( fBitDataCurrentPosInWord == 7 )
       return;
     fBitDataCurrentInput++;
     if ( fBitDataCurrentInput<fBitDataCurrentInputEnd )
       {
	 fBitDataCurrentWord = *fBitDataCurrentInput;
	 fBitDataCurrentPosInWord = 7;
       }
   }

bool AliHLTTPCCompDumpComponent::InputBytes( AliHLTUInt8_t* data, UInt_t const & byteCount )
   {
     // see header file for class documentation
     Pad8Bits();
     if ( fBitDataCurrentInput+byteCount>fBitDataCurrentInputEnd )
       return false;
     memcpy( data, fBitDataCurrentInput, byteCount );
     fBitDataCurrentInput += byteCount;
     if ( fBitDataCurrentInput<fBitDataCurrentInputEnd )
       {
	 fBitDataCurrentWord = *fBitDataCurrentInput;
	 fBitDataCurrentPosInWord = 7;
       }
     return true;
   }

int AliHLTTPCCompDumpComponent::DoInit( int argc, const char** argv )
    {
      // see header file for class documentation
      //char* cpErr;
      if ( argc )
	{
	  Logging( kHLTLogDebug, "HLT::TPCCompDump::DoInit", "Arguments", "argv[0] == %s", argv[0] );
	  Logging(kHLTLogError, "HLT::TPCCompDump::DoInit", "Unknown Option", "Unknown option '%s'", argv[0] );
	  return EINVAL;
	}
      return 0;
    }

int AliHLTTPCCompDumpComponent::DoDeinit()
    {
      // see header file for class documentation
      return 0;
    }

int AliHLTTPCCompDumpComponent::DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
						   AliHLTComponent_TriggerData& /*trigData*/, AliHLTUInt8_t* , 
				      AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>&  )
    {
      // see header file for class documentation
      // Process an event
      // Loop over all input blocks in the event
      for ( unsigned long n = 0; n < evtData.fBlockCnt; n++ )
	{
	  AliHLTUInt8_t slice, patch;
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkClustersDataType )
	    {
	      slice = AliHLTTPCDefinitions::GetMinSliceNr( blocks[n].fSpecification );
	      patch = AliHLTTPCDefinitions::GetMinPatchNr( blocks[n].fSpecification );
	      AliHLTTPCClusterData* clusters = (AliHLTTPCClusterData*)blocks[n].fPtr;
	      HLTInfo( "Cluster block slice %u - patch %u - %lu clusters", (unsigned)slice, (unsigned)patch, (unsigned long)clusters->fSpacePointCnt );
	      for ( unsigned long ii=0; ii<clusters->fSpacePointCnt; ii++ )
		{
		  HLTInfo( "  Cluster % 5lu: fX: %f - fY: %f - fZ: %f - fZ - fID: %u (0x%08X) - fPadRow: %u - fSigmaY2: %f - fSigmaZ2: %f - fCharge: %u - fUsed: %s - fTrackN: %d", 
			   ii, clusters->fSpacePoints[ii].fX, clusters->fSpacePoints[ii].fY, clusters->fSpacePoints[ii].fZ, clusters->fSpacePoints[ii].fID, clusters->fSpacePoints[ii].fID, (unsigned)clusters->fSpacePoints[ii].fPadRow, clusters->fSpacePoints[ii].fSigmaY2, clusters->fSpacePoints[ii].fSigmaZ2, clusters->fSpacePoints[ii].fCharge, (clusters->fSpacePoints[ii].IsUsed() ? "yes" : "no"), clusters->fSpacePoints[ii].GetTrackNumber() );
#if 0
		  Float_t xyzG[3] = { clusters->fSpacePoints[ii].fX, clusters->fSpacePoints[ii].fY, clusters->fSpacePoints[ii].fZ };
		  Float_t xyzR[3] = { clusters->fSpacePoints[ii].fX, clusters->fSpacePoints[ii].fY, clusters->fSpacePoints[ii].fZ };
		  AliHLTTPCTransform::LocHLT2Global( xyzG, slice, clusters->fSpacePoints[ii].fPadRow );
		  AliHLTTPCTransform::LocHLT2Raw( xyzG, slice, clusters->fSpacePoints[ii].fPadRow );
		  HLTInfo( "         Global: fX: %f - fY: %f - fZ: %f - fZ", 
			   xyzG[0], xyzG[1], xyzG[2] );
		  HLTInfo( "            Raw: fX: %f - fY: %f - fZ: %f - fZ", 
			   xyzR[0], xyzR[1], xyzR[2] );
#endif
		}
	      HLTInfo( "" );
	    }
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkTrackSegmentsDataType ||
	       blocks[n].fDataType == AliHLTTPCDefinitions::fgkTracksDataType)
	    {
	      //fConverter.SetInputTracks( (AliHLTTPCTrackletData*)blocks[n].fPtr );
	      AliHLTUInt8_t minSlice=0xFF, maxSlice=0xFF, minPatch=0xFF, maxPatch=0xFF;
	      minSlice = AliHLTTPCDefinitions::GetMinSliceNr( blocks[n].fSpecification );
	      maxSlice = AliHLTTPCDefinitions::GetMaxSliceNr( blocks[n].fSpecification );
	      minPatch = AliHLTTPCDefinitions::GetMinPatchNr( blocks[n].fSpecification );
	      maxPatch = AliHLTTPCDefinitions::GetMaxPatchNr( blocks[n].fSpecification );
	      AliHLTTPCTrackletData* tracks = (AliHLTTPCTrackletData*)blocks[n].fPtr;
	      AliHLTTPCTrackSegmentData* tracklet = tracks->fTracklets;
	      HLTInfo( "Track block slices %u-%u - patches %u-%u - %lu tracks", 
		       (unsigned)minSlice, (unsigned)maxSlice, (unsigned)minPatch, (unsigned)maxPatch,
		       (unsigned long)tracks->fTrackletCnt );
	      for ( unsigned long ii=0; ii<tracks->fTrackletCnt; ii++ )
		{
		  HLTInfo( "  Track % 5lu: fX: %f - fY: %f - fZ: %f - fLastX: %f - fLastY: %f - fLastZ: %f - fPt: %f - fPsi: %f - fTgl: %f - fPterr: %f - fPsierr: %f - fTglerr: %f - fCharge: %d - fNPoints: %u", 
			   ii, tracklet->fX, tracklet->fY, tracklet->fZ, tracklet->fLastX, tracklet->fLastY, tracklet->fLastZ, tracklet->fPt, tracklet->fPsi, tracklet->fTgl, tracklet->fPterr, tracklet->fPsierr, tracklet->fTglerr, tracklet->fCharge, tracklet->fNPoints );
		  for ( unsigned long jj=0; jj<tracklet->fNPoints; jj++ )
		    {
		      HLTInfo( "    Point % 5lu: %   8u / 0x%08X", jj, tracklet->fPointIDs[jj], tracklet->fPointIDs[jj] );
		    }
		  tracklet = (AliHLTTPCTrackSegmentData*) ( ((AliHLTUInt8_t*)tracklet)+sizeof(AliHLTTPCTrackSegmentData)+tracklet->fNPoints*sizeof(UInt_t) );
		}
	      
	    }
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkClusterTracksModelDataType )
	    {
	      AliHLTUInt8_t minSlice=0xFF, maxSlice=0xFF, minPatch=0xFF, maxPatch=0xFF;
	      minSlice = AliHLTTPCDefinitions::GetMinSliceNr( blocks[n].fSpecification );
	      maxSlice = AliHLTTPCDefinitions::GetMaxSliceNr( blocks[n].fSpecification );
	      minPatch = AliHLTTPCDefinitions::GetMinPatchNr( blocks[n].fSpecification );
	      maxPatch = AliHLTTPCDefinitions::GetMaxPatchNr( blocks[n].fSpecification );
	      unsigned long trackletCount = (blocks[n].fSize-sizeof(AliHLTUInt32_t)) / (sizeof(AliHLTTPCTrackModel)+AliHLTTPCTransform::GetNRows()*sizeof(AliHLTTPCClusterModel) );
	      HLTInfo( "Track model block version %u slices %u-%u - patches %u-%u - %lu tracks", 
		       (unsigned)*(AliHLTUInt32_t*)blocks[n].fPtr,
		       (unsigned)minSlice, (unsigned)maxSlice, (unsigned)minPatch, (unsigned)maxPatch,
		       (unsigned long)trackletCount );
	      AliHLTTPCTrackModel* trackModel = (AliHLTTPCTrackModel*)(((AliHLTUInt8_t*)blocks[n].fPtr)+sizeof(AliHLTUInt32_t));
	      for ( unsigned long ii=0; ii<trackletCount; ii++ )
		{
		  unsigned clusterCount=0;
		  AliHLTTPCClusterModel* clusters = (AliHLTTPCClusterModel*) ( ((AliHLTUInt8_t*)trackModel)+sizeof(AliHLTTPCTrackModel) );
		  for ( unsigned long jj=0; jj<(unsigned long)AliHLTTPCTransform::GetNRows(); jj++ )
		    {
		      if ( clusters[jj].fPresent )
			clusterCount++;
		    }
		  HLTInfo( "  Track Model % 5lu fKappa: %f - fPhi: %f - fD: %f - fZ0: %f - fTgl: %f - #clusters: %u",
			  ii, trackModel->fKappa, trackModel->fPhi, trackModel->fD, trackModel->fZ0, trackModel->fTgl, clusterCount );
		  clusterCount=0;
		  for ( unsigned long jj=0; jj<(unsigned long)AliHLTTPCTransform::GetNRows(); jj++ )
		    {
		      if ( clusters[jj].fPresent )
			{
#ifdef MODELDEBUG
			  HLTInfo( "    Cluster % 05u: fID: %u (0x%08X) - fDTime: %f - fDPad: %f - fDCharge: %f - fDSigmaY: %f - fDSigmaZ: %f - fNPads: %u - fSlice: %hd - padrow: %lu - fPresent: %u",
				   clusterCount, clusters[jj].fID, clusters[jj].fID, clusters[jj].fDTime, clusters[jj].fDPad, clusters[jj].fDCharge, clusters[jj].fDSigmaY, clusters[jj].fDSigmaZ, clusters[jj].fNPads, clusters[jj].fSlice, jj, (unsigned)clusters[jj].fPresent );
#else
			  HLTInfo( "    Cluster % 05u: fDTime: %f - fDPad: %f - fDCharge: %f - fDSigmaY: %f - fDSigmaZ: %f - fNPads: %u - fSlice: %hd - padrow: %lu - fPresent: %u",
				   clusterCount, clusters[jj].fDTime, clusters[jj].fDPad, clusters[jj].fDCharge, clusters[jj].fDSigmaY, clusters[jj].fDSigmaZ, clusters[jj].fNPads, clusters[jj].fSlice, jj, (unsigned)clusters[jj].fPresent );
#endif
			  clusterCount++;
			}
		    }
		  
		  trackModel = (AliHLTTPCTrackModel*) ( ((AliHLTUInt8_t*)trackModel)+sizeof(AliHLTTPCTrackModel)+AliHLTTPCTransform::GetNRows()*sizeof(AliHLTTPCClusterModel) );
		}
	      
	    }
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkRemainingClustersModelDataType )
	    {
	      AliHLTUInt8_t minSlice=0xFF, maxSlice=0xFF, minPatch=0xFF, maxPatch=0xFF;
	      minSlice = AliHLTTPCDefinitions::GetMinSliceNr( blocks[n].fSpecification );
	      maxSlice = AliHLTTPCDefinitions::GetMaxSliceNr( blocks[n].fSpecification );
	      minPatch = AliHLTTPCDefinitions::GetMinPatchNr( blocks[n].fSpecification );
	      maxPatch = AliHLTTPCDefinitions::GetMaxPatchNr( blocks[n].fSpecification );
	      unsigned long clusterCount=0;
	      HLTInfo( "Remaining cluster model block version %u", 
		       (unsigned)*(AliHLTUInt32_t*)blocks[n].fPtr );
	      AliHLTUInt8_t* readPtr = (AliHLTUInt8_t*)(((AliHLTUInt8_t*)blocks[n].fPtr)+sizeof(AliHLTUInt32_t));
	      for(Int_t lslice=0; lslice<36; lslice++)
		{
		  for(Int_t lpatch=0; lpatch < 6; lpatch++)
		    {
		      if ( !readPtr )
			{
			  readPtr++;
			  continue;
			}
		      unsigned rows = (unsigned)*readPtr;
		      readPtr++;
		      for ( unsigned ii=0; ii<rows; ii++ )
			{
			  AliHLTTPCRemainingRow* thisRow = (AliHLTTPCRemainingRow*)readPtr;
			  clusterCount += thisRow->fNClusters;
			  readPtr += sizeof(AliHLTTPCRemainingRow) + thisRow->fNClusters*sizeof(AliHLTTPCRemainingCluster);
			}
		    }
		}
	      HLTInfo( "Remaining cluster model block slices %u-%u - patches %u-%u - %lu clusters", 
		       (unsigned)minSlice, (unsigned)maxSlice, (unsigned)minPatch, (unsigned)maxPatch,
		       clusterCount );
	      readPtr = (AliHLTUInt8_t*)(((AliHLTUInt8_t*)blocks[n].fPtr)+sizeof(AliHLTUInt32_t));
	      clusterCount = 0;
	      for(Int_t lslice=0; lslice<36; lslice++)
		{
		  for(Int_t lpatch=0; lpatch < 6; lpatch++)
		    {
		      if ( !readPtr )
			{
			  readPtr++;
			  continue;
			}
		      unsigned rows = (unsigned)*readPtr;
		      readPtr++;
		      if ( rows )
			HLTInfo( "  Slice %d - Partition %d", lslice, lpatch );
		      for ( unsigned ii=0; ii<rows; ii++ )
			{
			  AliHLTTPCRemainingRow* thisRow = (AliHLTTPCRemainingRow*)readPtr;
			  for ( unsigned jj=0; jj<thisRow->fNClusters; jj++ )
			    {
#ifdef MODELDEBUG
			      HLTInfo( "    Cluster % 5lu: fID: %u (0x%08X) - fPadRow: %u - fPad: %f - fTime: %f - fSigmaY2: %f - fSigmaZ2: %f - fCharge: %hu",
				       clusterCount, thisRow->fClusters[jj].fID, thisRow->fClusters[jj].fID, (unsigned)thisRow->fPadRow, thisRow->fClusters[jj].fPad, thisRow->fClusters[jj].fTime, thisRow->fClusters[jj].fSigmaY2, thisRow->fClusters[jj].fSigmaZ2, thisRow->fClusters[jj].fCharge );
#else
			      HLTInfo( "    Cluster % 5lu: fPadRow: %u - fPad: %f - fTime: %f - fSigmaY2: %f - fSigmaZ2: %f - fCharge: %hu",
				       clusterCount, (unsigned)thisRow->fPadRow, thisRow->fClusters[jj].fPad, thisRow->fClusters[jj].fTime, thisRow->fClusters[jj].fSigmaY2, thisRow->fClusters[jj].fSigmaZ2, thisRow->fClusters[jj].fCharge );
#endif
			    }
			  readPtr += sizeof(AliHLTTPCRemainingRow) + thisRow->fNClusters*sizeof(AliHLTTPCRemainingCluster);
			}
		    }
		}
	    }
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkClusterTracksCompressedDataType )
	    {
	      AliHLTUInt8_t minSlice=0xFF, maxSlice=0xFF, minPatch=0xFF, maxPatch=0xFF;
	      minSlice = AliHLTTPCDefinitions::GetMinSliceNr( blocks[n].fSpecification );
	      maxSlice = AliHLTTPCDefinitions::GetMaxSliceNr( blocks[n].fSpecification );
	      minPatch = AliHLTTPCDefinitions::GetMinPatchNr( blocks[n].fSpecification );
	      maxPatch = AliHLTTPCDefinitions::GetMaxPatchNr( blocks[n].fSpecification );
	      HLTInfo( "Track model block slices %u-%u - patches %u-%u", 
		      (unsigned)minSlice, (unsigned)maxSlice, (unsigned)minPatch, (unsigned)maxPatch );
	      InitBitDataInput( (AliHLTUInt8_t*)blocks[n].fPtr, blocks[n].fSize );
	      HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
	      AliHLTUInt8_t version;
	      if ( !InputBits( version, 4 ) ) // Version information
		{
		  HLTError( "Corrupt input data. Cannot read data version number at position %lu / %u",
			    GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
		  continue;
		}
	      HLTInfo( "Data Format Version: %u", (unsigned)version );
	      HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
	      AliHLTUInt8_t readShape;
	      if ( !InputBit( readShape ) ) // Data format flag
		{
		  HLTError( "Corrupt input data. Cannot read shape flag at position %lu / %u",
			  GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
		  continue;
		}
	      HLTInfo( "Read shape: %s (%u)", (readShape ? "yes" : "no"), (unsigned)readShape );
	      HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
	      Pad8Bits();
	      HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
	      
	      
	      
	      bool inputError=false;
	      unsigned trackCount=0;
	      
	      while ( !EndOfBitInput() )
		{
		  AliHLTTPCTrackModel trackModel;
		  memset( &trackModel, 0, sizeof(trackModel) );
		  if ( !InputBytes( (AliHLTUInt8_t*)&trackModel, sizeof(AliHLTTPCTrackModel) ) )
		    {
		      HLTError( "Corrupt input data. Cannot read track model data at position %lu / %u",
				GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
		      inputError = true;
		      break;
		    }
		  HLTInfo( "  Track Model % 5lu fKappa: %f - fPhi: %f - fD: %f - fZ0: %f - fTgl: %f",
			   trackCount, trackModel.fKappa, trackModel.fPhi, trackModel.fD, trackModel.fZ0, trackModel.fTgl );
		  HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
		  
		  Int_t clustercount=0;
		  for(Int_t i=0; i<AliHLTTPCTransform::GetNRows(); i++)
		    {
		      AliHLTTPCClusterModel cluster;
		      memset( &cluster, 0, sizeof(cluster) );
		      AliHLTUInt8_t present;
		      if ( !InputBit( present ) )
			{
			  HLTError( "Corrupt input data. Cannot read  cluster presence bit at position %lu / %u",
				    GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
			  inputError = true;
			  break;
			}
		      HLTInfo( "Cluster %u present: %s (%u)", (unsigned)i, (present ? "yes" : "no"), (unsigned)present );
		      HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
		      cluster.fPresent = present;
		      if ( !present )
			continue;
		      if ( clustercount==0 )
			{
			  if ( !InputBits( slice,6 ) ) //Need 6 bits to encode slice number
			    {
			      HLTError( "Corrupt input data. Cannot read  cluster slice number at position %lu / %u",
					GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
			      inputError = true;
			      break;
			    }
			  HLTInfo( "First cluster slice: %u", (unsigned)slice );
			  HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
			}
		      else
			{
			AliHLTUInt8_t sliceChange;
			if ( !InputBit( sliceChange ) )
			    {
			      HLTError( "Corrupt input data. Cannot read cluster slice change bit at position %lu / %u",
				      GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
			      inputError = true;
			    break;
			    }
			HLTInfo( "Slice change: %s (%u)", (sliceChange ? "yes" : "no"), (unsigned)sliceChange );
			HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
			if ( sliceChange )
			  {  //Change of slice
			    if ( !InputBits( slice, 6 ) )
			      {
				HLTError( "Corrupt input data. Cannot read  cluster slice number at position %lu / %u",
					  GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
				inputError = true;
				break;
				}
			    HLTInfo( "Changed cluster slice: %u", (unsigned)slice );
			    HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
			  }
			}
		    HLTInfo( "Slice. %d", slice );
		    cluster.fSlice = slice;
		    if ( cluster.fSlice<0 || cluster.fSlice>35 )
			{
			  HLTError( "Inconsistent slice number %u (track %u, cluster %d)", cluster.fSlice, trackCount, i );
			  inputError = true;
			  break;
			}
		    AliHLTUInt8_t signBit;
		    Int_t sign;
		    AliHLTUInt64_t temp;
		    Int_t val;
		    //Read time information:
		    if ( !InputBit( signBit ) )
		      {
			HLTError( "Corrupt input data. Cannot read DTime sign bit at position %lu / %u",
				  GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
			inputError = true;
			break;
		      }
		    HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
		    sign = signBit;
		    sign = -1+sign*2;
		    if ( !InputBits( temp, AliHLTTPCCompDataCompressorHelper::GetNTimeBits()-1 ) )
		      {
			HLTError( "Corrupt input data. Cannot read DTime data at position %lu / %u",
				  GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
			inputError = true;
			break;
		      }
		    HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
		    val = (Int_t)temp;
		    cluster.fDTime = val*sign;
		    
		    
		    //Read pad information:
		    if ( !InputBit( signBit ) )
		      {
			HLTError( "Corrupt input data. Cannot read DPad sign bit at position %lu / %u",
				  GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
			inputError = true;
			break;
		      }
		    HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
		    sign = signBit;
		    sign = -1+sign*2;
		    if ( !InputBits( temp, AliHLTTPCCompDataCompressorHelper::GetNPadBits()-1 ) )
		      {
			HLTError( "Corrupt input data. Cannot read DPad data at position %lu / %u",
				  GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
			inputError = true;
			break;
		      }
		    HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
		    val = (Int_t)temp;
		    cluster.fDPad = val*sign;
		    
		    // Read charge information:
		    if ( !InputBits( temp, AliHLTTPCCompDataCompressorHelper::GetNChargeBits() ) )
		      {
			HLTError( "Corrupt input data. Cannot read charge data at position %lu / %u",
				  GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
			inputError = true;
			break;
		      }
		    HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
		    cluster.fDCharge = temp;
		    
		    if ( readShape )
		      {
			// Read shape information:
			if ( !InputBit( signBit ) )
			  {
			    HLTError( "Corrupt input data. Cannot read DSigmaY sign bit at position %lu / %u",
				      GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
			    inputError = true;
			    break;
			  }
			HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
			sign = signBit;
			sign = -1+sign*2;
			if ( !InputBits( temp, AliHLTTPCCompDataCompressorHelper::GetNShapeBits()-1 ) )
			  {
			    HLTError( "Corrupt input data. Cannot read DSigmaY data at position %lu / %u",
				      GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
			    inputError = true;
			    break;
			  }
			HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
			val = (Int_t)temp;
			cluster.fDSigmaY = val*sign;
			
			if ( !InputBit( signBit ) )
			  {
			    HLTError( "Corrupt input data. Cannot read DSigmaZ sign bit at position %lu / %u",
				      GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
			    inputError = true;
			    break;
			  }
			HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
			sign = signBit;
			sign = -1+sign*2;
			if ( !InputBits( temp, AliHLTTPCCompDataCompressorHelper::GetNShapeBits()-1 ) )
			  {
			    HLTError( "Corrupt input data. Cannot read DSigmaZ data at position %lu / %u",
				      GetCurrentByteInputPosition(), GetCurrentBitInputPosition() );
			    inputError = true;
			    break;
			  }
			HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
			val = (Int_t)temp;
			cluster.fDSigmaZ = val*sign;
		      }
		    
		    
		    HLTInfo( "    Cluster % 05u: fDTime: %f - fDPad: %f - fDCharge: %f - fDSigmaY: %f - fDSigmaZ: %f - fNPads: %u - fSlice: %hd - padrow: %lu - fPresent: %u",
			     clustercount, cluster.fDTime, cluster.fDPad, cluster.fDCharge, cluster.fDSigmaY, cluster.fDSigmaZ, cluster.fNPads, cluster.fSlice, (unsigned long)i, (unsigned)cluster.fPresent );
		    
		    
		    clustercount++;
		    
		    
		    }
		  if ( inputError )
		    break;
		  Pad8Bits();
		  HLTInfo( "Input position: %lu / %u (0x%02X)", GetCurrentByteInputPosition(), GetCurrentBitInputPosition(), (unsigned)GetCurrentInputByte() );
		  
		  trackCount++;
		}
	      if ( inputError )
		continue;
	      
	    }
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkRemainingClustersCompressedDataType )
	    {
	      AliHLTUInt8_t* inputPtr = (AliHLTUInt8_t*)blocks[n].fPtr;
	      
	      InitBitDataInput( inputPtr, blocks[n].fSize );
	      AliHLTUInt8_t version;
	      if ( !InputBits( version, 4 ) ) // Version information
		{
		  HLTError( "Corrupt input data. Cannot read data version number at position %u / %u",
			    (unsigned)GetCurrentByteInputPosition(), (unsigned)GetCurrentBitInputPosition() );
		  return EIO;
		}
	      HLTInfo( "Remaining cluster data version: %u", (unsigned)version );
	      if ( version != 0 )
		{
		  HLTError( "Unsupported version %hu. Only version 0 supported currently.", version );
		}
	    Pad8Bits();
	    
	    unsigned long clusterCount=0;
	    for(Int_t lslice=0; lslice<=35; lslice++)
	      {
		for(Int_t lpatch=0; lpatch < 6; lpatch++)
		  {
		    UInt_t i;
		    //Write number of padrows with clusters
		    UInt_t nRows;
		    if ( !InputBits( nRows, 8 ) )
		      {
			HLTError( "Corrupt input data. Cannot read padrow count at position %u / %u",
				  (unsigned)GetCurrentByteInputPosition(), (unsigned)GetCurrentBitInputPosition() );
			return EIO;
		      }
		    HLTInfo( "slice %u patch %u: %u padrows",
			     (unsigned)lslice, (unsigned)lpatch, (unsigned)nRows );
		    if ( !nRows )
		      {
			continue;
		      }
		    //HLTInfo( "  Slice %d - Partition %d", slice, patch );
		    for ( UInt_t jj=0; jj<nRows; jj++ )
		      {
			
			UInt_t padrow;
			if ( !InputBits(padrow,8) ) //Read padrow #
			  {
			    HLTError( "Corrupt input data. Cannot read padrow number at position %u / %u",
				      (unsigned)GetCurrentByteInputPosition(), (unsigned)GetCurrentBitInputPosition() );
			    return EIO;
			  }
			HLTInfo( "Padrow: %u", (unsigned)padrow );
			UInt_t nClusters;
			if ( !InputBits(nClusters,10) )//Read number of clusters on this padrow
			  {
			    HLTError( "Corrupt input data. Cannot read cluster count at position %u / %u",
				      (unsigned)GetCurrentByteInputPosition(), (unsigned)GetCurrentBitInputPosition() );
			    return EIO;
			  }
			HLTInfo( "  #Clusters: %u", (unsigned)nClusters );
			for ( i=0; i<nClusters; i++ )
			    {
			      //Read pad
			      AliHLTTPCRemainingCluster cl;
			      Int_t buff;
			      if ( !InputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNPadBitsRemaining()) )
				{
				  HLTError( "Corrupt input data. Cannot read cluster pad data at position %u / %u",
					    (unsigned)GetCurrentByteInputPosition(), (unsigned)GetCurrentBitInputPosition() );
				  return EIO;
				}
			      cl.fPad = (Float_t)( ((double)buff) / AliHLTTPCCompDataCompressorHelper::GetPadPrecisionFactor() );
			      
			      //Read time
			      if ( !InputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNTimeBitsRemaining()) )
				{
				  HLTError( "Corrupt input data. Cannot read cluster time data at position %u / %u",
					    (unsigned)GetCurrentByteInputPosition(), (unsigned)GetCurrentBitInputPosition() );
				  return EIO;
				}
			      cl.fTime = (Float_t)( ((double)buff) / AliHLTTPCCompDataCompressorHelper::GetTimePrecisionFactor() );
			      
			      //Read widths
			      if ( !InputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNShapeBitsRemaining()) )
				{
				  HLTError( "Corrupt input data. Cannot read cluster pad width data at position %u / %u",
					    (unsigned)GetCurrentByteInputPosition(), (unsigned)GetCurrentBitInputPosition() );
				  return EIO;
				}
			      Float_t padw = (Float_t)( ((double)buff) / AliHLTTPCCompDataCompressorHelper::GetPadPrecisionFactor() );
			    cl.fSigmaY2 = padw*padw;
			    
			    if ( !InputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNShapeBitsRemaining()) )
			      {
				HLTError( "Corrupt input data. Cannot read cluster time width data at position %u / %u",
					  (unsigned)GetCurrentByteInputPosition(), (unsigned)GetCurrentBitInputPosition() );
				return EIO;
			      }
			    Float_t timew = (Float_t)( ((double)buff) / AliHLTTPCCompDataCompressorHelper::GetTimePrecisionFactor() );
			    cl.fSigmaZ2 = timew*timew;
			    
			    //Read charge 
			    if ( !InputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNChargeBits()) )
			      {
				HLTError( "Corrupt input data. Cannot read cluster charge data at position %u / %u",
					  (unsigned)GetCurrentByteInputPosition(), (unsigned)GetCurrentBitInputPosition() );
				return EIO;
			      }
			    cl.fCharge = buff;
			    
			    HLTInfo( "  Cluster % 5lu (% 5u): fPadRow: %u - fPad: %f - fTime: %f - fSigmaY2: %f - fSigmaZ2: %f - fCharge: %hu",
				     clusterCount, (unsigned)i, (unsigned)padrow, cl.fPad, cl.fTime, cl.fSigmaY2, cl.fSigmaZ2, cl.fCharge );
			    clusterCount++;
			    }
		      }
		  }
		
	      }
	    
	    
	    }
	  
	}
      
      size = 0;
      return 0;
    }
