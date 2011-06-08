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

/** @file   AliHLTTPCCompModelConverterComponent.cxx
    @author Timm Steinbeck
    @author changed by J. Wagner
    @date   17-11-2007
    @brief  A copy processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCCompModelConverterComponent.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTTPCSpacePointContainer.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTComponentBenchmark.h"
//#include "AliHLTTPCCompModelAnalysis.h"
#include <errno.h>
#include <vector>

/** An implementiation of a converter component that
 * takes in clusters and tracks in the standard HLT format
 * and converts them into the Vestbo-format
 * such that the Vestbo compression can then be 
 * applied to these tracks and clusters
 */

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCCompModelConverterComponent)
    
AliHLTTPCCompModelConverterComponent::AliHLTTPCCompModelConverterComponent()
  : AliHLTProcessor()
  , fConverter(NULL)
  , fModelAnalysisInstance(NULL)
  , fDumpFileName()
  , fGraphFileName()
  , fModelAnalysis(0)
  , fTrackAnalysis(0)
  , fFillingFirstTrackArray(0)
  , fInputClusters(NULL)
  , fpBenchmark(NULL)
{
  // constructor
}

AliHLTTPCCompModelConverterComponent::~AliHLTTPCCompModelConverterComponent()
{
  // destructor
}

const char* AliHLTTPCCompModelConverterComponent::GetComponentID()
{
  // AliHLTComponent interface function: return component id
  return "TPCCompModelConverter";
}

void AliHLTTPCCompModelConverterComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
{
  // AliHLTComponent interface function: get list of input data types
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back( AliHLTTPCDefinitions::fgkClustersDataType );
  list.push_back( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC );
}

AliHLTComponent_DataType AliHLTTPCCompModelConverterComponent::GetOutputDataType()
{
  // AliHLTComponent interface function: get output data type
  return AliHLTTPCDefinitions::fgkClusterTracksModelDataType;
}

void AliHLTTPCCompModelConverterComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // AliHLTComponent interface function: output data size estimator
  constBase = 4+4+216; // Format versions + 1 byte per patch
  inputMultiplier = 4.;
}

AliHLTComponent* AliHLTTPCCompModelConverterComponent::Spawn()
{
  // AliHLTComponent interface function: return new instance of this class
  return new AliHLTTPCCompModelConverterComponent;
}

int AliHLTTPCCompModelConverterComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  Int_t i = 0;
  // get input argument either -modelanalysis or -trackanalysis
  while(i < argc)
    {
      if ( !strcmp( argv[i], "-modelanalysis" ) ) 
	{
	  fModelAnalysis = 1;
	  HLTInfo("Model analysis starting.");
	  ++i;
	  continue;
	}
 
	
      if ( !strcmp( argv[i], "-trackanalysis" ) ) 
	{
	  fTrackAnalysis = 1;
	  fFillingFirstTrackArray = 1;
	  HLTInfo("Tracking analysis starting.");
	  ++i;
	  continue;
	}

      if ( !strcmp( argv[i], "-dumptofile" ) ) 
	{

	  //check if any analysis has been specified (otherwise -dumptofile makes no sense!)
	  if(!fTrackAnalysis && !fModelAnalysis)
	    {
	      HLTError("Dump to file called without any model analysis specified.");
	      return EINVAL;
	    }

	  // read in filename (including path)
	   if ( argc <= i+1 ) 
	    {
	      HLTError("Missing filename to write analysis results");
	      return ENOTSUP;
	    }
	   
	   fDumpFileName = argv[i+1];
	   HLTInfo("File name of dump file for results set to %s.", fDumpFileName.Data());
	  ++i;
	  ++i;
	  continue;
	}

      // specify if graphical output is wanted (histograms, saved in a root file)
      if ( !strcmp( argv[i], "-graphs" ) ) 
	{

	  //check if any analysis has been specified (otherwise -graphs makes no sense!)
	  if(!fTrackAnalysis && !fModelAnalysis)
	    {
	      HLTError("Creation of histgrams called without any model analysis specified.");
	      return EINVAL;
	    }

	  // read in filename (including path like /afsuser/johndoe/TrackModelHistograms.root)
	   if ( argc <= i+1 ) 
	    {
	      HLTError("Missing filename to write histograms");
	      return ENOTSUP;
	    }
	   
	   fGraphFileName = argv[i+1];
	   HLTInfo("File name of file for graphical results set to %s.", fGraphFileName.Data());
	  ++i;
	  ++i;
	  continue;
	}

      HLTError("Unknown Option '%s'", argv[i] );
      return EINVAL;

    }
  
  // start new analysis by intialising respective arrays
  if(fModelAnalysis || fTrackAnalysis)
    {
      fModelAnalysisInstance = new AliHLTTPCCompModelAnalysis(fModelAnalysis, fTrackAnalysis, fDumpFileName, fGraphFileName);
      fModelAnalysisInstance->Init(); 
      fConverter = new AliHLTTPCCompModelConverter(fModelAnalysisInstance);

      if(fModelAnalysis)
	{
	  HLTInfo("Model Analysis initiated, calculating loss due to convertion to Vestbo-Model.");
	}
      else
	{
	  HLTInfo("Track Analysis initiated, showing influence of Vestbo-Model on tracking.");
	}
    }
  else
    {
      fConverter = new AliHLTTPCCompModelConverter();
    }
 
  /*if ( argc )
    {
    Logging( kHLTLogDebug, "HLT::TPCCompModelConverter::DoInit", "Arguments", "argv[0] == %s", argv[0] );
    Logging(kHLTLogError, "HLT::TPCCompModelConverter::DoInit", "Unknown Option", "Unknown option '%s'", argv[0] );
    return EINVAL;
    }*/
  
  fpBenchmark=new AliHLTComponentBenchmark;
  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->SetTimer(0,"total");
    GetBenchmarkInstance()->SetTimer(1,"clusterinput");
    GetBenchmarkInstance()->SetTimer(2,"trackinput");
  }

  return 0;
}

int AliHLTTPCCompModelConverterComponent::DoDeinit()
    {
      // see header file for class documentation
      // display results
      if(fModelAnalysisInstance)
	{
	  fModelAnalysisInstance->DisplayResults();

	  delete fModelAnalysisInstance;
	  fModelAnalysisInstance = NULL;
	};

      if(fConverter)
	{
	  delete fConverter;
	  fConverter = NULL;
	}

    return 0;
    }

int AliHLTTPCCompModelConverterComponent::DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
						   AliHLTComponent_TriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
				      AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks )
    {
      // see header file for class documentation
      int iResult=0;
      AliHLTUInt32_t capacity=size;
      size=0;
      if (!IsDataEvent()) return 0;

      if (GetBenchmarkInstance()) {
	GetBenchmarkInstance()->StartNewEvent();
	GetBenchmarkInstance()->Start(0);
      }

      fConverter->Init();
      // Process an event
      // Loop over all input blocks in the event
      AliHLTUInt8_t minSlice=0xFF, maxSlice=0xFF, minPatch=0xFF, maxPatch=0xFF;
      const AliHLTComponentBlockData* pDesc=NULL;

      /// input track array
      vector<AliHLTGlobalBarrelTrack> inputTrackArray;

      if (GetBenchmarkInstance()) {
	GetBenchmarkInstance()->Start(1);
      }
      for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType);
	   pDesc!=NULL; pDesc=GetNextInputBlock()) {
	if (GetBenchmarkInstance()) {
	  GetBenchmarkInstance()->AddInput(pDesc->fSize);
	}
	AliHLTUInt8_t slice = 0;
	AliHLTUInt8_t patch = 0;
	slice = AliHLTTPCDefinitions::GetMinSliceNr( pDesc->fSpecification );
	patch = AliHLTTPCDefinitions::GetMinPatchNr( pDesc->fSpecification );
	if ( minSlice==0xFF || slice<minSlice )	minSlice = slice;
	if ( maxSlice==0xFF || slice>maxSlice )	maxSlice = slice;
	if ( minPatch==0xFF || patch<minPatch )	minPatch = patch;
	if ( maxPatch==0xFF || patch>maxPatch )	maxPatch = patch;
	if (!fInputClusters) {
	  fInputClusters=new AliHLTTPCSpacePointContainer;	  
	  if (!fInputClusters) return -ENOMEM;
	}
	if (fInputClusters) {
	  fInputClusters->AddInputBlock(pDesc);
	}
      }
      if (GetBenchmarkInstance()) {
	GetBenchmarkInstance()->Stop(1);
	GetBenchmarkInstance()->Start(2);
      }

      for (pDesc=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
	   pDesc!=NULL; pDesc=GetNextInputBlock()) {
	if (GetBenchmarkInstance()) {
	  GetBenchmarkInstance()->AddInput(pDesc->fSize);
	}
	AliHLTUInt8_t slice = 0;
	AliHLTUInt8_t patch = 0;
	slice = AliHLTTPCDefinitions::GetMinSliceNr( pDesc->fSpecification );
	patch = AliHLTTPCDefinitions::GetMinPatchNr( pDesc->fSpecification );
	if ( minSlice==0xFF || slice<minSlice )	minSlice = slice;
	if ( maxSlice==0xFF || slice>maxSlice )	maxSlice = slice;
	if ( minPatch==0xFF || patch<minPatch )	minPatch = patch;
	if ( maxPatch==0xFF || patch>maxPatch )	maxPatch = patch;
	const AliHLTTracksData* pTracks=reinterpret_cast<const AliHLTTracksData*>(pDesc->fPtr);
	if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(pTracks, pDesc->fSize, inputTrackArray))<0) {
	  return iResult;
	}
	for (vector<AliHLTGlobalBarrelTrack>::const_iterator track=inputTrackArray.begin();
	     track!=inputTrackArray.end();
	     track++) {
	  if (!fInputClusters) continue;
	  int trackID=track->GetID();
	  if (trackID<0) {
	    // FIXME: error guard
	    HLTError("invalid track ID");
	    continue;
	  }
	  if ((iResult=fInputClusters->SetTrackID(trackID, track->GetPoints(), track->GetNumberOfPoints()))<0) {
	    HLTError("failed to set cluster id for track %d: error %d", trackID, iResult);
	    iResult=0;
	    continue;
	  }
	}
	// for (vector<AliHLTGlobalBarrelTrack>::const_iterator track=inputTrackArray.begin();
	//      track!=inputTrackArray.end();
	//      track++) {
	//   AliHLTSpacePointContainer* trackClusters=fInputClusters->SelectByTrack(track->GetID());
	//   if (!trackClusters) {
	//     HLTError("failed to select clusters assigned to track %d", track->GetID());
	//     continue;
	//   }
	//   track->Print();
	//   trackClusters->Print();
	// }
      }

      if (GetBenchmarkInstance()) {
	GetBenchmarkInstance()->Stop(2);
      }

      for ( unsigned long n = 0; n < evtData.fBlockCnt; n++ )
	{
	  AliHLTUInt8_t slice = 0;
	  AliHLTUInt8_t patch = 0;

	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkClustersDataType ||
	       blocks[n].fDataType == AliHLTTPCDefinitions::fgkTracksDataType )
	    {
	      slice = AliHLTTPCDefinitions::GetMinSliceNr( blocks[n].fSpecification );
	      patch = AliHLTTPCDefinitions::GetMinPatchNr( blocks[n].fSpecification );
	      if ( minSlice==0xFF || slice<minSlice )
		minSlice = slice;
	      if ( maxSlice==0xFF || slice>maxSlice )
		maxSlice = slice;
	      if ( minPatch==0xFF || patch<minPatch )
		minPatch = patch;
	      if ( maxPatch==0xFF || patch>maxPatch )
		maxPatch = patch;
	    }
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkClustersDataType )
	    {
	      fConverter->SetInputClusters( (AliHLTTPCClusterData*)blocks[n].fPtr, slice, patch );
	    }
	  if ( blocks[n].fDataType == kAliHLTDataTypeTrack )
	    {
	      fConverter->SetInputTracks( reinterpret_cast<AliHLTTracksData*>(blocks[n].fPtr), blocks[n].fSize );
	      
	      // if track analysis is desired, fill tracklets into track arrays of ModelAnalysis class to be compared
	      if(fTrackAnalysis)
		{
		  fModelAnalysisInstance->SetTracks((AliHLTTPCTrackletData*)blocks[n].fPtr, fFillingFirstTrackArray);
		  
		  // set flag for filling first array to zero --> second array is filled then
		  fFillingFirstTrackArray = 0;
		}
	    }
	}
      
      if(fTrackAnalysis)
	{
	  if(fModelAnalysis == 0) // stop processing if not required
	    return 0;
	};
      
      fConverter->Convert();
      
      unsigned long dataSize=0;
      unsigned long outputSize = fConverter->GetOutputModelDataSize();
      if ( outputSize+dataSize> capacity )
	{
	  HLTError( "Not enough output memory size for clusters&tracks model data. %lu needed",
		    outputSize );
	  return -ENOSPC;
	}

      fConverter->OutputModelData( outputPtr, outputSize );
      
      AliHLTComponent_BlockData ob;
      // Let the structure be filled with the default values.
      // This takes care of setting the shared memory and data type values to default values,
      // so that they can be filled in by the calling code.
      FillBlockData( ob );
      // This block's start (offset) is after all other blocks written so far
      ob.fOffset = dataSize;
      // the size of this block's data.
      ob.fSize = outputSize;
      // The specification of the data is copied from the input block.
      ob.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( minSlice, maxSlice, minPatch, maxPatch );
      // The type of the data is copied from the input block.
      ob.fDataType = AliHLTTPCDefinitions::fgkClusterTracksModelDataType;
      // Place this block into the list of output blocks
      outputBlocks.push_back( ob );
      
      outputPtr += ob.fSize;
      dataSize += ob.fSize;
      
      if ( dataSize+fConverter->GetRemainingClustersOutputDataSize()>capacity )
	{
	  HLTError( "Not enough output memory size for remaining clusters model data. %lu needed in total (clusters&tracks + rem. clusters)",
		    dataSize+fConverter->GetRemainingClustersOutputDataSize() );
	  return -ENOSPC;
	}
      fConverter->GetRemainingClusters( outputPtr, outputSize );
      printf( "clusterSize1: %lu\n", outputSize );
      
      FillBlockData( ob );
      // This block's start (offset) is after all other blocks written so far
      ob.fOffset = dataSize;
      // the size of this block's data.
      ob.fSize = outputSize;
      // The specification of the data is copied from the input block.
      ob.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( minSlice, maxSlice, minPatch, maxPatch );
      // The type of the data is copied from the input block.
      ob.fDataType = AliHLTTPCDefinitions::fgkRemainingClustersModelDataType;
      // Place this block into the list of output blocks
      outputBlocks.push_back( ob );
      
      outputPtr += ob.fSize;
      dataSize += ob.fSize;

      if (GetBenchmarkInstance()) {
	GetBenchmarkInstance()->Stop(0);
	HLTBenchmark(GetBenchmarkInstance()->GetStatistics());
      }
      
      // Finally we set the total size of output memory we consumed.
      size = dataSize;
      return 0;
    }
