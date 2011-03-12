// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de> *
//                  for The ALICE HLT Project.                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************


#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCTrackMCMarkerComponent.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCCADef.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCSpacePointData.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTTrackMCLabel.h"
#include "AliHLTTPCClusterFinder.h"
#include <climits>
#include <cstdlib>
#include <cerrno>


// ROOT macro for the implementation of ROOT specific class methods
ClassImp( AliHLTTPCTrackMCMarkerComponent )


AliHLTTPCTrackMCMarkerComponent::AliHLTTPCTrackMCMarkerComponent()
{
  // see header file for class documentation
  for( int i=0; i<36*6; i++ ){
    fClusterLabels[i] = 0;
    fNClusterLabels[i] = 0; 
  }
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char *AliHLTTPCTrackMCMarkerComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCTrackMCMarker";
}

void AliHLTTPCTrackMCMarkerComponent::GetInputDataTypes( AliHLTComponentDataTypeList &list )
{
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC );
  list.push_back(AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo);
}

AliHLTComponentDataType AliHLTTPCTrackMCMarkerComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeTrackMC|kAliHLTDataOriginTPC;
}

void AliHLTTPCTrackMCMarkerComponent::GetOutputDataSize( unsigned long &constBase, double &inputMultiplier )
{
  // see header file for class documentation
  constBase = 0;
  inputMultiplier = 1.0;
}

AliHLTComponent *AliHLTTPCTrackMCMarkerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCTrackMCMarkerComponent;
}




void AliHLTTPCTrackMCMarkerComponent::SetDefaultConfiguration()
{
  // Set default configuration for the CA merger component
  // Some parameters can be later overwritten from the OCDB
}

int AliHLTTPCTrackMCMarkerComponent::ReadConfigurationString(  const char* arguments )
{
  // Set configuration parameters for the CA merger component from the string

  int iResult = 0;
  if ( !arguments ) return iResult;

  TString allArgs = arguments;
  TString argument;
  int bMissingParam = 0;

  TObjArray* pTokens = allArgs.Tokenize( " " );
  
  int nArgs =  pTokens ? pTokens->GetEntries() : 0;
  
  for ( int i = 0; i < nArgs; i++ ) {
    argument = ( ( TObjString* )pTokens->At( i ) )->GetString();
    if ( argument.IsNull() ) continue;

    HLTError( "Unknown option \"%s\"", argument.Data() );
    iResult = -EINVAL;
  }
  delete pTokens;

  if ( bMissingParam ) {
    HLTError( "Specifier missed for parameter \"%s\"", argument.Data() );
    iResult = -EINVAL;
  }

  return iResult;
}


int AliHLTTPCTrackMCMarkerComponent::ReadCDBEntry( const char* /*cdbEntry*/, const char* /*chainId*/ )
{
  // see header file for class documentation

  // no settings for a moment, no CDB path, threfore return

  return 0;
/*
  const char* defaultNotify = "";

  if ( !cdbEntry ) {
    cdbEntry = "HLT/ConfigTPC/TPCTrackMCMarker";
    defaultNotify = " (default)";
    chainId = 0;
  }

  HLTInfo( "configure from entry \"%s\"%s, chain id %s", cdbEntry, defaultNotify, ( chainId != NULL && chainId[0] != 0 ) ? chainId : "<none>" );
  AliCDBEntry *pEntry = AliCDBManager::Instance()->Get( cdbEntry );//,GetRunNo());

  if ( !pEntry ) {
    HLTError( "cannot fetch object \"%s\" from CDB", cdbEntry );
    return -EINVAL;
  }

  TObjString* pString = dynamic_cast<TObjString*>( pEntry->GetObject() );

  if ( !pString ) {
    HLTError( "configuration object \"%s\" has wrong type, required TObjString", cdbEntry );
    return -EINVAL;
  }

  HLTInfo( "received configuration object string: \"%s\"", pString->GetString().Data() );

  return  ReadConfigurationString( pString->GetString().Data() );
*/
}



int AliHLTTPCTrackMCMarkerComponent::Configure( const char* cdbEntry, const char* chainId, const char *commandLine )
{
  // Configure the component
  // There are few levels of configuration,
  // parameters which are set on one step can be overwritten on the next step

  //* read hard-coded values

  SetDefaultConfiguration();

  //* read the default CDB entry

  int iResult1 = ReadCDBEntry( NULL, chainId );

  //* read the actual CDB entry if required

  int iResult2 = ( cdbEntry ) ? ReadCDBEntry( cdbEntry, chainId ) : 0;

  //* read extra parameters from input (if they are)

  int iResult3 = 0;

  if ( commandLine && commandLine[0] != '\0' ) {
    HLTInfo( "received configuration string from HLT framework: \"%s\"", commandLine );
    iResult3 = ReadConfigurationString( commandLine );
  }

  return iResult1 ? iResult1 : ( iResult2 ? iResult2 : iResult3 );
}




int AliHLTTPCTrackMCMarkerComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation

  TString arguments = "";
  for ( int i = 0; i < argc; i++ ) {
    if ( !arguments.IsNull() ) arguments += " ";
    arguments += argv[i];
  }

  return Configure( NULL, NULL, arguments.Data()  );
}


int AliHLTTPCTrackMCMarkerComponent::Reconfigure( const char* cdbEntry, const char* chainId )
{
  // Reconfigure the component from OCDB

  return Configure( cdbEntry, chainId, NULL );
}



int AliHLTTPCTrackMCMarkerComponent::DoDeinit()
{
  // see header file for class documentation

  return 0;
}

Int_t AliHLTTPCTrackMCMarkerComponent::GetTrackMCLabel( unsigned int *hits, int nHits )
{
  // get MC label for the track	 

  Int_t mcLabel = -1;
	    
  std::vector<int> labels;

  for( Int_t ih=0; ih<nHits; ih++){
    UInt_t id = hits[ih];
    int iSlice = AliHLTTPCSpacePointData::GetSlice(id);
    int iPatch = AliHLTTPCSpacePointData::GetPatch(id);
    int iCluster = AliHLTTPCSpacePointData::GetNumber(id);
    if( iSlice<0 || iSlice>=36 || iPatch<0 || iPatch>5 ){
      HLTError("Corrupted TPC cluster Id: slice %d, patch %d, cluster %d",
	       iSlice, iPatch,iCluster );
      continue;
    }
    AliHLTTPCClusterFinder::ClusterMCInfo *patchLabels = fClusterLabels[iSlice*6 + iPatch];
    if( !patchLabels ) continue;
    if( iCluster >= fNClusterLabels[iSlice*6 + iPatch] ){
      HLTError("TPC slice %d, patch %d: ClusterID==%d >= N MC labels==%d ",
	       iSlice, iPatch,iCluster, fNClusterLabels[iSlice*6 + iPatch] );
      continue;
    }
    AliHLTTPCClusterFinder::ClusterMCInfo &lab = patchLabels[iCluster];	    
    if ( lab.fClusterID[0].fMCID >= 0 ) labels.push_back( lab.fClusterID[0].fMCID );
    if ( lab.fClusterID[1].fMCID >= 0 ) labels.push_back( lab.fClusterID[1].fMCID );
    if ( lab.fClusterID[2].fMCID >= 0 ) labels.push_back( lab.fClusterID[2].fMCID );
  }
	  
  std::sort( labels.begin(), labels.end() );
	  
  labels.push_back( -1 ); // put -1 to the end
	  
  int labelMax = -1, labelCur = -1, nLabelsMax = 0, nLabelsCurr = 0;
  for ( unsigned int iLab = 0; iLab < labels.size(); iLab++ ) {
    if ( labels[iLab] != labelCur ) {	      
      if ( labelCur >= 0 && nLabelsMax< nLabelsCurr ) {
	nLabelsMax = nLabelsCurr;
	labelMax = labelCur;
      }
      labelCur = labels[iLab];
      nLabelsCurr = 0;
    }
    nLabelsCurr++;
  }
  
  if( labelMax>=0 && nLabelsMax < 0.9 * nHits ) labelMax = -labelMax;
  
  mcLabel = labelMax;

  return mcLabel;
}


int AliHLTTPCTrackMCMarkerComponent::DoEvent( const AliHLTComponentEventData &evtData,
    const AliHLTComponentBlockData *blocks, AliHLTComponentTriggerData &/*trigData*/,
    AliHLTUInt8_t *outputPtr, AliHLTUInt32_t &size, AliHLTComponentBlockDataList &outputBlocks )
{
  // see header file for class documentation

  int iResult = 0;
  unsigned int maxBufferSize = size;

  size = 0;

  if ( !outputPtr ) {
    return -ENOSPC;
  }
  if ( !IsDataEvent() ) {
    return 0;
  }

  for( int i=0; i<36*6; i++ ){
    fClusterLabels[i] = 0;
    fNClusterLabels[i] = 0;
  }

  int nBlocks = (int)evtData.fBlockCnt;

  int nInputMCLabels = 0;
  int nInputTracks = 0;

  // first read all the MC information
  for (int ndx=0; ndx<nBlocks && iResult>=0; ndx++) {
    const AliHLTComponentBlockData* iter = blocks+ndx;
    if(iter->fDataType == AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo ) {
      Int_t slice=AliHLTTPCDefinitions::GetMinSliceNr(iter->fSpecification);
      Int_t patch=AliHLTTPCDefinitions::GetMinPatchNr(iter->fSpecification);
      fClusterLabels[ slice*6 + patch] = (AliHLTTPCClusterFinder::ClusterMCInfo *)iter->fPtr;
      fNClusterLabels[ slice*6 + patch] = iter->fSize/sizeof(AliHLTTPCClusterFinder::ClusterMCInfo);
      nInputMCLabels+=fNClusterLabels[ slice*6 + patch];
    }
  }
      
  // read tracks and write output

  
  unsigned int mySize = 0;
  
  AliHLTTrackMCData* outPtr = ( AliHLTTrackMCData* )( outputPtr );
  AliHLTTrackMCLabel* currOutLabel = outPtr->fLabels;

  mySize =   ( ( AliHLTUInt8_t * )currOutLabel ) -  ( ( AliHLTUInt8_t * )outPtr );

  outPtr->fCount = 0;

  for (int ndx=0; ndx<nBlocks && iResult>=0; ndx++) {
    const AliHLTComponentBlockData* iter = blocks+ndx;
	
    if( iter->fDataType == ( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC ) ){	        
      AliHLTTracksData* dataPtr = ( AliHLTTracksData* ) iter->fPtr;
      int nTracks = dataPtr->fCount;
      AliHLTExternalTrackParam* currTrack = dataPtr->fTracklets;
      for( int itr=0; itr<nTracks; itr++ ){
	nInputTracks++;
	Int_t mcLabel = GetTrackMCLabel( currTrack->fPointIDs, currTrack->fNPoints );  

	currOutLabel->fTrackID = currTrack->fTrackID;
	currOutLabel->fMCLabel = mcLabel;

	if ( mySize + sizeof(AliHLTTrackMCLabel) > maxBufferSize ) {
	  HLTWarning( "Output buffer size exceed (buffer size %d, current size %d), %d mc labels are not stored", maxBufferSize, mySize, nTracks - itr + 1 );
	  iResult = -ENOSPC;
	  break;
	}
	mySize += sizeof(AliHLTTrackMCLabel);
	currOutLabel++;
	outPtr->fCount++;

	unsigned int dSize = sizeof( AliHLTExternalTrackParam ) + currTrack->fNPoints * sizeof( unsigned int );
	currTrack = ( AliHLTExternalTrackParam* )( (( Byte_t * )currTrack) + dSize );
      }
    }
  }
  

  AliHLTComponentBlockData resultData;
  FillBlockData( resultData );
  resultData.fOffset = 0;
  resultData.fSize = mySize;
  resultData.fDataType = kAliHLTDataTypeTrackMC|kAliHLTDataOriginTPC;
  resultData.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( 0, 35, 0, 5 );
  outputBlocks.push_back( resultData );
  size = resultData.fSize;

  HLTInfo( "TrackMCMarker:: input %d labels, %d tracks, output %d labels", nInputMCLabels, nInputTracks,outPtr->fCount );

  return iResult;
}

