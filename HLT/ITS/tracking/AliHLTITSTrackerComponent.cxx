// $Id: AliHLTITSTrackerComponent.cxx 32659 2009-06-02 16:08:40Z sgorbuno $
// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de> *
//                  Ivan Kisel <kisel@kip.uni-heidelberg.de>                *
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

///  @file   AliHLTITSTrackerComponent.cxx
///  @author Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de>
///  @date   June 2009
///  @brief  An ITS tracker processing component for the HLT


/////////////////////////////////////////////////////
//                                                 //
// a ITS tracker processing component for the HLT  //
//                                                 //
/////////////////////////////////////////////////////

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTITSTrackerComponent.h"
#include "AliHLTArray.h"
#include "AliExternalTrackParam.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "AliITStrackerHLT.h"
#include "AliHLTITSSpacePointData.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliHLTDataTypes.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliGeomManager.h"
#include "AliHLTTrackMCLabel.h"
#include <map>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp( AliHLTITSTrackerComponent )
AliHLTITSTrackerComponent::AliHLTITSTrackerComponent()
    :
    fSolenoidBz( 0 ),
    fBenchmark("ITSTracker"),
    fTracker(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTITSTrackerComponent::AliHLTITSTrackerComponent( const AliHLTITSTrackerComponent& )
    :
    AliHLTProcessor(),
    fSolenoidBz( 0 ),
    fBenchmark("ITSTracker"),
    fTracker(0)
{
  // see header file for class documentation
  HLTFatal( "copy constructor untested" );
}

AliHLTITSTrackerComponent& AliHLTITSTrackerComponent::operator=( const AliHLTITSTrackerComponent& )
{
  // see header file for class documentation
  HLTFatal( "assignment operator untested" );
  return *this;
}

AliHLTITSTrackerComponent::~AliHLTITSTrackerComponent()
{
  // see header file for class documentation
  delete fTracker;
}

//
// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process
//

const char* AliHLTITSTrackerComponent::GetComponentID()
{
  // see header file for class documentation
  return "ITSTracker";
}

void AliHLTITSTrackerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list )
{
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC );
  list.push_back( kAliHLTDataTypeTrackMC|kAliHLTDataOriginTPC );
  list.push_back( kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD );
  list.push_back( kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD );
  list.push_back( kAliHLTDataTypeClusters|kAliHLTDataOriginITSSDD );
  list.push_back( kAliHLTDataTypeClusters|kAliHLTDataOriginITS );
}

AliHLTComponentDataType AliHLTITSTrackerComponent::GetOutputDataType()
{
  // see header file for class documentation  
  return kAliHLTMultipleDataType;
}

int AliHLTITSTrackerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // see header file for class documentation  
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeTrack|kAliHLTDataOriginITS);
  tgtList.push_back(kAliHLTDataTypeTrack|kAliHLTDataOriginITSOut);
  tgtList.push_back( kAliHLTDataTypeTrackMC|kAliHLTDataOriginITS );
  return tgtList.size();
}

void AliHLTITSTrackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // define guess for the output data size
  constBase = 200;       // minimum size
  inputMultiplier = 2.; // size relative to input
}

AliHLTComponent* AliHLTITSTrackerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTITSTrackerComponent;
}

void AliHLTITSTrackerComponent::SetDefaultConfiguration()
{
  // Set default configuration for the CA tracker component
  // Some parameters can be later overwritten from the OCDB

  fSolenoidBz = -5.00668;
  
}

int AliHLTITSTrackerComponent::ReadConfigurationString(  const char* arguments )
{
  // Set configuration parameters for the CA tracker component from the string

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

    if ( argument.CompareTo( "-solenoidBz" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      HLTWarning("argument -solenoidBz is deprecated, magnetic field set up globally (%f)", GetBz());
      continue;
    }

    //if ( argument.CompareTo( "-minNClustersOnTrack" ) == 0 ) {
    //if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
    //fMinNTrackClusters = ( ( TObjString* )pTokens->At( i ) )->GetString().Atoi();
    //HLTInfo( "minNClustersOnTrack set to: %d", fMinNTrackClusters );
    //continue;
    //}

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


int AliHLTITSTrackerComponent::ReadCDBEntry( const char* cdbEntry, const char* chainId )
{
  // see header file for class documentation

  const char* defaultNotify = "";

  if ( !cdbEntry ) {
    return 0;// need to add the HLT/ConfigITS/ITSTracker directory to cdb SG!!!
    //cdbEntry = "HLT/ConfigITS/ITSTracker";
    //defaultNotify = " (default)";
    //chainId = 0;
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
}


int AliHLTITSTrackerComponent::Configure( const char* cdbEntry, const char* chainId, const char *commandLine )
{
  // Configure the component
  // There are few levels of configuration,
  // parameters which are set on one step can be overwritten on the next step

  //* read hard-coded values

  SetDefaultConfiguration();

  //* read the default CDB entry

  int iResult1 = ReadCDBEntry( NULL, chainId );

  //* read magnetic field

  int iResult2 = 0; //ReadCDBEntry( kAliHLTCDBSolenoidBz, chainId );
  fSolenoidBz = GetBz();

  //* read the actual CDB entry if required

  int iResult3 = ( cdbEntry ) ? ReadCDBEntry( cdbEntry, chainId ) : 0;

  //* read extra parameters from input (if they are)

  int iResult4 = 0;

  if ( commandLine && commandLine[0] != '\0' ) {
    HLTInfo( "received configuration string from HLT framework: \"%s\"", commandLine );
    iResult4 = ReadConfigurationString( commandLine );
  }

  // Initialise the tracker here

  return iResult1 ? iResult1 : ( iResult2 ? iResult2 : ( iResult3 ? iResult3 : iResult4 ) );
}



int AliHLTITSTrackerComponent::DoInit( int argc, const char** argv )
{
  // Configure the ITS tracker component

  if ( fTracker ) return -EINPROGRESS;

  if(AliGeomManager::GetGeometry()==NULL){
    AliGeomManager::LoadGeometry();
  }
  AliGeomManager::ApplyAlignObjsFromCDB("ITS");

  TString arguments = "";
  for ( int i = 0; i < argc; i++ ) {
    if ( !arguments.IsNull() ) arguments += " ";
    arguments += argv[i];
  }

  int ret = Configure( NULL, NULL, arguments.Data() );

  // Check field
  if (!TGeoGlobalMagField::Instance()) {
    HLTError("magnetic field not initialized, please set up TGeoGlobalMagField and AliMagF");
    return -ENODEV;
  }
  fSolenoidBz=GetBz();

  fTracker = new AliITStrackerHLT(0);
  fBenchmark.Reset();
  fBenchmark.SetTimer(0,"total");
  fBenchmark.SetTimer(1,"reco");
  return ret;
}


int AliHLTITSTrackerComponent::DoDeinit()
{
  // see header file for class documentation
  delete fTracker;
  fTracker = 0;
  return 0;
}



int AliHLTITSTrackerComponent::Reconfigure( const char* cdbEntry, const char* chainId )
{
  // Reconfigure the component from OCDB .

  return Configure( cdbEntry, chainId, NULL );
}



int AliHLTITSTrackerComponent::DoEvent
(
  const AliHLTComponentEventData& evtData,
  const AliHLTComponentBlockData* blocks,
  AliHLTComponentTriggerData& /*trigData*/,
  AliHLTUInt8_t* outputPtr,
  AliHLTUInt32_t& size,
  vector<AliHLTComponentBlockData>& outputBlocks )
{
  //* process event

  AliHLTUInt32_t maxBufferSize = size;
  size = 0; // output size
  
  if (!IsDataEvent()) return 0;

  if ( evtData.fBlockCnt <= 0 ) {
    HLTWarning( "no blocks in event" );
    return 0;
  }

  fBenchmark.StartNewEvent();
  fBenchmark.Start(0);

  // Event reconstruction in ITS

  int iResult=0;

  int nBlocks = evtData.fBlockCnt;

  
  vector< AliExternalTrackParam > tracksTPC;
  vector< int > tracksTPCLab;
  vector< int > tracksTPCId;

  int nClustersTotal = 0;

  for (int ndx=0; ndx<nBlocks && iResult>=0; ndx++) {

    const AliHLTComponentBlockData* iter = blocks+ndx;
 
    if ( (iter->fDataType == (kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD) ) || 
	 (iter->fDataType == (kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD) ) ||
	 (iter->fDataType == (kAliHLTDataTypeClusters|kAliHLTDataOriginITSSDD) ) ||
	 (iter->fDataType == (kAliHLTDataTypeClusters|kAliHLTDataOriginITS) ) 
	 ){      
      AliHLTITSClusterData *inPtr=reinterpret_cast<AliHLTITSClusterData*>( iter->fPtr );
      nClustersTotal+=inPtr->fSpacePointCnt;
    }         
  }


  fTracker->StartLoadClusters(nClustersTotal);

  // first read MC information (if present)
  
  std::map<int,int> mcLabels;

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrackMC|kAliHLTDataOriginTPC);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    
    fBenchmark.AddInput(pBlock->fSize);
    
    AliHLTTrackMCData* dataPtr = reinterpret_cast<AliHLTTrackMCData*>( pBlock->fPtr );
    if (sizeof(AliHLTTrackMCData)+dataPtr->fCount*sizeof(AliHLTTrackMCLabel)==pBlock->fSize) {
      for( unsigned int il=0; il<dataPtr->fCount; il++ ){
	AliHLTTrackMCLabel &lab = dataPtr->fLabels[il];
	mcLabels[lab.fTrackID] = lab.fMCLabel;
      }
    } else {
      HLTWarning("data mismatch in block %s (0x%08x): count %d, size %d -> ignoring track MC information", 
		 DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, 
		 dataPtr->fCount, pBlock->fSize);
    }
  }
  
  
  
  for (int ndx=0; ndx<nBlocks && iResult>=0; ndx++) {

    const AliHLTComponentBlockData* iter = blocks+ndx;
    
    // Read TPC tracks
    
    if( iter->fDataType == ( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC ) ){	  
      fBenchmark.AddInput(iter->fSize);
      AliHLTTracksData* dataPtr = ( AliHLTTracksData* ) iter->fPtr;
      int nTracks = dataPtr->fCount;
      AliHLTExternalTrackParam* currOutTrack = dataPtr->fTracklets;
      for( int itr=0; itr<nTracks; itr++ ){
	AliHLTGlobalBarrelTrack t(*currOutTrack);
	Int_t mcLabel = -1;
	if( mcLabels.find(currOutTrack->fTrackID)!=mcLabels.end() )
	  mcLabel = mcLabels[currOutTrack->fTrackID];
	
	tracksTPC.push_back( t );
	tracksTPCLab.push_back(mcLabel);
	tracksTPCId.push_back( currOutTrack->fTrackID );
	unsigned int dSize = sizeof( AliHLTExternalTrackParam ) + currOutTrack->fNPoints * sizeof( unsigned int );
	currOutTrack = ( AliHLTExternalTrackParam* )( (( Byte_t * )currOutTrack) + dSize );
      }
    }


    // Read ITS clusters

    if ( (iter->fDataType == (kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD) ) || 
	 (iter->fDataType == (kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD) ) ||
	 (iter->fDataType == (kAliHLTDataTypeClusters|kAliHLTDataOriginITSSDD) ) ||
	 (iter->fDataType == (kAliHLTDataTypeClusters|kAliHLTDataOriginITS) ) 
	 ){
      
      fBenchmark.AddInput(iter->fSize);

      AliHLTITSClusterData *inPtr=reinterpret_cast<AliHLTITSClusterData*>( iter->fPtr );
      int nClusters = inPtr->fSpacePointCnt;
      for( int icl=0; icl<nClusters; icl++ ){
	AliHLTITSSpacePointData &d = inPtr->fSpacePoints[icl];

	Int_t lab[4] = { d.fTracks[0], d.fTracks[1], d.fTracks[2], d.fIndex };
	Int_t info[3] = { d.fNy, d.fNz, d.fLayer };
	Float_t hit[6] = { d.fY, d.fZ, d.fSigmaY2, d.fSigmaZ2, d.fQ, d.fSigmaYZ };
	if( d.fLayer==4 ) hit[5] = -hit[5];
	fTracker->LoadCluster( AliITSRecPoint( lab, hit, info ) );
	//cout<<"SG "<<d.fLayer<<" "<<d.fTracks[0]<<endl;
      }   
    }
    
  }// end read input blocks
  
  // Reconstruct the event

    fBenchmark.Start(1);
    fTracker->Reconstruct( &(tracksTPC[0]), &(tracksTPCLab[0]), tracksTPC.size() );
    fBenchmark.Stop(1);

  
  // Fill output tracks
  int nITSUpdated = 0;
  {
    
    for( int iOut=0; iOut<=1; iOut++ ){

      unsigned int blockSize = 0;

      AliHLTTracksData* outPtr = ( AliHLTTracksData* )( outputPtr + size );
      AliHLTExternalTrackParam* currOutTrack = outPtr->fTracklets;

      blockSize =   ( ( AliHLTUInt8_t * )currOutTrack ) -  ( ( AliHLTUInt8_t * )outPtr );

      if ( size + blockSize  > maxBufferSize ) {
	HLTWarning( "Output buffer size exceed (buffer size %d, current size %d), tracks are not stored", maxBufferSize, size + blockSize );
	iResult = -ENOSPC;
	break;
      }

      outPtr->fCount = 0;
       AliHLTITSTrack *tracks=0;
      int nTracks = 0;
      if( iOut==0 ){
	tracks = fTracker->Tracks();
	nTracks = fTracker->NTracks();
      } else{
	tracks = fTracker->ITSOutTracks();
	nTracks = fTracker->NITSOutTracks();
      }
      
      for ( int itr = 0; itr < nTracks; itr++ ) {
	AliHLTITSTrack &t = tracks[itr];
	int id =  tracksTPCId[t.TPCtrackId()];      
	int nClusters = t.GetNumberOfClusters();
	if( iOut==0 && nClusters>0 ) nITSUpdated++;
	
	unsigned int dSize = sizeof( AliHLTExternalTrackParam ) + nClusters * sizeof( unsigned int );
	
	if ( size + blockSize + dSize > maxBufferSize ) {
	  HLTWarning( "Output buffer size exceed (buffer size %d, current size %d), %d tracks are not stored", maxBufferSize, size + blockSize + dSize, nTracks - itr + 1 );
	  iResult = -ENOSPC;
	  break;
	}
	
	currOutTrack->fAlpha = t.GetAlpha();
	currOutTrack->fX = t.GetX();
	currOutTrack->fY = t.GetY();
	currOutTrack->fZ = t.GetZ();            
	currOutTrack->fLastX = 0;
	currOutTrack->fLastY = 0;
	currOutTrack->fLastZ = 0;      
	currOutTrack->fq1Pt = t.GetSigned1Pt();
	currOutTrack->fSinPsi = t.GetSnp();
	currOutTrack->fTgl = t.GetTgl();
	for( int i=0; i<15; i++ ) currOutTrack->fC[i] = t.GetCovariance()[i];
	currOutTrack->fTrackID = id;
	currOutTrack->fFlags = 0;
	currOutTrack->fNPoints = nClusters;    
	for ( int i = 0; i < nClusters; i++ ) currOutTrack->fPointIDs[i] = t.GetClusterIndex( i );
	currOutTrack = ( AliHLTExternalTrackParam* )( (( Byte_t * )currOutTrack) + dSize );
	blockSize += dSize;
	outPtr->fCount++;
      }
  

      AliHLTComponentBlockData resultData;
      FillBlockData( resultData );
      resultData.fOffset = size;
      resultData.fSize = blockSize;
      if( iOut==0 ){
	resultData.fDataType = kAliHLTDataTypeTrack|kAliHLTDataOriginITS;
      } else {
	resultData.fDataType = kAliHLTDataTypeTrack|kAliHLTDataOriginITSOut;
      }
      fBenchmark.AddOutput(resultData.fSize);
      outputBlocks.push_back( resultData );
      size += resultData.fSize;       
    }
  }  

  {// fill MC labels

    unsigned int blockSize = 0;
      
    AliHLTTrackMCData* outPtr = ( AliHLTTrackMCData* )( outputPtr +size );
    AliHLTTrackMCLabel* currOutLabel = outPtr->fLabels;
    
    blockSize =   ( ( AliHLTUInt8_t * )currOutLabel ) -  ( ( AliHLTUInt8_t * )outPtr );
    
    outPtr->fCount = 0;
    
    AliHLTITSTrack *tracks= fTracker->Tracks();
    int nTracks = fTracker->NTracks();
    
    for ( int itr = 0; itr < nTracks; itr++ ) {
      AliHLTITSTrack &t = tracks[itr];
      //cout<<"SG out:"<<tracksTPCId[t.TPCtrackId()]<<" "<<t.GetLabel()<<endl;
      if( t.GetLabel()<0 ) continue;
      int id =  tracksTPCId[t.TPCtrackId()];
      
      if ( blockSize + sizeof(AliHLTTrackMCLabel) > maxBufferSize ) {
	HLTWarning( "Output buffer size exceed (buffer size %d, current size %d), %d mc labels are not stored", maxBufferSize, blockSize, nTracks - itr + 1 );
	iResult = -ENOSPC;
	break;
      }
      currOutLabel->fTrackID = id;
      currOutLabel->fMCLabel = t.GetLabel();
      blockSize += sizeof(AliHLTTrackMCLabel);
      currOutLabel++;
      outPtr->fCount++;
    }        
    if( iResult>=0 && outPtr->fCount>0 ){
      AliHLTComponentBlockData resultData;
      FillBlockData( resultData );
      resultData.fOffset = size;
      resultData.fSize = blockSize;
      resultData.fDataType = kAliHLTDataTypeTrackMC|kAliHLTDataOriginITS;
      outputBlocks.push_back( resultData );
      size+= resultData.fSize;
    }
  }
  
  fBenchmark.Stop(0);

  // Set log level to "Warning" for on-line system monitoring
  HLTInfo( "ITS Tracker: output %d tracks;  input %d clusters, %d tracks",
	   nITSUpdated, nClustersTotal, tracksTPC.size() );

  HLTInfo(fBenchmark.GetStatistics());
  return iResult;
}
