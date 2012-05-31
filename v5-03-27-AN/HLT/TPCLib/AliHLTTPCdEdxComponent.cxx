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

///  @file   AliHLTTPCdEdxComponent.cxx
///  @author Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de>
///  @date   June 2009
///  @brief  An ITS tracker processing component for the HLT


/////////////////////////////////////////////////////
//                                                 //
// dEdx calculation component for the HLT TPC       //
//                                                 //
/////////////////////////////////////////////////////

#include "AliHLTTPCdEdxComponent.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "AliHLTDataTypes.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCDefinitions.h"
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCRecoParam.h"
#include "AliTPCTransform.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp( AliHLTTPCdEdxComponent )
AliHLTTPCdEdxComponent::AliHLTTPCdEdxComponent()
    :
    fSolenoidBz( 0 ),
    fStatTime( 0 ),
    fStatNEvents( 0 )
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  for( int i=0; i<fkNPatches; i++ ){
    fPatchClusters[i] = 0;    
    fNPatchClusters[i] = 0;    
  }
}

AliHLTTPCdEdxComponent::AliHLTTPCdEdxComponent( const AliHLTTPCdEdxComponent& )
    :
    AliHLTProcessor(),
    fSolenoidBz( 0 ),
    fStatTime( 0 ),
    fStatNEvents( 0 )
{
  // dummy
  for( int i=0; i<fkNPatches; i++ ){
    fPatchClusters[i] = 0;    
    fNPatchClusters[i] = 0;    
  }
  HLTFatal( "copy constructor untested" );
}

AliHLTTPCdEdxComponent& AliHLTTPCdEdxComponent::operator=( const AliHLTTPCdEdxComponent& )
{
  // see header file for class documentation
  HLTFatal( "assignment operator untested" );
  for( int i=0; i<fkNPatches; i++ ){
    fPatchClusters[i] = 0;    
    fNPatchClusters[i] = 0;    
  }
  return *this;
}

AliHLTTPCdEdxComponent::~AliHLTTPCdEdxComponent()
{
  // see header file for class documentation
  for( int i=0; i<fkNPatches; i++ ){
    delete[] fPatchClusters[i];
  }
}

//
// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process
//

const char* AliHLTTPCdEdxComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCdEdx";
}

void AliHLTTPCdEdxComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list )
{
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC );
  list.push_back( AliHLTTPCDefinitions::fgkClustersDataType );
}

AliHLTComponentDataType AliHLTTPCdEdxComponent::GetOutputDataType()
{
  // see header file for class documentation  
  return kAliHLTDataTypedEdx|kAliHLTDataOriginTPC;
}

void AliHLTTPCdEdxComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // define guess for the output data size
  constBase = 200;       // minimum size
  inputMultiplier = 0.1; // size relative to input
}

AliHLTComponent* AliHLTTPCdEdxComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCdEdxComponent;
}

void AliHLTTPCdEdxComponent::SetDefaultConfiguration()
{
  // Set default configuration for the CA tracker component
  // Some parameters can be later overwritten from the OCDB

  fSolenoidBz = -5.00668;
  fStatTime = 0;
  fStatNEvents = 0;
}

int AliHLTTPCdEdxComponent::ReadConfigurationString(  const char* arguments )
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


int AliHLTTPCdEdxComponent::ReadCDBEntry( const char* cdbEntry, const char* chainId )
{
  // see header file for class documentation

  const char* defaultNotify = "";

  if ( !cdbEntry ) {
    return 0;// need to add the HLT/ConfigTPC/TPCdEdx directory to cdb SG!!!
    //cdbEntry = "HLT/ConfigTPC/TPCdEdx";
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


int AliHLTTPCdEdxComponent::Configure( const char* cdbEntry, const char* chainId, const char *commandLine )
{
  // Configure the component
  // There are few levels of configuration,
  // parameters which are set on one step can be overwritten on the next step

  //* read hard-coded values

  SetDefaultConfiguration();

  //* read the default CDB entry

  int iResult1 = ReadCDBEntry( NULL, chainId );

  //* read magnetic field

  fSolenoidBz=GetBz();

  //* read the actual CDB entry if required

  int iResult2 = ( cdbEntry ) ? ReadCDBEntry( cdbEntry, chainId ) : 0;

  //* read extra parameters from input (if they are)

  int iResult3 = 0;

  if ( commandLine && commandLine[0] != '\0' ) {
    HLTInfo( "received configuration string from HLT framework: \"%s\"", commandLine );
    iResult3 = ReadConfigurationString( commandLine );
  }

  // Initialise the tracker here

  return iResult1 ? iResult1 : ( iResult2 ? iResult2 : iResult3 );
}



int AliHLTTPCdEdxComponent::DoInit( int argc, const char** argv )
{
  // Configure the component

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

  AliTPCcalibDB::Instance()->SetRun(GetRunNo());
  AliTPCcalibDB::Instance()->UpdateRunInformations(GetRunNo());

  AliTPCcalibDB::Instance()->GetTransform()->SetCurrentRun(GetRunNo());
  AliTPCcalibDB::Instance()->GetTransform()->SetCurrentTimeStamp( GetTimeStamp() );
  AliTPCcalibDB::Instance()->GetTransform()->SetCurrentRecoParam(AliTPCRecoParam::GetHLTParam());

  //AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform();
  //if( transform ){
  //AliTPCRecoParam *reco = transform->GetCurrentRecoParam();
  //if( ! reco ){
  //reco = new AliTPCRecoParam;
  //reco->SetUseTotalCharge(0);
  //transform->SetCurrentRecoParam( reco );
  //}

  return ret;
}


int AliHLTTPCdEdxComponent::DoDeinit()
{
  // see header file for class documentation
  return 0;
}



int AliHLTTPCdEdxComponent::Reconfigure( const char* cdbEntry, const char* chainId )
{
  // Reconfigure the component from OCDB .

  return Configure( cdbEntry, chainId, NULL );
}



int AliHLTTPCdEdxComponent::DoEvent
(
  const AliHLTComponentEventData& evtData,
  const AliHLTComponentBlockData* blocks,
  AliHLTComponentTriggerData& /*trigData*/,
  AliHLTUInt8_t* outputPtr,
  AliHLTUInt32_t& size,
  vector<AliHLTComponentBlockData>& outputBlocks )
{
  //* process the event

  AliHLTUInt32_t maxBufferSize = size;
  size = 0; // output size

  if (!IsDataEvent()) return 0;

  if ( evtData.fBlockCnt <= 0 ) {
    HLTWarning( "no blocks in event" );
    return 0;
  }
  if ( !outputPtr ) {
    return -ENOSPC;
  }

  int iResult=0;

  TStopwatch timer;

  // Initialise the data pointers

  for( int i=0; i<fkNPatches; i++ ){
    delete[] fPatchClusters[i];    
    fPatchClusters[i] = 0;
    fNPatchClusters[i] = 0;    
  }

  int nBlocks = (int)evtData.fBlockCnt;

  int nInputClusters = 0;
  int nInputTracks = 0;

  // first read all the clusters

  for (int ndx=0; ndx<nBlocks && iResult>=0; ndx++) {
    const AliHLTComponentBlockData* iter = blocks+ndx;
    if ( iter->fDataType != AliHLTTPCDefinitions::fgkClustersDataType ) continue;
    Int_t slice=AliHLTTPCDefinitions::GetMinSliceNr(iter->fSpecification);
    Int_t patch=AliHLTTPCDefinitions::GetMinPatchNr(iter->fSpecification);
    Int_t slicepatch=slice*6+patch;
    if( slicepatch >= fkNPatches ){
      HLTWarning("Wrong header of TPC cluster data, slice %d, patch %d",
		 slice, patch );
      continue;
    }
    AliHLTTPCClusterData* inPtrSP = ( AliHLTTPCClusterData* )( iter->fPtr );
    nInputClusters += inPtrSP->fSpacePointCnt;

    delete[] fPatchClusters[slicepatch];
    fPatchClusters[slicepatch] = new AliTPCclusterMI[inPtrSP->fSpacePointCnt];
    fNPatchClusters[slicepatch] = inPtrSP->fSpacePointCnt;
    
    // create  off-line clusters out of the HLT clusters
    // todo: check which cluster information is really needed for the dEdx

    for ( unsigned int i = 0; i < inPtrSP->fSpacePointCnt; i++ ) {
      AliHLTTPCSpacePointData *chlt = &( inPtrSP->fSpacePoints[i] );
      AliTPCclusterMI *c = fPatchClusters[slicepatch]+i;
      c->SetX(chlt->fX);
      c->SetY(chlt->fY);
      c->SetZ(chlt->fZ);
      c->SetSigmaY2(chlt->fSigmaY2);
      c->SetSigmaYZ( 0 );
      c->SetSigmaZ2(chlt->fSigmaZ2);
      c->SetQ( chlt->fCharge );
      c->SetMax( chlt->fQMax );
      Int_t sector, row;
      Float_t padtime[3]={0,chlt->fY,chlt->fZ};
      AliHLTTPCTransform::Slice2Sector(slice,chlt->fPadRow, sector, row);
      AliHLTTPCTransform::Local2Raw( padtime, sector, row);
      c->SetDetector( sector );
      c->SetRow( row );
      c->SetPad( (Int_t) padtime[1] );
      c->SetTimeBin( (Int_t) padtime[2] );
    }
  }


  // loop over the input tracks: calculate dEdx and write output

  unsigned int outSize = 0;
  AliHLTFloat32_t *outPtr = ( AliHLTFloat32_t * )( outputPtr );

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    
    AliHLTComponentBlockData outBlock;
    FillBlockData( outBlock );
    outBlock.fOffset = outSize;
    outBlock.fSize = 0;
    outBlock.fDataType = kAliHLTDataTypedEdx|kAliHLTDataOriginTPC;
    outBlock.fSpecification = pBlock->fSpecification;
  
    AliHLTTracksData* dataPtr = ( AliHLTTracksData* ) pBlock->fPtr;
    int nTracks = dataPtr->fCount;
    AliHLTExternalTrackParam* currTrack = dataPtr->fTracklets;
    nInputTracks+=nTracks;

    for( int itr=0; 
	 itr<nTracks && ( (AliHLTUInt8_t *)currTrack < ((AliHLTUInt8_t *) pBlock->fPtr)+pBlock->fSize); 
	 itr++ ){

      // create an off-line track
      AliHLTGlobalBarrelTrack gb(*currTrack);
      AliTPCseed tTPC;
      tTPC.Set( gb.GetX(), gb.GetAlpha(),
		gb.GetParameter(), gb.GetCovariance() );

      // set the clusters
      
      for( UInt_t ic=0; ic<currTrack->fNPoints; ic++){	    
	UInt_t id = currTrack->fPointIDs[ic];
	int iSlice = AliHLTTPCSpacePointData::GetSlice(id);
	int iPatch = AliHLTTPCSpacePointData::GetPatch(id);
	int iCluster = AliHLTTPCSpacePointData::GetNumber(id);
	if( iSlice<0 || iSlice>36 || iPatch<0 || iPatch>5 ){
	  HLTError("Corrupted TPC cluster Id: slice %d, patch %d, cluster %d",
		   iSlice, iPatch,iCluster );
	  continue;
	}
	AliTPCclusterMI *patchClusters = fPatchClusters[iSlice*6 + iPatch];
	if( !patchClusters ){
	  HLTError("Clusters are missed for slice %d, patch %d",
		   iSlice, iPatch );
	  continue;
	}
	if( iCluster >= fNPatchClusters[iSlice*6 + iPatch] ){
	  HLTError("TPC slice %d, patch %d: ClusterID==%d >= N Cluaters==%d ",
		   iSlice, iPatch,iCluster, fNPatchClusters[iSlice*6 + iPatch] );
	  continue;
	}
	AliTPCclusterMI *c = &(patchClusters[iCluster]);	  
	
	int sec = c->GetDetector();
	int row = c->GetRow();
	if ( sec >= 36 ) {
	  row = row + AliHLTTPCTransform::GetNRowLow();
	}
	tTPC.SetClusterPointer( row, c );
	
	AliTPCTrackerPoint &point = *( tTPC.GetTrackPoint( row ) );
	tTPC.Propagate( TMath::DegToRad()*(sec%18*20.+10.), c->GetX(), fSolenoidBz );
	Double_t angle2 = tTPC.GetSnp()*tTPC.GetSnp();
	angle2 = (angle2<1) ?TMath::Sqrt(angle2/(1-angle2)) :10.; 
	point.SetAngleY( angle2 );
	point.SetAngleZ( tTPC.GetTgl() );
      }

      // Cook dEdx

      if( outSize+3*sizeof( AliHLTFloat32_t ) > maxBufferSize ){
        HLTWarning( "Output buffer size %d exceed", maxBufferSize );
        iResult = -ENOSPC;
        break;
      }

      //Old method
      /*
      tTPC.CookdEdx( 0.02, 0.6 );      
      outPtr[0] = tTPC.GetdEdx();
      outPtr[1] = tTPC.GetSDEDX(0);
      outPtr[2] = tTPC.GetNCDEDX(0);
      */

      //New method
      outPtr[0] = tTPC.CookdEdxAnalytical(0.02,0.6,1,0,159,0);
      outPtr[1] = 0.;
      outPtr[2] = 159;
      
      outPtr+=3;
      outSize+=3*sizeof( AliHLTFloat32_t );    
      outBlock.fSize+=3*sizeof( AliHLTFloat32_t );  

      unsigned int step = sizeof( AliHLTExternalTrackParam ) + currTrack->fNPoints * sizeof( unsigned int );
      currTrack = ( AliHLTExternalTrackParam* )( (( Byte_t * )currTrack) + step );  
    }
    if( iResult<0 ) break;
    outputBlocks.push_back( outBlock );    
    size = outSize;
  }

  timer.Stop();
  fStatTime += timer.RealTime();
  fStatNEvents++;
  
  // Set log level to "Warning" for on-line system monitoring
  int hz = ( int ) ( fStatTime > 1.e-10 ? fStatNEvents / fStatTime : 100000 );

  HLTInfo( "TPCdEdx : %d clusters and %d tracks processed; average time %d Hz",
	   nInputClusters, nInputTracks, hz );
  
  return iResult;
}
