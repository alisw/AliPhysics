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

///  @file   AliHLTITSSAPTrackerComponent.h
///  @author Ruben Shahoyan <ruben.shahoyan@cern.ch>
///  @date   August 2014
///  @brief  An ITS standalone primaries tracker/vertexer processing component for the HLT
///  Adapted from HLT/ITS/tracking/AliHLTITSSAPTrackerComponent component

/////////////////////////////////////////////////////
//                                                 //
// a ITS tracker processing component for the HLT  //
//                                                 //
/////////////////////////////////////////////////////

#include "AliHLTITSSAPTrackerComponent.h"
#include "AliHLTArray.h"
#include "AliExternalTrackParam.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "TObjString.h"
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include "AliITSSAPTracker.h"
#include "AliHLTITSSpacePointData.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliHLTDataTypes.h"
#include "AliHLTExternalTrackParam.h"
#include "AliGeomManager.h"
#include "AliHLTTrackMCLabel.h"
#include "AliITSRecPoint.h"
#include "AliHLTSAPTrackerData.h"
#include <map>

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp( AliHLTITSSAPTrackerComponent )
AliHLTITSSAPTrackerComponent::AliHLTITSSAPTrackerComponent()
: fSolenoidBz( 0 ),
  fBenchmark("ITSSAPTracker"),
  fTracker(0),
  fClusters(0)

{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTITSSAPTrackerComponent::AliHLTITSSAPTrackerComponent( const AliHLTITSSAPTrackerComponent& )
  :AliHLTProcessor(),
   fSolenoidBz( 0 ),
   fBenchmark("ITSSAPTracker"),
   fTracker(0),
   fClusters(0)
{
  // see header file for class documentation
  HLTFatal( "copy constructor untested" );
}

AliHLTITSSAPTrackerComponent& AliHLTITSSAPTrackerComponent::operator=( const AliHLTITSSAPTrackerComponent& )
{
  // see header file for class documentation
  HLTFatal( "assignment operator untested" );
  return *this;
}

AliHLTITSSAPTrackerComponent::~AliHLTITSSAPTrackerComponent()
{
  // see header file for class documentation
  delete fTracker;
  delete fClusters;
}

//
// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process
//

const char* AliHLTITSSAPTrackerComponent::GetComponentID()
{
  // see header file for class documentation
  return "ITSSAPTracker";
}

void AliHLTITSSAPTrackerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list )
{
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS );
  list.push_back( kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD );
  list.push_back( kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD );
  list.push_back( kAliHLTDataTypeClusters|kAliHLTDataOriginITSSDD );
  list.push_back( kAliHLTDataTypeClusters|kAliHLTDataOriginITS );
}

AliHLTComponentDataType AliHLTITSSAPTrackerComponent::GetOutputDataType()
{
  // see header file for class documentation  
  return kAliHLTMultipleDataType;
}

int AliHLTITSSAPTrackerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // see header file for class documentation  
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeTrack|kAliHLTDataOriginITS);
  tgtList.push_back(kAliHLTDataTypeTrack|kAliHLTDataOriginITSOut);
  tgtList.push_back(kAliHLTDataTypeTrackMC|kAliHLTDataOriginITS );
  tgtList.push_back(kAliHLTDataTypeESDVertex|kAliHLTDataOriginOut ); // RS??: is this correct?
  return tgtList.size();
}

void AliHLTITSSAPTrackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // define guess for the output data size
  constBase = 200;       // minimum size
  inputMultiplier = 2.; // size relative to input
}

AliHLTComponent* AliHLTITSSAPTrackerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTITSSAPTrackerComponent;
}

void AliHLTITSSAPTrackerComponent::SetDefaultConfiguration()
{
  // Set default configuration for the CA tracker component
  // Some parameters can be later overwritten from the OCDB

  fSolenoidBz = -5.00668;
  
}

int AliHLTITSSAPTrackerComponent::ReadConfigurationString(  const char* arguments )
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


int AliHLTITSSAPTrackerComponent::ReadCDBEntry( const char* cdbEntry, const char* chainId )
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


int AliHLTITSSAPTrackerComponent::Configure( const char* cdbEntry, const char* chainId, const char *commandLine )
{
  // Configure the component
  // There are few levels of configuration,
  // parameters which are set on one step can be overwritten on the next step

  //* read hard-coded values

  SetDefaultConfiguration();

  //* read the default CDB entry

  int iResult1 = ReadCDBEntry( NULL, chainId );

  //* read magnetic field

  fSolenoidBz = GetBz();

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



int AliHLTITSSAPTrackerComponent::DoInit( int argc, const char** argv )
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

  fTracker = new AliITSSAPTracker();
  fTracker->Init();
  fBenchmark.Reset();
  fBenchmark.SetTimer(0,"total");
  fBenchmark.SetTimer(1,"reco");
  return ret;
}


int AliHLTITSSAPTrackerComponent::DoDeinit()
{
  // see header file for class documentation
  delete fTracker;
  fTracker = 0;
  return 0;
}



int AliHLTITSSAPTrackerComponent::Reconfigure( const char* cdbEntry, const char* chainId )
{
  // Reconfigure the component from OCDB .

  return Configure( cdbEntry, chainId, NULL );
}



int AliHLTITSSAPTrackerComponent::DoEvent
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


  // Check if there is an SPD vertex
  const AliESDVertex *vertexSPD = 0;

  {
    const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITSSPD);
    if( iter != NULL  ) {
      if( !( vertexSPD = dynamic_cast<AliESDVertex*>(const_cast<TObject*>( iter ) ) ) ){    
	HLTError("ITS SPD vertex object is corrupted");
	iResult = -EINVAL;    
      }
    }
    else {
      HLTWarning("No SPD vertex, skip ITS standalone reconstruction");
      return 0;
    }
  }  

  int nBlocks = evtData.fBlockCnt;
  if (!fClusters) fClusters = new TClonesArray("AliITSRecPoint",1000);
  //
  int nClTotal = 0;
  for (int ndx=0; ndx<nBlocks && iResult>=0; ndx++) {

    const AliHLTComponentBlockData* iter = blocks+ndx;
    
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
	// tracker does not out the clusters, add them to transient array
	fTracker->AddCluster( new((*fClusters)[nClTotal++]) AliITSRecPoint(lab,hit,info) );
      }   
    }
    
  }// end read input blocks
  
  // Reconstruct the event
  
  fBenchmark.Start(1);
  fTracker->SetSPDVertex(vertexSPD);
  fTracker->ProcessEvent();
  fBenchmark.Stop(1);

  
  // Fill output tracks
  int nAddedTracks = 0;
  {  
    int nFoundTracks = fTracker->GetNTracks();
    AliHLTUInt32_t blockSize = sizeof(AliHLTITSSAPTrackerDataContainer) + nFoundTracks*sizeof(AliHLTITSSAPTrackerData);
    if( size + blockSize > maxBufferSize ){    	
      HLTWarning( "Output buffer size exceed (buffer size %d, current size %d), %d tracks are not stored", 
		  maxBufferSize, size + blockSize, nFoundTracks);
      iResult = -ENOSPC;
    }    
    if( iResult>=0 ){
      blockSize = sizeof(AliHLTITSSAPTrackerDataContainer);
      AliHLTITSSAPTrackerDataContainer *data = reinterpret_cast<AliHLTITSSAPTrackerDataContainer*>(outputPtr);
      data->fCount=0;
      for (int itr=0;itr<nFoundTracks;itr++) {
	const AliITSSAPTracker::ITStrack_t& track = fTracker->GetTrack(itr);
	// the track is just a struct of 2 AliExternalTrackParams (params at vertex and at the outer ITS layer)
	// + some extra info, see "struct ITStrack" in the AliITSSAPTracker.h
	if ( track.paramOut.TestBit(AliITSSAPTracker::kInvalidBit) || 
	     track.paramInw.TestBit(AliITSSAPTracker::kInvalidBit)) continue;
	// use only those tracks whose both inward and outward params are OK.
	AliHLTITSSAPTrackerData &trcHLT = data->fTracks[data->fCount];
	trcHLT.paramOut.SetExternalTrackParam(&track.paramOut);
	trcHLT.paramInw.SetExternalTrackParam(&track.paramInw);
	trcHLT.chi2 = track.chi2;
	trcHLT.ncl  = track.ncl;
	trcHLT.label = track.label;
	data->fCount++;
	blockSize += sizeof(AliHLTITSSAPTrackerData);
	nAddedTracks++;
      }
      
      AliHLTComponentBlockData resultData;
      FillBlockData( resultData );
      resultData.fOffset = size;
      resultData.fSize = blockSize;      
      resultData.fDataType = kAliHLTDataTypeITSSAPData|kAliHLTDataOriginITS;
      fBenchmark.AddOutput(resultData.fSize);
      outputBlocks.push_back( resultData );
      size += resultData.fSize;
    }
  }

  Bool_t vtxOK = kFALSE;
  { // Fill output vertexTracks  
    AliESDVertex& vtxTracks = fTracker->GetTrackVertex();
    if (vtxTracks.GetStatus()==1) {
      PushBack( (TObject*) &vtxTracks, kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS,0 );
      vtxOK = kTRUE;
    }
  }
  //
  fTracker->Clear();
  fClusters->Clear();
  //  
  fBenchmark.Stop(0);

  // Set log level to "Warning" for on-line system monitoring
  HLTInfo( "ITS SAP Tracker: output %d tracks;  input %d clusters, VertexTracks: %s",
	   nAddedTracks, nClTotal, vtxOK ? "OK" : "Found" );

  HLTInfo(fBenchmark.GetStatistics());
  return iResult;
}
