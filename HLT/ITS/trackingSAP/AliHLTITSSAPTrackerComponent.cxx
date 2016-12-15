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
#include <TObjString.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include "AliITSSAPTracker.h"
#include "AliITSSAPLayer.h"
#include "AliHLTITSSpacePointData.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliHLTDataTypes.h"
#include "AliHLTExternalTrackParam.h"
#include "AliGeomManager.h"
#include "AliHLTTrackMCLabel.h"
#include "AliITSRecPoint.h"
#include "AliITSRecoParam.h"
#include "AliHLTSAPTrackerData.h"
#include "AliHLTMessage.h"
#include "AliFlatESDVertex.h"
#include "AliHLTReadoutList.h"
#include "AliHLTCTPData.h"
#include "AliHLTITSTrackPoint.h"
#include <map>

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp( AliHLTITSSAPTrackerComponent )
AliHLTITSSAPTrackerComponent::AliHLTITSSAPTrackerComponent()
: fRecoParamType(AliRecoParam::kDefault),
  fSkipSDD(-1),
  fMaxMissL(1),
  fMaxTrackletsToRun(-1),
  fMaxVtxIter(-1),
  fStopScaleChange(-1),
  fMaxRSPDVtx(-1),
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
   fRecoParamType(AliRecoParam::kDefault),
   fSkipSDD(-1),
   fMaxMissL(1),
   fMaxTrackletsToRun(-1),
   fMaxVtxIter(-1),
   fStopScaleChange(-1),
   fMaxRSPDVtx(-1),
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
  list.push_back( kAliHLTDataTypeESDVertex|kAliHLTDataOriginITSSPD );
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
  tgtList.push_back(kAliHLTDataTypeITSSAPData|kAliHLTDataOriginITS);
  tgtList.push_back(kAliHLTDataTypeFlatESDVertex|kAliHLTDataOriginITS ); // RS??: is this correct?
  tgtList.push_back(kAliHLTDataTypeITSTrackPoint|kAliHLTDataOriginITS);
  return tgtList.size();
}

void AliHLTITSSAPTrackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // define guess for the output data size
  constBase = 200;       // minimum size
  inputMultiplier = 3.; // size relative to input
}

AliHLTComponent* AliHLTITSSAPTrackerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTITSSAPTrackerComponent;
}

void AliHLTITSSAPTrackerComponent::SetDefaultConfiguration()
{
  // Set default configuration for the ITSSAP tracker component
  // Some parameters can be later overwritten from the OCDB
    
}

int AliHLTITSSAPTrackerComponent::ReadConfigurationString(  const char* arguments )
{
  // Set configuration parameters for the CA tracker component from the string

  int iResult = 0;
  TString allArgs = arguments;
  if (allArgs.IsNull()) return iResult;

  TString argument;
  int bMissingParam = 0;

  TObjArray* pTokens = allArgs.Tokenize( " " );

  int nArgs =  pTokens ? pTokens->GetEntries() : 0;

  for ( int i = 0; i < nArgs; i++ ) {
    argument = ( ( TObjString* )pTokens->At( i ) )->GetString();
    argument.ToLower();
    if ( argument.IsNull() ) continue;
    
    if (argument.CompareTo("-lowflux")==0) {
      fRecoParamType = AliRecoParam::kLowMult;
      HLTInfo("Low flux reconstruction selected");
      continue;
    }
    if (argument.CompareTo("-highflux")==0) {
      fRecoParamType = AliRecoParam::kHighMult;
      HLTInfo("High flux reconstruction selected");
      continue;
    }
    if (argument.CompareTo("-cosmics")==0 ||
	     argument.CompareTo("-calib")==0) {
      HLTWarning("%s reconstruction selected: override to default",argument.Data());
      continue;
    }    
    //
    if (argument.CompareTo("-skipsdd")==0) {
      fSkipSDD = 1;
      HLTInfo("SDD will be ignored");
      continue;
    }    
    //
    if (argument.CompareTo("-maxmisslayers")==0) {
      if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
      fMaxMissL = ((TObjString*)pTokens->At(i))->GetString().Atoi();
      HLTInfo("Allow max active layers missed: %d", fMaxMissL);
      continue;
    } 
    //
    if (argument.CompareTo("-maxmult")==0) {
      if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
      fMaxTrackletsToRun = ((TObjString*)pTokens->At(i))->GetString().Atoi();
      HLTInfo("Skip tracking if N SPD tracklets > %d", fMaxTrackletsToRun);
      continue;
    } 
    //
    if (argument.CompareTo("-maxitervtx")==0) {
      if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
      fMaxVtxIter = ((TObjString*)pTokens->At(i))->GetString().Atoi();
      if (fMaxVtxIter<1) {
	HLTError("Incorrect maxitervtx %d supplied, ITSSAPTracker default will be used", fMaxVtxIter);
	fMaxVtxIter = -1;
      }
      else HLTInfo("Allow max %d iterations for vertexer", fMaxVtxIter);
      continue;
    } 
    //
    if (argument.CompareTo("-stopscalevtx")==0) {
      if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
      fStopScaleChange = ((TObjString*)pTokens->At(i))->GetString().Atof();
      if (fStopScaleChange<0.3) {
	HLTError("Incorrect stopscalevtx %.2f supplied, ITSSAPTracker default will be used", fStopScaleChange);
	fStopScaleChange = -1;
      }
      else HLTInfo("Stop vertexing on scale change below %.2f", fStopScaleChange);
      continue;
    } 
    //
    if (argument.CompareTo("-maxrspdv")==0) {
      if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
      fMaxRSPDVtx = ((TObjString*)pTokens->At(i))->GetString().Atof();
      HLTInfo("Accept SPD vertices with R<",fMaxRSPDVtx);
      continue;
    }
    //
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
  HLTInfo("cdbEnty: %s chaindId: %s commandLine: %s",cdbEntry,chainId,commandLine);
  //* read hard-coded values

  SetDefaultConfiguration();

  //* read the default CDB entry

  int iResult1 = ReadCDBEntry( NULL, chainId );

  //* read the actual CDB entry if required

  int iResult2 = ( cdbEntry ) ? ReadCDBEntry( cdbEntry, chainId ) : 0;

  //* read extra parameters from input (if they are)

  int iResult3 = 0;

  if ( commandLine && commandLine[0] != '\0' ) {
    //    HLTInfo( "received configuration string from HLT framework: \"%s\"", commandLine );
    iResult3 = ReadConfigurationString( commandLine );
  }

  // Initialise the tracker here

  return iResult1 ? iResult1 : ( iResult2 ? iResult2 : iResult3 );
}



int AliHLTITSSAPTrackerComponent::DoInit( int argc, const char** argv )
{
  // Configure the ITS tracker component

  if ( fTracker ) return -EINPROGRESS;

  SetupCTPData();

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

  fTracker = new AliITSSAPTracker();
  fTracker->SetBz(GetBz());
  fTracker->Init();  // init defaults
  //
  // check consistency of options (if provided)
  TObjArray* pArr = dynamic_cast<TObjArray*>(LoadAndExtractOCDBObject("ITS/Calib/RecoParam"));
  AliITSRecoParam* param = 0;   // fetch relevant recoparam
  if (pArr) {
    int np = pArr->GetEntriesFast();
    for (int ip=np;ip--;) {
      if (!(param=(AliITSRecoParam*)pArr->At(ip)) ||
	  param->GetEventSpecie()!=fRecoParamType) {
	param = 0;
	continue;
      }
    }
  }
  //
  if (fSkipSDD<0) {
    if (param && (param->GetLayersToSkip(2)||param->GetLayersToSkip(3))) {
      fSkipSDD = 1;
      HLTInfo("Force to skip SDD layers (recoparam)");
    }
    else fSkipSDD = 0;
  }
  //
  if (fMaxTrackletsToRun<0) {
    if (param) fMaxTrackletsToRun = param->GetMaxSPDcontrForSAToUseAllClusters();
    else fMaxTrackletsToRun = 99999;
  }
  HLTInfo("Max N SPD tracklets to run tracking: %d",fMaxTrackletsToRun);
  //
  if (fMaxMissL<0 && param) {
    fMaxMissL = 6 - param->GetMinNPointsSA();
  }
  if (fMaxMissL>3) fMaxMissL = 3;
  HLTInfo("Allow to skip at most %d layers",fMaxMissL);
  //
  if (fMaxRSPDVtx<0) fMaxRSPDVtx = 1.5;
  else HLTInfo("Accept SPD vertices with R<",fMaxRSPDVtx);
  fTracker->SetMaxRSPDVtx(fMaxRSPDVtx);

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
  AliHLTComponentTriggerData& trigData,
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

  // we don't use SDD if not in trigger or explicitly forbidden
  const AliHLTCTPData* ctpdata=CTPData();
  AliHLTReadoutList rd=ctpdata->ReadoutList(trigData);
  Bool_t skipSDD = fSkipSDD || !rd.DetectorEnabled(AliHLTReadoutList::kITSSDD);
  fTracker->SetSkipLayer(AliITSSAPTracker::kALrSDD1,skipSDD);
  fTracker->SetSkipLayer(AliITSSAPTracker::kALrSDD2,skipSDD);

  int maxMiss = fMaxMissL;
  if (skipSDD && maxMiss>1) maxMiss = 1; // there should be at least 1 hits above SPD
  fTracker->SetMaxMissedLayers(maxMiss);

  if (fMaxVtxIter>0)      fTracker->SetMaxVtxIter(fMaxVtxIter);
  if (fStopScaleChange>0) fTracker->SetStopScaleChange(fStopScaleChange);

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
      HLTInfo("No SPD vertex, skip ITS standalone reconstruction");
      return 0;
    }
  }  
  
  if (vertexSPD->GetNContributors()>fMaxTrackletsToRun) {
    HLTInfo("Skip tracking: HLT SPD vertex has %d>%d tracklets",
	    vertexSPD->GetNContributors(),fMaxTrackletsToRun );
    return 0;
  }
  fTracker->SetMaxTrackletsToRunTracking(fMaxTrackletsToRun);

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
  int nTrackClusters = 0;
  int trackIDMap[fTracker->GetNTracks()];
  for( int i=0; i<fTracker->GetNTracks(); i++) trackIDMap[i] = -1;
  AliHLTITSSAPTrackerDataContainer *outTrackData = reinterpret_cast<AliHLTITSSAPTrackerDataContainer*>(outputPtr);
  {  
    int nFoundTracks = fTracker->GetNTracks();
    AliHLTUInt32_t blockSize = sizeof(AliHLTITSSAPTrackerDataContainer) + nFoundTracks*sizeof(AliHLTITSSAPTrackerData);
    if( size + blockSize + AliFlatESDVertex::GetSize() > maxBufferSize ){    	
      HLTWarning( "Output buffer size exceed (buffer size %d, current size %d), %d tracks are not stored", 
		  maxBufferSize, size + blockSize, nFoundTracks);
      iResult = -ENOSPC;
    }    
    if( iResult>=0 ){
      blockSize = sizeof(AliHLTITSSAPTrackerDataContainer);
      outTrackData->fCount=0;
      outTrackData->fNSPDtracklets = fTracker->GetNTracklets();
      for (int i=AliITSSAPTracker::kNLrActive;i--;) outTrackData->fNclusters[i] = fTracker->GetLayer(i)->GetNClusters();
      for (int itr=0;itr<nFoundTracks;itr++) {
	const AliITSSAPTracker::ITStrack_t& track = fTracker->GetTrack(itr);
	// the track is just a struct of 2 AliExternalTrackParams (params at vertex and at the outer ITS layer)
	// + some extra info, see "struct ITStrack" in the AliITSSAPTracker.h
	if ( track.paramOut.TestBit(AliITSSAPTracker::kInvalidBit) || 
	     track.paramInw.TestBit(AliITSSAPTracker::kInvalidBit)) continue;
	// use only those tracks whose both inward and outward params are OK.
	AliHLTITSSAPTrackerData &trcHLT = outTrackData->fTracks[outTrackData->fCount];
	trcHLT.paramOut.SetExternalTrackParam(&track.paramOut);
	trcHLT.paramInw.SetExternalTrackParam(&track.paramInw);
	trcHLT.chi2 = track.chi2;
	trcHLT.ncl  = track.ncl;
	trcHLT.label = track.label;
	trackIDMap[outTrackData->fCount] = itr;
	outTrackData->fCount++;
	blockSize += sizeof(AliHLTITSSAPTrackerData);
	nAddedTracks++;
	nTrackClusters+=trcHLT.ncl;
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
    if ( iResult>=0 && vtxTracks.GetStatus()==1 ) {
      AliFlatESDVertex *flatVtx = reinterpret_cast<AliFlatESDVertex*>( outputPtr + size );
      flatVtx->SetFromESDVertex( vtxTracks );
      AliHLTComponentBlockData resultData;
      FillBlockData( resultData );
      resultData.fOffset = size;
      resultData.fSize = flatVtx->GetSize();      
      resultData.fDataType = kAliHLTDataTypeFlatESDVertex|kAliHLTDataOriginITS;
      fBenchmark.AddOutput(resultData.fSize);
      outputBlocks.push_back( resultData );
      size += resultData.fSize;
      vtxOK = kTRUE;
    }
  }

  if ( iResult>=0 ) do{ // Fill output track points

    AliHLTITSTrackPointData* outTrackPoints = reinterpret_cast<AliHLTITSTrackPointData*>( outputPtr + size );
    AliHLTUInt32_t blockSize = sizeof(AliHLTITSTrackPointData) + nTrackClusters*sizeof(AliHLTITSTrackPoint);
    if( size + blockSize > maxBufferSize ){    	
      HLTWarning( "Output buffer size exceed (buffer size %d, current size %d), %d track points are not stored", 
		  maxBufferSize, size + blockSize, nTrackClusters);
      iResult = -ENOSPC;
      break;
    }    
    outTrackPoints->fCount = 0; 
    blockSize = sizeof(AliHLTITSTrackPointData);

    for( int itr=0; itr<outTrackData->fCount; itr++){      
      AliHLTITSSAPTrackerData &trcHLT = outTrackData->fTracks[itr];
      const AliITSSAPTracker::ITStrack_t& track = fTracker->GetTrack(trackIDMap[itr]);
      int nCl=0;
      for( int iLayer=0; iLayer<6; iLayer++ ){
	int id = track.clID[iLayer];
	if( id<0 ) continue;
	nCl++;
	if( nCl > trcHLT.ncl ) break;
	if( fTracker->GetTrackPoint( iLayer, id, outTrackPoints->fPoints[outTrackPoints->fCount] )!=0 ){
	  HLTError( "wrong cluster pointer" );
	  outTrackPoints->fPoints[outTrackPoints->fCount].Reset();	
	}
	outTrackPoints->fCount++;
      }
      if( nCl != trcHLT.ncl ){
	HLTError( "ITS SAP Tracker: wrong n clusters in output track: %d instead of %d", trcHLT.ncl, nCl );
      }
      trcHLT.ncl = nCl;
    }

    AliHLTComponentBlockData resultData;
    FillBlockData( resultData );
    resultData.fOffset = size;
    resultData.fSize = blockSize;
    resultData.fDataType = kAliHLTDataTypeITSTrackPoint |kAliHLTDataOriginITS;
    fBenchmark.AddOutput(resultData.fSize);
    outputBlocks.push_back( resultData );
    size += resultData.fSize;   
  } while(0);

  //
  fTracker->Clear();
  fClusters->Clear();
  //  
  fBenchmark.Stop(0);

  // Set log level to "Warning" for on-line system monitoring
  HLTInfo( "ITS SAP Tracker: output %d tracks;  input %d clusters, VertexTracks: %s",
	   nAddedTracks, nClTotal, vtxOK ? "OK" : "Not Found" );

  HLTInfo(fBenchmark.GetStatistics());
  return iResult;
}
