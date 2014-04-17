// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

//  @file   AliHLTGlobalFlatEsdConverterComponent.cxx
//  @author Matthias Richter
//  @date   
//  @brief  Global ESD converter component.
// 

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include "AliHLTGlobalFlatEsdConverterComponent.h"
#include "AliFlatESDEvent.h"
#include "AliFlatESDTrack.h"
#include "AliFlatExternalTrackParam.h"
#include "AliExternalTrackParam.h"

#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTTrackMCLabel.h"
#include "AliHLTCTPData.h"
#include "AliHLTErrorGuard.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliPID.h"
#include "TTree.h"
#include "TList.h"
#include "TClonesArray.h"
//#include "AliHLTESDCaloClusterMaker.h"
//#include "AliHLTCaloClusterDataStruct.h"
//#include "AliHLTCaloClusterReader.h"
//#include "AliESDCaloCluster.h"
//#include "AliESDVZERO.h"
#include "AliHLTGlobalVertexerComponent.h"
#include "AliHLTVertexFinderBase.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCClusterMCData.h"
#include "AliHLTTPCTransform.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalFlatEsdConverterComponent)

AliHLTGlobalFlatEsdConverterComponent::AliHLTGlobalFlatEsdConverterComponent()
  : AliHLTProcessor()
  , fWriteClusters(0)
  , fVerbosity(0)  
  , fSolenoidBz(-5.00668)
  , fBenchmark("FlatEsdConverter")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTGlobalFlatEsdConverterComponent::~AliHLTGlobalFlatEsdConverterComponent()
{
  // see header file for class documentation
}

int AliHLTGlobalFlatEsdConverterComponent::Configure(const char* arguments)
{
  // see header file for class documentation
  int iResult=0;
  if (!arguments) return iResult;

  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->String();	
      if (argument.IsNull()) continue;      
      HLTError("unknown argument %s", argument.Data());
      iResult=-EINVAL;
      break;
    }  
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  return iResult;
}

int AliHLTGlobalFlatEsdConverterComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  const char* path=NULL;
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object string: \'%s\'", pString->String().Data());
	iResult=Configure(pString->String().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }
  
  return iResult;
}

void AliHLTGlobalFlatEsdConverterComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeTrack);
  list.push_back(kAliHLTDataTypeTrackMC);
  list.push_back(kAliHLTDataTypeCaloCluster);
  list.push_back(kAliHLTDataTypedEdx );
  list.push_back(kAliHLTDataTypeESDVertex );
  list.push_back(kAliHLTDataTypeESDObject);
  list.push_back(kAliHLTDataTypeTObject);
  list.push_back(kAliHLTDataTypeGlobalVertexer);
  list.push_back(kAliHLTDataTypeV0Finder); // array of track ids for V0s
  list.push_back(kAliHLTDataTypeKFVertex); // KFVertex object from vertexer
  list.push_back(kAliHLTDataTypePrimaryFinder); // array of track ids for prim vertex
  list.push_back(kAliHLTDataTypeESDContent);
  list.push_back(AliHLTTPCDefinitions::fgkClustersDataType| kAliHLTDataOriginTPC);
  list.push_back(AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo| kAliHLTDataOriginTPC);
}

AliHLTComponentDataType AliHLTGlobalFlatEsdConverterComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut;
}

void AliHLTGlobalFlatEsdConverterComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation
  constBase=2000000;
  inputMultiplier=10.0;
}

int AliHLTGlobalFlatEsdConverterComponent::DoInit(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  int bMissingParam=0;

  // default list of skiped ESD objects
  TString skipObjects=
    // "AliESDRun,"
    // "AliESDHeader,"
    // "AliESDZDC,"
    "AliESDFMD,"
    // "AliESDVZERO,"
    // "AliESDTZERO,"
    // "TPCVertex,"
    // "SPDVertex,"
    // "PrimaryVertex,"
    // "AliMultiplicity,"
    // "PHOSTrigger,"
    // "EMCALTrigger,"
    // "SPDPileupVertices,"
    // "TrkPileupVertices,"
    "Cascades,"
    "Kinks,"
    "AliRawDataErrorLogs,"
    "AliESDACORDE";

  iResult=Reconfigure(NULL, NULL);
  TString allArgs = "";
  for ( int i = 0; i < argc; i++ ) {
    if ( !allArgs.IsNull() ) allArgs += " ";
    allArgs += argv[i];
  }

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->String();	
      if (argument.IsNull()) continue;

      // -noclusters
      if (argument.CompareTo("-noclusters")==0) {
	fWriteClusters=0;	
	// -clusters
      } else if (argument.CompareTo("-clusters")==0) {
	fWriteClusters=1;
      } else if (argument.Contains("-skipobject=")) {
	argument.ReplaceAll("-skipobject=", "");
	skipObjects=argument;
      } else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  fSolenoidBz=GetBz();

  if (iResult>=0) {
    SetupCTPData();
  }

  fBenchmark.SetTimer(0,"total");

  return iResult;
}

int AliHLTGlobalFlatEsdConverterComponent::DoDeinit()
{
  // see header file for class documentation

  return 0;
}

int AliHLTGlobalFlatEsdConverterComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/,
						    const AliHLTComponentBlockData* /*blocks*/, 
						    AliHLTComponentTriggerData& /*trigData*/,
						    AliHLTUInt8_t* outputPtr, 
						    AliHLTUInt32_t& size,
						    AliHLTComponentBlockDataList& outputBlocks )
{
  // see header file for class documentation
  int iResult=0;

  if (!IsDataEvent()) return iResult;

  fBenchmark.StartNewEvent();
  fBenchmark.Start(0);

  size_t maxOutputSize = size;
  size = 0;

  AliFlatESDEvent *flatEsd = reinterpret_cast<AliFlatESDEvent*>(outputPtr); 
  new (flatEsd) AliFlatESDEvent;    

  /*
  pESD->Reset(); 
  pESD->SetMagneticField(fSolenoidBz);
  pESD->SetRunNumber(GetRunNo());
  pESD->SetPeriodNumber(GetPeriodNumber());
  pESD->SetOrbitNumber(GetOrbitNumber());
  pESD->SetBunchCrossNumber(GetBunchCrossNumber());
  pESD->SetTimeStamp(GetTimeStamp());
  
  const AliHLTCTPData* pCTPData=CTPData();
  if (pCTPData) {
    AliHLTUInt64_t mask=pCTPData->ActiveTriggers(trigData);
    for (int index=0; index<gkNCTPTriggerClasses; index++) {
      if ((mask&((AliHLTUInt64_t)0x1<<index)) == 0) continue;
      pESD->SetTriggerClass(pCTPData->Name(index), index);
    }
    pESD->SetTriggerMask(mask);
  }
  */

  // Barrel tracking
  // tracks are based on the TPC tracks, and only updated from the ITS information
  // Sequence:
  // 1) extract MC information for TPC and ITS from specific data blocks and store in
  //    intermediate vector arrays
  // 2) extract TPC tracks, update with MC labels if available, the track parameters
  //    are estimated at the first cluster position
  // 2.1) propagate to last cluster position and update kTPCout, sets also outer param (fOp)
  // 2.2) update kTPCin, sets also inner param (fIp) and TPC inner param (fTPCInner)
  // 2.3) update kTPCrefit using the same parameters at the first cluster position
  //      HLT has strictly spoking no refit, but we want the flag to be set
  //      can be changed to be done after all the individual barrel detector parameters
  //      have been updated by looping over the tracks again
  // 3) extract ITS tracks, the tracks are actually TPC tracks updated from the ITS
  //    tracking information
  // 3.1) TODO 2010-07-12: handle ITS standalone tracks by updating kITSout before kITSin
  // 3.2) update with kITSin
  //    TODO 2010-07-12 find out if the kITSrefit has to be set as well
  // 4) extract TRD tracks and add to ESD
  //    TODO 2010-07-12 at the moment there is no matching or merging of TPC and TRD tracks
  // 5) Add Trigger Detectors 
  //    VZERO, ZDC

  // 1) first read MC information (if present)

  std::map<int,int> mcLabelsTPC;
  std::map<int,int> mcLabelsITS;

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrackMC|kAliHLTDataOriginTPC);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    fBenchmark.AddInput(pBlock->fSize);
    AliHLTTrackMCData* dataPtr = reinterpret_cast<AliHLTTrackMCData*>( pBlock->fPtr );
    if (sizeof(AliHLTTrackMCData)+dataPtr->fCount*sizeof(AliHLTTrackMCLabel)==pBlock->fSize) {
      for( unsigned int il=0; il<dataPtr->fCount; il++ ){
	AliHLTTrackMCLabel &lab = dataPtr->fLabels[il];
	mcLabelsTPC[lab.fTrackID] = lab.fMCLabel;
      }
    } else {
      HLTWarning("data mismatch in block %s (0x%08x): count %d, size %d -> ignoring track MC information", 
		 DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, 
		 dataPtr->fCount, pBlock->fSize);
    }
  }
 
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrackMC|kAliHLTDataOriginITS);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    fBenchmark.AddInput(pBlock->fSize);
    AliHLTTrackMCData* dataPtr = reinterpret_cast<AliHLTTrackMCData*>( pBlock->fPtr );
    if (sizeof(AliHLTTrackMCData)+dataPtr->fCount*sizeof(AliHLTTrackMCLabel)==pBlock->fSize) {
      for( unsigned int il=0; il<dataPtr->fCount; il++ ){
	AliHLTTrackMCLabel &lab = dataPtr->fLabels[il];
	mcLabelsITS[lab.fTrackID] = lab.fMCLabel;
      }
    } else {
      HLTWarning("data mismatch in block %s (0x%08x): count %d, size %d -> ignoring track MC information", 
		 DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, 
		 dataPtr->fCount, pBlock->fSize);
    }
  }

  // 2) read dEdx information (if present)

  AliHLTFloat32_t *dEdxTPC = 0; 
  Int_t ndEdxTPC = 0;
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypedEdx|kAliHLTDataOriginTPC);
       pBlock!=NULL; pBlock=NULL/*GetNextInputBlock() there is only one block*/) {
    fBenchmark.AddInput(pBlock->fSize);
    dEdxTPC = reinterpret_cast<AliHLTFloat32_t*>( pBlock->fPtr );
    ndEdxTPC = pBlock->fSize / (3*sizeof(AliHLTFloat32_t));
  }

  // 3) read  TPC tracks 

  vector<AliHLTGlobalBarrelTrack> tracksTPC;

  { // there is only one block of TPC data expected
    const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
    if( pBlock ){
      fBenchmark.AddInput(pBlock->fSize);
      iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracksTPC);
    }
    if( iResult>=0 ){
	    HLTWarning("converted %d track(s) to AliESDtrack and added to ESD", tracksTPC.size());
    } else if (iResult<0) {
	    HLTError("can not extract tracks from data block of type %s (specification %08x) of size %d: error %d", 
		     DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize, iResult);
    }
  }
 
  // 4) read ITS refitted tracks

  vector<AliHLTGlobalBarrelTrack> tracksITS;
  vector<AliHLTGlobalBarrelTrack> tracksITSOut;

  if( iResult>=0 ) { // there is only one block of ITS tracks expected
     const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginITS);
     if( pBlock ){
       fBenchmark.AddInput(pBlock->fSize);
       iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracksITS);
     }
  } 

  if( iResult>=0 ) { // there is only one block of ITS tracks expected
    const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginITSOut);
    if( pBlock ){
      fBenchmark.AddInput(pBlock->fSize);
      iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracksITS);
    }
  } 

   
  // read TPC clusters
  const UInt_t kNSlices = 36;
  const UInt_t kNPatches = 6;

  const AliHLTTPCClusterData *clustersTPC[kNSlices][kNPatches];
  const AliHLTTPCClusterMCLabel *clustersTPCMC[kNSlices][kNPatches];
  for( UInt_t i=0; i<kNSlices; i++){
     for( UInt_t j=0; j<kNPatches; j++){
      clustersTPC[i][j] = 0;
      clustersTPCMC[i][j] =0;
     }
  }

  fWriteClusters = 1;

  if( fWriteClusters ){

    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType| kAliHLTDataOriginTPC);
	 pBlock!=NULL; pBlock=GetNextInputBlock()) {
      //fBenchmark.AddInput(pBlock->fSize);
      UInt_t slice     = AliHLTTPCDefinitions::GetMinSliceNr(*pBlock); 
      UInt_t patch  = AliHLTTPCDefinitions::GetMinPatchNr(*pBlock);
      if( slice >= kNSlices || patch>= kNPatches ){
	HLTWarning("Wrong slice / patch number of cluster block");
			  continue;
      } 
      clustersTPC[slice][patch] = reinterpret_cast <const AliHLTTPCClusterData*> ( pBlock->fPtr );
    }
    
    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo| kAliHLTDataOriginTPC);
	 pBlock!=NULL; pBlock=GetNextInputBlock()) {
      //fBenchmark.AddInput(pBlock->fSize);
      UInt_t slice     = AliHLTTPCDefinitions::GetMinSliceNr(*pBlock); 
      UInt_t patch  = AliHLTTPCDefinitions::GetMinPatchNr(*pBlock);
      if( slice >= kNSlices || patch>= kNPatches ){
	HLTWarning("Wrong slice / patch number of cluster MC block");
	continue;
      } 
      clustersTPCMC[slice][patch] = reinterpret_cast <const AliHLTTPCClusterMCLabel*> ( pBlock->fPtr );
    }

  }

  // Fill vertex information to the flat ESD

  const AliESDVertex *primaryVertex = 0;
  {
    const AliESDVertex *primaryVertexSPD = dynamic_cast<const AliESDVertex*>( GetFirstInputObject( kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS ) );
    const AliESDVertex *primaryVertexTracks = dynamic_cast<const AliESDVertex*>( GetFirstInputObject( kAliHLTDataTypeESDVertex|kAliHLTDataOriginOut ) );
    
    cout<<endl<<" Primary vertex Tracks: "<<primaryVertexTracks<<", SPD: "<< primaryVertexSPD <<endl<<endl;

    flatEsd->FillPrimaryVertices( primaryVertexSPD, primaryVertexTracks );
    
    primaryVertex = primaryVertexTracks;
    if( !primaryVertex ) primaryVertex = primaryVertexSPD;  
  }

  // Fill the track information to the flat ESD structure
  {
    UInt_t itsIter = 0;
    UInt_t itsOutIter = 0;
    for( UInt_t tpcIter=0; tpcIter < tracksTPC.size(); tpcIter++) {
       
       const AliHLTGlobalBarrelTrack *tpcTrack = &(tracksTPC[tpcIter]);
       /*
       Float_t points[4] = {
	  static_cast<Float_t>(tpcTrack->GetX()),
	  static_cast<Float_t>(tpcTrack->GetY()),
	  static_cast<Float_t>(tpcTrack->GetLastPointX()),
	  static_cast<Float_t>(tpcTrack->GetLastPointY())
	};
       */
	Int_t tpcLabel = -1;
	if( mcLabelsTPC.find(tpcTrack->TrackID())!=mcLabelsTPC.end() )
	  tpcLabel = mcLabelsTPC[tpcTrack->TrackID()];

	//tpcTrack->SetLabel( tpcLabel );
	// iotrack.SetID( tpcTrack->TrackID() );

	// set kTPCout - just propagate to the outermost TPC cluster
	
	AliHLTGlobalBarrelTrack outPar(*tpcTrack);
	{
	  //outPar.AliExternalTrackParam::PropagateTo( tpcTrack->GetLastPointX(), fSolenoidBz );
	  const Int_t N=10; // number of steps.
	  const Float_t xRange = tpcTrack->GetLastPointX() - tpcTrack->GetX();
	  const Float_t xStep = xRange / N ;
	  for(int i = 1; i <= N; ++i) {
	    if(!outPar.AliExternalTrackParam::PropagateTo(tpcTrack->GetX() + xStep * i, fSolenoidBz)) break;
	  }
	}
	
	//iotrack.SetTPCPoints(points);
	/*
	if( tpcTrack->TrackID()<ndEdxTPC ){
	  AliHLTFloat32_t *val = &(dEdxTPC[3*tpcTrack->TrackID()]);
	  iotrack.SetTPCsignal( val[0], val[1], (UChar_t) val[2] ); 
	  //AliTPCseed s;
	  //s.Set( tpcTrack->GetX(), tpcTrack->GetAlpha(),
	  //tpcTrack->GetParameter(), tpcTrack->GetCovariance() );
	  //s.SetdEdx( val[0] );
	  //s.CookPID();
	  //iotrack.SetTPCpid(s.TPCrPIDs() );
	} else {
	  if( dEdxTPC ) HLTWarning("Wrong number of dEdx TPC labels");
	}
	*/
	//iotrack.SetLabel(mcLabel);

	// ITS track

	AliHLTGlobalBarrelTrack *itsRefit=0;
	Int_t itsLabel = -1;
	
	for(; itsIter< tracksITS.size() && tracksITS[itsIter].TrackID()<(int) tpcIter; itsIter++ );

	if( itsIter< tracksITS.size() && tracksITS[itsIter].TrackID() == (int) tpcIter ){
	  itsRefit = &(tracksITS[itsIter]);
	  if( mcLabelsITS.find(tpcIter)!=mcLabelsITS.end() ) itsLabel = mcLabelsITS[tpcIter];
	  itsIter++;
	}

	// ITS Out track

	AliHLTGlobalBarrelTrack *itsOut=0;
	
	for(; itsOutIter< tracksITSOut.size() && tracksITSOut[itsOutIter].TrackID()<(int) tpcIter; itsOutIter++ );

	if( itsOutIter< tracksITSOut.size() && tracksITSOut[itsOutIter].TrackID() == (int) tpcIter ){
	  itsOut = &(tracksITSOut[itsOutIter]);
	  itsOutIter++;
	}

	// Fill DCA parameters for TPC tracks
	AliESDtrack cP;
	
	if( primaryVertex ){
	    cP.UpdateTrackParams( (itsRefit ?itsRefit :tpcTrack), AliESDtrack::kTPCin );
	    cP.RelateToVertex( primaryVertex, fSolenoidBz, 1000 );    	
	}
	
	AliFlatESDTrack *flatTrack = flatEsd->GetNextTrackPointer();

	UInt_t nClustersTPC = tpcTrack->GetNumberOfPoints();
	UInt_t nClustersITS = itsRefit ?itsRefit->GetNumberOfPoints() :0;

	flatTrack->SetNumberOfITSClusters( nClustersITS );

	if( flatEsd->GetSize() + flatTrack->EstimateSize( kTRUE, nClustersTPC ) >= maxOutputSize ){
		cout<<endl<<endl<<"NOT ENOUGH MEMORY!!!!"<<endl<<endl;
	   iResult=-ENOMEM;
	   break;
	}
	
	flatTrack->FillExternalTrackParam( itsRefit,  NULL, tpcTrack, &outPar, cP.GetConstrainedParam(), itsOut);
	
	if( fWriteClusters && tpcTrack->GetPoints() ){
	  const UInt_t* clusterIDs = tpcTrack->GetPoints();
	   for( UInt_t i=0; i<nClustersTPC; i++ ){
	      UInt_t id = clusterIDs[i];
	      UInt_t iSlice = AliHLTTPCSpacePointData::GetSlice(id);
	      UInt_t iPatch = AliHLTTPCSpacePointData::GetPatch(id);
	      UInt_t iCluster = AliHLTTPCSpacePointData::GetNumber(id);
	      if( iSlice >= kNSlices || iPatch>= kNPatches ){
		 HLTWarning("Wrong slice / patch number of TPC cluster");
		 continue;
	      } 
	      const AliHLTTPCClusterData *clusterBlock = clustersTPC[iSlice][iPatch];
	      if( !clusterBlock ){
		 HLTWarning("no cluster block found for slice %d, patch %d",iSlice,iPatch);
		 continue;
	      } 
	      if( iCluster>= clusterBlock->fSpacePointCnt ){
		      HLTWarning("no cluster block found for slice %d, patch %d, cluster %d",iSlice,iPatch,iCluster);
		 continue;
	      } 
	      const AliHLTTPCSpacePointData &cIn = clusterBlock->fSpacePoints[iCluster];
	      AliFlatTPCCluster *c= flatTrack->GetNextTPCClusterPointer();
	      c->fX = cIn.GetX();
	      c->fY = cIn.GetY();
	      c->fZ = cIn.GetZ();
	      c->fPadRow  = cIn.GetPadRow() + AliHLTTPCTransform::GetFirstRow(iPatch);
	      c->fSigmaY2 = cIn.GetSigmaY2();
	      c->fSigmaZ2 = cIn.GetSigmaZ2();
	      c->fCharge  = cIn.GetCharge();
	      c->fQMax    = cIn.GetQMax();
	      flatTrack->StoreLastTPCCluster();
	   }
	}

 	
	flatEsd->StoreLastTrack();
	
	if (fVerbosity>0) tpcTrack->Print();
    }    
  }

  // Fill v0's
  
  {    
    int nV0s =0;
    const AliHLTComponentBlockData* pP = GetFirstInputBlock(kAliHLTDataTypeGlobalVertexer|kAliHLTDataOriginOut);
    if (pP && pP->fSize && pP->fPtr) {
      const AliHLTGlobalVertexerComponent::AliHLTGlobalVertexerData *data = reinterpret_cast<AliHLTGlobalVertexerComponent::AliHLTGlobalVertexerData*>(pP->fPtr);
      const int* v0s = data->fTrackIndices + data->fNPrimTracks;
      nV0s = data->fNV0s;
      for (int i = 0; i < nV0s; ++i) {
	AliFlatESDV0 *v0 = flatEsd->GetNextV0Pointer();
	v0->fNegTrackID = v0s[2 * i];
	v0->fPosTrackID = v0s[2 * i + 1];
	flatEsd->StoreLastV0();
      }
    } else {
      HLTWarning("xxx No V0 data block");
    }
    cout<<"\nxxxx Found "<<nV0s<<" V0's\n"<<endl;
  }
  
  // Get ITS SPD vertex
  for( const AliHLTComponentBlockData *i= GetFirstInputBlock(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS); i!=NULL; i=GetNextInputBlock() ){
    fBenchmark.AddInput(i->fSize);
  }

  /*
  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS); iter != NULL; iter = GetNextInputObject() ) {
    AliESDVertex *vtx = dynamic_cast<AliESDVertex*>(const_cast<TObject*>( iter ) );
    pESD->SetPrimaryVertexSPD( vtx );
  }
  

  // update with  vertices and vertex-fitted tracks
  // output of the GlobalVertexerComponent
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeGlobalVertexer);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    fBenchmark.AddInput(pBlock->fSize);   
    AliHLTGlobalVertexerComponent::FillESD( pESD, reinterpret_cast<AliHLTGlobalVertexerComponent::AliHLTGlobalVertexerData* >(pBlock->fPtr) );
  }

  // update with  vertices and vertex-fitted tracks
  // output of PrimaryVertexer and V0Finder components
  TObject* pBase = (TObject*)GetFirstInputObject(kAliHLTDataTypeKFVertex | kAliHLTDataOriginOut);
  if (pBase) {
    AliKFVertex* kfVertex = dynamic_cast<AliKFVertex *>(pBase);
    if (kfVertex) {
      const AliHLTComponentBlockData* pP = GetFirstInputBlock(kAliHLTDataTypePrimaryFinder | kAliHLTDataOriginOut);
      if (pP && pP->fSize && pP->fPtr) {
	const AliHLTComponentBlockData* pV0 = GetFirstInputBlock(kAliHLTDataTypeV0Finder | kAliHLTDataOriginOut);
	if (pV0 && pV0->fPtr && pInputESD && pInputESD->GetNumberOfV0s()>0) {
	  const int* v0s = static_cast<const int*>(pV0->fPtr);
	  HLTWarning("V0 array already filled from the input esd block, additional filling from V0 block of %d entries might cause inconsistent content", v0s[0]);
	}
	AliHLTVertexFinderBase::FillESD(pESD, kfVertex, pP->fPtr, pV0?pV0->fPtr:NULL);
      } else
	HLTWarning("Problem with primary finder's data block");
    } else {
      HLTWarning("primary vertex block of wrong type, expecting AliKFVertex instead of %s", pBase->GetName());
    }
  } else {
    // throw an error if there is a V0 data block which can not be handled without
    // the AliKFVertex object
    if (GetFirstInputBlock(kAliHLTDataTypeV0Finder | kAliHLTDataOriginOut)!=NULL) {
      ALIHLTERRORGUARD(1, "missing AliKFVertex object ignoring V0 data block of type %s",
		       DataType2Text(kAliHLTDataTypeV0Finder|kAliHLTDataOriginOut).c_str());
    }
  }
  */


  // loop over all tracks and set the TPC refit flag by updating with the
  // original TPC inner parameter if not yet set
  // TODO: replace this by a proper refit
  // code is comented for the moment as it does not fully solve the problems with
  // the display
  // - would set the main parameters to the TPC inner wall again, or
  // - changes the inner param if the parameters are propagated, so we loose the track
  //   reference point for the display
  // with the current sequence we have the latter case as the DCA operations above
  // change the TPC inner parameters
  /*
  for (int i=0; i<pESD->GetNumberOfTracks(); i++) {
    if (!pESD->GetTrack(i) || 
	!pESD->GetTrack(i)->GetTPCInnerParam() ||
	pESD->GetTrack(i)->IsOn(AliESDtrack::kTPCrefit)) continue;
    AliESDtrack* tESD=pESD->GetTrack(i);
    AliHLTGlobalBarrelTrack inner(*tESD->GetTPCInnerParam());
    inner.SetLabel(tESD->GetLabel());
    tESD->UpdateTrackParams(&inner, AliESDtrack::kTPCrefit);
  }
  */
 
  if (iResult>=0) {            
 
    AliHLTComponentBlockData outBlock;
    FillBlockData( outBlock );
    outBlock.fOffset = size;
    outBlock.fSize = flatEsd->GetSize();
    outBlock.fDataType = kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut;
    outBlock.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( 0, 35, 0, 5 );

    outputBlocks.push_back( outBlock );

    fBenchmark.AddOutput(outBlock.fSize);
      
    size += outBlock.fSize;
  }

  fBenchmark.Stop(0);
  HLTWarning( fBenchmark.GetStatistics() );

  return iResult;
}

