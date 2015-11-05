// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Mikolaj Krzewicki <mikolaj.krzewicki@cern.ch>         *
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

//  @file   AliHLTGlobalPromptRecoQAComponent.cxx
//  @author Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//  @date   
//  @brief  Simple HLT reco QA/monitor
// 

#include <cassert>
#include "AliHLTGlobalPromptRecoQAComponent.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTTrackMCLabel.h"
#include "AliHLTCTPData.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliITStrackV2.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliHLTTPCGeometry.h"
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
#include "TTimeStamp.h"
#include "THnSparse.h"
#include "AliHLTESDCaloClusterMaker.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "AliESDCaloCluster.h"
#include "AliESDVZERO.h"
#include "AliHLTGlobalVertexerComponent.h"
#include "AliHLTVertexFinderBase.h"
#include "AliSysInfo.h"
#include "AliHLTSAPTrackerData.h"
#include "AliFlatESDVertex.h"
#include "tracking-ca/AliHLTTPCCADefinitions.h"
#include "tracking-ca/AliHLTTPCCACompressedInputData.h"
#include "tracking-ca/AliHLTTPCCASliceOutput.h"

#include "TH1I.h"
#include <string>
#include "AliHLTITSClusterDataFormat.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalPromptRecoQAComponent)

AliHLTGlobalPromptRecoQAComponent::AliHLTGlobalPromptRecoQAComponent()
  : AliHLTProcessor()
  , fVerbosity(0)
  , fBenchmark("PromptRecoQA")
  , fSkipEvents(0)
  , fPrintStats(0)
  , fEventsSinceSkip(0)
  , fHistSPDclusters_SPDrawSize(NULL)
  , fHistSSDclusters_SSDrawSize(NULL)
  , fHistSDDclusters_SDDrawSize(NULL)
  , fHistITSSAtracks_SPDclusters(NULL)
  , fHistSPDclusters_SSDclusters(NULL)
  , fHistTPCHLTclusters_TPCHLTclustersSize(NULL)
  , fHistTPCtracks_TPCtracklets(NULL)
  , fHistITStracks_ITSOutTracks(NULL)

{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
}

AliHLTGlobalPromptRecoQAComponent::~AliHLTGlobalPromptRecoQAComponent()
{
  // see header file for class documentation
  delete fHistSPDclusters_SPDrawSize;
  delete fHistSDDclusters_SDDrawSize;
  delete fHistSSDclusters_SSDrawSize;
  delete fHistITSSAtracks_SPDclusters;
  delete fHistSPDclusters_SSDclusters;
  delete fHistTPCHLTclusters_TPCHLTclustersSize;
  delete fHistTPCtracks_TPCtracklets;
  delete fHistITStracks_ITSOutTracks;
}

int AliHLTGlobalPromptRecoQAComponent::Configure(const char* arguments)
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
      
      if (argument.CompareTo("-skip-events")==0) {
	      argument.ReplaceAll("-skip-events=","");
        fSkipEvents = argument.Atoi();
      }	else if (argument.CompareTo("-print-stats")==0) {
        fPrintStats = 1;
      }	else if (!argument.Contains("pushback-period")) {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  return iResult;
}

int AliHLTGlobalPromptRecoQAComponent::Reconfigure(const char* cdbEntry, const char* chainId)
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

void AliHLTGlobalPromptRecoQAComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeClusters | kAliHLTDataOriginITSSPD);
  list.push_back(kAliHLTDataTypeClusters | kAliHLTDataOriginITSSDD);
  list.push_back(kAliHLTDataTypeClusters | kAliHLTDataOriginITSSSD);
  list.push_back(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSPD);
  list.push_back(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSDD);
  list.push_back(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSSD);
  list.push_back(kAliHLTDataTypeITSSAPData | kAliHLTDataOriginITS);
  list.push_back(kAliHLTDataTypeESDVertex | kAliHLTDataOriginITSSPD); //SPD Vertex
  list.push_back(AliHLTTPCDefinitions::fgkHWClustersDataType | kAliHLTDataOriginTPC); //HLT-TPC clusters from HWCF
  list.push_back(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC); //TPC DDL raw data
  list.push_back(AliHLTTPCDefinitions::fgkClustersDataType | kAliHLTDataOriginTPC); //Transformed HLT-TPC clusters
  list.push_back(AliHLTTPCCADefinitions::fgkCompressedInputDataType | kAliHLTDataOriginTPC); //Transformed HLT-TPC clusters (internally compressed form)
  list.push_back(AliHLTTPCDefinitions::fgkRawClustersDataType | kAliHLTDataOriginTPC); //Non-transformed HLT-TPC clusters
  list.push_back(AliHLTTPCCADefinitions::fgkTrackletsDataType); //HLT-TPC Tracklets (before TPC global merger)
  list.push_back(kAliHLTDataTypeTrack | kAliHLTDataOriginTPC); //HLT-TPC merged tracks
  list.push_back(kAliHLTDataTypeTrack | kAliHLTDataOriginITS); //TPC-ITS tracks
  list.push_back(kAliHLTDataTypeTrack | kAliHLTDataOriginITSOut); //ITS-Out merged tracks


  //All this is TPC Data compression
  list.push_back(AliHLTTPCDefinitions::DataCompressionDescriptorDataType());
  list.push_back(AliHLTTPCDefinitions::RawClustersDataType());
  list.push_back(AliHLTTPCDefinitions::RemainingClustersCompressedDataType());
  list.push_back(AliHLTTPCDefinitions::RemainingClusterIdsDataType());
  list.push_back(AliHLTTPCDefinitions::ClusterTracksCompressedDataType());
  list.push_back(AliHLTTPCDefinitions::ClusterIdTracksDataType());
}

AliHLTComponentDataType AliHLTGlobalPromptRecoQAComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTMultipleDataType;
}

int AliHLTGlobalPromptRecoQAComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList){ 
// see header file for class documentation

  tgtList.clear();
  tgtList.push_back( kAliHLTDataTypeHistogram|kAliHLTDataOriginOut );
  tgtList.push_back( kAliHLTDataTypeTObject|kAliHLTDataOriginOut );
  return tgtList.size();
}

void AliHLTGlobalPromptRecoQAComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation
  constBase=2000000;
  inputMultiplier=0.0;
}

int AliHLTGlobalPromptRecoQAComponent::DoInit(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;

  std::string argString = GetComponentArgs();
  if (Configure(argString.c_str())<0) return -EINVAL;

  fHistSPDclusters_SPDrawSize = new TH2I("SPDncls_SPDsize", "SPDncls vs SPD raw size", 50, 0., 3000., 50, 0., 10000.);
  fHistSSDclusters_SSDrawSize = new TH2I("SSDncls_SSDsize", "SSDncls vs SSD raw size", 50, 0., 2000., 50, 0., 40000.);
  fHistSDDclusters_SDDrawSize = new TH2I("SDDncls_SDDsize", "SDDncls vs SDD raw size", 50, 0., 1000., 50, 0., 20000.);
  fHistITSSAtracks_SPDclusters = new TH2I("ITSSAPntrk_SPDncls", "ITSSAP tracks vs SPD ncls", 50, 0., 1000., 50, 0., 10000.);
  fHistSPDclusters_SSDclusters = new TH2I("SPDncls_SSDncls", "SPDncls vs SSDncls", 50, 0., 3000., 50, 0., 2000.);
  fHistTPCHLTclusters_TPCHLTclustersSize = new TH2I("TPCHLTncls_TPCHLTsize", "TPCHLTSncls vs size", 50, 0., 5000000., 50, 0., 80000000.);
  fHistTPCtracks_TPCtracklets = new TH2I("TPCntrk_TPCntrl", "TPCntrk vs TPCnTracklets", 50, 0., 25000., 50, 0., 40000.);
  fHistITStracks_ITSOutTracks = new TH2I("ITSntrk_ITSOutntrk", "ITSntrk vs ITSOutntrk", 50, 0., 25000., 50, 0., 25000.);

  return iResult;
}

int AliHLTGlobalPromptRecoQAComponent::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

int AliHLTGlobalPromptRecoQAComponent::DoEvent( const AliHLTComponentEventData& evtData,
    const AliHLTComponentBlockData* blocks, 
    AliHLTComponentTriggerData& /*trigData*/,
    AliHLTUInt8_t* /*outputPtr*/, 
    AliHLTUInt32_t& /*size*/,
    AliHLTComponentBlockDataList& /*outputBlocks*/ )
{
  // see header file for class documentation
  int iResult=0;

  //perhaps downscale a bit
  if (fEventsSinceSkip++<fSkipEvents)
  {
    fEventsSinceSkip=0;
    return 0;
  }

  //what are we interested in?
  AliHLTUInt32_t nClustersSPD = 0;
  AliHLTUInt32_t rawSizeSPD = 0;
  AliHLTUInt32_t nClustersSDD = 0;
  AliHLTUInt32_t rawSizeSDD = 0;
  AliHLTUInt32_t nClustersSSD = 0;
  AliHLTUInt32_t rawSizeSSD = 0;
  AliHLTUInt32_t nClustersITS = 0;
  AliHLTUInt32_t rawSizeITS = 0;
  
  AliHLTUInt32_t nClustersTPC = 0;
  AliHLTUInt32_t rawSizeTPC = 0;
  AliHLTUInt32_t hwcfSizeTPC = 0;
  AliHLTUInt32_t clusterSizeTPC = 0;
  AliHLTUInt32_t compressedSizeTPC = 0;

  AliHLTUInt32_t nITSSAPtracks = 0;
  AliHLTUInt32_t nSPDtracklets =0;
  AliHLTUInt32_t nTPCtracklets = 0;
  AliHLTUInt32_t nTPCtracks = 0;
  AliHLTUInt32_t nITSTracks = 0;
  AliHLTUInt32_t nITSOutTracks = 0;
  
  Bool_t bITSSPDVertex = kFALSE;
  

  //loop over input blocks and extract basic stats
  int nBlocks = evtData.fBlockCnt;  
  for (int ndx=0; ndx<nBlocks; ndx++) {
    const AliHLTComponentBlockData* iter = blocks+ndx;
    
    //Vertex Found
    if (iter->fDataType == (kAliHLTDataTypeESDVertex | kAliHLTDataOriginITSSPD))
    {
      bITSSPDVertex = kTRUE;
    }

    //numbers of clusters
    if (iter->fDataType == (kAliHLTDataTypeClusters | kAliHLTDataOriginITSSPD))
    {
      AliHLTITSClusterData *inPtr=reinterpret_cast<AliHLTITSClusterData*>( iter->fPtr );
      nClustersSPD += inPtr->fSpacePointCnt;
    }
    if (iter->fDataType == (kAliHLTDataTypeClusters | kAliHLTDataOriginITSSDD))
    {
      AliHLTITSClusterData *inPtr=reinterpret_cast<AliHLTITSClusterData*>( iter->fPtr );
      nClustersSDD += inPtr->fSpacePointCnt;
    }
    if (iter->fDataType == (kAliHLTDataTypeClusters | kAliHLTDataOriginITSSSD))
    {
      AliHLTITSClusterData *inPtr=reinterpret_cast<AliHLTITSClusterData*>( iter->fPtr );
      nClustersSSD += inPtr->fSpacePointCnt;
    }
    if (iter->fDataType == (kAliHLTDataTypeClusters | kAliHLTDataOriginITS))
    {
      AliHLTITSClusterData *inPtr=reinterpret_cast<AliHLTITSClusterData*>( iter->fPtr );
      nClustersITS += inPtr->fSpacePointCnt;
    }

    if (iter->fDataType == AliHLTTPCDefinitions::fgkClustersDataType) //Transformed TPC clusters used in TPCCATracker
    {
      AliHLTTPCClusterData* inPtrSP = ( AliHLTTPCClusterData* )( iter->fPtr );
      nClustersTPC += inPtrSP->fSpacePointCnt;
    }
    else if (iter->fDataType == AliHLTTPCCADefinitions::fgkCompressedInputDataType) //Compressed (internally) form of transformed HLT TPC clusters (currently not used)
    {
      const AliHLTUInt8_t * inPtr = (const AliHLTUInt8_t *)iter->fPtr;
      while(inPtr< ((const AliHLTUInt8_t *) iter->fPtr) + iter->fSize)
      {
        AliHLTTPCCACompressedClusterRow *row = (AliHLTTPCCACompressedClusterRow*) inPtr;
        nClustersTPC+= row->fNClusters;
        inPtr = (const AliHLTUInt8_t *)(row->fClusters+row->fNClusters);
      }
    }

    if (iter->fDataType == AliHLTTPCDefinitions::DataCompressionDescriptorDataType() || //Used
      iter->fDataType == AliHLTTPCDefinitions::RawClustersDataType() ||
      iter->fDataType == AliHLTTPCDefinitions::RemainingClustersCompressedDataType() ||
      iter->fDataType == AliHLTTPCDefinitions::RemainingClusterIdsDataType() ||
      iter->fDataType == AliHLTTPCDefinitions::ClusterTracksCompressedDataType() || //Used
      iter->fDataType == AliHLTTPCDefinitions::ClusterIdTracksDataType()) //Used
    {
      compressedSizeTPC += iter->fSize;
    }

    //RAW sizes
    if (iter->fDataType == (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSPD))
    {
      rawSizeSPD += iter->fSize;
    }
    if (iter->fDataType == (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSDD))
    {
      rawSizeSDD += iter->fSize;
    }
    if (iter->fDataType == (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSSD))
    {
      rawSizeSSD += iter->fSize;
    }
    if (iter->fDataType == (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITS))
    {
      rawSizeITS += iter->fSize;
    }
    if (iter->fDataType == (AliHLTTPCDefinitions::fgkHWClustersDataType | kAliHLTDataOriginTPC)) //Size of HLT-TPC hardware clusters
    {
      hwcfSizeTPC += iter->fSize;
    }
    
    if (iter->fDataType == (AliHLTTPCDefinitions::fgkRawClustersDataType | kAliHLTDataOriginTPC)) //Size of HLT-TPC clusters
    {
      clusterSizeTPC += iter->fSize;
    }
    if (iter->fDataType == (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC)) //TPC RAW DDL Size
    {
      rawSizeTPC += iter->fSize;
    }

    //numbers of tracks
    if (iter->fDataType == AliHLTTPCCADefinitions::fgkTrackletsDataType) //HLT-TPC CA-trackets (before TPC global merger)
    {
      AliHLTTPCCASliceOutput* out = reinterpret_cast<AliHLTTPCCASliceOutput*>(iter->fPtr);
      nTPCtracklets += out->NTracks();
    }

    if (iter->fDataType == (kAliHLTDataTypeTrack | kAliHLTDataOriginTPC))
    {
      nTPCtracks += ((AliHLTTracksData*) iter->fPtr)->fCount;
    }

    if (iter->fDataType == (kAliHLTDataTypeTrack | kAliHLTDataOriginITS))
    {
      nITSTracks += ((AliHLTTracksData*) iter->fPtr)->fCount;
    }
    if (iter->fDataType == (kAliHLTDataTypeTrack | kAliHLTDataOriginITSOut))
    {
      nITSOutTracks += ((AliHLTTracksData*) iter->fPtr)->fCount;
    }

    if (iter->fDataType == (kAliHLTDataTypeITSSAPData | kAliHLTDataOriginITS))
    {
      AliHLTITSSAPTrackerDataContainer* inPtr = reinterpret_cast<AliHLTITSSAPTrackerDataContainer*>(iter->fPtr);
      nITSSAPtracks += inPtr->fCount;
    }
  }// end read input blocks
  

  //fill histograms
  fHistSPDclusters_SPDrawSize->Fill(nClustersSPD, rawSizeSPD);
  if (fHistSPDclusters_SPDrawSize->GetEntries() && PushBack(fHistSPDclusters_SPDrawSize, kAliHLTDataTypeHistogram|kAliHLTDataOriginOut) > 0)
    fHistSPDclusters_SPDrawSize->Reset();

  fHistSDDclusters_SDDrawSize->Fill(nClustersSDD, rawSizeSDD);
  if (fHistSDDclusters_SDDrawSize->GetEntries() && PushBack(fHistSDDclusters_SDDrawSize, kAliHLTDataTypeHistogram|kAliHLTDataOriginOut) > 0)
    fHistSDDclusters_SDDrawSize->Reset();

  fHistSSDclusters_SSDrawSize->Fill(nClustersSSD, rawSizeSSD);
  if (fHistSSDclusters_SSDrawSize->GetEntries() && PushBack(fHistSSDclusters_SSDrawSize, kAliHLTDataTypeHistogram|kAliHLTDataOriginOut) > 0)
    fHistSSDclusters_SSDrawSize->Reset();

  fHistITSSAtracks_SPDclusters->Fill(nITSSAPtracks, nClustersSPD);
  if (fHistITSSAtracks_SPDclusters->GetEntries() && PushBack(fHistITSSAtracks_SPDclusters, kAliHLTDataTypeHistogram|kAliHLTDataOriginOut) > 0)
    fHistITSSAtracks_SPDclusters->Reset();

  fHistSPDclusters_SSDclusters->Fill(nClustersSPD, nClustersSSD);
  if (fHistSPDclusters_SSDclusters->GetEntries() && PushBack(fHistSPDclusters_SSDclusters, kAliHLTDataTypeHistogram|kAliHLTDataOriginOut) > 0)
    fHistSPDclusters_SSDclusters->Reset();

  fHistTPCHLTclusters_TPCHLTclustersSize->Fill(nClustersTPC, clusterSizeTPC);
  if (fHistTPCHLTclusters_TPCHLTclustersSize->GetEntries() && PushBack(fHistTPCHLTclusters_TPCHLTclustersSize, kAliHLTDataTypeHistogram|kAliHLTDataOriginOut) > 0)
    fHistTPCHLTclusters_TPCHLTclustersSize->Reset();

  fHistTPCtracks_TPCtracklets->Fill(nTPCtracks, nTPCtracklets);
  if (fHistTPCtracks_TPCtracklets->GetEntries() && PushBack(fHistTPCtracks_TPCtracklets, kAliHLTDataTypeHistogram|kAliHLTDataOriginOut) > 0)
    fHistTPCtracks_TPCtracklets->Reset();


  int pushed_something = 0;
  fHistITStracks_ITSOutTracks->Fill(nITSTracks, nITSOutTracks);
  if (fHistITStracks_ITSOutTracks->GetEntries() && PushBack(fHistITStracks_ITSOutTracks, kAliHLTDataTypeHistogram|kAliHLTDataOriginOut) > 0)
  {
    fHistITStracks_ITSOutTracks->Reset();
    pushed_something = 1;
  }

  if (fPrintStats && pushed_something) //Don't print this for every event if we use a pushback period
  {
    HLTImportant("Blocks %d: HLT Reco QA Stats: SPD-Cl %d (%d), SDD-Cl %d (%d), SSD-Cl %d (%d) ITS-Cl %d (%d) TPC-Cl %d (%d / %d / %d), TPC-Comp (%d), ITSSAP-Tr %d, SPD-Tr %d, TPC-Tr %d / %d, ITS-Tr %d / %d, SPD-Ver %d",
      nBlocks, nClustersSPD, rawSizeSPD, nClustersSDD, rawSizeSDD, nClustersSSD, rawSizeSSD, nClustersITS, rawSizeITS, nClustersTPC, rawSizeTPC, hwcfSizeTPC, clusterSizeTPC, compressedSizeTPC, nITSSAPtracks, nSPDtracklets, nTPCtracklets, nTPCtracks, nITSTracks, nITSOutTracks, (int) bITSSPDVertex);
  }

  return iResult;
}
