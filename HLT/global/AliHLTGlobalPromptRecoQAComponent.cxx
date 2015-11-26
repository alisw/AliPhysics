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
#include "AliESDZDC.h"
#include "AliHLTGlobalVertexerComponent.h"
#include "AliHLTVertexFinderBase.h"
#include "AliSysInfo.h"
#include "AliHLTSAPTrackerData.h"
#include "AliFlatESDVertex.h"
#include "tracking-ca/AliHLTTPCCADefinitions.h"
#include "tracking-ca/AliHLTTPCCACompressedInputData.h"
#include "tracking-ca/AliHLTTPCCASliceOutput.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliHLTTPCHWCFData.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"

#include "TH1I.h"
#include "TH2F.h"
#include <string>
#include "AliHLTITSClusterDataFormat.h"
#include "AliHLTCDHWrapper.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalPromptRecoQAComponent)

AliHLTGlobalPromptRecoQAComponent::AliHLTGlobalPromptRecoQAComponent()
  : AliHLTProcessor()
  , fVerbosity(0)
  , fBenchmark("PromptRecoQA")
  , fpHWCFData(NULL)
  , fSkipEvents(0)
  , fPrintStats(0)
  , fEventsSinceSkip(0)
  , fHistSPDclusters_SPDrawSize(NULL)
  , fHistSSDclusters_SSDrawSize(NULL)
  , fHistSDDclusters_SDDrawSize(NULL)
  , fHistITSSAtracks_SPDclusters(NULL)
  , fHistSPDclusters_SSDclusters(NULL)
  , fHistSPDclusters_SDDclusters(NULL)
  , fHistSSDclusters_SDDclusters(NULL)
  , fHistTPCHLTclusters_TPCCompressionRatio(NULL)
  , fHistTPCHLTclusters_TPCFullCompressionRatio(NULL)
  , fHistHLTSize_HLTInOutRatio(NULL)
  , fHistTPCtracks_TPCtracklets(NULL)
  , fHistITStracks_ITSOutTracks(NULL)
  , fHistTPCClusterSize_TPCCompressedSize(NULL)
  , fHistTPCRawSize_TPCCompressedSize(NULL)
  , fHistHLTInSize_HLTOutSize(NULL)
  , fHistZNA_VZEROTrigChargeA(NULL)
  , fHistZNC_VZEROTrigChargeC(NULL)
  , fHistZNT_VZEROTrigChargeT(NULL)
  , fHistVZERO_SPDClusters(NULL)
  , fHistVZERO_ITSSAPTracks(NULL)
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
  delete fHistSPDclusters_SDDclusters;
  delete fHistSSDclusters_SDDclusters;
  delete fHistTPCHLTclusters_TPCCompressionRatio;
  delete fHistTPCHLTclusters_TPCFullCompressionRatio;
  delete fHistHLTSize_HLTInOutRatio;
  delete fHistTPCtracks_TPCtracklets;
  delete fHistITStracks_ITSOutTracks;
  delete fHistTPCClusterSize_TPCCompressedSize;
  delete fHistHLTInSize_HLTOutSize;
  delete fHistTPCRawSize_TPCCompressedSize;
  delete fHistZNA_VZEROTrigChargeA;
  delete fHistZNC_VZEROTrigChargeC;
  delete fHistZNT_VZEROTrigChargeT;
  delete fHistVZERO_SPDClusters;
  delete fHistVZERO_ITSSAPTracks;
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
      }	else if (argument.CompareTo("-print-stats-verbose")==0) {
        fPrintStats = 2;
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
  
  list.push_back(kAliHLTDataTypeESDContent | kAliHLTDataOriginVZERO); //VZERO-RECO output
  list.push_back(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginVZERO);
  list.push_back(kAliHLTDataTypeESDContent | kAliHLTDataOriginZDC); //ZDC-RECO output
  list.push_back(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginZDC);
  
  list.push_back(AliHLTEMCALDefinitions::fgkTriggerPatchDataType); //EMCAL-RECO output
  list.push_back(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginEMCAL);
  list.push_back(AliHLTEMCALDefinitions::fgkTriggerSTUDataType); //STU
  list.push_back(AliHLTEMCALDefinitions::fgkTriggerRawDigitDataType); //TRU
  
  list.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginOut);
  list.push_back(kAliHLTDataTypeESDfriendObject|kAliHLTDataOriginOut);
  list.push_back(kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut);
  list.push_back(kAliHLTDataTypeFlatESDFriend|kAliHLTDataOriginOut);

  //All this is TPC Data compression
  list.push_back(AliHLTTPCDefinitions::DataCompressionDescriptorDataType());
  list.push_back(AliHLTTPCDefinitions::RawClustersDataTypeNotCompressed());
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

  AliGRPManager mgr;
  mgr.ReadGRPEntry();
  float scaleSPDClusters, scaleSSDClusters, scaleSDDClusters, scaleTPCClusters, scaleITSSAPTracks, scaleTPCTracks, scaleTPCTracklets, scaleHLTIn, scaleHLTOut, scaleSPDSize,
    scaleSDDSize, scaleSSDSize, scaleTPCSize, scaleVZEROChargeA, scaleVZEROChargeC, scaleZDCCharge,
    scaleTPCCompressionRatio, scaleFullTPCCompressionRatio, scaleTPCCompressedSize;
  if (mgr.GetGRPData()->GetBeamType() == "Pb-Pb" || mgr.GetGRPData()->GetBeamType() == "PbPb" || mgr.GetGRPData()->GetBeamType() == "A-A" || mgr.GetGRPData()->GetBeamType() == "AA")
  {
    scaleSPDClusters = 5000.;
    scaleSSDClusters = 3500.;
    scaleSDDClusters = 2000.;
    scaleSPDSize = 25000.;
    scaleSSDSize = 80000.;
    scaleSDDSize = 80000.;
    scaleITSSAPTracks = 1000.;
    scaleTPCCompressionRatio = 8;
    scaleFullTPCCompressionRatio = 12;
    scaleTPCSize = 40000000;
    scaleTPCTracks = 30000;
    scaleTPCTracklets = 40000;
    scaleTPCCompressedSize = 6000000;
    scaleHLTIn = 60000000;
    scaleHLTOut = 7000000;
    scaleVZEROChargeA = 25000.;
    scaleVZEROChargeC = 50000.;
    scaleZDCCharge = 1500.;
    scaleTPCClusters = 5000000;
  }
  else
  {
    scaleSPDClusters = 500.;
    scaleSSDClusters = 600.;
    scaleSDDClusters = 500.;
    scaleSPDSize = 10000.;
    scaleSSDSize = 80000.;
    scaleSDDSize = 50000.;
    scaleITSSAPTracks = 100.;
    scaleTPCCompressionRatio = 7;
    scaleFullTPCCompressionRatio = 10;
    scaleTPCSize = 5000000;
    scaleTPCTracks = 700;
    scaleTPCTracklets = 1000;
    scaleTPCCompressedSize = 1500000;
    scaleHLTIn = 6000000;
    scaleHLTOut = 2000000;
    scaleVZEROChargeA = 2000.;
    scaleVZEROChargeC = 4000.;
    scaleZDCCharge = 500.;
    scaleTPCClusters = 20000;
  }

  fHistSPDclusters_SPDrawSize = new TH2I("SPDncls_SPDsize", "SPD clusters vs SPD raw size", 50, 4000., scaleSPDSize, 50, 0., scaleSPDClusters);
  fHistSSDclusters_SSDrawSize = new TH2I("SSDncls_SSDsize", "SSD clusters vs SSD raw size", 50, 20000., scaleSSDSize, 50, 0., scaleSSDClusters);
  fHistSDDclusters_SDDrawSize = new TH2I("SDDncls_SDDsize", "SDD clusters vs SDD raw size", 50, 0., scaleSDDSize, 50, 0., scaleSDDClusters);
  fHistITSSAtracks_SPDclusters = new TH2I("ITSSAPntrk_SPDncls", "ITSSAP tracks vs SPD clusters", 50, 0., scaleSPDClusters, 50, 0., scaleITSSAPTracks);
  fHistSPDclusters_SSDclusters = new TH2I("SSDncls_SPDncls", "SSD clusters vs SPD clusters", 50, 0., scaleSPDClusters, 50, 0., scaleSSDClusters);
  fHistSPDclusters_SDDclusters = new TH2I("SDDncls_SPDncls", "SDD clusters vs SPD clusters", 50, 0., scaleSPDClusters, 50, 0., scaleSDDClusters);
  fHistSSDclusters_SDDclusters = new TH2I("SDDncls_SSDncls", "SDD clusters vs SSD clusters", 50, 0., scaleSSDClusters, 50, 0., scaleSDDClusters);
  fHistTPCHLTclusters_TPCCompressionRatio = new TH2F("TPCHLTncls_TPCCompRatio", "Huffman compression ratio vs TPC HLT clusters", 50, 0., scaleTPCClusters, 50, 0., scaleTPCCompressionRatio);
  fHistTPCHLTclusters_TPCFullCompressionRatio = new TH2F("TPCHLTncls_TPCCompRatioFull", "Full compression ratio vs TPC HLT clusters", 50, 0., scaleTPCClusters, 50, 0., scaleFullTPCCompressionRatio);
  fHistHLTSize_HLTInOutRatio = new TH2F("TPCHLTSize_HLTInOutRatio", "HLT Out/In Size Ratio vs HLT Input Size", 50, 0., scaleHLTIn, 50, 0., 1.);
  fHistTPCClusterSize_TPCCompressedSize = new TH2I("TPCHWCFSize_TPCCompSize", "TPC compressed size vs TPC HWCF Size", 50, 0., scaleTPCSize, 50, 0., scaleTPCCompressedSize);
  fHistTPCRawSize_TPCCompressedSize = new TH2I("TPCRawSize_TPCCompSize", "TPC compressed size vs TPC Raw Size", 50, 0., scaleTPCSize, 50, 0., scaleTPCCompressedSize);
  fHistHLTInSize_HLTOutSize = new TH2I("HLTInSize_HLTOutSize", "HLT Out Size vs HLT In Size", 50, 0., scaleHLTIn, 50, 0., scaleHLTOut);
  fHistTPCtracks_TPCtracklets = new TH2I("TPCntrk_TPCntrl", "TPC Tracks vs TPC Tracklets", 50, 0., scaleTPCTracklets, 50, 0., scaleTPCTracks);
  fHistITStracks_ITSOutTracks = new TH2I("ITSntrk_ITSOutntrk", "ITS Tracks vs ITS Out Tracks", 50, 0., scaleTPCTracks, 50, 0., scaleTPCTracks);
  fHistZNA_VZEROTrigChargeA = new TH2F("VZEROTrigChargeA_ZNA", "ZNA vs. VZERO Trigger Charge A", 100, 0., scaleVZEROChargeA, 100, 0., scaleZDCCharge);
  fHistZNC_VZEROTrigChargeC = new TH2F("VZEROTrigChargeC_ZNC", "ZNC vs. VZERO Trigger Charge C", 100, 0., scaleVZEROChargeC, 100, 0., scaleZDCCharge);
  fHistZNT_VZEROTrigChargeT = new TH2F("VZEROTrigChargeT_ZNT", "ZN (A+C) vs. VZERO Trigger Charge (A+C)", 100, 0., scaleVZEROChargeA + scaleVZEROChargeC, 100, 0., 2 * scaleZDCCharge);
  fHistVZERO_SPDClusters = new TH2F("VZERO_SPDClusters", "SPD Clusters vs VZERO Multiplicity", 100, 0., scaleVZEROChargeA + scaleVZEROChargeC, 100, 0., scaleSPDClusters);
  fHistVZERO_ITSSAPTracks = new TH2F("VZERO_ITSSAPTracks", "ITS SAP Tracks vs VZERO Multiplicity", 100, 0., scaleVZEROChargeA + scaleVZEROChargeC, 100, 0., scaleITSSAPTracks);
  
  fpHWCFData = new AliHLTTPCHWCFData;

  return iResult;
}

int AliHLTGlobalPromptRecoQAComponent::DoDeinit()
{
  // see header file for class documentation
  
  delete fpHWCFData;
  return 0;
}

//Root cannot do templates...
//template <class T, class S, class U> void AliHLTGlobalPromptRecoQAComponent::FillHist(int check, T* hist, S val1, U val2, int& flag)

#define FillHist(check, hist, val1, val2, flag) \
{ \
  if (check) hist->Fill(val1, val2); \
  if (hist->GetEntries() && PushBack(hist, kAliHLTDataTypeHistogram|kAliHLTDataOriginOut) > 0) \
  { \
    flag = 1; \
    hist->Reset(); \
  } \
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

  if (!IsDataEvent()) {return iResult;}

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
  AliHLTUInt32_t rawSizeVZERO = 0;
  AliHLTUInt32_t rawSizeEMCAL = 0;
  AliHLTUInt32_t rawSizeZDC = 0;
  
  AliHLTUInt32_t nClustersTPC = 0;
  AliHLTUInt32_t rawSizeTPC = 0;
  AliHLTUInt32_t hwcfSizeTPC = 0;
  AliHLTUInt32_t clusterSizeTPCtransformed = 0;
  AliHLTUInt32_t clusterSizeTPC = 0;
  AliHLTUInt32_t compressedSizeTPC = 0;

  AliHLTUInt32_t nITSSAPtracks = 0;
  AliHLTUInt32_t nTPCtracklets = 0;
  AliHLTUInt32_t nTPCtracks = 0;
  AliHLTUInt32_t nITSTracks = 0;
  AliHLTUInt32_t nITSOutTracks = 0;

  float vZEROMultiplicity = 0.;
  UShort_t vZEROTriggerChargeA = 0;
  UShort_t vZEROTriggerChargeC = 0;

  Double_t zdcZNC = 0.;
  Double_t zdcZNA = 0.;

  int zdcRecoSize = 0;
  int emcalRecoSize = 0;
  int emcalTRU = 0;
  int emcalSTU = 0;
  
  static int nEvents = 0;
  
  Bool_t bITSSPDVertex = kFALSE;
  
  float compressionRatio, compressionRatioFull;
  
  AliHLTUInt32_t nESDSize = 0;
  AliHLTUInt32_t nESDFriendSize = 0;
  AliHLTUInt32_t nFlatESDSize = 0;
  AliHLTUInt32_t nFlatESDFriendSize = 0;
  
  AliHLTUInt32_t nHLTInSize = 0;
  AliHLTUInt32_t nHLTOutSize = 0;
  float hltRatio;

  nEvents++;

  //loop over input blocks and extract basic stats
  int nBlocks = evtData.fBlockCnt;  
  for (int ndx=0; ndx<nBlocks; ndx++) {
    const AliHLTComponentBlockData* iter = blocks+ndx;
    
    //SPD Vertex Found
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

    if (iter->fDataType == AliHLTTPCDefinitions::DataCompressionDescriptorDataType() ||
      iter->fDataType == AliHLTTPCDefinitions::RawClustersDataTypeNotCompressed() ||
      iter->fDataType == AliHLTTPCDefinitions::RemainingClustersCompressedDataType() ||
      iter->fDataType == AliHLTTPCDefinitions::RemainingClusterIdsDataType() ||
      iter->fDataType == AliHLTTPCDefinitions::ClusterTracksCompressedDataType() ||
      iter->fDataType == AliHLTTPCDefinitions::ClusterIdTracksDataType())
    {
      compressedSizeTPC += iter->fSize;
    }
    
    //VZERO Multiplicity
    if (iter->fDataType == (kAliHLTDataTypeESDContent | kAliHLTDataOriginVZERO))
    {
      const TObject* o = GetInputObjectFromIndex(ndx);
      if (o)
      {
        const AliESDVZERO* esdVZERO = dynamic_cast<const AliESDVZERO*>(o);
        if (esdVZERO)
        {
          for (int i = 0;i < 64;i++) vZEROMultiplicity += esdVZERO->GetMultiplicity(i);
          vZEROTriggerChargeA = esdVZERO->GetTriggerChargeA();
          vZEROTriggerChargeC = esdVZERO->GetTriggerChargeC();
        }
      }
    }

    //EMCAL Reco Size
    if (iter->fDataType == AliHLTEMCALDefinitions::fgkTriggerPatchDataType)
    {
      emcalRecoSize += iter->fSize;
    }
    if (iter->fDataType == AliHLTEMCALDefinitions::fgkTriggerSTUDataType)
    {
      emcalSTU += iter->fSize;
    }
    if (iter->fDataType == AliHLTEMCALDefinitions::fgkTriggerRawDigitDataType)
    {
      emcalTRU += iter->fSize;
    }

    //ZDC tower energies
    //ZDC Reco size
    if (iter->fDataType == (kAliHLTDataTypeESDContent | kAliHLTDataOriginZDC))
    {
      zdcRecoSize += iter->fSize;
      const TObject* o = GetInputObjectFromIndex(ndx);
      const AliESDZDC* esdZDC = dynamic_cast<const AliESDZDC*>(o);
      if (esdZDC)
      {
        zdcZNC = esdZDC->GetZNCTowerEnergy()[0];
        zdcZNA = esdZDC->GetZNATowerEnergy()[0];
      }
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
      AliHLTCDHWrapper header(iter->fPtr);
      AliHLTUInt8_t* pData = reinterpret_cast<AliHLTUInt8_t*>(iter->fPtr);
      pData+=header.GetHeaderSize();
      fpHWCFData->Init(pData, iter->fSize-header.GetHeaderSize());
      const AliHLTUInt32_t* pRCUTrailer = reinterpret_cast<const AliHLTUInt32_t*>(fpHWCFData->GetRCUTrailer());
      AliHLTUInt32_t payloadSize = (*pRCUTrailer) & 0x00ffffff;
      rawSizeTPC += header.GetHeaderSize() + payloadSize * sizeof(AliHLTUInt32_t) + fpHWCFData->GetRCUTrailerSize();
    }
    if (iter->fDataType == (AliHLTTPCDefinitions::fgkClustersDataType | kAliHLTDataOriginTPC)) //Size of transformed HLT-TPC clusters
    {
      clusterSizeTPCtransformed += iter->fSize;
    }
    if (iter->fDataType == (AliHLTTPCDefinitions::fgkRawClustersDataType | kAliHLTDataOriginTPC)) //Size of HLT-TPC clusters (uncompressed from HWCF, not yet transformed)
    {
      clusterSizeTPC += iter->fSize;
    }
    if (iter->fDataType == (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC)) //TPC RAW DDL Size
    {
      rawSizeTPC += iter->fSize;
    }
    if (iter->fDataType == (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginVZERO))
    {
      rawSizeVZERO += iter->fSize;
    }
    if (iter->fDataType == (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginZDC))
    {
      rawSizeZDC += iter->fSize;
    }
    if (iter->fDataType == (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginEMCAL))
    {
      rawSizeEMCAL += iter->fSize;
    }

    //esd size
    if (iter->fDataType == (kAliHLTDataTypeESDObject|kAliHLTDataOriginOut))
    {
      nESDSize += iter->fSize;
    }
    if (iter->fDataType == (kAliHLTDataTypeESDfriendObject|kAliHLTDataOriginOut))
    {
      nESDFriendSize += iter->fSize;
    }
    if (iter->fDataType == (kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut))
    {
      nFlatESDSize += iter->fSize;
    }
    if (iter->fDataType == (kAliHLTDataTypeFlatESDFriend|kAliHLTDataOriginOut))
    {
      nFlatESDFriendSize += iter->fSize;
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

  compressionRatio = compressedSizeTPC > 0 ? ((float) hwcfSizeTPC / (float) compressedSizeTPC) : 0.f;
  compressionRatioFull = compressedSizeTPC > 0 ? ((float) rawSizeTPC / (float) compressedSizeTPC) : 0.f;

  nHLTInSize = rawSizeSPD + rawSizeSSD + rawSizeSDD + rawSizeTPC + rawSizeVZERO + rawSizeEMCAL + rawSizeZDC;
  nHLTOutSize = nESDSize + nESDFriendSize + nFlatESDSize + nFlatESDFriendSize + compressedSizeTPC;
  hltRatio = nHLTInSize > 0 ? ((float) nHLTOutSize / (float) nHLTInSize) : 0.f;

  int pushed_something = 0;
  //fill histograms
  FillHist(nClustersSPD, fHistSPDclusters_SPDrawSize, rawSizeSPD, nClustersSPD, pushed_something);
  FillHist(nClustersSDD, fHistSDDclusters_SDDrawSize, rawSizeSDD, nClustersSDD, pushed_something);
  FillHist(nClustersSSD, fHistSSDclusters_SSDrawSize, rawSizeSSD, nClustersSSD, pushed_something);
  FillHist(nITSSAPtracks, fHistITSSAtracks_SPDclusters, nClustersSPD, nITSSAPtracks, pushed_something);
  FillHist(nClustersSPD, fHistSPDclusters_SSDclusters, nClustersSPD, nClustersSSD, pushed_something);
  FillHist(nClustersSPD, fHistSPDclusters_SDDclusters, nClustersSPD, nClustersSDD, pushed_something);
  FillHist(nClustersSSD, fHistSSDclusters_SDDclusters, nClustersSSD, nClustersSDD, pushed_something);
  FillHist(nTPCtracks, fHistTPCtracks_TPCtracklets, nTPCtracklets, nTPCtracks, pushed_something);

  //TPC Compression Plots
  FillHist(compressedSizeTPC, fHistTPCClusterSize_TPCCompressedSize, hwcfSizeTPC, compressedSizeTPC, pushed_something);
  FillHist(compressedSizeTPC, fHistTPCRawSize_TPCCompressedSize, rawSizeTPC, compressedSizeTPC, pushed_something);
  FillHist(nClustersTPC, fHistTPCHLTclusters_TPCCompressionRatio, nClustersTPC, compressionRatio, pushed_something);
  FillHist(nClustersTPC, fHistTPCHLTclusters_TPCFullCompressionRatio, nClustersTPC, compressionRatioFull, pushed_something);

  //HLT In Out Plots
  FillHist(nHLTInSize, fHistHLTInSize_HLTOutSize, nHLTInSize, nHLTOutSize, pushed_something);
  FillHist(nHLTInSize, fHistHLTSize_HLTInOutRatio, nHLTInSize, hltRatio, pushed_something);

  //VZERO vs ZDC
  FillHist(1, fHistZNA_VZEROTrigChargeA, vZEROTriggerChargeA, zdcZNA, pushed_something);
  FillHist(1, fHistZNC_VZEROTrigChargeC, vZEROTriggerChargeC, zdcZNC, pushed_something);
  FillHist(1, fHistZNT_VZEROTrigChargeT, vZEROTriggerChargeA+vZEROTriggerChargeC, zdcZNA+zdcZNC, pushed_something);

  //VZERO vs SPD
  FillHist(nClustersSPD, fHistVZERO_SPDClusters, vZEROTriggerChargeA+vZEROTriggerChargeC, nClustersSPD, pushed_something);
  FillHist(nITSSAPtracks, fHistVZERO_ITSSAPTracks, vZEROTriggerChargeA+vZEROTriggerChargeC, nITSSAPtracks, pushed_something);

  if (fPrintStats && (fPrintStats == 2 || pushed_something)) //Don't print this for every event if we use a pushback period
  {
    HLTImportant("Events %d Blocks %4d: HLT Reco QA Stats: HLTInOut %d / %d / %4.1f%%, SPD-Cl %d (%d), SDD-Cl %d (%d), SSD-Cl %d (%d) TPC-Cl %d (%d / %d / %d / %d), TPC-Comp %5.3fx / %5.3fx (%d)"
      ", ITSSAP-Tr %d, TPC-Tr %d / %d, ITS-Tr %d / %d, SPD-Ver %d, V0 %6.2f (%d), EMCAL %d (%d / %d / %d), ZDC %d (%d), ESD %d / %d (%d / %d)",
      nEvents, nBlocks, nHLTInSize, nHLTOutSize, hltRatio * 100, nClustersSPD, rawSizeSPD, nClustersSDD, rawSizeSDD, nClustersSSD, rawSizeSSD, nClustersTPC, rawSizeTPC, hwcfSizeTPC, clusterSizeTPC, clusterSizeTPCtransformed, compressionRatio, compressionRatioFull, compressedSizeTPC,
      nITSSAPtracks, nTPCtracklets, nTPCtracks, nITSTracks, nITSOutTracks, (int) bITSSPDVertex, vZEROMultiplicity, rawSizeVZERO, emcalRecoSize, emcalTRU, emcalSTU, rawSizeEMCAL, zdcRecoSize, rawSizeZDC, nESDSize, nFlatESDSize, nESDFriendSize, nFlatESDFriendSize);
  }

  return iResult;
}
