// $Id$

//**************************************************************************
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

//  @file   AliHLTTPCAgent.cxx
//  @author Matthias Richter
//  @date   
//  @brief  Agent of the libAliHLTTPC library
//  @note   

#include "AliHLTTPCAgent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTOUT.h"
#include "AliHLTOUTHandlerChain.h"
#include "AliHLTMisc.h"
#include "AliRunLoader.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTPCParam.h"
#include "AliTPCRecoParam.h"
#include "AliDAQ.h"
#include "TObject.h"
#include "AliHLTPluginBase.h"
#include "AliHLTSystem.h"

/** global instance for agent registration */
AliHLTTPCAgent gAliHLTTPCAgent;

// component headers
#include "AliHLTTPCCAInputDataCompressorComponent.h"
#include "AliHLTTPCCATrackerComponent.h"
#include "AliHLTTPCCATrackerOutputConverter.h"
#include "AliHLTTPCTrackMCMarkerComponent.h"
#include "AliHLTTPCCAGlobalMergerComponent.h"
#include "AliHLTTPCdEdxComponent.h"
#include "AliHLTTPCdEdxMonitoringComponent.h"
#include "AliHLTTPCClusterFinderComponent.h"
#include "AliHLTTPCDigitPublisherComponent.h"
#include "AliHLTTPCDigitDumpComponent.h"
#include "AliHLTTPCClusterDumpComponent.h"
#include "AliHLTTPCClusterHistoComponent.h"
#include "AliHLTTPCHistogramHandlerComponent.h"
#include "AliHLTTPCTrackHistoComponent.h"
#include "AliHLTTPCHWCFDataReverterComponent.h"
#include "AliHLTTPCCFComparisonComponent.h"
#include "AliHLTTPCDataCheckerComponent.h"
#include "AliHLTTPCHWCFEmulatorComponent.h"
#include "AliHLTTPCHWCFConsistencyControlComponent.h"
#include "AliHLTTPCClusterStatComponent.h"
#include "AliHLTTPCDataCompressionComponent.h"
#include "AliHLTTPCDataCompressionMonitorComponent.h"
#include "AliHLTTPCDataCompressionUnpackerComponent.h"
#include "AliHLTTPCDataCompressionFilterComponent.h"
#include "AliHLTTPCDataPublisherComponent.h"
#include "AliHLTTPCHWClusterDecoderComponent.h"
#include "AliHLTTPCClusterTransformationComponent.h"
#include "AliHLTTPCClusterTransformationPrepareComponent.h"
#include "AliHLTTPCOfflinePreprocessorWrapperComponent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCAgent)

AliHLTTPCAgent::AliHLTTPCAgent()
  : AliHLTModuleAgent("TPC")
  , fRawDataHandler(NULL)
  , fTracksegsDataHandler(NULL)
  , fClustersDataHandler(NULL)
  , fCompressionMonitorHandler(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCAgent::~AliHLTTPCAgent()
{
  // see header file for class documentation
}

UInt_t AliHLTTPCAgent::GetDetectorMask() const
{
  return AliDAQ::kTPC;
}

int AliHLTTPCAgent::CreateConfigurations(AliHLTConfigurationHandler* handler,
					 AliRawReader* rawReader,
					 AliRunLoader* runloader) const
{
  // see header file for class documentation
  if (handler) {
    // AliSimulation: use the AliRawReaderPublisher if the raw reader is available
    // AliReconstruction: indicated by runloader==NULL, run always on raw data
    bool bPublishRaw=rawReader!=NULL || runloader==NULL;

    AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

    bool isRawHLTOUT = 1;
    int tpcInputMode = 0;
    if( pHLT ){
      TString hltoptions = pHLT->GetConfigurationString();
      TObjArray* pTokens=hltoptions.Tokenize(" ");
      if( pTokens ){
	int iEntries=pTokens->GetEntriesFast();
	for (int i=0; i<iEntries; i++) {
	  if (!pTokens->At(i)) continue;
	  TString token = pTokens->At(i)->GetName();
	  if (token.CompareTo("ignore-hltout")==0) {
	    isRawHLTOUT = 0;
	  } else if (token.CompareTo("run-online-config")==0) {
	    isRawHLTOUT = 0;
	  } else if (token.Contains("TPC-input=")) {
	    TString param=token.ReplaceAll("TPC-input=", "");
	    if (param == "default") {
	      tpcInputMode = 0;
	    } else if (param == "raw") {
	      tpcInputMode = 1;
	    } else if (param == "compressed") {
	      tpcInputMode = 2;
	    } else {
	      HLTWarning("wrong parameter \'%s\' for option \'TPC-input=\', expected \'default\'/\'raw\'/\'compressed\'",param.Data() );
	    }
	  }
	}
	delete pTokens;
      }
    }

    // This the tracking configuration for the full TPC
    // - 216 clusterfinders (1 per partition)
    // - 36 slice trackers
    // - one global merger
    // - the esd converter
    // The ESD is shipped embedded into a TTree
    int iMinSlice=0; 
    int iMaxSlice=35;
    int iMinPart=0;
    int iMaxPart=5;
    TString arg;
    TString mergerInput;
    TString sinkRawData;
    TString sinkHWClusterInput;
    TString dEdXInput;
    TString hwcfemuInput;
    TString hwclustOutput;
    TString compressorInput;
    TString trackerInput;

    // Default TPC input: raw data or compressed clusters 

    if( isRawHLTOUT ){
      // Compressed clusters are already copied from the raw file to the HLTOUT.
      // The system should only reconstruct raw data for the sectors where the compressed clusters are missing
      arg.Form("-publish-clusters off -publish-raw filtered");
    } else {
      // HLTOUT is initially empty, therefore the HLT system reconstructs everything
      // The prefered input data is compressed clusters. TPC raw data is only published when the clusters are not present.
      arg.Form("-publish-clusters all -publish-raw filtered");
    }

    // Overwrite the default input when "TPC-input=.." option is set by user

    if( tpcInputMode==1 ){
      arg.Form("-publish-clusters off -publish-raw all");
    } else if ( tpcInputMode==2 ){
      arg.Form("-publish-clusters all -publish-raw off");
    }

    handler->CreateConfiguration("TPC-DP", "TPCDataPublisher", NULL , arg.Data());

    if (bPublishRaw) {
      hwcfemuInput = "TPC-DP";
    }

    for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
      for (int part=iMinPart; part<=iMaxPart; part++) {
	TString publisher;
	// digit publisher components
	publisher.Form("TPC-DP_%02d_%d", slice, part);
	if (bPublishRaw) {
	  // AliSimulation: use the AliRawReaderPublisher if the raw reader is available
	  // AliReconstruction: indicated by runloader==NULL, run always on raw data
	  int ddlno=768;
	  if (part>1) ddlno+=72+4*slice+(part-2);
	  else ddlno+=2*slice+part;

	  // arg.Form("-datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x", slice, slice, part, part);
	  arg.Form("-dataspec 0x%02x%02x%02x%02x", slice, slice, part, part); // publish also unpacked clusters, not only raw data

	  handler->CreateConfiguration(publisher.Data(), "BlockFilter", "TPC-DP" , arg.Data());
	  if (sinkRawData.Length()>0) sinkRawData+=" ";
	  sinkRawData+=publisher;
	} else {
	  arg.Form("-slice %d -partition %d", slice, part);
	  handler->CreateConfiguration(publisher.Data(), "TPCDigitPublisher", NULL , arg.Data());
	  if (hwcfemuInput.Length()>0) hwcfemuInput+=" ";
	  hwcfemuInput+=publisher;
	}
      }
    }

    // Hardware CF emulator
    TString hwcfemu;
    hwcfemu.Form("TPC-HWCFEmu");
    arg="";
    if (!bPublishRaw) arg+=" -do-mc 1";
    handler->CreateConfiguration(hwcfemu.Data(), "TPCHWClusterFinderEmulator", hwcfemuInput.Data(), arg.Data());
    if (hwclustOutput.Length()>0) hwclustOutput+=" ";
    hwclustOutput+=hwcfemu;
 
    TString hwcfDecoder = "TPC-HWCFDecoder";
    handler->CreateConfiguration(hwcfDecoder.Data(), "TPCHWClusterDecoder",hwclustOutput.Data(), "");

    TString clusterTransformation = "TPC-ClusterTransformation";
    handler->CreateConfiguration(clusterTransformation.Data(), "TPCClusterTransformation",hwcfDecoder.Data(), "-offline-mode");

    if (trackerInput.Length()>0) trackerInput+=" ";
    trackerInput+=clusterTransformation;
    if (dEdXInput.Length()>0) dEdXInput+=" ";
    dEdXInput+=clusterTransformation;
    if (sinkHWClusterInput.Length()>0) sinkHWClusterInput+=" ";
    sinkHWClusterInput+=clusterTransformation;

    // tracker finder component
    // 2012-01-05 changing the configuration according to online setup
    // the tracking strategy has been changed in the online system in Sep 2011
    // the tracker now processes all clusters, and several of this 'global' trackers
    // run in parallel. The GlobalMerger is still in the chain as it produces the final
    // fit.
    TString tracker;
    tracker.Form("TPC-TR");
    handler->CreateConfiguration(tracker.Data(), "TPCCATracker", trackerInput.Data(), "-GlobalTracking");

    if (mergerInput.Length()>0) mergerInput+=" ";
    mergerInput+=tracker;

    // GlobalMerger component
    handler->CreateConfiguration("TPC-globalmerger","TPCCAGlobalMerger",mergerInput.Data(),"");

    // dEdx component
    if (dEdXInput.Length()>0) dEdXInput+=" ";
    dEdXInput+="TPC-globalmerger";

    handler->CreateConfiguration("TPC-dEdx","TPCdEdx",dEdXInput.Data(),"");

    // compression component
    if (compressorInput.Length()>0) compressorInput+=" ";
    //compressorInput+=hwclustOutput;
    compressorInput+=hwcfDecoder;

    // special configuration to run the emulation automatically if the compressed clusters
    // of a particular partition is missing. This configuration is announced for reconstruction
    // of raw data if the HLT mode of the TPC reconstruction is enabled. Compression component
    // always needs to run in mode 1. Even if the recorded data is mode 3 (optimized partition
    // clusters), 2 (track model compression), or 4. The emulation can not be in mode 2 or 4,
    // since the track model block can not be identified with a partition. Have to duplicate the
    // configuration of the compression component
    handler->CreateConfiguration("TPC-compression-emulation", "TPCDataCompressor", compressorInput.Data(), "-mode 1");

    if (compressorInput.Length()>0) compressorInput+=" ";
    compressorInput+="TPC-globalmerger";
    handler->CreateConfiguration("TPC-compression", "TPCDataCompressor", compressorInput.Data(), "");
    handler->CreateConfiguration("TPC-compression-huffman-trainer", "TPCDataCompressor", compressorInput.Data(),"-deflater-mode 3");
    handler->CreateConfiguration("TPC-compression-monitoring-component", "TPCDataCompressorMonitor", "TPC-compression TPC-hwcfdata","-pushback-period=30");
    handler->CreateConfiguration("TPC-compression-monitoring", "ROOTFileWriter", "TPC-compression-monitoring-component","-concatenate-events -overwrite -datafile HLT.TPCDataCompression-statistics.root");

    // the esd converter configuration
    TString converterInput="TPC-globalmerger";
    if (!bPublishRaw) {
      // propagate cluster info to the esd converter in order to fill the MC information
      handler->CreateConfiguration("TPC-clustermc-info", "BlockFilter"   , sinkHWClusterInput.Data(), "-datatype 'CLMCINFO' 'TPC '");  
      handler->CreateConfiguration("TPC-mcTrackMarker","TPCTrackMCMarker","TPC-globalmerger TPC-clustermc-info","" );
      converterInput+=" ";
      converterInput+="TPC-mcTrackMarker";
    }
    handler->CreateConfiguration("TPC-esd-converter", "TPCEsdConverter"   , converterInput.Data(), "");

    // cluster dump collection
    handler->CreateConfiguration("TPC-clusters", "BlockFilter"   , sinkHWClusterInput.Data(), "-datatype 'CLUSTERS' 'TPC ' -datatype 'CLMCINFO' 'TPC '");
    handler->CreateConfiguration("TPC-raw-clusters", "BlockFilter"   , sinkHWClusterInput.Data(), "-datatype 'CLUSTRAW' 'TPC ' -datatype 'CLMCINFO' 'TPC '");
    handler->CreateConfiguration("TPC-hwclusters", "BlockFilter"   , sinkHWClusterInput.Data(), "-datatype 'CLUSTERS' 'TPC ' -datatype 'CLMCINFO' 'TPC '");
    handler->CreateConfiguration("TPC-raw-hwclusters", "BlockFilter"   , sinkHWClusterInput.Data(), "-datatype 'CLUSTRAW' 'TPC ' -datatype 'CLMCINFO' 'TPC '");

    // raw data
    handler->CreateConfiguration("TPC-raw-data", "BlockFilter"   , sinkRawData.Data(), "");

    handler->CreateConfiguration("TPC-hwcfdata", "BlockFilter"   , hwclustOutput.Data(), "-datatype 'HWCLUST1' 'TPC '");

    /////////////////////////////////////////////////////////////////////////////////////
    //
    // dumps on the ALTRO digit level
    //
    // selected channel dump
    arg.Form("-datafile selected-channel.dump -specfmt=_0x%%08x -subdir -blcknofmt= -idfmt=");
    handler->CreateConfiguration("TPC-selected-altro-digits", "TPCDigitDump", "RCU-channelselect", arg.Data());

    // raw channel dump
    arg.Form("-datafile channel.dump -specfmt=_0x%%08x -subdir -blcknofmt= -idfmt=");
    handler->CreateConfiguration("TPC-raw-altro-digits", "TPCDigitDump", "TPC-raw-data", arg.Data());

    /////////////////////////////////////////////////////////////////////////////////////
    //
    // a kChain HLTOUT configuration for processing of {'TRAKSEGS':'TPC '} data blocks
    // collects the data blocks, merges the tracks and produces an ESD object

    // publisher component
    handler->CreateConfiguration("TPC-hltout-tracksegs-publisher", "AliHLTOUTPublisher"   , NULL, "");

    // GlobalMerger component
    handler->CreateConfiguration("TPC-hltout-tracksegs-merger", "TPCGlobalMerger", "TPC-hltout-tracksegs-publisher", "");

    // the esd converter configuration
    handler->CreateConfiguration("TPC-hltout-tracksegs-esd-converter", "TPCEsdConverter", "TPC-hltout-tracksegs-merger", "");

    /////////////////////////////////////////////////////////////////////////////////////
    //
    // a kChain HLTOUT configuration for processing of {'TRACKS  ':'TPC '} data blocks
    // produces an ESD object from the track structure

    // publisher component
    handler->CreateConfiguration("TPC-hltout-tracks-publisher", "AliHLTOUTPublisher"   , NULL, "");

    // the esd converter configuration
    handler->CreateConfiguration("TPC-hltout-tracks-esd-converter", "TPCEsdConverter", "TPC-hltout-tracks-publisher", "");

    /////////////////////////////////////////////////////////////////////////////////////
    //
    // a kChain HLTOUT configuration for processing of {'CLUSTERS':'TPC '} data blocks
    // stores the blocks in file HLT.TPC.Clusters.root in HOMER format

    // publisher component
    handler->CreateConfiguration("TPC-hltout-cluster-publisher", "AliHLTOUTPublisher"   , NULL, "");

    // the HLTOUT component collects the blocks and stores the file
    handler->CreateConfiguration("TPC-hltout-cluster-dump", "HLTOUT", "TPC-hltout-cluster-publisher", "-digitfile HLT.TPC.Clusters.root -rawout=off -links 2");

    /////////////////////////////////////////////////////////////////////////////////////
    //
    // monitoring of compressed TPC data {CLUSTRAW:TPC }, {REMCLSCM,TPC }, {CLSTRKCM,TPC }
    // 

    // publisher component
    handler->CreateConfiguration("TPC-hltout-compressionmonitor-publisher", "AliHLTOUTPublisher"   , NULL,
				 "-datatype HWCLUST1 'TPC ' "
				 "-datatype CLUSTRAW 'TPC ' "
				 "-datatype REMCLSCM 'TPC ' "
				 "-datatype CLSTRKCM 'TPC ' "
				 "-datatype REMCLIDS 'TPC ' "
				 "-datatype CLIDSTRK 'TPC ' "
				 );

    // the HLTOUT component collects the blocks and stores the file
    handler->CreateConfiguration("TPC-hltout-compressionmonitor", "TPCDataCompressorMonitor", "TPC-hltout-compressionmonitor-publisher", "-histogram-file HLT.TPC-compression-statistics.root -publishing-mode off");
  }

  return 0;
}

const char* AliHLTTPCAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
						    AliRunLoader* runloader) const
{
  // see header file for class documentation
  if (runloader) {
    // reconstruction chains for AliRoot simulation
    return "TPC-compression";
  } else {
    return "TPC-compression-emulation";
  }
  return NULL;
}

const char* AliHLTTPCAgent::GetRequiredComponentLibraries() const
{
  // see header file for class documentation

  // actually, the TPC library has dependencies to Util and RCU
  // so the two has to be loaded anyhow before we get here
  //return "libAliHLTUtil.so libAliHLTRCU.so";
  return "libAliHLTUtil.so";
}

int AliHLTTPCAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
  // see header file for class documentation
  if (!pHandler) return -EINVAL;

  pHandler->AddComponent(new AliHLTTPCCAInputDataCompressorComponent);
  pHandler->AddComponent(new AliHLTTPCCATrackerComponent);
  pHandler->AddComponent(new AliHLTTPCCATrackerOutputConverter);
  pHandler->AddComponent(new AliHLTTPCCAGlobalMergerComponent);
  pHandler->AddComponent(new AliHLTTPCTrackMCMarkerComponent);
  pHandler->AddComponent(new AliHLTTPCdEdxComponent);
  pHandler->AddComponent(new AliHLTTPCdEdxMonitoringComponent);
  pHandler->AddComponent(new AliHLTTPCClusterFinderComponent(AliHLTTPCClusterFinderComponent::kClusterFinderPacked));
  pHandler->AddComponent(new AliHLTTPCClusterFinderComponent(AliHLTTPCClusterFinderComponent::kClusterFinderUnpacked));
  pHandler->AddComponent(new AliHLTTPCClusterFinderComponent(AliHLTTPCClusterFinderComponent::kClusterFinderDecoder));
  pHandler->AddComponent(new AliHLTTPCClusterFinderComponent(AliHLTTPCClusterFinderComponent::kClusterFinder32Bit));
  pHandler->AddComponent(new AliHLTTPCDigitPublisherComponent);
  pHandler->AddComponent(new AliHLTTPCDigitDumpComponent);
  pHandler->AddComponent(new AliHLTTPCClusterDumpComponent);
  pHandler->AddComponent(new AliHLTTPCClusterHistoComponent);
  pHandler->AddComponent(new AliHLTTPCHistogramHandlerComponent);
  pHandler->AddComponent(new AliHLTTPCTrackHistoComponent);
  pHandler->AddComponent(new AliHLTTPCHWCFDataReverterComponent);
  pHandler->AddComponent(new AliHLTTPCCFComparisonComponent);
  pHandler->AddComponent(new AliHLTTPCDataCheckerComponent);
  pHandler->AddComponent(new AliHLTTPCHWCFEmulatorComponent);
//  pHandler->AddComponent(new AliHLTTPCHWCFConsistencyControlComponent);  //FIXME: Causes crash: https://savannah.cern.ch/bugs/?83677
  pHandler->AddComponent(new AliHLTTPCClusterStatComponent);
  pHandler->AddComponent(new AliHLTTPCDataCompressionComponent);
  pHandler->AddComponent(new AliHLTTPCDataCompressionMonitorComponent);
  pHandler->AddComponent(new AliHLTTPCDataCompressionUnpackerComponent);
  pHandler->AddComponent(new AliHLTTPCDataCompressionFilterComponent);
  pHandler->AddComponent(new AliHLTTPCDataPublisherComponent);
  pHandler->AddComponent(new AliHLTTPCHWClusterDecoderComponent);
  pHandler->AddComponent(new AliHLTTPCClusterTransformationComponent);
  pHandler->AddComponent(new AliHLTTPCClusterTransformationPrepareComponent);
  pHandler->AddComponent(new AliHLTTPCOfflinePreprocessorWrapperComponent);
  return 0;
}

int AliHLTTPCAgent::GetHandlerDescription(AliHLTComponentDataType dt,
					  AliHLTUInt32_t spec,
					  AliHLTOUTHandlerDesc& desc) const
{
  // see header file for class documentation

  // raw data blocks to be fed into offline reconstruction
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC)) {
    int slice=AliHLTTPCDefinitions::GetMinSliceNr(spec);
    int part=AliHLTTPCDefinitions::GetMinPatchNr(spec);
    if (slice==AliHLTTPCDefinitions::GetMaxSliceNr(spec) &&
	part==AliHLTTPCDefinitions::GetMaxPatchNr(spec)) {
      desc=AliHLTOUTHandlerDesc(kRawReader, dt, GetModuleId());
      return 1;
    } else {
      HLTWarning("handler can not process merged data from multiple ddls:"
		 " min slice %d, max slice %d, min part %d, max part %d",
		 slice, AliHLTTPCDefinitions::GetMaxSliceNr(spec),
		 part, AliHLTTPCDefinitions::GetMaxPatchNr(spec));
      return 0;
    }
  }

  // dump for {'CLUSTERS':'TPC '} blocks stored in a 'digit' file
  if (dt==AliHLTTPCDefinitions::fgkClustersDataType) {
      desc=AliHLTOUTHandlerDesc(kChain, dt, GetModuleId());
      return 1;
  }

  // define handlers for all blocks related to compression, flag if the
  // cluster id blocks are existing, this will be used to decide
  // whether to create the handler or not
  // {'CLUSTRAW':'TPC '}
  // {'HWCLUST1':'TPC '}
  // {'REMCLSCM':'TPC '}
  // {'CLSTRKCM':'TPC '}
  // {'REMCLIDS':'TPC '}
  // {'CLIDSTRK':'TPC '}
  if (dt==AliHLTTPCDefinitions::RawClustersDataType() ||
      dt==AliHLTTPCDefinitions::RawClustersDataTypeNotCompressed() ||
      dt==AliHLTTPCDefinitions::HWClustersDataType() ||
      dt==AliHLTTPCDefinitions::RemainingClustersCompressedDataType() ||
      dt==AliHLTTPCDefinitions::ClusterTracksCompressedDataType()) {
      desc=AliHLTOUTHandlerDesc(kProprietary, dt, GetModuleId());
      return 1;
  }
  if (dt==AliHLTTPCDefinitions::RemainingClusterIdsDataType() ||
      dt==AliHLTTPCDefinitions::ClusterIdTracksDataType()) {
      desc=AliHLTOUTHandlerDesc(kProprietary, dt, GetModuleId());
      const_cast<AliHLTTPCAgent*>(this)->SetBit(kHaveCompressedClusterIdDataBlock);
      return 1;
  }

  // {'CLMCINFO':'TPC '} 
  if (dt==AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo) {
      desc=AliHLTOUTHandlerDesc(kProprietary, dt, GetModuleId());
      return 1;
  }

  // afterburner for {'TRAKSEGS':'TPC '} blocks to be converted to ESD format
  if (dt==AliHLTTPCDefinitions::fgkTrackSegmentsDataType) {
      desc=AliHLTOUTHandlerDesc(kChain, dt, GetModuleId());
      return 1;
  }

  // afterburner for {'TRACKS  ':'TPC '} block to be converted to ESD format
  // there is only one data block
  if (dt==AliHLTTPCDefinitions::fgkTracksDataType) {
      desc=AliHLTOUTHandlerDesc(kChain, dt, GetModuleId());
      return 1;
  }
  return 0;
}

AliHLTOUTHandler* AliHLTTPCAgent::GetOutputHandler(AliHLTComponentDataType dt,
						   AliHLTUInt32_t /*spec*/)
{
  // see header file for class documentation

  // raw data blocks to be fed into offline reconstruction
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC)) {
    if (!fRawDataHandler) {
      fRawDataHandler=new AliHLTTPCAgent::AliHLTTPCRawDataHandler;
    }
    return fRawDataHandler;
  }

  // dump for {'CLUSTERS':'TPC '}, stored in a file HLT.TPC.Clusters.root in HOMER format
  if (dt==AliHLTTPCDefinitions::fgkClustersDataType) {
    if (fClustersDataHandler==NULL)
      fClustersDataHandler=new AliHLTOUTHandlerChain("chains=TPC-hltout-cluster-dump libHLTsim.so libAliHLTUtil.so");
    return fClustersDataHandler;
  }

  // afterburner for {'TRAKSEGS':'TPC '} blocks to be converted to ESD format
  // in a kChain HLTOUT handler
  if (dt==AliHLTTPCDefinitions::fgkTrackSegmentsDataType) {
    if (fTracksegsDataHandler==NULL)
      fTracksegsDataHandler=new AliHLTOUTHandlerChain("chains=TPC-hltout-tracksegs-esd-converter");
    return fTracksegsDataHandler;
  }

  // afterburner for {'TRACKS  ':'TPC '} block to be converted to ESD format
  // there is only one data block
  if (dt==AliHLTTPCDefinitions::fgkTracksDataType) {
    return new AliHLTOUTHandlerChain("chains=TPC-hltout-tracks-esd-converter");
  }

  // monitoring of compressed data if cluster verification blocks exist
  // {'REMCLIDS':'TPC '}
  // {'CLIDSTRK':'TPC '}
  // FIXME: needs to be commissioned
  // if (dt==AliHLTTPCDefinitions::RawClustersDataType() ||
  //     dt==AliHLTTPCDefinitions::RawClustersDataTypeNotCompressed() ||
  //     dt==AliHLTTPCDefinitions::HWClustersDataType() ||
  //     dt==AliHLTTPCDefinitions::RemainingClustersCompressedDataType() ||
  //     dt==AliHLTTPCDefinitions::ClusterTracksCompressedDataType() ||
  //     dt==AliHLTTPCDefinitions::RemainingClusterIdsDataType() ||
  //     dt==AliHLTTPCDefinitions::ClusterIdTracksDataType()) {
  //   const char* arg="chains=TPC-hltout-compressionmonitor";
  //   if (!TestBit(kHaveCompressedClusterIdDataBlock))
  //     arg="chains=TPC-hltout-compressionmonitorpublisher";
  //   if (!fCompressionMonitorHandler)
  //     fCompressionMonitorHandler=new AliHLTOUTHandlerChain(arg);
  //   return fCompressionMonitorHandler;
  // }

  return NULL;
}

int AliHLTTPCAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance)
{
  // see header file for class documentation
  if (pInstance==NULL) return -EINVAL;

  if (pInstance==fRawDataHandler) {
    delete fRawDataHandler;
    fRawDataHandler=NULL;
  }

  if (pInstance==fTracksegsDataHandler) {
    delete fTracksegsDataHandler;
    fTracksegsDataHandler=NULL;
  }

  if (pInstance==fClustersDataHandler) {
    delete fClustersDataHandler;
    fClustersDataHandler=NULL;
  }

  if (pInstance==fCompressionMonitorHandler) {
    delete fCompressionMonitorHandler;
    fCompressionMonitorHandler=NULL;
  }

  return 0;
}

AliHLTTPCAgent::AliHLTTPCRawDataHandler::AliHLTTPCRawDataHandler()
{
  // see header file for class documentation
}

AliHLTTPCAgent::AliHLTTPCRawDataHandler::~AliHLTTPCRawDataHandler()
{
  // see header file for class documentation
}

int AliHLTTPCAgent::AliHLTTPCRawDataHandler::ProcessData(AliHLTOUT* pData)
{
  // see header file for class documentation
  if (!pData) return -EINVAL;
  AliHLTComponentDataType dt=kAliHLTVoidDataType;
  AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
  int iResult=pData->GetDataBlockDescription(dt, spec);
  if (iResult>=0) {
    int slice=AliHLTTPCDefinitions::GetMinSliceNr(spec);
    int part=AliHLTTPCDefinitions::GetMinPatchNr(spec);
    if (slice==AliHLTTPCDefinitions::GetMaxSliceNr(spec) &&
	part==AliHLTTPCDefinitions::GetMaxPatchNr(spec)) {
      iResult=768;
      if (part>1) iResult+=72+4*slice+(part-2);
      else iResult+=2*slice+part;
    } else {
      HLTError("handler can not process merged data from multiple ddls:"
	       " min slice %d, max slice %d, min part %d, max part %d",
	       slice, AliHLTTPCDefinitions::GetMaxSliceNr(spec),
	       part, AliHLTTPCDefinitions::GetMaxPatchNr(spec));
      iResult=-EBADMSG;
    }
  }
  return iResult;
}
