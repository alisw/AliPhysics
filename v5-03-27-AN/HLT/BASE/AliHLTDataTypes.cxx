// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
 *                  Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTDataTypes.cxx
    @author Matthias Richter, Timm Steinbeck, Jochen Thaeder
    @date   
    @brief  Implementation of data types. */

// those types can not be implemented in the header files as rootcint
// can not cope with the type id and origin defines.
//
// change Aug 01 2008
// some compilers can not cope with the fomerly used initialization of the
// default data type variable by using the operator | like e.g
//   const AliHLTComponentDataType kAliHLTDataTypeComponentTable = (AliHLTComponentDataType) {
//     sizeof(AliHLTComponentDataType),
//     kAliHLTComponentTableDataTypeID,
//     kAliHLTDataOriginAny
//   }|kAliHLTDataOriginPrivate;
// Mainly the compined type cast and utilization of the operator| is the problem.
// An initializer function has been defined in order to work around this issue.

#include "AliHLTDataTypes.h"

/** multiple output data types */
const char kAliHLTMultipleDataTypeIDstring[8] = {'M','U','L','T','I','P','L','E'};
const AliHLTComponentDataType kAliHLTMultipleDataType =  AliHLTComponentDataTypeInitializer(kAliHLTMultipleDataTypeIDstring, kAliHLTDataOriginPrivate);

/** data to file exchange subscriber */
const char kAliHLTFXSCalibDataTypeIDstring[8] = kAliHLTFXSCalibDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeFXSCalib = AliHLTComponentDataTypeInitializer(kAliHLTFXSCalibDataTypeIDstring, kAliHLTDataOriginOut);

/** DDL list data type */
const char kAliHLTDDLDataTypeIDstring[8] = kAliHLTDDLDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeDDL = AliHLTComponentDataTypeInitializer(kAliHLTDDLDataTypeIDstring, kAliHLTDataOriginOut);

/** DAQ readout list */
const AliHLTComponentDataType kAliHLTDataTypeDAQRDOUT = AliHLTComponentDataTypeInitializer(kAliHLTDAQRDOUTDataTypeID, kAliHLTDataOriginAny);

/** SOR data type */
const char kAliHLTSORDataTypeIDstring[8] = kAliHLTSORDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeSOR = AliHLTComponentDataTypeInitializer(kAliHLTSORDataTypeIDstring, kAliHLTDataOriginPrivate);

/** EOR data type */
const char kAliHLTEORDataTypeIDstring[8] = kAliHLTEORDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeEOR = AliHLTComponentDataTypeInitializer(kAliHLTEORDataTypeIDstring, kAliHLTDataOriginPrivate);

/** run type data block */
const char kAliHLTRunTypeDataTypeIDstring[8] = kAliHLTRunTypeDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeRunType = AliHLTComponentDataTypeInitializer(kAliHLTRunTypeDataTypeIDstring, kAliHLTDataOriginPrivate);

/** Event type specification */
const char kAliHLTEventDataTypeIDstring[8] = kAliHLTEventDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeEvent = AliHLTComponentDataTypeInitializer(kAliHLTEventDataTypeIDstring, kAliHLTDataOriginPrivate);

/** ECS parameter string */
const char kAliHLTECSParamDataTypeIDstring[8] = kAliHLTECSParamDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeECSParam = AliHLTComponentDataTypeInitializer(kAliHLTECSParamDataTypeIDstring, kAliHLTDataOriginPrivate);

/** Configuration event data type */
const char kAliHLTComConfDataTypeIDstring[8] = kAliHLTComConfDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeComConf = AliHLTComponentDataTypeInitializer(kAliHLTComConfDataTypeIDstring, kAliHLTDataOriginPrivate);

/** DCS value update event */
const char kAliHLTUpdtDCSDataTypeIDstring[8] = kAliHLTUpdtDCSDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeUpdtDCS = AliHLTComponentDataTypeInitializer(kAliHLTUpdtDCSDataTypeIDstring, kAliHLTDataOriginPrivate);

/** RAW DDL data specification, data publisher will set type id and origin correctly */
const char kAliHLTDDLRawDataTypeIDstring[8] = kAliHLTDDLRawDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeDDLRaw = AliHLTComponentDataTypeInitializer(kAliHLTDDLRawDataTypeIDstring, kAliHLTDataOriginAny);

/** CLUSTERS data type */
const char kAliHLTClustersDataTypeIDstring[8] = kAliHLTClustersDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeClusters = AliHLTComponentDataTypeInitializer(kAliHLTClustersDataTypeIDstring, kAliHLTDataOriginAny);

/** MC data specification */
const char kAliHLTMCObjectDataTypeIDstring[8] = kAliHLTMCObjectDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeMCObject = AliHLTComponentDataTypeInitializer(kAliHLTMCObjectDataTypeIDstring, kAliHLTDataOriginOffline);

/** ESD vertex data specification */
const char kAliHLTESDVertexDataTypeIDstring[8] = kAliHLTESDVertexDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeESDVertex = AliHLTComponentDataTypeInitializer(kAliHLTESDVertexDataTypeIDstring, kAliHLTDataOriginAny);

/** KF vertex data specification */
const char kAliHLTKFVertexDataTypeIDstring[8] = kAliHLTKFVertexDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeKFVertex = AliHLTComponentDataTypeInitializer(kAliHLTKFVertexDataTypeIDstring, kAliHLTDataOriginAny);

/** Global Vertexer data specification */
const char kAliHLTDataTypeGlobalVertexerIDstring[8] = kAliHLTDataTypeGlobalVertexerID;
const AliHLTComponentDataType kAliHLTDataTypeGlobalVertexer = AliHLTComponentDataTypeInitializer(kAliHLTDataTypeGlobalVertexerIDstring, kAliHLTDataOriginAny);

/** Primary finder data specification */
const char kAliHLTPrimaryFinderDataTypeIDstring[8] = kAliHLTDataTypePrimaryFinderID;
const AliHLTComponentDataType kAliHLTDataTypePrimaryFinder = AliHLTComponentDataTypeInitializer(kAliHLTPrimaryFinderDataTypeIDstring, kAliHLTDataOriginAny);

/** V0 finder data specification */
const char kAliHLTV0FinderDataTypeIDstring[8] = kAliHLTDataTypeV0FinderID;
const AliHLTComponentDataType kAliHLTDataTypeV0Finder = AliHLTComponentDataTypeInitializer(kAliHLTV0FinderDataTypeIDstring, kAliHLTDataOriginAny);

/** ESD data specification */
const char kAliHLTESDObjectDataTypeIDstring[8] = kAliHLTESDObjectDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeESDObject = AliHLTComponentDataTypeInitializer(kAliHLTESDObjectDataTypeIDstring, kAliHLTDataOriginAny);

/** ESD content specification */
const char kAliHLTESDContentDataTypeIDstring[8] = kAliHLTESDContentDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeESDContent = AliHLTComponentDataTypeInitializer(kAliHLTESDContentDataTypeIDstring, kAliHLTDataOriginAny);

/** ESD tree data specification */
const char kAliHLTESDTreeDataTypeIDstring[8] = kAliHLTESDTreeDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeESDTree = AliHLTComponentDataTypeInitializer(kAliHLTESDTreeDataTypeIDstring, kAliHLTDataOriginAny);

/** AliRoot TreeD data specification */
const char kAliHLTTreeDDataTypeIDstring[8] = kAliHLTTreeDDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeAliTreeD = AliHLTComponentDataTypeInitializer(kAliHLTTreeDDataTypeIDstring, kAliHLTDataOriginAny);

/** AliRoot TreeR data specification */
const char kAliHLTTreeRDataTypeIDstring[8] = kAliHLTTreeRDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeAliTreeR = AliHLTComponentDataTypeInitializer(kAliHLTTreeRDataTypeIDstring, kAliHLTDataOriginAny);

/** 16 bit Hardware address selection data specification, origin is 'any' */
const char kAliHLTHwAddr16DataTypeIDstring[8] = kAliHLTHwAddr16DataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeHwAddr16 = AliHLTComponentDataTypeInitializer(kAliHLTHwAddr16DataTypeIDstring, kAliHLTDataOriginAny);

/** Event statistics */
const char kAliHLTEventStatisticsDataTypeIDstring[8] = kAliHLTEventStatisticsDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeEventStatistics = AliHLTComponentDataTypeInitializer(kAliHLTEventStatisticsDataTypeIDstring, kAliHLTDataOriginAny);

/** Event summary */
const char kAliHLTEventSummaryDataTypeIDstring[8] = kAliHLTEventSummaryDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeEventSummary = AliHLTComponentDataTypeInitializer(kAliHLTEventSummaryDataTypeIDstring, kAliHLTDataOriginOut);

/** Run statistics */
const char kAliHLTRunStatisticsDataTypeIDstring[8] = kAliHLTRunStatisticsDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeRunStatistics = AliHLTComponentDataTypeInitializer(kAliHLTRunStatisticsDataTypeIDstring, kAliHLTDataOriginAny);

/** Run summary */
const char kAliHLTRunSummaryDataTypeIDstring[8] = kAliHLTRunSummaryDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeRunSummary = AliHLTComponentDataTypeInitializer(kAliHLTRunSummaryDataTypeIDstring, kAliHLTDataOriginOut);

/** Trigger decision */
const char kAliHLTTriggerDecisionDataTypeIDstring[8] = kAliHLTTriggerDecisionDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeTriggerDecision = AliHLTComponentDataTypeInitializer(kAliHLTTriggerDecisionDataTypeIDstring, kAliHLTDataOriginOut);

/** HLT readout list from trigger component */
const char kAliHLTReadoutListDataTypeIDstring[8] = kAliHLTReadoutListDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeReadoutList = AliHLTComponentDataTypeInitializer(kAliHLTReadoutListDataTypeIDstring, kAliHLTDataOriginOut);

/** Global trigger decision */
const char kAliHLTGlobalTriggerDataTypeIDstring[8] = kAliHLTGlobalTriggerDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeGlobalTrigger = AliHLTComponentDataTypeInitializer(kAliHLTGlobalTriggerDataTypeIDstring, kAliHLTDataOriginOut);

/** Component statistics */
const char  kAliHLTComponentStatisticsDataTypeIDstring[8] = kAliHLTComponentStatisticsDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeComponentStatistics = AliHLTComponentDataTypeInitializer(kAliHLTComponentStatisticsDataTypeIDstring, kAliHLTDataOriginPrivate);

/** Component table */
const char kAliHLTComponentTableDataTypeIDstring[8] = kAliHLTComponentTableDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeComponentTable = AliHLTComponentDataTypeInitializer(kAliHLTComponentTableDataTypeIDstring, kAliHLTDataOriginPrivate);

/** Forwarded component table */
const char kAliHLTComponentFwdTableDataTypeIDstring[8] = kAliHLTComponentFwdTableDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeComponentFwdTable = AliHLTComponentDataTypeInitializer(kAliHLTComponentFwdTableDataTypeIDstring, kAliHLTDataOriginPrivate);

/** general ROOT TObject */
const char kAliHLTTObjectDataTypeIDstring[8] = kAliHLTTObjectDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeTObject = AliHLTComponentDataTypeInitializer(kAliHLTTObjectDataTypeIDstring, kAliHLTDataOriginAny);

/** ROOT streamer info */
const char kAliHLTStreamerInfoDataTypeIDstring[8] = kAliHLTStreamerInfoDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeStreamerInfo = AliHLTComponentDataTypeInitializer(kAliHLTStreamerInfoDataTypeIDstring, kAliHLTDataOriginHLT);

/** ROOT TObjArray */
const char kAliHLTTObjArrayDataTypeIDstring[8] = kAliHLTTObjArrayDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeTObjArray = AliHLTComponentDataTypeInitializer(kAliHLTTObjArrayDataTypeIDstring, kAliHLTDataOriginAny);

/** ROOT TTree */
const char kAliHLTTTreeDataTypeIDstring[8] = kAliHLTTTreeDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeTTree = AliHLTComponentDataTypeInitializer(kAliHLTTTreeDataTypeIDstring, kAliHLTDataOriginAny);

/** ROOT TH1 (can be used for all histograms, they derive from TH1) */
const char kAliHLTHistogramDataTypeIDstring[8] = kAliHLTHistogramDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeHistogram = AliHLTComponentDataTypeInitializer(kAliHLTHistogramDataTypeIDstring, kAliHLTDataOriginAny);

/** ROOT TNtuple */
const char kAliHLTTNtupleDataTypeIDstring[8] = kAliHLTTNtupleDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeTNtuple = AliHLTComponentDataTypeInitializer(kAliHLTTNtupleDataTypeIDstring, kAliHLTDataOriginAny);

/** Array of HLT Tracks (AliHLTTracksData) */
const char kAliHLTTrackDataTypeIDstring[8] = kAliHLTTrackDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeTrack = AliHLTComponentDataTypeInitializer(kAliHLTTrackDataTypeIDstring, kAliHLTDataOriginAny);

/** Array of Track MC ids */
const char kAliHLTTrackMCDataTypeIDstring[8] = kAliHLTTrackMCDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeTrackMC = AliHLTComponentDataTypeInitializer(kAliHLTTrackMCDataTypeIDstring, kAliHLTDataOriginAny);

/** TClonesArray of AliExternalTrackParam */
const char kAliHLTExternalTrackParamDataTypeIDstring[8] = kAliHLTExternalTrackParamDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeExternalTrackParam = AliHLTComponentDataTypeInitializer(kAliHLTExternalTrackParamDataTypeIDstring, kAliHLTDataOriginAny);

/** Container of HLT Jets (AliHLTJETJets) */
const char kAliHLTJetDataTypeIDstring[8] = kAliHLTJetDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeJet = AliHLTComponentDataTypeInitializer(kAliHLTJetDataTypeIDstring, kAliHLTDataOriginAny);

/** Container of HLT ITS tracks */
const AliHLTComponentDataType fgkITSTracksDataType = AliHLTComponentDataTypeInitializer( "ITSTRACK", kAliHLTDataOriginITS );

/** Container of HLT calorimeter clusters */
const AliHLTComponentDataType kAliHLTDataTypeCaloCluster = AliHLTComponentDataTypeInitializer( "CALOCLUS", kAliHLTDataOriginAny );

/** Container of dEdx */
const AliHLTComponentDataType kAliHLTDataTypedEdx = AliHLTComponentDataTypeInitializer( "DEDX    ", kAliHLTDataOriginAny );

/** Container of dNdPt */
const AliHLTComponentDataType kAliHLTDataTypedNdPt = AliHLTComponentDataTypeInitializer( "DNDPT   ", kAliHLTDataOriginAny );

/** Input trigger counters */
const char kAliHLTInputTriggerCountersDataTypeIDstring[8] = kAliHLTInputTriggerCountersDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeInputTriggerCounters = AliHLTComponentDataTypeInitializer(kAliHLTInputTriggerCountersDataTypeIDstring, kAliHLTDataOriginHLT);

/** Input trigger counters */
const char kAliHLTOutputTriggerCountersDataTypeIDstring[8] = kAliHLTOutputTriggerCountersDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeOutputTriggerCounters = AliHLTComponentDataTypeInitializer(kAliHLTOutputTriggerCountersDataTypeIDstring, kAliHLTDataOriginHLT);

/** Meta data block for the Common Data Header (CDH) and readout list forwarded by TCPDumpSubscriber. */
const char kAliHLTMetaDataTypeIDstring[8] = kAliHLTMetaDataTypeID;
const AliHLTComponentDataType kAliHLTDataTypeTriggerMetaBlock = AliHLTComponentDataTypeInitializer(kAliHLTMetaDataTypeIDstring, kAliHLTDataOriginPrivate);

//////////////////////////////////////////////////////////////////////////
//
// Data origin variables, to be used with the operator|
//
// AliHLTComponentDataType dt;
// dt = kAliHLTDataTypeDDLRaw | gkAliHLTDataOriginTPC;
//
//////////////////////////////////////////////////////////////////////////

/** HLT out */
const char kAliHLTDataOriginOut[kAliHLTComponentDataTypefOriginSize]     = {'H','L','T',' '};

/** HLT */
const char kAliHLTDataOriginHLT[kAliHLTComponentDataTypefOriginSize]     = {'H','L','T',' '};

/** Offline */
const char kAliHLTDataOriginOffline[kAliHLTComponentDataTypefOriginSize] = {'O','F','F','L'};

/** HLT/PubSub private internal */
const char kAliHLTDataOriginPrivate[kAliHLTComponentDataTypefOriginSize] = {'P','R','I','V'};

/** TPC */
const char kAliHLTDataOriginTPC[kAliHLTComponentDataTypefOriginSize]     = {'T','P','C',' '};

/** PHOS */
const char kAliHLTDataOriginPHOS[kAliHLTComponentDataTypefOriginSize]    = {'P','H','O','S'};

/** FMD */
const char kAliHLTDataOriginFMD[kAliHLTComponentDataTypefOriginSize]     = {'F','M','D',' '};

/** MUON */
const char kAliHLTDataOriginMUON[kAliHLTComponentDataTypefOriginSize]    = {'M','U','O','N'};

/** TRD */
const char kAliHLTDataOriginTRD[kAliHLTComponentDataTypefOriginSize]     = {'T','R','D',' '};

/** ITS */
const char kAliHLTDataOriginITS[kAliHLTComponentDataTypefOriginSize]     = {'I','T','S',' '};

/** ITSOut */
const char kAliHLTDataOriginITSOut[kAliHLTComponentDataTypefOriginSize]     = {'I','T','S','O'};

/** ITS-SPD */
const char kAliHLTDataOriginITSSPD[kAliHLTComponentDataTypefOriginSize]  = {'I','S','P','D'};

/** ITS-SDD */
const char kAliHLTDataOriginITSSDD[kAliHLTComponentDataTypefOriginSize]  = {'I','S','D','D'};

/** ITS-SSD */
const char kAliHLTDataOriginITSSSD[kAliHLTComponentDataTypefOriginSize]  = {'I','S','S','D'};

/** Sample */
const char kAliHLTDataOriginSample[kAliHLTComponentDataTypefOriginSize]  = {'S','M','P','L'};

/** EMCAL */
const char kAliHLTDataOriginEMCAL[kAliHLTComponentDataTypefOriginSize]   = {'E','M','C','A'};

/** TOF */
const char kAliHLTDataOriginTOF[kAliHLTComponentDataTypefOriginSize]   = {'T','O','F',' '};

/** HMPID */
const char kAliHLTDataOriginHMPID[kAliHLTComponentDataTypefOriginSize]   = {'H','M','P','I'};

/** CPV */
const char kAliHLTDataOriginCPV[kAliHLTComponentDataTypefOriginSize]   = {'C','P','V',' '};

/** PMD */
const char kAliHLTDataOriginPMD[kAliHLTComponentDataTypefOriginSize]   = {'P','M','D',' '};

/** T0 */
const char kAliHLTDataOriginT0[kAliHLTComponentDataTypefOriginSize]   = {'T','Z','R','O'};

/** VZERO */
const char kAliHLTDataOriginVZERO[kAliHLTComponentDataTypefOriginSize]   = {'V','Z','R','O'};

/** ZDC */
const char kAliHLTDataOriginZDC[kAliHLTComponentDataTypefOriginSize]   = {'Z','D','C',' '};

/** ACORDE */
const char kAliHLTDataOriginACORDE[kAliHLTComponentDataTypefOriginSize]   = {'A','C','O','R'};

/** TRG */
const char kAliHLTDataOriginTRG[kAliHLTComponentDataTypefOriginSize]   = {'T','R','G',' '};
