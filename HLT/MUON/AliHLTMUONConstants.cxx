/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/**
 * @file   AliHLTMUONConstants.cxx
 * @author Indranil Das <indra.das@saha.ac.in>,
 *         Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Definitions of the various dimuon HLT constants.
 */

#include "AliHLTMUONConstants.h"
#include "AliHLTMUONTriggerRecordsBlockStruct.h"
#include "AliHLTMUONTrigRecsDebugBlockStruct.h"
#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTMUONClustersBlockStruct.h"
#include "AliHLTMUONChannelsBlockStruct.h"
#include "AliHLTMUONMansoTracksBlockStruct.h"
#include "AliHLTMUONMansoCandidatesBlockStruct.h"
#include "AliHLTMUONSinglesDecisionBlockStruct.h"
#include "AliHLTMUONPairsDecisionBlockStruct.h"

ClassImp(AliHLTMUONConstants);

const AliHLTMUONTriggerRecordStruct
AliHLTMUONConstants::fgkNilTriggerRecordStruct = {
	0, 0, 0, 0, 0,
	{
	 AliHLTMUONConstants::fgkNilRecHitStruct,
	 AliHLTMUONConstants::fgkNilRecHitStruct,
	 AliHLTMUONConstants::fgkNilRecHitStruct,
	 AliHLTMUONConstants::fgkNilRecHitStruct
	}
};

const AliHLTMUONTrigRecInfoStruct
AliHLTMUONConstants::fgkNilTrigRecInfoStruct = {0, {0, 0, 0, 0}, 0, 0, {0, 0, 0, 0, 0}};
	
const AliHLTMUONRecHitStruct
AliHLTMUONConstants::fgkNilRecHitStruct = {0, 0, 0, 0};

const AliHLTMUONClusterStruct
AliHLTMUONConstants::fgkNilClusterStruct = {
	0, AliHLTMUONConstants::fgkNilRecHitStruct, 0, 0, 0
};

const AliHLTMUONChannelStruct
AliHLTMUONConstants::fgkNilChannelStruct = {0, 0, 0, 0, 0, 0};

const AliHLTMUONMansoTrackStruct
AliHLTMUONConstants::fgkNilMansoTrackStruct = {
	0, 0, 0, 0, 0, 0, 0,
	{
	 AliHLTMUONConstants::fgkNilRecHitStruct,
	 AliHLTMUONConstants::fgkNilRecHitStruct,
	 AliHLTMUONConstants::fgkNilRecHitStruct,
	 AliHLTMUONConstants::fgkNilRecHitStruct
	}
};
	
const AliHLTMUONMansoRoIStruct
AliHLTMUONConstants::fgkNilMansoRoIStruct = {0, 0, 0, 0};

const AliHLTMUONMansoCandidateStruct
AliHLTMUONConstants::fgkNilMansoCandidateStruct = {
	AliHLTMUONConstants::fgkNilMansoTrackStruct,
	{
	 AliHLTMUONConstants::fgkNilMansoRoIStruct,
	 AliHLTMUONConstants::fgkNilMansoRoIStruct,
	 AliHLTMUONConstants::fgkNilMansoRoIStruct,
	 AliHLTMUONConstants::fgkNilMansoRoIStruct
	},
	0, 0
};

const AliHLTMUONTrackDecisionStruct
AliHLTMUONConstants::fgkNilTrackDecisionStruct = {0, 0, 0};

const AliHLTMUONPairDecisionStruct
AliHLTMUONConstants::fgkNilPairDecisionStruct = {0, 0, 0, 0};


const AliHLTComponentDataType
AliHLTMUONConstants::fgkDDLRawDataType = AliHLTComponentDataTypeInitializer(kAliHLTDataTypeDDLRaw.fID, kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkTriggerRecordsBlockDataType = AliHLTComponentDataTypeInitializer("TRIGRECS", kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkTrigRecsDebugBlockDataType = AliHLTComponentDataTypeInitializer("TRIGRDBG", kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkRecHitsBlockDataType = AliHLTComponentDataTypeInitializer("RECHITS ", kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkClusterBlockDataType = AliHLTComponentDataTypeInitializer("CLUSTERS", kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkChannelBlockDataType = AliHLTComponentDataTypeInitializer("CHANNELS", kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkMansoTracksBlockDataType = AliHLTComponentDataTypeInitializer("MANTRACK", kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkMansoCandidatesBlockDataType = AliHLTComponentDataTypeInitializer("MNCANDID", kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkSinglesDecisionBlockDataType = AliHLTComponentDataTypeInitializer("DECIDSIN", kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkPairsDecisionBlockDataType = AliHLTComponentDataTypeInitializer("DECIDPAR", kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkRootifiedEventDataType = AliHLTComponentDataTypeInitializer("ROOTEVNT", kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkESDDataType = AliHLTComponentDataTypeInitializer(kAliHLTDataTypeESDObject.fID, kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkClusterStoreDataType = AliHLTComponentDataTypeInitializer("CLUSTORE", kAliHLTDataOriginMUON);

const AliHLTComponentDataType
AliHLTMUONConstants::fgkHistogramDataType = AliHLTComponentDataTypeInitializer("ROOTHIST", kAliHLTDataOriginMUON);

const char* AliHLTMUONConstants::fgkRecHitsSourceId = "MUONRecHitsSource";
const char* AliHLTMUONConstants::fgkTriggerRecordsSourceId = "MUONTriggerRecordsSource";
const char* AliHLTMUONConstants::fgkTracksSourceId = "MUONTracksSource";
const char* AliHLTMUONConstants::fgkDigitPublisherId = "MUONDigitPublisher";
const char* AliHLTMUONConstants::fgkTriggerReconstructorId = "MUONTriggerReconstructor";
const char* AliHLTMUONConstants::fgkHitReconstructorId = "MUONHitReconstructor";
const char* AliHLTMUONConstants::fgkMansoTrackerFSMId = "MUONMansoTrackerFSM";
const char* AliHLTMUONConstants::fgkDecisionComponentId = "MUONDecisionComponent";
const char* AliHLTMUONConstants::fgkESDMakerId = "MUONESDMaker";
const char* AliHLTMUONConstants::fgkRootifierComponentId = "MUONRootifier";
const char* AliHLTMUONConstants::fgkEmptyEventFilterComponentId = "MUONEmptyEventFilter";
const char* AliHLTMUONConstants::fgkDataCheckerComponentId = "MUONDataChecker";
const char* AliHLTMUONConstants::fgkClusterFinderId = "MUONClusterFinder";
const char* AliHLTMUONConstants::fgkRawDataHistogrammerId = "MUONRawDataHistogrammer";

const char* AliHLTMUONConstants::fgkTriggerReconstructorCDBPath = "HLT/ConfigMUON/TriggerReconstructor";
const char* AliHLTMUONConstants::fgkHitReconstructorCDBPath = "HLT/ConfigMUON/HitReconstructor";
const char* AliHLTMUONConstants::fgkMansoTrackerFSMCDBPath = "HLT/ConfigMUON/MansoTrackerFSM";
const char* AliHLTMUONConstants::fgkDecisionComponentCDBPath = "HLT/ConfigMUON/DecisionComponent";

