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

// $Id$

/**
 * @file   AliHLTMUONConstants.cxx
 * @author Indranil Das <indra.das@saha.ac.in>,
 *         Artur Szostak <artursz@iafrica.com>
 * @date   17 May 2007
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
#include "AliHLTMUONTracksBlockStruct.h"
#include "AliHLTMUONSinglesDecisionBlockStruct.h"
#include "AliHLTMUONPairsDecisionBlockStruct.h"

ClassImp(AliHLTMUONConstants);


const AliHLTMUONTriggerRecordStruct& AliHLTMUONConstants::NilTriggerRecordStruct()
{
	// Returns a nil trigger record structure.
	static const AliHLTMUONTriggerRecordStruct nilTriggerRecordStruct = {
		0, 0, 0, 0, 0,
		{
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct()
		}
	};
	return nilTriggerRecordStruct;
}


const AliHLTMUONTrigRecInfoStruct& AliHLTMUONConstants::NilTrigRecInfoStruct()
{
	// Returns a nil trigger record debug information structure.
	static const AliHLTMUONTrigRecInfoStruct nilTrigRecInfoStruct = {
		0, {0, 0, 0, 0}, 0, 0,
		{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}
	};
	return nilTrigRecInfoStruct;
}


const AliHLTMUONRecHitStruct& AliHLTMUONConstants::NilRecHitStruct()
{
	// Returns a nil reconstructed hit structure.
	static const AliHLTMUONRecHitStruct nilRecHitStruct = {0, 0, 0, 0};
	return nilRecHitStruct;
}


const AliHLTMUONChannelStruct& AliHLTMUONConstants::NilChannelStruct()
{
	// Returns a nil tracking chamber channel structure.
	static const AliHLTMUONChannelStruct nilChannelStruct = {0, 0, 0, 0, 0, 0};
	return nilChannelStruct;
}


const AliHLTMUONClusterStruct& AliHLTMUONConstants::NilClusterStruct()
{
	// Returns a nil tracking chamber cluster.
	static const AliHLTMUONClusterStruct nilClusterStruct = {
		0, AliHLTMUONConstants::NilRecHitStruct(), 0, 0, 0, 0, 0
	};
	return nilClusterStruct;
}


const AliHLTMUONMansoTrackStruct& AliHLTMUONConstants::NilMansoTrackStruct()
{
	// Returns a nil Manso track structure.
	static const AliHLTMUONMansoTrackStruct nilMansoTrackStruct = {
		0, 0, 0, 0, 0, 0, 0,
		{
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct()
		}
	};
	return nilMansoTrackStruct;
}


const AliHLTMUONMansoRoIStruct& AliHLTMUONConstants::NilMansoRoIStruct()
{
	// Returns a nil Manso region of interest structure.
	static const AliHLTMUONMansoRoIStruct nilMansoRoIStruct = {0, 0, 0, 0};
	return nilMansoRoIStruct;
}


const AliHLTMUONMansoCandidateStruct& AliHLTMUONConstants::NilMansoCandidateStruct()
{
	// Returns a nil Manso candidate track structure.
	static const AliHLTMUONMansoCandidateStruct nilMansoCandidateStruct = {
		AliHLTMUONConstants::NilMansoTrackStruct(),
		{
			AliHLTMUONConstants::NilMansoRoIStruct(),
			AliHLTMUONConstants::NilMansoRoIStruct(),
			AliHLTMUONConstants::NilMansoRoIStruct(),
			AliHLTMUONConstants::NilMansoRoIStruct()
		},
		0, 0
	};
	return nilMansoCandidateStruct;
}


const AliHLTMUONTrackStruct& AliHLTMUONConstants::NilTrackStruct()
{
	// Returns a nil track structure.
	static const AliHLTMUONTrackStruct nilTrackStruct = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		{
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct(),
			AliHLTMUONConstants::NilRecHitStruct()
		}
	};
	return nilTrackStruct;
}


const AliHLTMUONTrackDecisionStruct& AliHLTMUONConstants::NilTrackDecisionStruct()
{
	// Returns a nil decision structure for single track.
	static const AliHLTMUONTrackDecisionStruct nilTrackDecisionStruct = {0, 0, 0};
	return nilTrackDecisionStruct;
}


const AliHLTMUONPairDecisionStruct& AliHLTMUONConstants::NilPairDecisionStruct()
{
	// Returns a nil decision structure for track pair.
	static const AliHLTMUONPairDecisionStruct nilPairDecisionStruct = {0, 0, 0, 0};
	return nilPairDecisionStruct;
}


const AliHLTComponentDataType& AliHLTMUONConstants::DDLRawDataType()
{
	// Returns the raw data type for MUON DDLs.
	static const AliHLTComponentDataType ddlRawDataType = AliHLTComponentDataTypeInitializer(kAliHLTDataTypeDDLRaw.fID, kAliHLTDataOriginMUON);
	return ddlRawDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::TriggerRecordsBlockDataType()
{
	// Returns a trigger records block type generated by trigger DDL translation components.
	static const AliHLTComponentDataType triggerRecordsBlockDataType = AliHLTComponentDataTypeInitializer("TRIGRECS", kAliHLTDataOriginMUON);
	return triggerRecordsBlockDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::TrigRecsDebugBlockDataType()
{
	// Returns a debugging information block type generated by trigger DDL translation components.
	static const AliHLTComponentDataType trigRecsDebugBlockDataType = AliHLTComponentDataTypeInitializer("TRIGRDBG", kAliHLTDataOriginMUON);
	return trigRecsDebugBlockDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::RecHitsBlockDataType()
{
	// Returns a reconstructed hits block type generated by hit reconstruction components.
	static const AliHLTComponentDataType recHitsBlockDataType = AliHLTComponentDataTypeInitializer("RECHITS ", kAliHLTDataOriginMUON);
	return recHitsBlockDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::ClusterBlockDataType()
{
	// Returns a debugging information block type for reconstructed hit clusters.
	static const AliHLTComponentDataType clusterBlockDataType = AliHLTComponentDataTypeInitializer("CLUSTERS", kAliHLTDataOriginMUON);
	return clusterBlockDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::ChannelBlockDataType()
{
	// Returns a debugging information block type for channels corresponding to clusters.
	static const AliHLTComponentDataType channelBlockDataType = AliHLTComponentDataTypeInitializer("CHANNELS", kAliHLTDataOriginMUON);
	return channelBlockDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::MansoTracksBlockDataType()
{
	// Returns a Manso tracks block type generated by Manso tracker components.
	static const AliHLTComponentDataType mansoTracksBlockDataType = AliHLTComponentDataTypeInitializer("MANTRACK", kAliHLTDataOriginMUON);
	return mansoTracksBlockDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::MansoCandidatesBlockDataType()
{
	// Returns a data type for debugging information data blocks about track candidates generated by the Manso algorithm.
	static const AliHLTComponentDataType mansoCandidatesBlockDataType = AliHLTComponentDataTypeInitializer("MNCANDID", kAliHLTDataOriginMUON);
	return mansoCandidatesBlockDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::TracksBlockDataType()
{
	// Returns a full tracks block type generated by the tracker components.
	static const AliHLTComponentDataType tracksBlockDataType = AliHLTComponentDataTypeInitializer("TRACKS  ", kAliHLTDataOriginMUON);
	return tracksBlockDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::SinglesDecisionBlockDataType()
{
	// Returns a trigger decision block type for single track decisions.
	static const AliHLTComponentDataType singlesDecisionBlockDataType = AliHLTComponentDataTypeInitializer("DECIDSIN", kAliHLTDataOriginMUON);
	return singlesDecisionBlockDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::PairsDecisionBlockDataType()
{
	// Returns a trigger decision block type for pairs of particles.
	static const AliHLTComponentDataType pairsDecisionBlockDataType = AliHLTComponentDataTypeInitializer("DECIDPAR", kAliHLTDataOriginMUON);
	return pairsDecisionBlockDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::RootifiedEventDataType()
{
	// Returns an AliHLTMUONEvent ROOT object data type.
	static const AliHLTComponentDataType rootifiedEventDataType = AliHLTComponentDataTypeInitializer("ROOTEVNT", kAliHLTDataOriginMUON);
	return rootifiedEventDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::ESDDataType()
{
	// Returns the ESD data type with origin equal to MUON.
	static const AliHLTComponentDataType esdDataType = AliHLTComponentDataTypeInitializer(kAliHLTDataTypeESDObject.fID, kAliHLTDataOriginMUON);
	return esdDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::ClusterStoreDataType()
{
	// Returns the offline algorithm's cluster store object type.
	static const AliHLTComponentDataType clusterStoreDataType = AliHLTComponentDataTypeInitializer("CLUSTORE", kAliHLTDataOriginMUON);
	return clusterStoreDataType;
}


const AliHLTComponentDataType& AliHLTMUONConstants::HistogramDataType()
{
	// Returns the TH1/2/3 histogram data type.
	static const AliHLTComponentDataType histogramDataType = AliHLTComponentDataTypeInitializer("ROOTHIST", kAliHLTDataOriginMUON);
	return histogramDataType;
}


const char* AliHLTMUONConstants::RecHitsSourceId()
{
	// Returns the name of source component for reconstructed hits for debugging.
	static const char* recHitsSourceId = "MUONRecHitsSource";
	return recHitsSourceId;
}


const char* AliHLTMUONConstants::TriggerRecordsSourceId()
{
	// Returns the name of source component for trigger records for debugging.
	static const char* triggerRecordsSourceId = "MUONTriggerRecordsSource";
	return triggerRecordsSourceId;
}


const char* AliHLTMUONConstants::TracksSourceId()
{
	// Returns the name of source component for tracks for debugging.
	static const char* tracksSourceId = "MUONTracksSource";
	return tracksSourceId;
}


const char* AliHLTMUONConstants::DigitPublisherId()
{
	// Returns the component name for publishing DDL streams from digits.
	static const char* digitPublisherId = "MUONDigitPublisher";
	return digitPublisherId;
}

const char* AliHLTMUONConstants::TriggerReconstructorId()
{
	// Returns the trigger record reconstructor component name.
	static const char* triggerReconstructorId = "MUONTriggerReconstructor";
	return triggerReconstructorId;
}


const char* AliHLTMUONConstants::HitReconstructorId()
{
	// Returns the centre of gravity cluster finder component name.
	static const char* hitReconstructorId = "MUONHitReconstructor";
	return hitReconstructorId;
}


const char* AliHLTMUONConstants::MansoTrackerFSMId()
{
	// Returns the Manso tracker FSM implementation component name.
	static const char* mansoTrackerFSMId = "MUONMansoTrackerFSM";
	return mansoTrackerFSMId;
}


const char* AliHLTMUONConstants::FullTrackerId()
{
	// Returns the full tracker implementation component name.
	static const char* fullTrackerId = "MUONFullTracker";
	return fullTrackerId;
}


const char* AliHLTMUONConstants::DecisionComponentId()
{
	// Returns the dHLT decision component name.
	static const char* decisionComponentId = "MUONDecisionComponent";
	return decisionComponentId;
}


const char* AliHLTMUONConstants::ESDMakerId()
{
	// Returns the name of ESD maker component which converts dHLT data to AliESDEvent classes. 
	static const char* esdMakerId = "MUONESDMaker";
	return esdMakerId;
}


const char* AliHLTMUONConstants::RootifierComponentId()
{
	// Returns the name of the event filter debugging component.
	static const char* rootifierComponentId = "MUONRootifier";
	return rootifierComponentId;
}


const char* AliHLTMUONConstants::EmptyEventFilterComponentId()
{
	// Returns the name of the event filter debugging component.
	static const char* emptyEventFilterComponentId = "MUONEmptyEventFilter";
	return emptyEventFilterComponentId;
}


const char* AliHLTMUONConstants::DataCheckerComponentId()
{
	// Returns the name of data checking component for debugging.
	static const char* dataCheckerComponentId = "MUONDataChecker";
	return dataCheckerComponentId;
}


const char* AliHLTMUONConstants::ClusterFinderId()
{
	// Returns the name of cluster finder implementing offline algorithms.
	static const char* clusterFinderId = "MUONClusterFinder";
	return clusterFinderId;
}


const char* AliHLTMUONConstants::RawDataHistogrammerId()
{
	// Returns the raw data histogrammer component name.
	static const char* rawDataHistogrammerId = "MUONRawDataHistogrammer";
	return rawDataHistogrammerId;
}


const char* AliHLTMUONConstants::ClusterHistogrammerId()
{
	// Returns the cluster data histogrammer component name.
	static const char* clusterHistogrammerId = "MUONClusterHistogrammer";
	return clusterHistogrammerId;
}


const char* AliHLTMUONConstants::TriggerReconstructorCDBPath()
{
	// Returns the path to the CDB entry for the trigger reconstruction component.
	static const char* triggerReconstructorCDBPath = "HLT/ConfigMUON/TriggerReconstructor";
	return triggerReconstructorCDBPath;
}


const char* AliHLTMUONConstants::HitReconstructorCDBPath()
{
	// Returns the path to the CDB entry for the hit reconstruction component.
	static const char* hitReconstructorCDBPath = "HLT/ConfigMUON/HitReconstructor";
	return hitReconstructorCDBPath;
}


const char* AliHLTMUONConstants::MansoTrackerFSMCDBPath()
{
	// Returns the path to the CDB entry for the Manso FSM tracker component.
	static const char* mansoTrackerFSMCDBPath = "HLT/ConfigMUON/MansoTrackerFSM";
	return mansoTrackerFSMCDBPath;
}


const char* AliHLTMUONConstants::DecisionComponentCDBPath()
{
	// Returns the path to the CDB entry for trigger decision component.
	static const char* decisionComponentCDBPath = "HLT/ConfigMUON/DecisionComponent";
	return decisionComponentCDBPath;
}


const char* AliHLTMUONConstants::FieldIntegralsCDBPath()
{
	// Returns the path to the CDB entry for magnetic field integrals.
	static const char* fieldIntegralsCDBPath = "HLT/ConfigMUON/FieldIntegrals";
	return fieldIntegralsCDBPath;
}


double AliHLTMUONConstants::DefaultNonBendingReso()
{
	// Returns the typical X (non-bending plane) resolution of the hit reconstruction (units = cm).
	static const double resolution = 0.144;
	return resolution;
}


double AliHLTMUONConstants::DefaultBendingReso()
{
	// Returns the typical Y (bending plane) resolution of the hit reconstruction (units = cm).
	static const double resolution = 0.01;
	return resolution;
}

