#ifndef ALIHLTMUONCONSTANTS_H
#define ALIHLTMUONCONSTANTS_H
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
 * @file   AliHLTMUONConstants.h
 * @author Indranil Das <indra.das@saha.ac.in>,
 *         Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Class containing various dimuon HLT constants used in the system.
 */

// Forward declare structures.
extern "C" {
struct AliHLTComponentDataType;
struct AliHLTMUONTriggerRecordStruct;
struct AliHLTMUONTrigRecInfoStruct;
struct AliHLTMUONTriggerChannelStruct;
struct AliHLTMUONRecHitStruct;
struct AliHLTMUONChannelStruct;
struct AliHLTMUONClusterStruct;
struct AliHLTMUONMansoTrackStruct;
struct AliHLTMUONMansoRoIStruct;
struct AliHLTMUONMansoCandidateStruct;
struct AliHLTMUONTrackDecisionStruct;
struct AliHLTMUONPairDecisionStruct;
} // extern "C"

/**
 * AliHLTMUONConstants contains a list of global dimuon HLT specific constants
 * and constant structures used in the system.
 * Static methods are provided to access these values.
 */
class AliHLTMUONConstants
{
public:

	static const AliHLTMUONTriggerRecordStruct& NilTriggerRecordStruct()
	{
		return fgkNilTriggerRecordStruct;
	}

	static const AliHLTMUONTrigRecInfoStruct& NilTrigRecInfoStruct()
	{
		return fgkNilTrigRecInfoStruct;
	}

	static const AliHLTMUONTriggerChannelStruct& NilTriggerChannelStruct()
	{
		return fgkNilTriggerChannelStruct;
	}

	static const AliHLTMUONRecHitStruct& NilRecHitStruct()
	{
		return fgkNilRecHitStruct;
	}

	static const AliHLTMUONChannelStruct& NilChannelStruct()
	{
		return fgkNilChannelStruct;
	}

	static const AliHLTMUONClusterStruct& NilClusterStruct()
	{
		return fgkNilClusterStruct;
	}

	static const AliHLTMUONMansoTrackStruct& NilMansoTrackStruct()
	{
		return fgkNilMansoTrackStruct;
	}

	static const AliHLTMUONMansoRoIStruct& NilMansoRoIStruct()
	{
		return fgkNilMansoRoIStruct;
	}

	static const AliHLTMUONMansoCandidateStruct& NilMansoCandidateStruct()
	{
		return fgkNilMansoCandidateStruct;
	}

	static const AliHLTMUONTrackDecisionStruct& NilTrackDecisionStruct()
	{
		return fgkNilTrackDecisionStruct;
	}

	static const AliHLTMUONPairDecisionStruct& NilPairDecisionStruct()
	{
		return fgkNilPairDecisionStruct;
	}

	static const AliHLTComponentDataType& TriggerDDLRawDataType()
	{
		return fgkTriggerDDLRawDataType;
	}

	static const AliHLTComponentDataType& TrackingDDLRawDataType()
	{
		return fgkTrackingDDLRawDataType;
	}

	static const AliHLTComponentDataType& TriggerRecordsBlockDataType()
	{
		return fgkTriggerRecordsBlockDataType;
	}

	static const AliHLTComponentDataType& TrigRecsDebugBlockDataType()
	{
		return fgkTrigRecsDebugBlockDataType;
	}

	static const AliHLTComponentDataType& TriggerChannelBlockDataType()
	{
		return fgkTriggerChannelBlockDataType;
	}

	static const AliHLTComponentDataType& RecHitsBlockDataType()
	{
		return fgkRecHitsBlockDataType;
	}

	static const AliHLTComponentDataType& ClusterBlockDataType()
	{
		return fgkClusterBlockDataType;
	}

	static const AliHLTComponentDataType& ChannelBlockDataType()
	{
		return fgkChannelBlockDataType;
	}

	static const AliHLTComponentDataType& MansoTracksBlockDataType()
	{
		return fgkMansoTracksBlockDataType;
	}

	static const AliHLTComponentDataType& MansoCandidatesBlockDataType()
	{
		return fgkMansoCandidatesBlockDataType;
	}

	static const AliHLTComponentDataType& SinglesDecisionBlockDataType()
	{
		return fgkSinglesDecisionBlockDataType;
	}

	static const AliHLTComponentDataType& PairsDecisionBlockDataType()
	{
		return fgkPairsDecisionBlockDataType;
	}
	
	static const char* RecHitsSourceId()
	{
		return fgkRecHitsSourceId;
	}
	
	static const char* TriggerRecordsSourceId()
	{
		return fgkTriggerRecordsSourceId;
	}
	
	static const char* MansoTracksSourceId()
	{
		return fgkMansoTracksSourceId;
	}
	
	static const char* TriggerReconstructorId()
	{
		return fgkTriggerReconstructorId;
	}
	
	static const char* HitReconstructorId()
	{
		return fgkHitReconstructorId;
	}
	
	static const char* MansoTrackerFSMId()
	{
		return fgkMansoTrackerFSMId;
	}
	
	static const char* DecisionComponentId()
	{
		return fgkDecisionComponentId;
	}

private:

	// Should never have to create, destroy or copy this object.
	AliHLTMUONConstants();
	AliHLTMUONConstants(const AliHLTMUONConstants& obj);
	~AliHLTMUONConstants();
	AliHLTMUONConstants& operator = (const AliHLTMUONConstants& obj);
	
	// The following are null/nil structures that can also be used as sentinels:
	static const AliHLTMUONTriggerRecordStruct fgkNilTriggerRecordStruct; // Nil trigger record.
	static const AliHLTMUONTrigRecInfoStruct fgkNilTrigRecInfoStruct; // Nil trigger record debug information.
	static const AliHLTMUONTriggerChannelStruct fgkNilTriggerChannelStruct; // Nil trigger chamber channel.
	static const AliHLTMUONRecHitStruct fgkNilRecHitStruct; // Nil reconstructed hit.
	static const AliHLTMUONChannelStruct fgkNilChannelStruct; // Nil tracking chamber channel.
	static const AliHLTMUONClusterStruct fgkNilClusterStruct; // Nil tracking chamber cluster.
	static const AliHLTMUONMansoTrackStruct fgkNilMansoTrackStruct; // Nil manso track.
	static const AliHLTMUONMansoRoIStruct fgkNilMansoRoIStruct; // Nil manso region of interest.
	static const AliHLTMUONMansoCandidateStruct fgkNilMansoCandidateStruct; // Nil manso candidate track.
	static const AliHLTMUONTrackDecisionStruct fgkNilTrackDecisionStruct; // Nil decision for single track.
	static const AliHLTMUONPairDecisionStruct fgkNilPairDecisionStruct; // Nil decision for track pair.

	// HLT component input and output data block types:
	static const AliHLTComponentDataType fgkTriggerDDLRawDataType; // DDL packed data block type from dimuon trigger stations.
	static const AliHLTComponentDataType fgkTrackingDDLRawDataType; // DDL packed data block type from dimuon tracking stations.
	static const AliHLTComponentDataType fgkTriggerRecordsBlockDataType; // Trigger records block type generated by trigger DDL translation components.
	static const AliHLTComponentDataType fgkTrigRecsDebugBlockDataType; // Debugging information block type generated by trigger DDL translation components.
	static const AliHLTComponentDataType fgkTriggerChannelBlockDataType; // Debugging information about the channels from the hardware trigger.
	static const AliHLTComponentDataType fgkRecHitsBlockDataType; // Reconstructed hits block type generated by hit reconstruction components.
	static const AliHLTComponentDataType fgkClusterBlockDataType; // Debugging information block type for reconstructed hit clusters.
	static const AliHLTComponentDataType fgkChannelBlockDataType; // Debugging information block type for channels corresponding to clusters.
	static const AliHLTComponentDataType fgkMansoTracksBlockDataType; // Manso tracks block type generated by Manso tracker components.
	static const AliHLTComponentDataType fgkMansoCandidatesBlockDataType; // Debugging information about a track candidate generated by the Manso algorithm.
	static const AliHLTComponentDataType fgkSinglesDecisionBlockDataType; // Trigger decision block type for single track decisions.
	static const AliHLTComponentDataType fgkPairsDecisionBlockDataType; // Trigger decision block type for pairs of particles.
	
	// Component ID names:
	static const char* fgkRecHitsSourceId; // Name of source component for reconstructed hits for debugging.
	static const char* fgkTriggerRecordsSourceId; // Name of source component for trigger records for debugging.
	static const char* fgkMansoTracksSourceId; // Name of source component for Manso tracks for debugging.
	static const char* fgkTriggerReconstructorId; // Trigger record reconstructor component name.
	static const char* fgkHitReconstructorId; // Centre of gravity cluster finder component name.
	static const char* fgkMansoTrackerFSMId; // Manso tracker FSM implementation component name.
	static const char* fgkDecisionComponentId; // dHLT decision component name.
};

#endif // ALIHLTMUONCONSTANTS_H
