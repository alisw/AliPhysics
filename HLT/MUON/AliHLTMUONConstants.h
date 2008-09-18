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

#include "TObject.h"

// Forward declare structures.
extern "C" {
struct AliHLTComponentDataType;
struct AliHLTMUONTriggerRecordStruct;
struct AliHLTMUONTrigRecInfoStruct;
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

	// Returns the raw data type for MUON DDLs. To figure out if the DDL format
	// will be for a tracking DDL or trigger DDL one needs to also check the
	// sepcification word of the input data block. If one of the first 20 least
	// significant bits are set then it is a tracker DDL otherwise if it is
	// the 21st or 22nd bit then it is from the muon trigger.
	static const AliHLTComponentDataType& DDLRawDataType()
	{
		return fgkDDLRawDataType;
	}

	static const AliHLTComponentDataType& TriggerRecordsBlockDataType()
	{
		return fgkTriggerRecordsBlockDataType;
	}

	static const AliHLTComponentDataType& TrigRecsDebugBlockDataType()
	{
		return fgkTrigRecsDebugBlockDataType;
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

	static const AliHLTComponentDataType& ESDDataType()
	{
		return fgkESDDataType;
	}

	static const AliHLTComponentDataType& ClusterStoreDataType()
	{
		return fgkClusterStoreDataType;
	}

	static const AliHLTComponentDataType& HistogramDataType()
	{
		return fgkHistogramDataType;
	}
	
	static const char* RecHitsSourceId()
	{
		return fgkRecHitsSourceId;
	}
	
	static const char* TriggerRecordsSourceId()
	{
		return fgkTriggerRecordsSourceId;
	}
	
	static const char* TracksSourceId()
	{
		return fgkTracksSourceId;
	}
	
	static const char* DigitPublisherId()
	{
		return fgkDigitPublisherId;
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
	
	static const char* ESDMakerId()
	{
		return fgkESDMakerId;
	}
	
	static const char* RootifierComponentId()
	{
		return fgkRootifierComponentId;
	}
	
	static const char* EmptyEventFilterComponentId()
	{
		return fgkEmptyEventFilterComponentId;
	}
	
	static const char* DataCheckerComponentId()
	{
		return fgkDataCheckerComponentId;
	}
	
	static const char* ClusterFinderId()
	{
		return fgkClusterFinderId;
	}
	
	static const char* RawDataHistogrammerId()
	{
		return fgkRawDataHistogrammerId;
	}
	
	static const char* TriggerReconstructorCDBPath()
	{
		return fgkTriggerReconstructorCDBPath;
	}
	
	static const char* HitReconstructorCDBPath()
	{
		return fgkHitReconstructorCDBPath;
	}
	
	static const char* MansoTrackerFSMCDBPath()
	{
		return fgkMansoTrackerFSMCDBPath;
	}
	
	static const char* DecisionComponentCDBPath()
	{
		return fgkDecisionComponentCDBPath;
	}
	
	static double DefaultNonBendingReso() { return 0.01; }
	
	static double DefaultBendingReso() { return 0.144; }

private:

	// Should never have to create, destroy or copy this object.
	AliHLTMUONConstants();
	AliHLTMUONConstants(const AliHLTMUONConstants& obj);
	~AliHLTMUONConstants();
	AliHLTMUONConstants& operator = (const AliHLTMUONConstants& obj);
	
	// The following are null/nil structures that can also be used as sentinels:
	static const AliHLTMUONTriggerRecordStruct fgkNilTriggerRecordStruct; // Nil trigger record.
	static const AliHLTMUONTrigRecInfoStruct fgkNilTrigRecInfoStruct; // Nil trigger record debug information.
	static const AliHLTMUONRecHitStruct fgkNilRecHitStruct; // Nil reconstructed hit.
	static const AliHLTMUONChannelStruct fgkNilChannelStruct; // Nil tracking chamber channel.
	static const AliHLTMUONClusterStruct fgkNilClusterStruct; // Nil tracking chamber cluster.
	static const AliHLTMUONMansoTrackStruct fgkNilMansoTrackStruct; // Nil manso track.
	static const AliHLTMUONMansoRoIStruct fgkNilMansoRoIStruct; // Nil manso region of interest.
	static const AliHLTMUONMansoCandidateStruct fgkNilMansoCandidateStruct; // Nil manso candidate track.
	static const AliHLTMUONTrackDecisionStruct fgkNilTrackDecisionStruct; // Nil decision for single track.
	static const AliHLTMUONPairDecisionStruct fgkNilPairDecisionStruct; // Nil decision for track pair.

	// HLT component input and output data block types:
	static const AliHLTComponentDataType fgkDDLRawDataType; // DDL packed data block type from dimuon spectrometer.
	static const AliHLTComponentDataType fgkTriggerRecordsBlockDataType; // Trigger records block type generated by trigger DDL translation components.
	static const AliHLTComponentDataType fgkTrigRecsDebugBlockDataType; // Debugging information block type generated by trigger DDL translation components.
	static const AliHLTComponentDataType fgkRecHitsBlockDataType; // Reconstructed hits block type generated by hit reconstruction components.
	static const AliHLTComponentDataType fgkClusterBlockDataType; // Debugging information block type for reconstructed hit clusters.
	static const AliHLTComponentDataType fgkChannelBlockDataType; // Debugging information block type for channels corresponding to clusters.
	static const AliHLTComponentDataType fgkMansoTracksBlockDataType; // Manso tracks block type generated by Manso tracker components.
	static const AliHLTComponentDataType fgkMansoCandidatesBlockDataType; // Debugging information about a track candidate generated by the Manso algorithm.
	static const AliHLTComponentDataType fgkSinglesDecisionBlockDataType; // Trigger decision block type for single track decisions.
	static const AliHLTComponentDataType fgkPairsDecisionBlockDataType; // Trigger decision block type for pairs of particles.
	static const AliHLTComponentDataType fgkESDDataType; // The ESD data type with origin equal to MUON.
	static const AliHLTComponentDataType fgkClusterStoreDataType; // Offline algorithm cluster store object.
	static const AliHLTComponentDataType fgkHistogramDataType; // TH1/2/3 histogram type.
	
	// Component ID names:
	static const char* fgkRecHitsSourceId; // Name of source component for reconstructed hits for debugging.
	static const char* fgkTriggerRecordsSourceId; // Name of source component for trigger records for debugging.
	static const char* fgkTracksSourceId; // Name of source component for tracks for debugging.
	static const char* fgkDigitPublisherId; // Component name for publishing DDL streams from digits.
	static const char* fgkTriggerReconstructorId; // Trigger record reconstructor component name.
	static const char* fgkHitReconstructorId; // Centre of gravity cluster finder component name.
	static const char* fgkMansoTrackerFSMId; // Manso tracker FSM implementation component name.
	static const char* fgkDecisionComponentId; // dHLT decision component name.
	static const char* fgkESDMakerId; // Name of ESD maker component which converts dHLT data to AliESDEvent classes.
	static const char* fgkRootifierComponentId; // The name of the event filter debugging component.
	static const char* fgkEmptyEventFilterComponentId; // The name of the event filter debugging component.
	static const char* fgkDataCheckerComponentId; // Name of data checking component for debugging.
	static const char* fgkClusterFinderId; // Name of cluster finder implementing offline algorithms.
	static const char* fgkRawDataHistogrammerId; // Raw data histogrammer component name.
	
	// CDB path entries to configuration information.
	static const char* fgkTriggerReconstructorCDBPath; // Path to CDB entry for the trigger reconstruction component.
	static const char* fgkHitReconstructorCDBPath; // Path to CDB entry for the hit reconstruction component.
	static const char* fgkMansoTrackerFSMCDBPath; // Path to CDB entry for the Manso FSM tracker component.
	static const char* fgkDecisionComponentCDBPath; // Path to CDB entry for trigger decision component.
	
	ClassDef(AliHLTMUONConstants, 0);  // Interface class to dHLT constants.
};

#endif // ALIHLTMUONCONSTANTS_H
