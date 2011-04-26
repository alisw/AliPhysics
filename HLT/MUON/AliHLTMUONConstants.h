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

// $Id$

/**
 * @file   AliHLTMUONConstants.h
 * @author Indranil Das <indra.das@saha.ac.in>,
 *         Artur Szostak <artursz@iafrica.com>
 * @date   17 May 2007
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
struct AliHLTMUONTrackStruct;
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

	// The following methods return null/nil structures that can also be used as sentinels:
	static const AliHLTMUONTriggerRecordStruct& NilTriggerRecordStruct();
	static const AliHLTMUONTrigRecInfoStruct& NilTrigRecInfoStruct();
	static const AliHLTMUONRecHitStruct& NilRecHitStruct();
	static const AliHLTMUONChannelStruct& NilChannelStruct();
	static const AliHLTMUONClusterStruct& NilClusterStruct();
	static const AliHLTMUONMansoTrackStruct& NilMansoTrackStruct();
	static const AliHLTMUONMansoRoIStruct& NilMansoRoIStruct();
	static const AliHLTMUONMansoCandidateStruct& NilMansoCandidateStruct();
	static const AliHLTMUONTrackStruct& NilTrackStruct();
	static const AliHLTMUONTrackDecisionStruct& NilTrackDecisionStruct();
	static const AliHLTMUONPairDecisionStruct& NilPairDecisionStruct();

	// Methods returning HLT component input and output data block types:
	
	/**
	 * Returns the raw data type for MUON DDLs. To figure out if the DDL format
	 * will be for a tracking DDL or trigger DDL one needs to also check the
	 * sepcification word of the input data block. If one of the first 20 least
	 * significant bits are set then it is a tracker DDL otherwise if it is
	 * the 21st or 22nd bit then it is from the muon trigger.
	 */
	static const AliHLTComponentDataType& DDLRawDataType();

	static const AliHLTComponentDataType& TriggerRecordsBlockDataType();
	static const AliHLTComponentDataType& TrigRecsDebugBlockDataType();
	static const AliHLTComponentDataType& RecHitsBlockDataType();
	static const AliHLTComponentDataType& ClusterBlockDataType();
	static const AliHLTComponentDataType& ChannelBlockDataType();
	static const AliHLTComponentDataType& MansoTracksBlockDataType();
	static const AliHLTComponentDataType& MansoCandidatesBlockDataType();
	static const AliHLTComponentDataType& TracksBlockDataType();
	static const AliHLTComponentDataType& SinglesDecisionBlockDataType();
	static const AliHLTComponentDataType& PairsDecisionBlockDataType();
	static const AliHLTComponentDataType& RootifiedEventDataType();
	static const AliHLTComponentDataType& ESDDataType();
	static const AliHLTComponentDataType& ClusterStoreDataType();
	static const AliHLTComponentDataType& HistogramDataType();
	
	
	// Methods to return component ID names:
	static const char* RecHitsSourceId();
	static const char* TriggerRecordsSourceId();
	static const char* TracksSourceId();
	static const char* DigitPublisherId();
	static const char* TriggerReconstructorId();
	static const char* HitReconstructorId();
	static const char* MansoTrackerFSMId();
	static const char* FullTrackerId();
	static const char* DecisionComponentId();
	static const char* ESDMakerId();
	static const char* RootifierComponentId();
	static const char* EmptyEventFilterComponentId();
	static const char* DataCheckerComponentId();
	static const char* ClusterFinderId();
	static const char* RawDataHistogrammerId();
	static const char* ClusterHistogrammerId();
	
	// Methods for returning the CDB path entries to configuration information:
	static const char* TriggerReconstructorCDBPath();
	static const char* HitReconstructorCDBPath();
	static const char* MansoTrackerFSMCDBPath();
	static const char* DecisionComponentCDBPath();
	static const char* FieldIntegralsCDBPath();
	
	/// Returns the typical X (non-bending plane) resolution of the hit reconstruction (units = cm).
	static double DefaultNonBendingReso();
	
	/// Returns the typical Y (bending plane) resolution of the hit reconstruction (units = cm).
	static double DefaultBendingReso();

private:

	// Should never have to create, destroy or copy this object.
	AliHLTMUONConstants() {}
	AliHLTMUONConstants(const AliHLTMUONConstants& obj);
	virtual ~AliHLTMUONConstants() {}
	AliHLTMUONConstants& operator = (const AliHLTMUONConstants& obj);
	
	ClassDef(AliHLTMUONConstants, 0);  // Interface class to dHLT constants.
};

#endif // ALIHLTMUONCONSTANTS_H
