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
#include "AliHLTMUONTriggerChannelsBlockStruct.h"
#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTMUONClustersBlockStruct.h"
#include "AliHLTMUONChannelsBlockStruct.h"
#include "AliHLTMUONMansoTracksBlockStruct.h"
#include "AliHLTMUONMansoCandidatesBlockStruct.h"
#include "AliHLTMUONSinglesDecisionBlockStruct.h"
#include "AliHLTMUONPairsDecisionBlockStruct.h"


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
AliHLTMUONConstants::fgkNilTrigRecInfoStruct = {0, 0, 0, 0};

const AliHLTMUONTriggerChannelStruct
AliHLTMUONConstants::fgkNilTriggerChannelStruct = {0, 0, 0, 0};
	
const AliHLTMUONRecHitStruct
AliHLTMUONConstants::fgkNilRecHitStruct = {0, 0, 0};

const AliHLTMUONClusterStruct
AliHLTMUONConstants::fgkNilClusterStruct = {
	0, AliHLTMUONConstants::fgkNilRecHitStruct, 0, 0
};

const AliHLTMUONChannelStruct
AliHLTMUONConstants::fgkNilChannelStruct = {0, 0, 0, 0, 0};

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
	}
};

const AliHLTMUONTrackDecisionStruct
AliHLTMUONConstants::fgkNilTrackDecisionStruct = {0, 0};

const AliHLTMUONPairDecisionStruct
AliHLTMUONConstants::fgkNilPairDecisionStruct = {0, 0, 0, 0};


const AliHLTComponentDataType
AliHLTMUONConstants::fgkTriggerDDLRawDataType = {
	sizeof(AliHLTComponentDataType),
	{'D','D','L','T','R','I','G','R'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkTrackingDDLRawDataType = {
	sizeof(AliHLTComponentDataType),
	{'D','D','L','T','R','A','C','K'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkTriggerRecordsBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'T','R','I','G','R','E','C','S'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkTrigRecsDebugBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'T','R','I','G','R','D','B','G'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkTriggerChannelBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'T','R','I','G','C','H','N','L'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkRecHitsBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'R','E','C','H','I','T','S',' '},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkClusterBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'C','L','U','S','T','E','R','S'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkChannelBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'C','H','A','N','N','E','L','S'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkMansoTracksBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'M','A','N','T','R','A','C','K'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkMansoCandidatesBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'M','N','C','A','N','D','I','D'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkSinglesDecisionBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'D','E','C','I','D','S','G','L'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkPairsDecisionBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'D','E','C','I','D','P','A','R'},
	{'D','I','M','U'}
};

const char* AliHLTMUONConstants::fgkRecHitsSourceId = "DimuRecHitsSource";
const char* AliHLTMUONConstants::fgkTriggerRecordsSourceId = "DimuTriggerRecordsSource";
const char* AliHLTMUONConstants::fgkMansoTracksSourceId = "DimuMansoTracksSource";
const char* AliHLTMUONConstants::fgkTriggerReconstructorId = "DimuTriggerReconstructor";
const char* AliHLTMUONConstants::fgkHitReconstructorId = "DimuHitReconstructor";
const char* AliHLTMUONConstants::fgkMansoTrackerFSMId = "DimuMansoTrackerFSM";
const char* AliHLTMUONConstants::fgkDecisionComponentId = "DimuDecisionComponent";
