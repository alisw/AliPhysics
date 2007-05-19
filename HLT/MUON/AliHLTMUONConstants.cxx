/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
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


const AliHLTMUONTriggerRecordStruct
AliHLTMUONConstants::fgkNilTriggerRecordStruct = {
	0, 0, 0, 0, 0,
	{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
};

const AliHLTMUONTrigRecInfoStruct
AliHLTMUONConstants::fgkNilTrigRecInfoStruct = {0, 0, 0, 0};

const AliHLTMUONTriggerChannelStruct
AliHLTMUONConstants::fgkNilTriggerChannelStruct = {0, 0, 0, 0};
	
const AliHLTMUONRecHitStruct
AliHLTMUONConstants::fgkNilRecHitStruct = {0, 0, 0};

const AliHLTMUONClusterStruct
AliHLTMUONConstants::fgkNilClusterStruct = {0, {0, 0, 0}, 0, 0};

const AliHLTMUONChannelStruct
AliHLTMUONConstants::fgkNilChannelStruct = {0, 0, 0, 0, 0};


const AliHLTComponentDataType
AliHLTMUONConstants::fgkTriggerDDLStreamDataType = {
	sizeof(AliHLTComponentDataType),
	{'D','D','L','T','R','I','G','R'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkTrackingDDLStreamDataType = {
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
AliHLTMUONConstants::fgkMansoRoIBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'M','A','N','S','O','R','O','I'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkMansoTrialsBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'M','A','N','T','R','I','A','L'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkDecisionBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'D','E','C','I','S','I','O','N'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkDecisionDebugBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'D','E','C','I','S','D','B','G'},
	{'D','I','M','U'}
};
