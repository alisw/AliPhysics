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

const AliHLTMUONRecHitStruct
AliHLTMUONConstants::fgkNilRecHitStruct = {0, 0, 0};

const AliHLTMUONChannelInfoStruct
AliHLTMUONConstants::fgkNilChannelInfoStruct = {0, 0, 0, 0};

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
AliHLTMUONConstants::fgkRecHitsBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'R','E','C','H','I','T','S',' '},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkRecHitsDebugBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'R','H','I','T','S','d','b','g'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkTriggerRecordsBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'T','R','I','G','R','E','C','S'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkTriggerRecordsDebugBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'T','R','I','G','R','d','b','g'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkMansoTracksBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'M','A','N','T','R','A','C','K'},
	{'D','I','M','U'}
};

const AliHLTComponentDataType
AliHLTMUONConstants::fgkMansoTracksDebugBlockDataType = {
	sizeof(AliHLTComponentDataType),
	{'M','N','T','R','K','d','b','g'},
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
	{'D','E','C','I','S','d','b','g'},
	{'D','I','M','U'}
};
