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
 * @file   AliHLTMUONUtils.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Implementation of AliHLTMUONUtils utility routines.
 */

#include "AliHLTMUONUtils.h"


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONRecHitsDebugBlockStruct& block)
{
	// Perform some sanity checks:
	// We expect the start and end offsets to correspond to the number
	// of elements in the array. The offsets should not overlap and
	// the arrays should be consecutive.
	if (block.fClustersEndOffset < block.fClustersStartOffset) return false;
	if (block.fChannelsEndOffset < block.fChannelsStartOffset) return false;

	if ((block.fClustersEndOffset - block.fClustersStartOffset)
		/ sizeof(AliHLTMUONClusterInfoStruct) != block.fNclusters
	   )
		return false;

	if ((block.fChannelsEndOffset - block.fChannelsStartOffset)
		/ sizeof(AliHLTMUONChannelInfoStruct) != block.fNchannels
	   )
		return false;

	if (block.fClustersEndOffset != block.fChannelsStartOffset) return false;

	if (block.fClustersStartOffset != sizeof(AliHLTMUONRecHitsDebugBlockStruct))
		return false;

	// TODO: check the channel ID's etc...

	return true;
}
