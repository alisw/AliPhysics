#ifndef ALIHLTMUONMANSOTRACKERFSMCALLBACK_H
#define ALIHLTMUONMANSOTRACKERFSMCALLBACK_H
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
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
 *  @file   AliHLTMUONMansoTrackerFSMCallback.h
 *  @author Artur Szostak <artursz@iafrica.com>
 *  @date   
 *  @brief  This is the abstract callback interface required by the Finite
 *          State Machine (FSM) implementation of the Manso algorithm.
 */

#include "AliHLTMUONDataTypes.h"
#include <cassert>

class AliHLTMUONMansoTrackerFSM;

class AliHLTMUONMansoTrackerFSMCallback
{
public:

	virtual ~AliHLTMUONMansoTrackerFSMCallback() {};
	
	/* All clusters that fall within the specified boundary box on the specified
	   chamber should be returned to the tracker, by calling the ReturnClusters
	   method of the given tracker. The same tag parameter must be passed on the 
	   ReturnClusters method's parameter list.
	   @param left    The left border of the boundary box (x direction).
	   @param right   The right border of the boundary box (x direction).
	   @param bottom  The bottom border of the boundary box (y direction).
	   @param top     The top border of the boundary box (y direction).
	 */
	virtual void RequestClusters(
			AliHLTMUONMansoTrackerFSM* tracker,
			AliHLTFloat32_t left, AliHLTFloat32_t right,
			AliHLTFloat32_t bottom, AliHLTFloat32_t top,
			AliHLTMUONChamberName chamber, const void* tag
		) = 0;

	/* When this method is called then one knows no more RequestClusters method
	   calls are expected.
	 */
	virtual void EndOfClusterRequests(AliHLTMUONMansoTrackerFSM* tracker) = 0;

	/* This method is called when the tracker has found a track. The FillTrackData
	   method of the given tracker should be called to receive the track data.
	   At this point all cluster blocks can be released.
	 */
	virtual void FoundTrack(AliHLTMUONMansoTrackerFSM* tracker) = 0;
	
	/* When the tracker is finished with its work but no track was found then
	   this method is called. At this point no more work should be performed by
	   the tracker and all cluster blocks can be released.
	 */
	virtual void NoTrackFound(AliHLTMUONMansoTrackerFSM* tracker) = 0;
};

#endif // ALIHLTMUONMANSOTRACKERFSMCALLBACK_H
