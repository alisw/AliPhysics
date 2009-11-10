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

// $Id$

/**
 * @file   AliHLTMUONSinglesDecisionBlockStruct.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   21 May 2007
 * @brief  Implementation of useful stream and comparison operators.
 */

#include "AliHLTMUONSinglesDecisionBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include <cassert>


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONSinglesDecisionBlockStruct& block
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(block) );

	const AliHLTMUONTrackDecisionStruct* decision =
		reinterpret_cast<const AliHLTMUONTrackDecisionStruct*>(&block + 1);
	
	stream 	<< "{fHeader = " << block.fHeader
		<< ", fNlowPt = " << block.fNlowPt
		<< ", fNhighPt = " << block.fNhighPt
		<< ", fDecision[] = [";
	if (block.fHeader.fNrecords > 0) stream << decision[0];
	for (AliHLTUInt32_t i = 1; i < block.fHeader.fNrecords; i++)
		stream << ", " << decision[i];
	stream << "]}";
	return stream;
}


bool operator == (
		const AliHLTMUONSinglesDecisionBlockStruct& a,
		const AliHLTMUONSinglesDecisionBlockStruct& b
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(a) );
	assert( AliHLTMUONUtils::IntegrityOk(b) );

	const AliHLTMUONTrackDecisionStruct* decisionA =
		reinterpret_cast<const AliHLTMUONTrackDecisionStruct*>(&a + 1);
	const AliHLTMUONTrackDecisionStruct* decisionB =
		reinterpret_cast<const AliHLTMUONTrackDecisionStruct*>(&b + 1);
	
	// First check if the blocks have the same header. If they do then
	// check if every track decision is the same. In either case if we find
	// a difference return false.
	if (a.fHeader != b.fHeader) return false;
	if (a.fNlowPt != b.fNlowPt) return false;
	if (a.fNhighPt != b.fNhighPt) return false;
	for (AliHLTUInt32_t i = 0; i < a.fHeader.fNrecords; i++)
		if (decisionA[i] != decisionB[i]) return false;
	return true;
}
