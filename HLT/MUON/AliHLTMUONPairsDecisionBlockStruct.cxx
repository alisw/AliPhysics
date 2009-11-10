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
 * @file   AliHLTMUONPairsDecisionBlockStruct.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   21 May 2007
 * @brief  Implementation of useful stream and comparison operators.
 */

#include "AliHLTMUONPairsDecisionBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include <cassert>


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONPairsDecisionBlockStruct& block
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(block) );
	
	const AliHLTMUONPairDecisionStruct* decision =
		reinterpret_cast<const AliHLTMUONPairDecisionStruct*>(&block + 1);

	stream 	<< "{fHeader = " << block.fHeader
		<< ", fNunlikeAnyPt = " << block.fNunlikeAnyPt
		<< ", fNunlikeLowPt = " << block.fNunlikeLowPt
		<< ", fNunlikeHighPt = " << block.fNunlikeHighPt
		<< ", fNlikeAnyPt = " << block.fNlikeAnyPt
		<< ", fNlikeLowPt = " << block.fNlikeLowPt
		<< ", fNlikeHighPt = " << block.fNlikeHighPt
		<< ", fNmassAny = " << block.fNmassAny
		<< ", fNmassLow = " << block.fNmassLow
		<< ", fNmassHigh = " << block.fNmassHigh
		<< ", fDecision[] = [";
	if (block.fHeader.fNrecords > 0) stream << decision[0];
	for (AliHLTUInt32_t i = 1; i < block.fHeader.fNrecords; i++)
		stream << ", " << decision[i];
	stream << "]}";
	return stream;
}


bool operator == (
		const AliHLTMUONPairsDecisionBlockStruct& a,
		const AliHLTMUONPairsDecisionBlockStruct& b
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(a) );
	assert( AliHLTMUONUtils::IntegrityOk(b) );
	
	const AliHLTMUONPairDecisionStruct* decisionA =
		reinterpret_cast<const AliHLTMUONPairDecisionStruct*>(&a + 1);
	const AliHLTMUONPairDecisionStruct* decisionB =
		reinterpret_cast<const AliHLTMUONPairDecisionStruct*>(&b + 1);

	// First check if the blocks have the same header. If they do then
	// check if every track pair's decision is the same. In either case if
	// we find a difference return false.
	if (a.fHeader != b.fHeader) return false;
	if (a.fNunlikeAnyPt != b.fNunlikeAnyPt) return false;
	if (a.fNunlikeLowPt != b.fNunlikeLowPt) return false;
	if (a.fNunlikeHighPt != b.fNunlikeHighPt) return false;
	if (a.fNlikeAnyPt != b.fNlikeAnyPt) return false;
	if (a.fNlikeLowPt != b.fNlikeLowPt) return false;
	if (a.fNlikeHighPt != b.fNlikeHighPt) return false;
	if (a.fNmassAny != b.fNmassAny) return false;
	if (a.fNmassLow != b.fNmassLow) return false;
	if (a.fNmassHigh != b.fNmassHigh) return false;
	for (AliHLTUInt32_t i = 0; i < a.fHeader.fNrecords; i++)
		if (decisionA[i] != decisionB[i]) return false;
	return true;
}
