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
 * @file   AliHLTMUONRecHitsBlockStruct.h
 * @author Indranil Das <indra.das@saha.ac.in>,
 *         Artur Szostak <artursz@iafrica.com>
 * @date   17 May 2007
 * @brief  Implementation of useful stream and comparison operators.
 */

#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include <cassert>

std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONRecHitsBlockStruct& block
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(block) );

	const AliHLTMUONRecHitStruct* hit =
		reinterpret_cast<const AliHLTMUONRecHitStruct*>(&block + 1);
	stream 	<< "{fHeader = " << block.fHeader << ", fHit[] = [";
	if (block.fHeader.fNrecords > 0) stream << hit[0];
	for (AliHLTUInt32_t i = 1; i < block.fHeader.fNrecords; i++)
		stream << ", " << hit[i];
	stream << "]}";
	return stream;
}


bool operator == (
		const AliHLTMUONRecHitsBlockStruct& a,
		const AliHLTMUONRecHitsBlockStruct& b
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(a) );
	assert( AliHLTMUONUtils::IntegrityOk(b) );
	
	const AliHLTMUONRecHitStruct* hitA =
		reinterpret_cast<const AliHLTMUONRecHitStruct*>(&a + 1);
	const AliHLTMUONRecHitStruct* hitB =
		reinterpret_cast<const AliHLTMUONRecHitStruct*>(&b + 1);

	// First check if the blocks have the same header. If they do then check if
	// every hit is the same. In either case if we find a difference return false.
	if (a.fHeader != b.fHeader) return false;
	for (AliHLTUInt32_t i = 0; i < a.fHeader.fNrecords; i++)
		if (hitA[i] != hitB[i]) return false;
	return true;
}
