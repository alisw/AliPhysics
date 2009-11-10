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
 * @file   AliHLTMUONMansoCandidatesBlockStruct.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   21 May 2007
 * @brief  Implementation of useful stream and comparison operators.
 */

#include "AliHLTMUONMansoCandidatesBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include <cassert>


std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONMansoRoIStruct& roi
	)
{
	stream	<< "{fX = " << roi.fX
		<< ", fY = " << roi.fY
		<< ", fZ = " << roi.fZ
		<< ", fRadius = " << roi.fRadius
		<< "}";
	return stream;
}

std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONMansoCandidateStruct& candidate
	)
{
	stream	<< "{fTrack = " << candidate.fTrack
		<< ", fRoI[0] = " << candidate.fRoI[0]
		<< ", fRoI[1] = " << candidate.fRoI[1]
		<< ", fRoI[2] = " << candidate.fRoI[2]
		<< ", fRoI[3] = " << candidate.fRoI[3]
		<< ", fZmiddle = " << candidate.fZmiddle
		<< ", fBl = " << candidate.fBl
		<< "}";
	return stream;
}


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONMansoCandidatesBlockStruct& block
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(block) );

	const AliHLTMUONMansoCandidateStruct* candidate =
		reinterpret_cast<const AliHLTMUONMansoCandidateStruct*>(&block + 1);
	stream 	<< "{fHeader = " << block.fHeader << ", fCandidate[] = [";
	if (block.fHeader.fNrecords > 0) stream << candidate[0];
	for (AliHLTUInt32_t i = 1; i < block.fHeader.fNrecords; i++)
		stream << ", " << candidate[i];
	stream << "]}";
	return stream;
}


bool operator == (
		const AliHLTMUONMansoCandidatesBlockStruct& a,
		const AliHLTMUONMansoCandidatesBlockStruct& b
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(a) );
	assert( AliHLTMUONUtils::IntegrityOk(b) );
	
	const AliHLTMUONMansoCandidateStruct* candidateA =
		reinterpret_cast<const AliHLTMUONMansoCandidateStruct*>(&a + 1);
	const AliHLTMUONMansoCandidateStruct* candidateB =
		reinterpret_cast<const AliHLTMUONMansoCandidateStruct*>(&b + 1);

	// First check if the blocks have the same header. If they do then check
	// if every track candidate is the same. In either case if we find a
	// difference return false.
	if (a.fHeader != b.fHeader) return false;
	for (AliHLTUInt32_t i = 0; i < a.fHeader.fNrecords; i++)
		if (candidateA[i] != candidateB[i]) return false;
	return true;
}
