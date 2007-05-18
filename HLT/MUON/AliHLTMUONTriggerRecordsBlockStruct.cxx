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
 * @file   AliHLTMUONTriggerRecordsBlockStruct.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Implementation of useful stream and comparison operators.
 */

#include "AliHLTMUONTriggerRecordsBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include <cassert>


std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONTriggerRecordStruct& trigrec
	)
{
	stream	<< "{fId = " << trigrec.fFlags
		<< ", fFlags = " << std::showbase << std::hex
		<< trigrec.fFlags << std::dec
		<< ", fPx = " << trigrec.fPx
		<< ", fPy = " << trigrec.fPy
		<< ", fPz = " << trigrec.fPz
		<< ", fHit[0] = " << trigrec.fHit[0]
		<< ", fHit[1] = " << trigrec.fHit[1]
		<< ", fHit[2] = " << trigrec.fHit[2]
		<< ", fHit[3] = " << trigrec.fHit[3]
		<< "}";
	return stream;
}


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONTriggerRecordsBlockStruct& block
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(block) );

	stream 	<< "{fHeader = " << block.fHeader << ", fTriggerRecord[] = [";
	if (block.fHeader.fNrecords > 0) stream << block.fTriggerRecord[0];
	for (AliHLTUInt32_t i = 1; i < block.fHeader.fNrecords; i++)
		stream << ", " << block.fTriggerRecord[i];
	stream << "]}";
	return stream;
}


bool operator == (
		const AliHLTMUONTriggerRecordsBlockStruct& a,
		const AliHLTMUONTriggerRecordsBlockStruct& b
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(a) );
	assert( AliHLTMUONUtils::IntegrityOk(b) );

	// First check if the blocks have the same header. If they do then check
	// if every trigger record is the same. In either case if we find a
	// difference return false.
	if (a.fHeader != b.fHeader) return false;
	for (AliHLTUInt32_t i = 0; i < a.fHeader.fNrecords; i++)
		if (a.fTriggerRecord[i] != b.fTriggerRecord[i]) return false;
	return true;
}
