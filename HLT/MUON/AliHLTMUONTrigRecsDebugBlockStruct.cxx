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
 * @file   AliHLTMUONTrigRecsDebugBlockStruct.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   19 May 2007
 * @brief  Implementation of useful stream and comparison operators.
 */

#include "AliHLTMUONTrigRecsDebugBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include <cassert>


std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONTrigRecInfoStruct& info
	)
{
	stream	<< "{fTrigRecId = " << info.fTrigRecId
		<< ", fDetElemId = [" << info.fDetElemId[0]
		<< ", " << info.fDetElemId[1]
		<< ", " << info.fDetElemId[2]
		<< ", " << info.fDetElemId[3]
		<< "], fZmiddle = " << info.fZmiddle
		<< ", fBl = " << info.fBl
		<< ", fL0Struct = {fX2X1 = " << std::showbase << std::hex
		<< info.fL0Struct.fX2X1
		<< ", fX4X3 = " << info.fL0Struct.fX4X3
		<< ", fY2Y1 = " << info.fL0Struct.fY2Y1
		<< ", fY4Y3 = " << info.fL0Struct.fY4Y3
		<< ", fTriggerBits = " << info.fL0Struct.fTriggerBits
		<< std::dec << "}"
		<< ", fL0StructPrev = {fX2X1 = " << std::showbase << std::hex
		<< info.fL0StructPrev.fX2X1
		<< ", fX4X3 = " << info.fL0StructPrev.fX4X3
		<< ", fY2Y1 = " << info.fL0StructPrev.fY2Y1
		<< ", fY4Y3 = " << info.fL0StructPrev.fY4Y3
		<< ", fTriggerBits = " << info.fL0StructPrev.fTriggerBits
		<< std::dec << "}"
		<< ", fL0StructNext = {fX2X1 = " << std::showbase << std::hex
		<< info.fL0StructNext.fX2X1
		<< ", fX4X3 = " << info.fL0StructNext.fX4X3
		<< ", fY2Y1 = " << info.fL0StructNext.fY2Y1
		<< ", fY4Y3 = " << info.fL0StructNext.fY4Y3
		<< ", fTriggerBits = " << info.fL0StructNext.fTriggerBits
		<< std::dec << "}}";
	return stream;
}


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONTrigRecsDebugBlockStruct& block
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(block) );

	const AliHLTMUONTrigRecInfoStruct* trigRecInfo =
		reinterpret_cast<const AliHLTMUONTrigRecInfoStruct*>(&block + 1);
	stream 	<< "{fHeader = " << block.fHeader << ", fTrigRecInfo[] = [";
	if (block.fHeader.fNrecords > 0) stream << trigRecInfo[0];
	for (AliHLTUInt32_t i = 1; i < block.fHeader.fNrecords; i++)
		stream << ", " << trigRecInfo[i];
	stream << "]}";
	return stream;
}


bool operator == (
		const AliHLTMUONTrigRecInfoStruct& a,
		const AliHLTMUONTrigRecInfoStruct& b
	)
{
	return a.fTrigRecId == b.fTrigRecId and a.fDetElemId == b.fDetElemId
		and a.fZmiddle == b.fZmiddle and a.fBl == b.fBl
		and a.fL0Struct.fX2X1 == b.fL0Struct.fX2X1
		and a.fL0Struct.fX4X3 == b.fL0Struct.fX4X3
		and a.fL0Struct.fY2Y1 == b.fL0Struct.fY2Y1
		and a.fL0Struct.fY4Y3 == b.fL0Struct.fY4Y3
		and a.fL0Struct.fTriggerBits == b.fL0Struct.fTriggerBits
		and a.fL0StructPrev.fX2X1 == b.fL0StructPrev.fX2X1
		and a.fL0StructPrev.fX4X3 == b.fL0StructPrev.fX4X3
		and a.fL0StructPrev.fY2Y1 == b.fL0StructPrev.fY2Y1
		and a.fL0StructPrev.fY4Y3 == b.fL0StructPrev.fY4Y3
		and a.fL0StructPrev.fTriggerBits == b.fL0StructPrev.fTriggerBits
		and a.fL0StructNext.fX2X1 == b.fL0StructNext.fX2X1
		and a.fL0StructNext.fX4X3 == b.fL0StructNext.fX4X3
		and a.fL0StructNext.fY2Y1 == b.fL0StructNext.fY2Y1
		and a.fL0StructNext.fY4Y3 == b.fL0StructNext.fY4Y3
		and a.fL0StructNext.fTriggerBits == b.fL0StructNext.fTriggerBits;
}


bool operator == (
		const AliHLTMUONTrigRecsDebugBlockStruct& a,
		const AliHLTMUONTrigRecsDebugBlockStruct& b
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(a) );
	assert( AliHLTMUONUtils::IntegrityOk(b) );

	const AliHLTMUONTrigRecInfoStruct* trigRecInfoA =
		reinterpret_cast<const AliHLTMUONTrigRecInfoStruct*>(&a + 1);
	const AliHLTMUONTrigRecInfoStruct* trigRecInfoB =
		reinterpret_cast<const AliHLTMUONTrigRecInfoStruct*>(&b + 1);
	
	// First check if the blocks have the same header. If they do then check
	// if all debug information is the same. In either case if we find a
	// difference return false.
	if (a.fHeader != b.fHeader) return false;
	for (AliHLTUInt32_t i = 0; i < a.fHeader.fNrecords; i++)
		if (trigRecInfoA[i] != trigRecInfoB[i]) return false;
	return true;
}
