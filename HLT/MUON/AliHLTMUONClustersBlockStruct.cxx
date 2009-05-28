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
 * @file   AliHLTMUONClustersBlockStruct.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   18 May 2007
 * @brief  Implementation of the stream and comparison operators.
 */

#include "AliHLTMUONClustersBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include <cassert>


std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONClusterStruct& cluster
	)
{
	stream	<< "{fId = " << cluster.fId
		<< ", fHit = " << cluster.fHit
		<< ", fDetElemId = " << cluster.fDetElemId
		<< ", fNchannels = " << cluster.fNchannels
		<< ", fCharge = " << cluster.fCharge << "}";
	return stream;
}


std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONClustersBlockStruct& block
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(block) );

	const AliHLTMUONClusterStruct* cluster =
		reinterpret_cast<const AliHLTMUONClusterStruct*>(&block + 1);
	stream 	<< "{fHeader = " << block.fHeader << ", fCluster[] = [";
	if (block.fHeader.fNrecords > 0) stream << cluster[0];
	for (AliHLTUInt32_t i = 1; i < block.fHeader.fNrecords; i++)
		stream << ", " << cluster[i];
	stream << "]}";
	return stream;
}


bool operator == (
		const AliHLTMUONClustersBlockStruct& a,
		const AliHLTMUONClustersBlockStruct& b
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(a) );
	assert( AliHLTMUONUtils::IntegrityOk(b) );
	
	const AliHLTMUONClusterStruct* clusterA =
		reinterpret_cast<const AliHLTMUONClusterStruct*>(&a + 1);
	const AliHLTMUONClusterStruct* clusterB =
		reinterpret_cast<const AliHLTMUONClusterStruct*>(&b + 1);

	// First check if the blocks have the same header. If they do then check
	// if every cluster is the same. In either case if we find a difference
	// return false.
	if (a.fHeader != b.fHeader) return false;
	for (AliHLTUInt32_t i = 0; i < a.fHeader.fNrecords; i++)
		if (clusterA[i] != clusterB[i]) return false;
	return true;
}
