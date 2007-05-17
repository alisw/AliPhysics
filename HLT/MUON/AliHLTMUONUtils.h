#ifndef ALIHLTMUONUTILS_H
#define ALIHLTMUONUTILS_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/**
 * @file   AliHLTMUONUtils.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Class containing various dimuon HLT utility routines.
 */

#include "AliHLTMUONRecHitsDebugBlockStruct.h"

/**
 * AliHLTMUONUtils contains arbitrary utility methods to be used in various
 * parts of the dimuon HLT system.
 * These include methods to perform basic sanity checks on the integrity of
 * data blocks.
 */
class AliHLTMUONUtils
{
public:

	/**
	 * Checks if the integrity of the block is Ok and returns true in that case.
	 */
	static bool IntegrityOk(const AliHLTMUONRecHitsDebugBlockStruct& block);

private:
	// Should never have to create or destroy this object.
	AliHLTMUONUtils();
	~AliHLTMUONUtils();
};

#endif // ALIHLTMUONUTILS_H
