#ifndef AliHLTMUONTRIGGERRECONSTRUCTOR_H
#define AliHLTMUONTRIGGERRECONSTRUCTOR_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/**********************************************************************
 Created on : 16/05/2007
 Purpose    : This class is supposed to read the trigger DDL files and 
              give the output AliHLTMUONTriggerRecordStruct
 Author     : Indranil Das, HEP Division, SINP
 Email      : indra.das@saha.ac.in | indra.ehep@gmail.com

 Artur Szostak <artursz@iafrica.com>:
  Completely reimplemented the lookup table to a simplified format.
**********************************************************************/

#include "AliHLTLogging.h"
#include "AliHLTMUONTriggerRecordsBlockStruct.h"
#include "AliHLTMUONHitReconstructor.h"


class AliHLTMUONTriggerReconstructor : public AliHLTLogging
{
public:

	AliHLTMUONTriggerReconstructor();
	virtual ~AliHLTMUONTriggerReconstructor();

	bool Run(
			const AliHLTUInt32_t* rawData,
			AliHLTUInt32_t rawDataSize,
			AliHLTMUONTriggerRecordStruct* trigRecord,
			AliHLTUInt32_t& nofTrigRec,
			bool suppressPartialTrigs = false
		);

	void* LookupTableBuffer() { return &fLookupTable; }
	size_t LookupTableSize() const { return sizeof(fLookupTable); }

private:

	AliHLTMUONTriggerReconstructor(const AliHLTMUONTriggerReconstructor& rhs); // copy constructor
	AliHLTMUONTriggerReconstructor& operator=(const AliHLTMUONTriggerReconstructor& rhs); // assignment operator

	struct LookupTableRow
	{
		float fX, fY, fZ;
	};

	AliHLTUInt32_t fMaxRecPointsCount;   // max nof reconstructed hit
	AliHLTInt32_t fTrigRecId;  // A running counter for the trigger record ID.
	
	// [regional header index][local board ID][chamber][cathode - X/Y][bit set in bit pattern]
	LookupTableRow fLookupTable[8][16][4][2][16];  // pointer to the array of Lookuptable data
};

#endif // AliHLTMUONTRIGGERRECONSTRUCTOR_H
