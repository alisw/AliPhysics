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
#include "AliHLTMUONDataTypes.h"

extern "C" struct AliHLTMUONTriggerRecordStruct;


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

	AliHLTMUONTriggerRecoLookupTable* LookupTableBuffer() { return &fLookupTable; }
	size_t LookupTableSize() const { return sizeof(fLookupTable); }

private:

	AliHLTMUONTriggerReconstructor(const AliHLTMUONTriggerReconstructor& rhs); // copy constructor
	AliHLTMUONTriggerReconstructor& operator=(const AliHLTMUONTriggerReconstructor& rhs); // assignment operator

	AliHLTUInt32_t fMaxRecPointsCount;   // max nof reconstructed hit
	AliHLTInt32_t fTrigRecId;  // A running counter for the trigger record ID.
	AliHLTMUONTriggerRecoLookupTable fLookupTable;  // The lookup table used for mapping between channel addresses and geometrical information.
};

#endif // AliHLTMUONTRIGGERRECONSTRUCTOR_H
