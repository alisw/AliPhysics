#ifndef ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
#define ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
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

///
///  @file   AliHLTMUONHitReconstructorComponent.h
///  @author Indranil Das <indra.das@saha.ac.in> | <indra.ehep@gmail.com>
///  @date   
///  @brief  Hit Reconstruction processing component for the dimuon HLT. 
///

#include "AliHLTProcessor.h"
#include <TString.h>
#include "AliHLTMUONHitReconstructor.h"

#if __GNUC__ && __GNUC__ < 3
#define std
#endif


class AliHLTMUONHitReconstructorComponent : public AliHLTProcessor
{
    public:
	AliHLTMUONHitReconstructorComponent();
	virtual ~AliHLTMUONHitReconstructorComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	virtual const char* GetComponentID();
	virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
	virtual AliHLTComponent* Spawn();
	
    protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	virtual int DoInit(int argc, const char** argv);
	virtual int DoDeinit();
	virtual int DoEvent(
			const AliHLTComponentEventData& evtData,
			const AliHLTComponentBlockData* blocks,
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr, 
			AliHLTUInt32_t& size,
			std::vector<AliHLTComponentBlockData>& outputBlocks
		);
	
	using AliHLTProcessor::DoEvent;
	
    private:

	// Do not allow copying of this class.
	AliHLTMUONHitReconstructorComponent(const AliHLTMUONHitReconstructorComponent& /*obj*/);
	AliHLTMUONHitReconstructorComponent& operator = (const AliHLTMUONHitReconstructorComponent& /*obj*/);
	
	AliHLTMUONHitReconstructor* fHitRec;   // Internal class instance implementing the hit reconstruction algorithm.
	bool ReadLookUpTable(AliHLTMUONHitReconstructor::DHLTLut* lookupTable, const char* lutpath);
	bool ReadBusPatchToDetElemFile(BusToDetElem& busToDetElem, BusToDDL& busToDDL, const char* buspatchmappath);

	TString fDDLDir;   //TODO: Deprecated, should be removed.
	Int_t fDDL;        // DDL number in the range [12..19].
	bool fReaderType;  //TODO: Deprecated, should be removed.
	bool fWarnForUnexpecedBlock;  // Flag indicating if we should log a warning if we got a block of an unexpected type.

	ClassDef(AliHLTMUONHitReconstructorComponent, 0)
};

#endif // ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
