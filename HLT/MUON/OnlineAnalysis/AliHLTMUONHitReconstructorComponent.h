#ifndef ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
#define ALIHLTMUONHITREONSTRUCTORCCOMPONENT_H
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

///*  @file   AliHLTMUONHitReconstructorComponent.h
// *  @author Indranil Das <indra.das@saha.ac.in> | <indra.ehep@gmail.com>
// *  @date   
// *  @brief  Hit Reconstruction processing component for the dimuon HLT. 
// */

#include "AliHLTProcessor.h"
#include <TString.h>
#include "AliHLTMUONHitReconstructor.h"

#if __GNUC__ < 3
#define std
#endif


class AliHLTMUONHitReconstructorComponent : public AliHLTProcessor
{
    public:
	AliHLTMUONHitReconstructorComponent();
	virtual ~AliHLTMUONHitReconstructorComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	const char* GetComponentID();
	void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);
	AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
	AliHLTComponent* Spawn();
	
    protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	int DoInit(int argc, const char** argv);
	int DoDeinit();
	int DoEvent(
			const AliHLTComponentEventData& evtData,
			const AliHLTComponentBlockData* blocks,
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr, 
			AliHLTUInt32_t& size,
			std::vector<AliHLTComponentBlockData>& outputBlocks
		);
	
    private:
	
	AliHLTMUONHitReconstructor* fHitRec;
	bool ReadLookUpTable(AliHLTMUONHitReconstructor::DHLTLut* lookupTable, const char* lutpath);
	bool ReadBusPatchToDetElemFile(BusToDetElem& busToDetElem, BusToDDL& busToDDL, const char* buspatchmappath);

	TString fDDLDir;
	Int_t fDDL;
	bool fReaderType;

	ClassDef(AliHLTMUONHitReconstructorComponent, 0)
};

#endif // ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
