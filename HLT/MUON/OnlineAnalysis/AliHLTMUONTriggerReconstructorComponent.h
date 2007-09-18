#ifndef AliHLTMUONTRIGGERRECONSTRUCTORCOMPONENT_H
#define AliHLTMUONTRIGGERRECONSTRUCTORCOMPONENT_H
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

/** @file   AliHLTMUONTriggerReconstructorComponent.h
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  Declaration of a dummy component. */


#include "AliHLTProcessor.h"
#include "AliHLTMUONConstants.h"

#include "AliHLTMUONTriggerReconstructor.h"


/**
 * @class AliHLTMUONTriggerReconstructorComponent
 * @brief A dummy HLT processing component. 
 *
 * An implementiation of a dummy component that just copies its input data
 * as a test, demonstration, and example of the HLT component scheme.
 * @ingroup alihlt_tutorial
 */
class AliHLTMUONTriggerReconstructorComponent : public AliHLTProcessor
    {
    public:
	AliHLTMUONTriggerReconstructorComponent();
	virtual ~AliHLTMUONTriggerReconstructorComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	const char* GetComponentID();
	void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
	AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
	AliHLTComponent* Spawn();
	
    protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	int DoInit( int argc, const char** argv );
	int DoDeinit();
	int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
	
    private:

	// The size of the output data produced, as a percentage of the input data's size.
	// Can be greater than 100 (%)
	unsigned fOutputPercentage;
	
	AliHLTMUONTriggerReconstructor* fTrigRec;
	bool ReadLookUpTable(AliHLTMUONHitReconstructor::DHLTLut* lookupTable, const char* lutpath);

	TString fDDLDir;
	Int_t fDDL;

	ClassDef(AliHLTMUONTriggerReconstructorComponent, 0)

    };

#endif // AliHLTMUONTRIGGERRECONSTRUCTORCOMPONENT_H
