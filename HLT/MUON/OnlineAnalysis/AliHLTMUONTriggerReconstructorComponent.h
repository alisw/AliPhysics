#ifndef AliHLTMUONTRIGGERRECONSTRUCTORCOMPONENT_H
#define AliHLTMUONTRIGGERRECONSTRUCTORCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/**
 * @file   AliHLTMUONTriggerReconstructorComponent.h
 * @author Indranil Das <indra.das@saha.ac.in>, Artur Szostak <artursz@iafrica.com>
 * @date
 * @brief  A processing component for the dHLT trigger DDL reconstruction.
 */

#include "AliHLTProcessor.h"
#include "AliHLTMUONDataTypes.h"

#if __GNUC__ && __GNUC__ < 3
#define std
#endif

class AliHLTMUONTriggerReconstructor;

/**
 * @class AliHLTMUONTriggerReconstructorComponent
 * @brief A processing component for the dHLT trigger DDL reconstruction.
 */
class AliHLTMUONTriggerReconstructorComponent : public AliHLTProcessor
{
public:
	AliHLTMUONTriggerReconstructorComponent();
	virtual ~AliHLTMUONTriggerReconstructorComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	const char* GetComponentID();
	void GetInputDataTypes( std::vector<AliHLTComponentDataType>& list);
	AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
	AliHLTComponent* Spawn();

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

private:

	bool ReadLookUpTable(const char* lutpath);
	
	AliHLTMUONTriggerReconstructor* fTrigRec; // The trigger reconstructor class implementing the algorithm.
	AliHLTInt32_t fDDL;   // The DDL number in the range 20..21 from which to expect input.
	bool fWarnForUnexpecedBlock;  // Flag indicating if we should log a warning if we got a block of an unexpected type.
	bool fSuppressPartialTrigs;   // Flag indicating if we should suppress triggers that did not trigger the L0

	ClassDef(AliHLTMUONTriggerReconstructorComponent, 0)

};

#endif // AliHLTMUONTRIGGERRECONSTRUCTORCOMPONENT_H
