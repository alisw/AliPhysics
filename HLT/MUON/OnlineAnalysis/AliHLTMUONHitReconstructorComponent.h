#ifndef ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
#define ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
///  @file   AliHLTMUONHitReconstructorComponent.h
///  @author Indranil Das <indra.das@saha.ac.in> | <indra.ehep@gmail.com>
///  @date   
///  @brief  Hit Reconstruction processing component for the dimuon HLT.
///

#include "AliHLTProcessor.h"

#if __GNUC__ && __GNUC__ < 3
#define std
#endif

// TODO: see if we can remove the following header somehow.
#include "AliHLTMUONHitReconstructor.h"
//class AliHLTMUONHitReconstructor;

extern "C" struct AliHLTMUONHitRecoLutRow;


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
	
	void FreeMemory();
	bool ReadLookUpTable(AliHLTMUONHitRecoLutRow* lookupTable, const char* lutpath);
	bool ReadCDB(AliHLTMUONHitRecoLutRow*& lookupTable, const char* cdbpath, Int_t run);
	bool GetLutLine(const char* lutPath, AliHLTUInt32_t& iLine); // To count the nof lines in lookuptable
 
	AliHLTMUONHitReconstructor* fHitRec;  // Internal class instance implementing the hit reconstruction algorithm.
	AliHLTInt32_t fDDL;  // DDL number in the range [13..20]. Set to -1 for invalid/unspecified value.
	IdManuChannelToEntry fIdToEntry; // id to line mapping.
	bool fWarnForUnexpecedBlock;  // Flag indicating if we should log a warning if we got a block of an unexpected type.
	
	ClassDef(AliHLTMUONHitReconstructorComponent, 0)
};

#endif // ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
