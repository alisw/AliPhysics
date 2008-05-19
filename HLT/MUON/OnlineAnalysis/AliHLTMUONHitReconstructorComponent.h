#ifndef ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
#define ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
///  @file   AliHLTMUONHitReconstructorComponent.h
///  @author Indranil Das <indra.das@saha.ac.in> | <indra.ehep@gmail.com>
///  @date   28 May 2007
///  @brief  Hit Reconstruction processing component for the dimuon HLT.
///

#include "AliHLTMUONProcessor.h"
#include "AliHLTMUONHitReconstructor.h"

#if __GNUC__ && __GNUC__ < 3
#define std
#endif


extern "C" struct AliHLTMUONHitRecoLutRow;

/**
 * @class AliHLTMUONHitReconstructorComponent
 * @brief A processing component for the dHLT tracker DDL reconstruction.
 */
class AliHLTMUONHitReconstructorComponent : public AliHLTMUONProcessor
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
	
	/**
	 * Generates a ASCII text file containing the lookup table (LUT) from
	 * the CDB, which can be used for the hit reconstructor component later.
	 * @param ddl  Must be the DDL for which to generate the DDL,
	 *             in the range [12..19].
	 * @param filename  The name of the LUT file to generate.
	 * @param cdbPath  The CDB path to use.
	 * @param run  The run number to use for the CDB.
	 * @return  True if the generation of the LUT file succeeded.
	 */
	static bool GenerateLookupTable(
			AliHLTInt32_t ddl, const char* filename,
			const char* cdbPath, Int_t run
		);
	
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
	int ReadLookUpTable(const char* lutpath);
	int ReadCDB(const char* cdbpath, Int_t run);
	
	AliHLTMUONHitReconstructor* fHitRec;  ///< Internal class instance implementing the hit reconstruction algorithm.
	AliHLTInt32_t fDDL;  ///< DDL number in the range [12..19]. Set to -1 for invalid/unspecified value.
	AliHLTUInt32_t fLutSize;  ///< The number of rows / entries in the LUT.
	AliHLTMUONHitRecoLutRow* fLut;  ///< The lookup table used by the hit reconstruction algorithm (Owned by this component however).
	IdManuChannelToEntry fIdToEntry; ///< id to line mapping.
	bool fWarnForUnexpecedBlock;  ///< Flag indicating if we should log a warning if we got a block of an unexpected type.
	
	ClassDef(AliHLTMUONHitReconstructorComponent, 0) // Hit reconstructor component for dHLT tracker DDL raw data.
};

#endif // ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
