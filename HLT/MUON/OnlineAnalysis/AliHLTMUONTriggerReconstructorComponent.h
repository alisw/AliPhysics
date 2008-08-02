#ifndef AliHLTMUONTRIGGERRECONSTRUCTORCOMPONENT_H
#define AliHLTMUONTRIGGERRECONSTRUCTORCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// @file   AliHLTMUONTriggerReconstructorComponent.h
/// @author Indranil Das <indra.das@saha.ac.in>, Artur Szostak <artursz@iafrica.com>
/// @date   18 Sep 2007
/// @brief  A processing component for the dHLT trigger DDL reconstruction.
///

#include "AliHLTMUONProcessor.h"
#include "AliHLTMUONDataTypes.h"

#if __GNUC__ && __GNUC__ < 3
#define std
#endif

class AliHLTMUONTriggerReconstructor;

/**
 * @class AliHLTMUONTriggerReconstructorComponent
 * @brief A processing component for the dHLT trigger DDL reconstruction.
 */
class AliHLTMUONTriggerReconstructorComponent : public AliHLTMUONProcessor
{
public:
	AliHLTMUONTriggerReconstructorComponent();
	virtual ~AliHLTMUONTriggerReconstructorComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	virtual const char* GetComponentID();
	virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
	virtual AliHLTComponent* Spawn();
	
	/**
	 * Generates a binary file containing the lookup table (LUT) from the
	 * CDB, which can be used for the trigger reconstructor component later.
	 * @param ddl  Must be the DDL for which to generate the DDL,
	 *             in the range [20..21].
	 * @param filename  The name of the LUT file to generate.
	 * @param cdbPath  The CDB path to use.
	 * @param run  The run number to use for the CDB.
	 * @param useCrateId  Indicates that the crate ID should be used rather
	 *             than a sequencial number (default is true).
	 * @return  True if the generation of the LUT file succeeded.
	 */
	static bool GenerateLookupTable(
			AliHLTInt32_t ddl, const char* filename,
			const char* cdbPath, Int_t run, bool useCrateId = true
		);

protected:

	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component.
	
	virtual int DoInit(int argc, const char** argv);
	virtual int Reconfigure(const char* cdbEntry, const char* componentId);
	virtual int ReadPreprocessorValues(const char* modules);
	virtual int DoDeinit();

	virtual int DoEvent(
			const AliHLTComponentEventData& evtData,
			const AliHLTComponentBlockData* blocks,
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr,
			AliHLTUInt32_t& size,
			AliHLTComponentBlockDataList& outputBlocks
		);
	
	using AliHLTProcessor::DoEvent;

private:

	// Do not allow copying of this class.
	/// Not implemented.
	AliHLTMUONTriggerReconstructorComponent(const AliHLTMUONTriggerReconstructorComponent& /*obj*/);
	/// Not implemented.
	AliHLTMUONTriggerReconstructorComponent& operator = (const AliHLTMUONTriggerReconstructorComponent& /*obj*/);

	int ReadLookUpTable(const char* lutpath);
	int ReadLutFromCDB();
	
	/**
	 * Reads this component's configuration parameters from the CDB.
	 * These include the middle of the dipole Z coordinate (zmiddle) and the
	 * integrated magnetic field of the dipole.
	 * \param setZmiddle Indicates if the zmiddle parameter should be set
	 *       (default true).
	 * \param setBL Indicates if the integrated magnetic field parameter should
	 *       be set (default true).
	 * \return 0 if no errors occured and negative error code compatible with
	 *       the HLT framework on errors.
	 */
	int ReadConfigFromCDB(bool setZmiddle = true, bool setBL = true);
	
	AliHLTMUONTriggerReconstructor* fTrigRec; ///< The trigger reconstructor class implementing the algorithm.
	AliHLTInt32_t fDDL;   ///< The DDL number in the range 20..21 from which to expect input. Set to -1 for invalid/unspecified value.
	bool fWarnForUnexpecedBlock;  ///< Flag indicating if we should log a warning if we got a block of an unexpected type.
	bool fStopOnOverflow;  ///< Flag indicating if we should fail in the DoEvent method if the output buffer was overflowed.
	bool fUseCrateId;  ///< Flag to indicate if the crate ID as found in the regional header structures should be used or not.
	bool fZmiddleSpecified;  ///< Indicates if the zmiddle parameter was specified on the command line.
	bool fBLSpecified;  ///< Indicates if the bfieldintegral parameter was specified on the command line.
	bool fLutInitialised;  ///< Flag to indicate if the LUT was loaded yet or not.

	ClassDef(AliHLTMUONTriggerReconstructorComponent, 0) // Trigger reconstructor component for dHLT trigger DDL raw data.

};

#endif // AliHLTMUONTRIGGERRECONSTRUCTORCOMPONENT_H
