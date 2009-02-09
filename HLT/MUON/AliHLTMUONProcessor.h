#ifndef ALIHLTMUONPROCESSOR_H
#define ALIHLTMUONPROCESSOR_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///
/// @file   AliHLTMUONProcessor.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 May 2008
/// @brief  Declaration of a common processor component abstract interface for dHLT components.
///

#include "AliHLTProcessor.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTMUONUtils.h"

class TMap;
class AliMUONRecoParam;

/**
 * @class AliHLTMUONProcessor
 * This component class is an abstract base class for dHLT components.
 * Some common methods useful to all dHLT specific components are implemented
 * by this class.
 *
 * @ingroup alihlt_dimuon_component
 */
class AliHLTMUONProcessor : public AliHLTProcessor
{
public:
	/// Default constructor.
	AliHLTMUONProcessor();
	
	/// Default destructor.
	virtual ~AliHLTMUONProcessor() {}

protected:

	/**
	 * This method parses the common arguments for dHLT processing components
	 * and initialises the common internal state.
	 * Deriving classes can use the ArgumentAlreadyHandled method to check if
	 * the parent class has processed a particular argument. The following is
	 * an example of this:
	 *
	 * \code
	 * int DerivedClass::DoInit(int argc, const char** argv)
	 * {
	 *   int result = AliHLTMUONProcessor::DoInit(argc, argv);
	 *   if (result != 0) return result;
	 *   for (int i = 0; i < argc; i++)
	 *   {
	 *     if (ArgumentAlreadyHandled(i, argv[i])) continue;
	 *     // ... handle custom arguments here ...
	 *   }
	 * }
	 * \endcode
	 */
	virtual int DoInit(int argc, const char** argv);

	/**
	 * This method can be used by the derivind child class to check if a particular
	 * argument in argv was already processed.
	 * \note This assumes that the deriving class called the DoInit method of the
	 * parent class in its own DoInit method.
	 */
	virtual bool ArgumentAlreadyHandled(int& i, const char* argi) const;
	
	/**
	 * This method returns the command line arguments that should not be parsed
	 * by this class. This method can be used by child classes that derive from
	 * AliHLTMUONProcessor, to indicate which arguments should not be handled by
	 * the AliHLTMUONProcessor::DoInit method. Default return value is false.
	 */
	virtual bool IgnoreArgument(const char* /*arg*/) const { return false; }
	
	/**
	 * Returns true if the component was told to delay initialisation from
	 * CDB until the first start of run event. This gets set by the -delaysetup
	 * flag which is processed in AliHLTMUONProcessor::DoInit.
	 */
	bool DelaySetup() const { return fDelaySetup; }

	/**
	 * This method should be called when a derived component has handled a
	 * delayed setup requested on the command line with -delaysetup and indicated
	 * by the flag returned by the DelaySetup method.
	 */
	void DoneDelayedSetup() { fDelaySetup = false; }
	
	/**
	 * Returns true if the component has the flag set indicating to dump raw
	 * data when an error occurs. The DumpEvent method should be used by the
	 * deriving components to actually dump data at the appropriate point.
	 * \note This facility is intended for debugging.
	 */
	bool DumpDataOnError() const { return fDumpDataOnError; }

	/**
	 * Returns the path where the dump files will be written to by the Dump*
	 * methods. Defaults to the current working directory.
	 * \note This facility is intended for debugging.
	 */
	const char* DumpPath() const { return fDumpPath; }

	/**
	 * Method to check the block structure and log appropriate error messages.
	 * If a problem is found with the data block then an appropriate HLT error
	 * message is logged and the method returns false.
	 * \param block  The lightweight block reader whose data block should be checked.
	 * \param name  A string containing a descriptive name of the data block
	 *          type. This name is used in the logged error messages.
	 * \param checkHeader  Indicates if the common data block header should be checked.
	 * \returns  true if the structure of the block looks OK and false otherwise.
	 * \note  The BlockType should be a class deriving from AliHLTMUONDataBlockReader.
	 */
	template <class BlockType>
	bool BlockStructureOk(
			const BlockType& block, const char* name,
			bool checkHeader = true
		) const;
	
	/// Checks the structure of a trigger records data block.
	bool BlockStructureOk(
			const AliHLTMUONTriggerRecordsBlockReader& block,
			bool checkHeader = true
		) const
	{
		return BlockStructureOk(block, "trigger records", checkHeader);
	}
	
	/// Checks the structure of a trigger records debug information data block.
	bool BlockStructureOk(
			const AliHLTMUONTrigRecsDebugBlockReader& block,
			bool checkHeader = true
		) const
	{
		return BlockStructureOk(block, "trigger records debug information", checkHeader);
	}

	/// Checks the structure of a reconstructed hits data block.
	bool BlockStructureOk(
			const AliHLTMUONRecHitsBlockReader& block,
			bool checkHeader = true
		) const
	{
		return BlockStructureOk(block, "reconstructed hits", checkHeader);
	}
	
	/// Checks the structure of a clusters data block.
	bool BlockStructureOk(
			const AliHLTMUONClustersBlockReader& block,
			bool checkHeader = true
		) const
	{
		return BlockStructureOk(block, "clusters", checkHeader);
	}
	
	/// Checks the structure of a ADC channels data block.
	bool BlockStructureOk(
			const AliHLTMUONChannelsBlockReader& block,
			bool checkHeader = true
		) const
	{
		return BlockStructureOk(block, "channels", checkHeader);
	}

	/// Checks the structure of a Manso tracks data block.
	bool BlockStructureOk(
			const AliHLTMUONMansoTracksBlockReader& block,
			bool checkHeader = true
		) const
	{
		return BlockStructureOk(block, "Manso tracks", checkHeader);
	}
	
	/// Checks the structure of a Manso track candidates data block.
	bool BlockStructureOk(
			const AliHLTMUONMansoCandidatesBlockReader& block,
			bool checkHeader = true
		) const
	{
		return BlockStructureOk(block, "Manso track candidates", checkHeader);
	}

	/// Checks the structure of a single track trigger decision data block.
	bool BlockStructureOk(
			const AliHLTMUONSinglesDecisionBlockReader& block,
			bool checkHeader = true
		) const
	{
		return BlockStructureOk(block, "singles decision", checkHeader);
	}

	/// Checks the structure of a track pairs trigger decision data block.
	bool BlockStructureOk(
			const AliHLTMUONPairsDecisionBlockReader& block,
			bool checkHeader = true
		) const
	{
		return BlockStructureOk(block, "pairs decision", checkHeader);
	}
	
	/**
	 * Sets the CDB path and run number to read from.
	 * \param cdbPath  The CDB path to use. If set to NULL and the path has
	 *      not been set in the CDB manager then the default path
	 *      "local://$ALICE_ROOT/OCDB" is used if the 'useDefault' flag is also true.
	 * \param run  The run number to use. If set to -1 and the run number has
	 *      not been set in the CDB manager then a value of zero is used if
	 *      the 'useDefault' flag is also true.
	 * \param useDefault  If set to true then a default CDB path and/or run number
	 *      is used if they have not been set and 'cdbPath' == NULL or
	 *      'run' == -1. (false by default).
	 * \return Zero if the object could be loaded. Otherwise an error code,
	 *      compatible with the HLT framework, is returned.
	 */
	int SetCDBPathAndRunNo(
			const char* cdbPath, Int_t run, bool useDefault = false
		) const;
	
	/**
	 * Fetches the DDL and detector element store objects for MUON mapping.
	 * \return Zero if the objects could be loaded. Otherwise an error code,
	 *      which is compatible with the HLT framework, is returned.
	 * \note AliMpDDLStore::Instance() and AliMpDEStore::Instance() must be used
	 *      to fetch the objects after this method returns a code equal to zero.
	 */
	int FetchMappingStores() const;
	
	/**
	 * Fetches a TMap object from the CDB.
	 * [in] \param pathToEntry  The relative path to the entry in the CDB to fetch.
	 * [out] \param map  This will be filled with the TMap object found if
	 *      a successful status code is returned. Otherwise it will be unchanged.
	 * \return Zero if the object could be found. Otherwise an error code,
	 *      which is compatible with the HLT framework, is returned.
	 */
	int FetchTMapFromCDB(const char* pathToEntry, TMap*& map) const;
	
	/**
	 * Tries to find the string value associated with a certain parameter in a TMap.
	 * [in] \param map  The TMap object to search in.
	 * [in] \param paramName  The name of the parameter to search for.
	 * [out] \param value  Will be filled with the object found.
	 * [in] \param prettyName  Should be the name of the parameter which will
	 *      be used when printing error messages. If this is set to NULL then
	 *      the paramName will be used instead (default is NULL).
	 * \return Zero if the object could be found. Otherwise an error code,
	 *      which is compatible with the HLT framework, is returned.
	 */
	int GetValueFromTMap(
			TMap* map, const char* paramName, TString& value,
			const char* pathToEntry = NULL, const char* prettyName = NULL
		) const;
	
	/**
	 * Tries to find a certain parameter in the TMap object and convert it to
	 * an integer value.
	 * [in] \param map  The TMap object to search in.
	 * [in] \param paramName  The name of the parameter to search for.
	 * [out] \param value  Will be filled with the integer value for the parameter,
	 *       if it was found and it was an integer value.
	 * [in] \param prettyName  Should be the name of the parameter which will
	 *      be used when printing error messages. If this is set to NULL then
	 *      the paramName will be used instead (default is NULL).
	 * \return Zero if the object could be found and is valid. Otherwise an
	 *       error code, which is compatible with the HLT framework, is returned.
	 */
	int GetIntFromTMap(
			TMap* map, const char* paramName, Int_t& value,
			const char* pathToEntry = NULL, const char* prettyName = NULL
		) const;
	
	/**
	 * Tries to find a certain parameter in the TMap object and convert it to
	 * a positive integer value.
	 * [in] \param map  The TMap object to search in.
	 * [in] \param paramName  The name of the parameter to search for.
	 * [out] \param value  Will be filled with the integer value for the parameter,
	 *       if it was found and it was a positive integer value.
	 * [in] \param prettyName  Should be the name of the parameter which will
	 *      be used when printing error messages. If this is set to NULL then
	 *      the paramName will be used instead (default is NULL).
	 * \return Zero if the object could be found and is valid. Otherwise an
	 *       error code, which is compatible with the HLT framework, is returned.
	 */
	int GetPositiveIntFromTMap(
			TMap* map, const char* paramName, Int_t& value,
			const char* pathToEntry = NULL, const char* prettyName = NULL
		) const;
	
	/**
	 * Tries to find a certain parameter in the TMap object and convert it to
	 * an floating point value.
	 * [in] \param map  The TMap object to search in.
	 * [in] \param paramName  The name of the parameter to search for.
	 * [out] \param value  Will be filled with the floating point value for the
	 *       parameter, if it was found and it was a floating point value.
	 * [in] \param prettyName  Should be the name of the parameter which will
	 *      be used when printing error messages. If this is set to NULL then
	 *      the paramName will be used instead (default is NULL).
	 * \return Zero if the object could be found and is valid. Otherwise an
	 *       error code, which is compatible with the HLT framework, is returned.
	 */
	int GetFloatFromTMap(
			TMap* map, const char* paramName, Double_t& value,
			const char* pathToEntry = NULL, const char* prettyName = NULL
		) const;
	
	/**
	 * Tries to find a certain parameter in the TMap object and convert it to
	 * an positive floating point value.
	 * [in] \param map  The TMap object to search in.
	 * [in] \param paramName  The name of the parameter to search for.
	 * [out] \param value  Will be filled with the floating point value for the
	 *       parameter, if it was found and it was a positive floating point value.
	 * [in] \param prettyName  Should be the name of the parameter which will
	 *      be used when printing error messages. If this is set to NULL then
	 *      the paramName will be used instead (default is NULL).
	 * \return Zero if the object could be found and is valid. Otherwise an
	 *       error code, which is compatible with the HLT framework, is returned.
	 */
	int GetPositiveFloatFromTMap(
			TMap* map, const char* paramName, Double_t& value,
			const char* pathToEntry = NULL, const char* prettyName = NULL
		) const;
	
	/**
	 * Fetches the reconstruction parameters object from the CDB for MUON.
	 * [out] \param params  This will be filled with the reconstruction
	 *      parameters object found if a successful status code is returned.
	 *      Otherwise it will be unchanged.
	 * \return Zero if the object could be found. Otherwise an error code,
	 *      which is compatible with the HLT framework, is returned.
	 */
	int LoadRecoParamsFromCDB(AliMUONRecoParam*& params) const;

	/**
	 * Dumps the data contained in a buffer to file as is.
	 */
	void DumpBuffer(
			const void* buffer, AliHLTUInt32_t size,
			const char* filename
		) const;

	/**
	 * Dumps the data block to file.
	 */
	void DumpBlock(
			const AliHLTComponentBlockData* block,
			const char* fileNamePrefix
		) const;
	
	/**
	 * Dumps the event information to files in the dump path given by the
	 * method DumpPath, which can be set by the command line argument -dumppath.
	 */
	void DumpEvent(
			const AliHLTComponentEventData& evtData,
			const AliHLTComponentBlockData* blocks,
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr,
			AliHLTUInt32_t& size,
			AliHLTComponentBlockDataList& outputBlocks
		) const;
	
	/**
	 * Dumps the event information to files in the dump path given by the
	 * method DumpPath, which can be set by the command line argument -dumppath.
	 */
	void DumpEvent(
			const AliHLTComponentEventData& evtData,
			AliHLTComponentTriggerData& trigData
		) const;

private:

	// Do not allow copying of this class.
	/// Not implemented.
	AliHLTMUONProcessor(const AliHLTMUONProcessor& /*obj*/);
	/// Not implemented.
	AliHLTMUONProcessor& operator = (const AliHLTMUONProcessor& /*obj*/);

	bool fWarnForUnexpecedBlock;  ///< Flag indicating if we should log a warning if we got a block of an unexpected type.
	bool fDelaySetup;  ///< Indicates if the component should delay loading and initialising from the CDB to the start of run event.
	bool fDumpDataOnError; ///< Flag indicating if we should dump data when an error occurs in the reconstruction class.
	const char* fDumpPath; ///< This is the path prefix to use to dump event data too when an error occurs.
	
	ClassDef(AliHLTMUONProcessor, 0)  // Abstract base class for dHLT specific components.
};

//______________________________________________________________________________

template <class BlockType>
bool AliHLTMUONProcessor::BlockStructureOk(
		const BlockType& block, const char* name, bool checkHeader
	) const
{
	/// Performs basic checks to see if the input data block structure is OK,
	/// that it is not corrupt, too short etc...
	
	if (not block.BufferSizeOk())
	{
		size_t headerSize = sizeof(typename BlockType::HeaderType);
		if (block.BufferSize() < headerSize)
		{
			HLTError("Received a %s data block with a size of %d bytes,"
				" which is smaller than the minimum valid header size of %d bytes."
				" The block must be corrupt.",
				name, block.BufferSize(), headerSize
			);
			return false;
		}
		
		size_t expectedWidth = sizeof(typename BlockType::ElementType);
		if (block.CommonBlockHeader().fRecordWidth != expectedWidth)
		{
			HLTError("Received a %s data block with a record"
				" width of %d bytes, but the expected value is %d bytes."
				" The block might be corrupt.",
				name,
				block.CommonBlockHeader().fRecordWidth,
				expectedWidth
			);
			return false;
		}
		
		HLTError("Received a %s data block with a size of %d bytes,"
			" but the block header claims the block should be %d bytes."
			" The block might be corrupt.",
			name, block.BufferSize(), block.BytesUsed()
		);
		return false;
	}
	
	if (checkHeader)
	{
		AliHLTMUONUtils::WhyNotValid reason;
		if (not AliHLTMUONUtils::HeaderOk(block.BlockHeader(), &reason))
		{
			HLTError("Received a %s data block which might be corrupt. %s",
				name, AliHLTMUONUtils::FailureReasonToMessage(reason)
			);
			return false;
		}
	}
	
	return true;
}

#endif // ALIHLTMUONPROCESSOR_H
