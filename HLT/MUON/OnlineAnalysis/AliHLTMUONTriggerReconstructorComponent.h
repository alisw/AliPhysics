#ifndef AliHLTMUONTRIGGERRECONSTRUCTORCOMPONENT_H
#define AliHLTMUONTRIGGERRECONSTRUCTORCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

// $Id$

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
 *
 * The trigger reconstructor component is used to decode the raw data coming
 * from the trigger chambers and electronics of the muon spectrometer.
 * The local trigger decisions are converted into trigger records which is a
 * usable format by the tracking stage.
 * No full cluster finding is performed, rather just the fired strip information
 * as received from the trigger electronics is converted into global coordinates
 * to be used by the tracker as track seeds.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b MUONTriggerReconstructor <br>
 * Library: \b libAliHLTMUON.so   <br>
 * Input Data Types: AliHLTMUONConstants::DDLRawDataType() = "DDL_RAW :MUON" <br>
 * Output Data Types: \li AliHLTMUONConstants::TriggerRecordsBlockDataType() = "TRIGRECS:MUON"
 *                    \li AliHLTMUONConstants::TrigRecsDebugBlockDataType() = "TRIGRDBG:MUON" <br>
 *
 * <h2>Mandatory arguments:</h2>
 * \li -ddl <i>number</i> <br>
 *      This indicates the DDL from which the component is expect to receive data
 *      and for which it should load the appropriate electronics mapping lookup
 *      table.
 *      The <i>number</i> should be in the range [21..22], following local dimuon
 *      spectrometer DDL numbering. If either the -ddlid, -lut or -delaysetup
 *      arguments are used, then -ddl becomes optional. <br>
 * \li -ddlid <i>number</i> <br>
 *      This indicates the DDL by equipment ID, from which the component is expect
 *      to receive data and for which it should load the appropriate electronics
 *      mapping lookup table.
 *      The <i>number</i> should be in the range [2816..2817].
 *      If either the -ddl, -lut or -delaysetup arguments are used, then -ddlid
 *      becomes optional. <br>
 * \li -delaysetup <br>
 *      Specifying this option causes the component to initialise the lookup table
 *      and magnetic field parameters from CDB only after receiving the first event
 *      to process in DoEvent.
 *      If -ddl or -ddlid were not used, then the DDL number will be taken from
 *      the first block's specification during runtime from the first
 *      event (i.e. Start-of-Run event).
 *      Using the -lut, -zmiddle or -bfieldintegral arguments will override loading
 *      from CDB for a delayed setup. <br>
 *
 * <h2>Optional arguments:</h2>
 * \li -lut <i>filename</i> <br>
 *      A pree-built lookup table for the electronics mapping and calibration
 *      information can be loaded with this argument. The file should have been
 *      generated with the GenerateLookupTable method. The location of the file
 *      is given by the parameter <i>filename</i> <br>
 * \li -cdb <br>
 *      Indicates that the component should load from CDB. This option is implicit
 *      if the -cdbpath is given. It will also override the -lut option.<br>
 * \li -cdbpath <i>path</i> <br>
 *      Specifies the CDB path to use, given by <i>path</i>. This option will override
 *      the CDB path automatically set by the HLT framework. <br>
 * \li -run <i>number</i> <br>
 *      Specifies the run number to use, given by <i>number</i>. This option will
 *      override the current run number automatically set by the HLT framework. <br>
 * \li -zmiddle <i>position</i> <br>
 *      This indicates the Z coordinate position of the middle of the dipole magnetic
 *      field. <i>position</i> is a floating point value in centimeters. Specifying
 *      this argument on the will override the value loaded from CDB. <br>
 * \li -bfieldintegral <i>field</i> <br>
 *      This indicates the magnetic field integral for the dipole magnetic field.
 *      <i>field</i> must be a floating point value in Tesla meters (T.m).
 *      The sign of the value will indicate the polarity setting of the dipole magnet.
 *      Specifying this argument on the will override the value loaded from CDB. <br>
 * \li -warn_on_unexpected_block <br>
 *      This will cause the component to generate warnings when it receives data block
 *      types it does not know how to handle. Without this option the component only
 *      generates debug messages when they are compiled in. <br>
 * \li -suppress_partial_triggers <br>
 *      This option forces all trigger records that have less than 3 hits in them
 *      to be removed from the output. This is the default setting. <br>
 * \li -generate_partial_triggers <br>
 *      With this option all trigger records, even partial ones with just one or two
 *      hits is written to the output. <br>
 * \li -stop_on_buffer_overflow <br>
 *      If this option is specified then the component will stop processing and generate
 *      an error code in the DoEvent method as soon as the output buffer has been filled.
 *      Otherwise the component normally just keeps processing but some data might be lost
 *      due to full buffers. <br>
 * \li -tryrecover <br>
 *      This is a special option to the raw data decoder to turn on logic which will
 *      try and recover from corrupt raw DDL data. This is off by default. <br>
 * \li -dont_use_crateid <br>
 *      This option indicates that the crate ID values found in the regional structures
 *      in the raw DDL data should not be used to identify the channels in the offline
 *      mapping. Rather the position of the raw data structure instead. <br>
 * \li -dont_use_localid <br>
 *      This option indicates that the local structure ID values found in the raw DDL
 *      data should not be used to identify the channels in the offline mapping, but
 *      rather the position of the local structure in the DDL should be used instead. <br>
 * \li -dumponerror <br>
 *      This flag will cause the component to dump the data blocks it received if
 *      an error occurs during the processing of an event. <br>
 * \li -dumppath <i>path</i> <br>
 *      Allows one to specify the path in which to dump the received data blocks
 *      if an error occurs. <br>
 * \li -makedebuginfo <br>
 *      If specified then the trigger record debug informaiton data blocks are generated. <br>
 * \li -dontprintwrongeventerror <br>
 *      If specified the error message about an incorrect event type found in the DDL DARC
 *      header is not generated or logged. <br>
 *
 * <h2>Standard configuration:</h2>
 * The configuration is taken from the CDB by default. It can be overridden with
 * the command line arguments.
 *
 * <h2>Default CDB entries:</h2>
 * HLT/ConfigHLT/HLTGlobalTrigger - Contains the global trigger menu.
 *
 * <h2>Performance:</h2>
 * This is a linear function of the number of input triggers (AliHLTTrigger) that
 * need to be processed.
 * For a modest trigger menu configurations the processing time per event should
 * be on the order of a few milliseconds.
 *
 * <h2>Memory consumption:</h2>
 * This is a linear function of the input data size, but only a fraction. Thus the
 * memory usage is minimal. It should be under 1 MBytes.
 *
 * <h2>Output size:</h2>
 * This will depend linearly on the number of tracks found. But for nominal multiplicities
 * this should be less than 16 kBytes.
 *
 * @ingroup alihlt_dimuon_component
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
	virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& list);
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

	/**
	 * Read in the lookup table from file.
	 * \param lutpath  The file to read the lookup table from.
	 */
	int ReadLookUpTable(const char* lutpath);
	
	/**
	 * Loads the lookup table containing channel and geometrical position
	 * information about trigger strips from CDB.
	 *
	 * \note To override the default CDB path and/or run number the
	 * SetCDBPathAndRunNo(cdbPath, run) method should be called before this
	 * method.
	 *
	 * \return 0 on success and non zero codes for errors.
	 */
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
