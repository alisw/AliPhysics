#ifndef ALIHLTMUONESDMAKER_H
#define ALIHLTMUONESDMAKER_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

// $Id: $

///
/// @file   AliHLTMUONESDMaker.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   30 June 2008
/// @brief  Component for converting dHLT raw data into AliESDEvent objects.
///

#include "AliHLTMUONProcessor.h"
#include <vector>

extern "C" class AliHLTMUONTriggerRecordStruct;
class AliESDMuonTrack;

/**
 * @class AliHLTMUONESDMaker
 * \brief Component for converting dHLT results into ESD format.
 *
 * The component is used to convert dHLT reconstructed data into AliESDEvent
 * objects which can be stored in ROOT files during offline reconstruction.
 * Only the dHLT track and trigger record data is converted, then filled in the ESD.
 * These should then be merged together with ESDs from all the other parts of
 * HLT (eg. TPC HLT).<br>
 * This component can also be run online to have ESDs directly in the raw
 * HLT output data stream.<br>
 *
 * Component ID: \b MUONESDMaker <br>
 * Library: \b libAliHLTMUON.so  <br>
 *
 * Optional arguments:<br>
 * \li -make_minimal_esd <br>
 *       Indicates that AliESDEvent objects should be created with only the TClonesArray
 *       for the muon tracks created. (default is to generate all standard ESD objects)<br>
 * \li -add_rootified_objects <br>
 *       If specified then the any rootified dHLT event data that is found is added to the
 *       ESD list of objects as a custom data object.<br>
 * \li -makeclonesarray <br>
 *      This option will cause the component to generate a TClonesArray of MUON ESD tracks
 *      and send it as a kAliHLTDataTypeTObject data block type. <br>
 * \li -makeonlyclonesarray <br>
 *      Same as the -makeclonesarray option, however the data block with the AliESDEvent
 *      object is not generated at all. <br>
 * \li -warn_on_unexpected_block <br>
 *       If set, then warning messages are generated for any data block types that
 *       were not expected. (default is to generate only debug messages)<br>
 * \li -dumponerror <br>
 *      This flag will cause the component to dump the data blocks it received if
 *      an error occurs during the processing of an event. <br>
 * \li -dumppath <i>path</i> <br>
 *      Allows one to specify the path in which to dump the received data blocks
 *      if an error occurs. <br>
 *
 * @ingroup alihlt_dimuon_component
 */
class AliHLTMUONESDMaker : public AliHLTMUONProcessor
{
public:

	AliHLTMUONESDMaker();
	virtual ~AliHLTMUONESDMaker();
	
	virtual const char* GetComponentID();

	virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& list);
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

	virtual AliHLTComponent* Spawn();

protected:

	virtual int DoInit(int argc, const char** argv);
	virtual int DoDeinit();
	virtual int DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
	virtual bool IgnoreArgument(const char* arg) const;
	using AliHLTProcessor::DoEvent;
	
private:

	// Prevent copying of these objects.
	AliHLTMUONESDMaker(const AliHLTMUONESDMaker& /*object*/);
	AliHLTMUONESDMaker& operator = (const AliHLTMUONESDMaker& /*object*/);
	
	typedef std::vector<const AliHLTMUONTriggerRecordStruct*> AliTriggerRecordList;
	
	/**
	 * Finds the trigger record with the specified ID in the list of trigger records
	 * and then fills the ESD muon track structure with the information.
	 * \param [in] triggerRecords  The list of trigger records to search in.
	 * \param [in] trigRecId  The trigger record ID to seach for.
	 * \param [in] trackId  The track ID of the track structure where the trigger
	 *                     record ID comes from.
	 * \param [out] muTrack  The track structure to fill.
	 * \param [in,out] nHits  The number of hits added. Will increment this value
	 *                        for every new hit added.
	 */
	void FillTriggerInfo(
			const AliTriggerRecordList& triggerRecords,
			AliHLTInt32_t trigRecId, AliHLTInt32_t trackId,
			AliESDMuonTrack& muTrack, Int_t& nHits
		);
	
	bool fWarnForUnexpecedBlock;  /// Flag indicating if we should log a warning if we got a block of an unexpected type.
	bool fMakeMinimalESD;  /// Flag to indicate if a minimal ESD object should be created.
	bool fAddCustomData;  /// Flag to turn on adding of all dHLT rootified objects to the ESD.
	bool fMakeClonesArray;  /// Flag indicating if a data block of TClonesArray with AliESDMuonTrack objects should be generated.
	bool fMakeESDDataBlock;  /// Flag indicating if the ESD data block should generated.
	AliHLTUInt32_t fClusterIndex;  /// Running counter for the unique cluster index number.

	ClassDef(AliHLTMUONESDMaker, 0); // Component for converting dHLT reconstructed data into the ESD format.
};

#endif // ALIHLTMUONESDMAKER_H
