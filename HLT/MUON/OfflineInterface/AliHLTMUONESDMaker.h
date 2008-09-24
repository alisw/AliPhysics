#ifndef ALIHLTMUONESDMAKER_H
#define ALIHLTMUONESDMAKER_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///
/// @file   AliHLTMUONESDMaker.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   30 June 2008
/// @brief  Component for converting dHLT raw data into AliESDEvent objects.
///

#include "AliHLTMUONProcessor.h"

/**
 * @class AliHLTMUONESDMaker
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
 * \li -warn_on_unexpected_block <br>
 *       If set, then warning messages are generated for any data block types that
 *       were not expected. (default is to generate only debug messages)<br>
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
	
	bool fWarnForUnexpecedBlock;  /// Flag indicating if we should log a warning if we got a block of an unexpected type.
	bool fMakeMinimalESD;  /// Flag to indicate if a minimal ESD object should be created.
	bool fAddCustomData;  /// Flag to turn on adding of all dHLT rootified objects to the ESD.

	ClassDef(AliHLTMUONESDMaker, 0); // Component for converting dHLT reconstructed data into the ESD format.
};

#endif // ALIHLTMUONESDMAKER_H
