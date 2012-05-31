#ifndef ALIHLTMUONROOTIFIERCOMPONENT_H
#define ALIHLTMUONROOTIFIERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

// $Id$

///
/// @file   AliHLTMUONRootifierComponent.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 Sep 2007
/// @brief  Component for converting dHLT raw data into ROOT objects.
///

#include "AliHLTMUONProcessor.h"

class AliHLTMUONEvent;
class AliHLTMUONMansoTrack;
class AliHLTMUONTrack;
extern "C" struct AliHLTMUONMansoTrackStruct;

/**
 * \class AliHLTMUONRootifierComponent
 * \brief Converts dHLT raw data blocks into ROOT objects.
 *
 * This component class is used to convert all internal raw dHLT data blocks into
 * ROOT object that can be stored in '.root' files in a platform independant manner.
 * This can also make some of the analysis easier because the dHLT internal data
 * will be available in TObjects.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b MUONRootifier <br>
 * Library: \b libAliHLTMUON.so <br>
 * Input Data Types:  kAliHLTAnyDataType = "*******:***" <br>
 * Output Data Types: AliHLTMUONConstants::RootifiedEventDataType() = "ROOTEVNT:MUON" <br>
 *
 * <h2>Mandatory arguments:</h2>
 * None.
 *
 * <h2>Optional arguments:</h2>
 * \li -warn_on_unexpected_block <br>
 *      This will cause the component to generate warnings when it receives data block
 *      types it does not know how to handle. Without this option the component only
 *      generates debug messages when they are compiled in. <br>
 * \li -dumponerror <br>
 *      This flag will cause the component to dump the data blocks it received if
 *      an error occurs during the processing of an event. <br>
 * \li -dumppath <i>path</i> <br>
 *      Allows one to specify the path in which to dump the received data blocks
 *      if an error occurs. <br>
 *
 * <h2>Standard configuration:</h2>
 * There is no special configuration for this component.
 *
 * <h2>Default CDB entries:</h2>
 * None.
 *
 * <h2>Performance:</h2>
 * A few milliseconds per event.
 *
 * <h2>Memory consumption:</h2>
 * A few MBytes.
 *
 * <h2>Output size:</h2>
 * A few kBytes.
 *
 * @ingroup alihlt_dimuon_component
 */
class AliHLTMUONRootifierComponent : public AliHLTMUONProcessor
{
public:

	AliHLTMUONRootifierComponent();
	virtual ~AliHLTMUONRootifierComponent();
	
	virtual const char* GetComponentID();

	virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
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
	AliHLTMUONRootifierComponent(const AliHLTMUONRootifierComponent& /*object*/);
	AliHLTMUONRootifierComponent& operator = (const AliHLTMUONRootifierComponent& /*object*/);
	
	/**
	 * This method creates a AliHLTMUONMansoTrack object from the track structure
	 * and adds it to the dHLT event object.
	 * \param event  The dHLT event object.
	 * \param track  The track structure to convert and add to the event.
	 */
	AliHLTMUONMansoTrack* AddTrack(AliHLTMUONEvent& event, const AliHLTMUONMansoTrackStruct& track);
	
	/**
	 * This method creates a AliHLTMUONTrack object from the given track structure
	 * and adds it to the dHLT event object.
	 * \param event  The dHLT event object.
	 * \param track  The track structure to convert and add to the event.
	 */
	AliHLTMUONTrack* AddTrack(AliHLTMUONEvent& event, const AliHLTMUONTrackStruct& track);
	
	bool fWarnForUnexpecedBlock;  /// Flag indicating if we should log a warning if we got a block of an unexpected type.

	ClassDef(AliHLTMUONRootifierComponent, 0); // Converter component of dHLT raw data.
};

#endif // ALIHLTMUONROOTIFIERCOMPONENT_H
