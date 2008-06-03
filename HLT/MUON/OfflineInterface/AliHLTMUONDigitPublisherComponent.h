#ifndef ALIHLTMUONDIGITPUBLISHERCOMPONENT_H
#define ALIHLTMUONDIGITPUBLISHERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id: AliHLTMUONDigitPublisherComponent.h 26179 2008-05-29 22:27:27Z aszostak $ */

///
/// @file   AliHLTMUONDigitPublisherComponent.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 May 2008
/// @brief  Declaration of a component to publish MUON data from digits.
///

#include "AliHLTOfflineDataSource.h"

#if __GNUC__ && __GNUC__ < 3
#define std
#endif

/**
 * @class AliHLTMUONDigitPublisherComponent
 * The component is used to convert simulated or reconstructed digits into DDL
 * raw data streams on the fly and publish them. This is useful for running
 * dHLT simulations on digit data where the raw data files are not available.
 * The component is also used during dHLT simulation under AliSimulation where
 * we do not have access to a rawReader.<br>
 *
 * Component ID: \b MUONDigitPublisher <br>
 * Library: \b libAliHLTMUON.so  <br>
 *
 * Manditory arguments:<br>
 * \li -ddl <number> <br>
 *       Specifies which DDL to generate raw data for. <number> should be in the
 *       range [1..22], where 1..20 are for tracker DDLs and 21 or 22 is for the
 *       dimuon trigger.<br>
 *
 * Optional arguments:<br>
 * \li -simdata <br>
 *       Indicates that the simulated digits should be used. (default option)<br>
 * \li -recdata <br>
 *       Indicates that the reconstructed digits tree should be used.<br>
 * \li -firstevent <number> <br>
 *      Indicates the first event number to fetch from AliRoot. The default is to
 *      start from zero and increment the event number after every GetEvent call.
 *      This mode causes the component to ignore the event number passed to it by
 *      the system and rather use an internal counter. This mode can be overriden
 *      with the -event_number_literal flag. <br>
 * \li -event_number_literal <br>
 *      This flag indicates to use the event numbers as literal indices into the
 *      AliRoot trees. This option will cause the component to ignore the -firstevent
 *      flag. <br>
 *
 * @ingroup alihlt_dimuon_component
 */
class AliHLTMUONDigitPublisherComponent : public AliHLTOfflineDataSource
{
public:
	AliHLTMUONDigitPublisherComponent();
	virtual ~AliHLTMUONDigitPublisherComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the component registration process.

	virtual const char* GetComponentID();
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
	virtual AliHLTComponent* Spawn();
	
protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component.
	
	virtual int DoInit(int argc, const char** argv);
	virtual int DoDeinit();
	
	virtual int GetEvent(
			const AliHLTComponentEventData& evtData,
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr,
			AliHLTUInt32_t& size,
			vector<AliHLTComponentBlockData>& outputBlocks
		);
	
	using AliHLTOfflineDataSource::GetEvent;
	
private:

	// Do not allow copying of this class.
	AliHLTMUONDigitPublisherComponent(const AliHLTMUONDigitPublisherComponent& /*obj*/);
	AliHLTMUONDigitPublisherComponent& operator = (const AliHLTMUONDigitPublisherComponent& /*obj*/);
	
	AliHLTInt32_t fDDL;  ///< DDL number in the range [0..21]. Set to -1 for invalid/unspecified value.
	
	Int_t fCurrentEventIndex;  ///< The current event index that is to be loaded.
	                           //  -1 indicates that we should rather use the event
	                           // numbers as given by the system.
	
	ClassDef(AliHLTMUONDigitPublisherComponent, 0)  // dHLT component for publishing DDL streams from digits on the fly.
};

#endif // ALIHLTMUONDIGITPUBLISHERCOMPONENT_H
