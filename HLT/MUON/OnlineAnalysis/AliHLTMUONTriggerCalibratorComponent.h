#ifndef ALIHLTMUONTRIGGERCALIBRATIONCOMPONENT_H
#define ALIHLTMUONTRIGGERCALIBRATIONCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// @file   AliHLTMUONTriggerCalibratorComponent.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   
/// @brief  The calibration component for the muon trigger.
///

#include "AliHLTCalibrationProcessor.h"

/**
 * @class AliHLTMUONTriggerCalibratorComponent
 *
 * This class addapts the code found in MUON/MUONTRGda.cxx
 */
class AliHLTMUONTriggerCalibratorComponent : public AliHLTCalibrationProcessor
{
public:

	AliHLTMUONTriggerCalibratorComponent();
	virtual ~AliHLTMUONTriggerCalibratorComponent();
	
	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process
	
	virtual const char* GetComponentID();
	virtual void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
	virtual AliHLTComponent* Spawn();
	
protected:
	
	using AliHLTCalibrationProcessor::ProcessCalibration;
	using AliHLTCalibrationProcessor::ShipDataToFXS;
	
	// Protected functions to implement AliHLTCalibrationProcessor's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component.
	
	virtual Int_t InitCalibration();
	virtual Int_t DeinitCalibration();
	virtual Int_t ScanArgument(int argc, const char** argv);
	virtual Int_t ProcessCalibration(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
	virtual Int_t ShipDataToFXS(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
	
private:

	// Do not allow copying of this class.
	AliHLTMUONTriggerCalibratorComponent(const AliHLTMUONTriggerCalibratorComponent& /*obj*/);
	AliHLTMUONTriggerCalibratorComponent& operator = (const AliHLTMUONTriggerCalibratorComponent& /*obj*/);
	
	Int_t fSkipEvents;  //! Number of events to skip initially.
	Int_t fMaxEvents;   //! Maximum number of events to process.
	
	ClassDef(AliHLTMUONTriggerCalibratorComponent, 0); // The online calibration processing component for the muon trigger.

};

#endif // ALIHLTMUONTRIGGERCALIBRATIONCOMPONENT_H
