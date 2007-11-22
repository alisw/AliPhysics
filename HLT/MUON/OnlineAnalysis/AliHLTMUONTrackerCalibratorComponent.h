#ifndef ALIHLTMUONTRACKERCALIBRATIONCOMPONENT_H
#define ALIHLTMUONTRACKERCALIBRATIONCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// @file   AliHLTMUONTrackerCalibratorComponent.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   
/// @brief  The calibration component for the muon tracking stations.
///

#include "AliHLTCalibrationProcessor.h"

class TFitter;

/**
 * @class AliHLTMUONTrackerCalibratorComponent
 *
 * This class addapts the code found in MUON/MUONTRKda.cxx
 */
class AliHLTMUONTrackerCalibratorComponent : public AliHLTCalibrationProcessor
{
public:

	AliHLTMUONTrackerCalibratorComponent();
	virtual ~AliHLTMUONTrackerCalibratorComponent();
	
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
	AliHLTMUONTrackerCalibratorComponent(const AliHLTMUONTrackerCalibratorComponent& /*obj*/);
	AliHLTMUONTrackerCalibratorComponent& operator = (const AliHLTMUONTrackerCalibratorComponent& /*obj*/);
	
	TFitter* fFitter;  //! Fitter used for gain calibration.
	Int_t fSkipEvents;  //! Number of events to skip initially.
	Int_t fMaxEvents;   //! Maximum number of events to process.
	
	TString fFlatOutputFile;   //! Flat file name.
	TString fCrocusOutputFile; //! Crocus command file name.
	TString fCrocusConfigFile; //! Crocus config file name.
	
	Int_t fInjCharge;  //! Injection charge.
	Double_t fNoSigma; //! Number of sigmas.
	Int_t fThreshold;  //! Threshold.
	
	ClassDef(AliHLTMUONTrackerCalibratorComponent, 0); // The online calibration processing component for the muon trigger.

};

#endif // ALIHLTMUONTRACKERCALIBRATIONCOMPONENT_H
