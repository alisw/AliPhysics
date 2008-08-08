#ifndef AliHLTMUONRAWDATAHISTOCOMPONENT_H
#define AliHLTMUONRAWDATAHISTOCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///
///  @file   AliHLTMUONRawDataHistoComponent.h
///  @author Artur Szostak <artursz@iafrica.com>
///  @date   30 April 2008
///  @brief  Declaration of a component to generate basic monitoring histograms of raw data.
///

#include "AliHLTMUONProcessor.h"
#include "AliHLTMUONDataTypes.h"

class TH1D;

/**
 * @class AliHLTMUONRawDataHistoComponent
 * @brief Dimuon HLT component for generating basic monitoring histograms for raw data.
 */
class AliHLTMUONRawDataHistoComponent : public AliHLTMUONProcessor
{
public:
	AliHLTMUONRawDataHistoComponent();
	virtual ~AliHLTMUONRawDataHistoComponent();

	// Public functions to implement the AliHLTProcessor interface.
	// These functions are required for the registration process.
	virtual const char* GetComponentID();
	virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
	virtual AliHLTComponent* Spawn();

protected:

	// Protected functions to implement the AliHLTProcessor interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component.
	virtual int DoInit(int argc, const char** argv);
	virtual bool IgnoreArgument(const char* arg) const;
	virtual int DoDeinit();
	virtual int DoEvent(
			const AliHLTComponentEventData& evtData,
			AliHLTComponentTriggerData& trigData
		);
	
	using AliHLTProcessor::DoEvent;

private:

	// Do not allow copying of this class.
	AliHLTMUONRawDataHistoComponent(const AliHLTMUONRawDataHistoComponent& /*obj*/);
	AliHLTMUONRawDataHistoComponent& operator = (const AliHLTMUONRawDataHistoComponent& /*obj*/);

	void ProcessTrackerDDL(const AliHLTComponentBlockData* block);
	void ProcessTriggerDDL(const AliHLTComponentBlockData* block);
	void FreeObjects();
	
	double fLastPublishTime;  /// Timestamp for the last time we published data (seconds).
	double fCurrentEventTime;  /// Timestamp for the current event being processed (seconds).
	double fPublishDelay;  /// Delay in second to wait between publishing data.
	TH1D* fErrorHist[22]; /// Histograms for error codes per DDL.
	TH1D* fManuHist[20]; /// Histograms for MANU distributions per DDL.
	TH1D* fSignalHist[20]; /// Histograms for signal distributions per DDL.
	
	ClassDef(AliHLTMUONRawDataHistoComponent, 0);  // Trigger decision component for the dimuon HLT.
};

#endif // AliHLTMUONRAWDATAHISTOCOMPONENT_H
