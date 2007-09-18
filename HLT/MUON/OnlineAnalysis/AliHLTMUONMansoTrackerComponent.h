#ifndef AliHLTMUONMANSOTRACKERCOMPONENT_H
#define AliHLTMUONMANSOTRACKERCOMPONENT_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/**
 *  @file   AliHLTMUONMansoTrackerComponent.h
 *  @author Artur Szostak <artursz@iafrica.com>
 *  @date   
 *  @brief  Tracker component implementing the Manso algorithm for the dimuon HLT.
 */

#include "AliHLTProcessor.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONMansoTracker.h"

/**
 * @class AliHLTMUONMansoTrackerComponent
 * @brief Manso tracker component
 */
class AliHLTMUONMansoTrackerComponent : public AliHLTProcessor
    {
    public:
	AliHLTMUONMansoTrackerComponent();
	virtual ~AliHLTMUONMansoTrackerComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	const char* GetComponentID();
	void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
	AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
	AliHLTComponent* Spawn();
	
    protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	int DoInit( int argc, const char** argv );
	int DoDeinit();
	int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
	
    private:

	// The size of the output data produced, as a percentage of the input data's size.
	// Can be greater than 100 (%)
	unsigned fOutputPercentage;

	// For Tracking

	//TrackerCallback callback;
	AliHLTMUONMansoTracker* fTracker;

/* 	struct TriggerBlock */
/* 	{ */
/* 	  const AliHLTMUONCoreTriggerRecord *data; */
/* 	}; */

	AliHLTMUONCoreTriggerRecord *fTrigData;
	//std::vector<TriggerBlock> fTriggerBlocks;  // array of trigger record blocks
	
	ClassDef(AliHLTMUONMansoTrackerComponent, 0);

    };

#endif // AliHLTMUONMANSOTRACKERCOMPONENT_H
