#ifndef AliHLTMUONMANSOTRACKERFSMCOMPONENT_H
#define AliHLTMUONMANSOTRACKERFSMCOMPONENT_H
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///
///  @file   AliHLTMUONMansoTrackerFSMComponent.h
///  @author Artur Szostak <artursz@iafrica.com>,
///          Indranil Das <indra.das@saha.ac.in>
///  @date   
///  @brief  Tracker component for the dimuon HLT using the Manso algorithm
///          implemented as a finite state machine.
///

#include "AliHLTProcessor.h"
#include "AliHLTMUONDataTypes.h"
#include "AliHLTMUONMansoTrackerFSMCallback.h"
#include <vector>

#if __GNUC__ && __GNUC__ < 3
#define std
#endif

class AliHLTMUONMansoTrackerFSM;
//class AliHLTMUONMansoTracksBlockWriter;
extern "C" {
struct AliHLTMUONRecHitStruct;
} // extern C


/**
 * @class AliHLTMUONMansoTrackerFSMComponent
 * @brief Dimuon HLT tracker component using the Manso tracking algorithm
 *        implemented as a finite state machine.
 */
class AliHLTMUONMansoTrackerFSMComponent
	: public AliHLTProcessor, public AliHLTMUONMansoTrackerFSMCallback
{
public:
	AliHLTMUONMansoTrackerFSMComponent();
	virtual ~AliHLTMUONMansoTrackerFSMComponent();

	// Public functions to implement the AliHLTProcessor interface.
	// These functions are required for the registration process.
	virtual const char* GetComponentID();
	virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
	virtual AliHLTComponent* Spawn();
	
	// Inherited from AliHLTMUONMansoTrackerFSMCallback:
	virtual void RequestClusters(
			AliHLTMUONMansoTrackerFSM* tracker,
			AliHLTFloat32_t left, AliHLTFloat32_t right,
			AliHLTFloat32_t bottom, AliHLTFloat32_t top,
			AliHLTMUONChamberName chamber, const void* tag
		);
	virtual void EndOfClusterRequests(AliHLTMUONMansoTrackerFSM* tracker);
	virtual void FoundTrack(AliHLTMUONMansoTrackerFSM* tracker);
	virtual void NoTrackFound(AliHLTMUONMansoTrackerFSM* tracker);

protected:

	// Protected functions to implement the AliHLTProcessor interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 
	virtual int DoInit(int argc, const char** argv);
	virtual int DoDeinit();
	virtual int DoEvent(
			const AliHLTComponentEventData& evtData,
			const AliHLTComponentBlockData* blocks,
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr,
			AliHLTUInt32_t& size,
			std::vector<AliHLTComponentBlockData>& outputBlocks
		);
	
	using AliHLTProcessor::DoEvent;

private:

	// Do not allow copying of this class.
	AliHLTMUONMansoTrackerFSMComponent(const AliHLTMUONMansoTrackerFSMComponent& /*obj*/);
	AliHLTMUONMansoTrackerFSMComponent& operator = (const AliHLTMUONMansoTrackerFSMComponent& /*obj*/);

	void Reset();
	
	void AddRecHits(
			AliHLTUInt32_t specification,	
			const AliHLTMUONRecHitStruct* recHits,
			AliHLTUInt32_t count
		);

	AliHLTMUONMansoTrackerFSM* fTracker;  //! Tracker to do the actual work.
	
	AliHLTUInt32_t fTrackCount;  //! Number of tracks currently found.
	/*AliHLTMUONMansoTracksBlockWriter*/void* fBlock;  //! The current data block we are writing.
	
	struct RecHitBlockInfo
	{
		AliHLTUInt32_t fCount;  // Number of elements in fData.
		const AliHLTMUONRecHitStruct* fData; // Pointer to the array of rec hits.
	};
	
	std::vector<RecHitBlockInfo> fRecHitBlock[4];  //! Arrays of rec hit block data.

	bool fWarnForUnexpecedBlock;  // Flag indicating if we should log a warning if we got a block of an unexpected type.
	
	ClassDef(AliHLTMUONMansoTrackerFSMComponent, 0);  // Manso tracker component implemented as a finite state machine (FSM).
};

#endif // AliHLTMUONMANSOTRACKERFSMCOMPONENT_H
