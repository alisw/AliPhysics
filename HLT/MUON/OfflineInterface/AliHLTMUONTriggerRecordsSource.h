#ifndef ALIHLTMUONTRIGGERRECORDSSOURCE_H
#define ALIHLTMUONTRIGGERRECORDSSOURCE_H
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
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

/**
 * @file   AliHLTMUONTriggerRecordsSource.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Class for generating trigger record data blocks from AliRoot data.
 */

#include "AliHLTOfflineDataSource.h"

class AliMUONMCDataInterface;
class AliMUONDataInterface;

/**
 * AliHLTMUONTriggerRecordsSource is a HLT-AliRoot data source object which generates
 * and serves AliHLTMUONTriggerRecordsBlockStruct type data blocks to the HLT system.
 * This is meant as a debugging utility which can optionally generate the data
 * blocks from simulated GEANT hits, simulated local trigger objects or MUON
 * offline reconstructed local trigger objects.
 *
 * Command line flags:
 *  -hitdata
 *      Specify this option to publish trigger records constructed from GEANT hits.
 *  -simdata
 *      Specify this option to publish trigger records constructed from simulated
 *      local trigger objects.
 *  -recdata
 *      Specify this option to publish trigger records constructed from offline
 *      reconstructed local trigger objects.
 *  -plane left|right|all
 *      Specifies if data from the left (x < 0), right (x >= 0) or the whole XY
 *      plane should be published.
 *  -firstevent <number>
 *      Indicates the first event number to fetch from AliRoot. The default is to
 *      start from zero and increment the event number after every GetEvent call.
 *      This mode causes the component to ignore the event number passed to it by
 *      the system and rather use an internal counter. This mode can be overriden
 *      with the -event_number_literal flag.
 *  -event_number_literal
 *      This flag indicates to use the event numbers as literal indices into the
 *      AliRoot trees. This option will cause the component to ignore the -firstevent
 *      flag.
 */
class AliHLTMUONTriggerRecordsSource : public AliHLTOfflineDataSource
{
public:

	AliHLTMUONTriggerRecordsSource();
	virtual ~AliHLTMUONTriggerRecordsSource();
	
	virtual int GetEvent(
			const AliHLTComponentEventData& evtData,
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr, 
			AliHLTUInt32_t& size,
			vector<AliHLTComponentBlockData>& outputBlocks
		);
	
	virtual const char* GetComponentID();

	virtual AliHLTComponentDataType GetOutputDataType();

	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

	virtual AliHLTComponent* Spawn();

protected:

	virtual int DoInit(int argc, const char** argv);
	virtual int DoDeinit();
	
private:

	// Prevent copying of these objects.
	AliHLTMUONTriggerRecordsSource(const AliHLTMUONTriggerRecordsSource& /*object*/);
	AliHLTMUONTriggerRecordsSource& operator = (const AliHLTMUONTriggerRecordsSource& /*object*/);
	
	AliMUONMCDataInterface* fMCDataInterface; //! access to MUON MC-related data
	AliMUONDataInterface* fDataInterface; //! access to MUON data
	bool fBuildFromHits;  //! Flag indicating if trigger records should be built from GEANT hits.

	enum SelectionType
	{
		kLeftPlane,  // everything from x < 0
		kRightPlane, // everything from x >= 0
		kWholePlane  // for all x
	};

	SelectionType fSelection; //! Indicates if we should publish from the left, right or whole XY plane.
	
	Int_t fCurrentEvent;  //! The current event index that is loaded.
	                      //  -1 indicates that we should rather use the event
	                      // numbers as given by the system.
	
	ClassDef(AliHLTMUONTriggerRecordsSource, 0); // dHLT data source for trigger record data blocks.
};

#endif // ALIHLTMUONTRIGGERRECORDSSOURCE_H
