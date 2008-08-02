#ifndef ALIHLTMUONRECHITSSOURCE_H
#define ALIHLTMUONRECHITSSOURCE_H
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

///
/// @file   AliHLTMUONRecHitsSource.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   
/// @brief  Class for generating reconstructed hits data blocks from AliRoot data.
///

#include "AliHLTOfflineDataSource.h"

class AliMUONMCDataInterface;
class AliMUONDataInterface;

/**
 * AliHLTMUONRecHitsSource is a HLT-AliRoot data source object which generates
 * and serves AliHLTMUONRecHitsBlockStruct type data blocks to the HLT system.
 * This is meant as a debugging utility which can optionally generate the data
 * blocks from simulated GEANT hits or MUON offline reconstructed hits.
 *
 * Command line flags:
 *  -simdata
 *      Specify this option to publish GEANT hits.
 *  -recdata
 *      Specify this option to publish offline reconstructed raw clusters.
 *  -plane left|right|all
 *      Specifies if data from the left (x < 0), right (x >= 0) or the whole XY
 *      plane should be published.
 *  -chamber <number>|<number>-<number>[,<number>|<number>-<number>]...
 *      Selects the chambers from which to publish data. Valid chamber numbers
 *      in the range [1..10]. The string after '-chamber' is a comma separated
 *      list of numbers or ranges. Some examples of strings:
 *      1  1-2  1,2,3  1,3-5,7 etc...
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
class AliHLTMUONRecHitsSource : public AliHLTOfflineDataSource
{
public:

	AliHLTMUONRecHitsSource();
	virtual ~AliHLTMUONRecHitsSource();
	
	virtual const char* GetComponentID();

	virtual AliHLTComponentDataType GetOutputDataType();

	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

	virtual AliHLTComponent* Spawn();

protected:

	virtual int DoInit(int argc, const char** argv);
	virtual int DoDeinit();
	
	virtual int GetEvent(
			const AliHLTComponentEventData& evtData,
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr, 
			AliHLTUInt32_t& size,
			AliHLTComponentBlockDataList& outputBlocks
		);
	
	using AliHLTOfflineDataSource::GetEvent;
	
private:

	// Prevent copying of these objects.
	AliHLTMUONRecHitsSource(const AliHLTMUONRecHitsSource& /*object*/);
	AliHLTMUONRecHitsSource& operator = (const AliHLTMUONRecHitsSource& /*object*/);
	
	AliMUONMCDataInterface* fMCDataInterface; // access to MUON MC-related data
	AliMUONDataInterface* fDataInterface; // access to MUON data

	enum SelectionType
	{
		kLeftPlane,  // everything from x < 0
		kRightPlane, // everything from x >= 0
		kWholePlane  // for all x
	};
	
	/**
	 * Parses a string with the following format:
	 *   <number>|<number>-<number>[,<number>|<number>-<number>]...
	 * For example: 1  1,2,3  1-2   1,2-4,5  etc...
	 * Flags in the fServeChamber will be set to 'true' for all appropriate
	 * values parsed.
	 * @param str  The string to parse.
	 * @return  Zero on success and EINVAL if there is a parse error.
	 */
	int ParseChamberString(const char* str);

	SelectionType fSelection; //! Indicates if we should publish from the left, right or whole XY plane.
	bool fServeChamber[10]; //! Flag to indicate if hits from a given chamber should be published.
	
	Int_t fCurrentEventIndex;  //! The current event index that is loaded.
	                           //  -1 indicates that we should rather use the event
	                           // numbers as given by the system.

	ClassDef(AliHLTMUONRecHitsSource, 0); // dHLT data source for reconstructed hit data blocks.
};

#endif // ALIHLTMUONRECHITSSOURCE_H
