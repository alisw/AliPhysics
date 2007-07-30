#ifndef ALIHLTMUONRECHITSSOURCE_H
#define ALIHLTMUONRECHITSSOURCE_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/**
 * @file   AliHLTMUONRecHitsSource.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Class for generating reconstructed hits data blocks from AliRoot data.
 */

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
 */
class AliHLTMUONRecHitsSource : public AliHLTOfflineDataSource
{
public:

	AliHLTMUONRecHitsSource();
	virtual ~AliHLTMUONRecHitsSource();
	
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

	ClassDef(AliHLTMUONRecHitsSource, 0); // dHLT data source for reconstructed hit data blocks.
};

#endif // ALIHLTMUONRECHITSSOURCE_H
