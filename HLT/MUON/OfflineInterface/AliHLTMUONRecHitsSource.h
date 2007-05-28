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

class AliMUONSimData;
class AliMUONRecData;
class AliRunLoader;
class AliLoader;

/**
 * AliHLTMUONRecHitsSource is a HLT-AliRoot data source object which generates
 * and serves AliHLTMUONRecHitsBlockStruct type data blocks to the HLT system.
 * This is meant as a debugging utility which can optionally generate the data
 * blocks from simulate GEANT hits or MUON offline reconstructed hits.
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

	AliMUONSimData* fSimData; //! MUON module interface to simulated data.
	AliMUONRecData* fRecData; //! MUON module interface to reconstructed data.
	AliRunLoader* fRunLoader; //! A pointer to the AliRunLoader instance.
	AliLoader* fLoader; //! Pointer to the MUON loader instance.

	ClassDef(AliHLTMUONRecHitsSource, 0); // dHLT data source for reconstructed hit data blocks.
};

#endif // ALIHLTMUONRECHITSSOURCE_H
