#ifndef ALIHLTTPCRAWCLUSTERDUMPCOMPONENT_H
#define ALIHLTTPCRAWCLUSTERDUMPCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#include "AliHLTComponentBenchmark.h"
#include "AliHLTProcessor.h"
#include "AliOptionParser.h"

class AliHLTTPCRawClusterDumpComponent : public AliHLTProcessor, public AliOptionParser
{
public:
	/** standard constructor */
	AliHLTTPCRawClusterDumpComponent();
	/** destructor */
	virtual ~AliHLTTPCRawClusterDumpComponent();

	struct RawClusterDumpHeader
	{
		int fSector;
		int fPatch;
		int fNClusters;
	};

	const char* GetComponentID() {return "RawClusterDump";};
	void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	AliHLTComponentDataType GetOutputDataType();
	void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
	AliHLTComponent* Spawn() {return new AliHLTTPCRawClusterDumpComponent;}

protected:
	    // interface methods of base class
	int DoInit(int argc, const char **argv);
	int DoDeinit();
	int DoEvent(const AliHLTComponentEventData &evtData,
	            const AliHLTComponentBlockData *blocks,
	            AliHLTComponentTriggerData &trigData,
	            AliHLTUInt8_t *outputPtr,
	            AliHLTUInt32_t &size,
	            AliHLTComponentBlockDataList &outputBlocks);

	using AliHLTProcessor::DoEvent;
	int ProcessOption(TString option, TString value);

private:
	/** copy constructor prohibited */
	AliHLTTPCRawClusterDumpComponent(const AliHLTTPCRawClusterDumpComponent &);
	/** assignment operator prohibited */
	AliHLTTPCRawClusterDumpComponent &operator=(const AliHLTTPCRawClusterDumpComponent &);

protected:
	ClassDef(AliHLTTPCRawClusterDumpComponent, 0)
};
#endif
