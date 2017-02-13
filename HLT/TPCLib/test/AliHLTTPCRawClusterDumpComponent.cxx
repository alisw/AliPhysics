// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Mikolaj Krzewicki <mikolaj.krzewicki@cern.ch>         *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

//  @file   AliHLTTPCRawClusterDumpComponent.cxx
//  @author David Rohr <drohr@cern.ch>
//

#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCRawClusterDumpComponent.h"

ClassImp(AliHLTTPCRawClusterDumpComponent)

AliHLTTPCRawClusterDumpComponent::AliHLTTPCRawClusterDumpComponent() : AliHLTProcessor()
{
}

AliHLTTPCRawClusterDumpComponent::~AliHLTTPCRawClusterDumpComponent()
{
}

void AliHLTTPCRawClusterDumpComponent::GetInputDataTypes(AliHLTComponentDataTypeList &list)
{
	list.push_back(AliHLTTPCDefinitions::fgkRawClustersDataType | kAliHLTDataOriginTPC);
}

AliHLTComponentDataType AliHLTTPCRawClusterDumpComponent::GetOutputDataType()
{
	return kAliHLTDataTypeHistogram | kAliHLTDataOriginOut;
}

void AliHLTTPCRawClusterDumpComponent::GetOutputDataSize(unsigned long &constBase, double &inputMultiplier)
{
	constBase = 1000;
	inputMultiplier = 0.0;
}

int AliHLTTPCRawClusterDumpComponent::ProcessOption(TString option, TString value)
{
	int iResult = 0;

	if (0) {}
	else
	{
		HLTError("invalid option: %s", value.Data());
		return -EINVAL;
	}
	return iResult;
}

int AliHLTTPCRawClusterDumpComponent::DoInit(int argc, const char **argv)
{
	int iResult = 0;

	if (ProcessOptionString(GetComponentArgs()) < 0)
	{
		HLTFatal("wrong config string! %s", GetComponentArgs().c_str());
		return -EINVAL;
	}

	return iResult;
}

int AliHLTTPCRawClusterDumpComponent::DoDeinit()
{
	return 0;
}

int AliHLTTPCRawClusterDumpComponent::DoEvent(const AliHLTComponentEventData &evtData, const AliHLTComponentBlockData *blocks, AliHLTComponentTriggerData & /*trigData*/, AliHLTUInt8_t * /*outputPtr*/, AliHLTUInt32_t & /*size*/, AliHLTComponentBlockDataList & /*outputBlocks*/)
{
	int iResult = 0;

	if (!IsDataEvent()) {
		return iResult;
	}
	
	int nBlocks = evtData.fBlockCnt;
	static int nEvent = 0;
	char filename[1024];
	sprintf(filename, "event_tpc_dump.%d", nEvent++);
	FILE* fp = fopen(filename, "w+b");
	if (fp == NULL) return(-EINVAL);
	
	int nTotalClusters = 0;
	int nClusterBlocks = 0;

	for (int ndx = 0; ndx < nBlocks; ndx++)
	{
		const AliHLTComponentBlockData *iter = blocks + ndx;

		if (iter->fDataType == (AliHLTTPCDefinitions::fgkRawClustersDataType | kAliHLTDataOriginTPC))
		{
			int slice = AliHLTTPCDefinitions::GetMinSliceNr(iter->fSpecification);
			int patch = AliHLTTPCDefinitions::GetMinPatchNr(iter->fSpecification);

			AliHLTTPCRawClusterData* clusters = (AliHLTTPCRawClusterData *) (iter->fPtr);
			RawClusterDumpHeader header;
			header.fSector = slice;
			header.fPatch = patch;
			header.fNClusters = clusters->fCount;
			fwrite(&header, sizeof(header), 1, fp);
			fwrite(clusters->fClusters, sizeof(clusters->fClusters[0]), header.fNClusters, fp);
			nTotalClusters += header.fNClusters;
			nClusterBlocks++;
		}

	}
	
	fclose(fp);
	
	HLTImportant("Dumped %d clusters to %s (%d blocks)", nTotalClusters, filename, nClusterBlocks);

	return iResult;
}
