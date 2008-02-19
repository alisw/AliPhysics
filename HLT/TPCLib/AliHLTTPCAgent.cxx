// @(#) $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCAgent.cxx
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTTPC library
*/

#include "AliHLTTPCAgent.h"
#include "AliHLTConfiguration.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTOUT.h"

/** global instance for agent registration */
AliHLTTPCAgent gAliHLTTPCAgent;

// component headers
#include "AliHLTTPCDigitDumpComponent.h"
#include "AliHLTTPCEsdWriterComponent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCAgent)

AliHLTTPCAgent::AliHLTTPCAgent()
  :
  AliHLTModuleAgent("TPC"),
  fRawDataHandler(NULL),
  fNofRawDataHandler(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCAgent::~AliHLTTPCAgent()
{
  // see header file for class documentation
}

int AliHLTTPCAgent::CreateConfigurations(AliHLTConfigurationHandler* handler,
					 AliRawReader* /*rawReader*/,
					 AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  if (handler) {
    int iMinSlice=0; 
    int iMaxSlice=1;
    int iMinPart=0;
    int iMaxPart=1;
    TString fileWriterInput;
    TString esdWriterInput;
    for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
      TString trackerInput;
      for (int part=iMinPart; part<=iMaxPart; part++) {
	TString arg, publisher, cf;

	// digit publisher components
	arg.Form("-slice %d -partition %d", slice, part);
	publisher.Form("DP_%02d_%d", slice, part);
	handler->CreateConfiguration(publisher.Data(), "TPCDigitPublisher", NULL , arg.Data());

	// cluster finder components
	cf.Form("CF_%02d_%d", slice, part);
	handler->CreateConfiguration(cf.Data(), "TPCClusterFinderUnpacked", publisher.Data(), "pp-run timebins 446");
	if (trackerInput.Length()>0) trackerInput+=" ";
	trackerInput+=cf;
      }
      TString tracker;
      // tracker finder components
      tracker.Form("TR_%02d", slice);
      handler->CreateConfiguration(tracker.Data(), "TPCSliceTracker", trackerInput.Data(), "pp-run bfield 0.5");

      // input for the global file writer
      if (fileWriterInput.Length()>0) fileWriterInput+=" ";
      fileWriterInput+=trackerInput;

      // input for the esd writer
      if (esdWriterInput.Length()>0) esdWriterInput+=" ";
      esdWriterInput+=tracker;
    }

    // the writer configuration
    handler->CreateConfiguration("sink1", "FileWriter"   , fileWriterInput.Data(), "-specfmt -subdir=test_%d -blcknofmt=_0x%x -idfmt=_0x%08x");
    // the esd writer configuration
    handler->CreateConfiguration("esd-writer", "TPCEsdWriter"   , esdWriterInput.Data(), "-datafile AliESDs.root");
  }
  return 0;
}

const char* AliHLTTPCAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
						    AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  return NULL;
  //return "sink1";
  //return "esd-writer";
}

const char* AliHLTTPCAgent::GetRequiredComponentLibraries() const
{
  // see header file for class documentation
  return NULL;
}

int AliHLTTPCAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
  // see header file for class documentation
  if (!pHandler) return -EINVAL;
  pHandler->AddComponent(new AliHLTTPCDigitDumpComponent);
  pHandler->AddComponent(new AliHLTTPCEsdWriterComponent::AliWriter);
  pHandler->AddComponent(new AliHLTTPCEsdWriterComponent::AliConverter);

  return 0;
}

int AliHLTTPCAgent::GetHandlerDescription(AliHLTComponentDataType dt,
					  AliHLTUInt32_t spec,
					  AliHLTOUTHandlerDesc& desc) const
{
  // see header file for class documentation
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC)) {
    int slice=AliHLTTPCDefinitions::GetMinSliceNr(spec);
    int part=AliHLTTPCDefinitions::GetMinPatchNr(spec);
    if (slice==AliHLTTPCDefinitions::GetMaxSliceNr(spec) &&
	part==AliHLTTPCDefinitions::GetMaxPatchNr(spec)) {
      desc=AliHLTOUTHandlerDesc(kRawReader, dt, GetModuleId());
      return 1;
    } else {
      HLTWarning("handler can not process merged data from multiple ddls:"
		 " min slice %d, max slice %d, min part %d, max part %d",
		 slice, AliHLTTPCDefinitions::GetMaxSliceNr(spec),
		 part, AliHLTTPCDefinitions::GetMaxPatchNr(spec));
      return 0;
    }
  }
  return 0;
}

AliHLTOUTHandler* AliHLTTPCAgent::GetOutputHandler(AliHLTComponentDataType dt,
						   AliHLTUInt32_t /*spec*/)
{
  // see header file for class documentation
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC)) {
    if (!fRawDataHandler) {
      fRawDataHandler=new AliHLTTPCAgent::AliHLTTPCRawDataHandler;
    }
    fNofRawDataHandler++;
    return fRawDataHandler;
  }
  return NULL;
}

int AliHLTTPCAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance)
{
  // see header file for class documentation
  if (pInstance==NULL) return -EINVAL;

  if (pInstance==fRawDataHandler) {
    if (--fNofRawDataHandler<=0) {
      delete fRawDataHandler;
      fRawDataHandler=NULL;
    }
  }
  return 0;
}

AliHLTTPCAgent::AliHLTTPCRawDataHandler::AliHLTTPCRawDataHandler()
{
  // see header file for class documentation
}

AliHLTTPCAgent::AliHLTTPCRawDataHandler::~AliHLTTPCRawDataHandler()
{
  // see header file for class documentation
}

int AliHLTTPCAgent::AliHLTTPCRawDataHandler::ProcessData(AliHLTOUT* pData)
{
  // see header file for class documentation
  if (!pData) return -EINVAL;
  AliHLTComponentDataType dt=kAliHLTVoidDataType;
  AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
  int iResult=pData->GetDataBlockDescription(dt, spec);
  if (iResult>=0) {
    int slice=AliHLTTPCDefinitions::GetMinSliceNr(spec);
    int part=AliHLTTPCDefinitions::GetMinPatchNr(spec);
    if (slice==AliHLTTPCDefinitions::GetMaxSliceNr(spec) &&
	part==AliHLTTPCDefinitions::GetMaxPatchNr(spec)) {
      iResult=768;
      if (part>1) iResult+=72+4*slice+(part-2);
      else iResult+=2*slice+part;
    } else {
      HLTError("handler can not process merged data from multiple ddls:"
	       " min slice %d, max slice %d, min part %d, max part %d",
	       slice, AliHLTTPCDefinitions::GetMaxSliceNr(spec),
	       part, AliHLTTPCDefinitions::GetMaxPatchNr(spec));
      iResult=-EBADMSG;
    }
  }
  return iResult;
}
