// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/// @file   AliHLTRCUAgent.cxx
/// @author Matthias Richter
/// @date   
/// @brief  Agent of the libAliHLTRCU library
///

#include <cassert>
#include "AliHLTRCUAgent.h"
#include "AliHLTDAQ.h"

// header files of library components
#include "AliHLTAltroChannelSelectorComponent.h"
#include "AliHLTAltroTimebinAverageComponent.h"

/** global instance for agent registration */
AliHLTRCUAgent gAliHLTRCUAgent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTRCUAgent)

AliHLTRCUAgent::AliHLTRCUAgent()
  : AliHLTModuleAgent("RCU")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTRCUAgent::~AliHLTRCUAgent()
{
  // see header file for class documentation
}

int AliHLTRCUAgent::CreateConfigurations(AliHLTConfigurationHandler* handler,
					 AliRawReader* rawReader,
					 AliRunLoader* runloader) const
{
  // add configurations for the RCU library
  if (handler) {
    if (rawReader) {
      // AliSimulation: use the AliRawReaderPublisher if the raw reader is available
      // Alireconstruction: indicated by runloader==NULL, run always on raw data
      int iMinSlice=0; 
      int iMaxSlice=35;
      int iMinPart=0;
      int iMaxPart=5;
      TString sinkChannelSelectors;
      for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
	for (int part=iMinPart; part<=iMaxPart; part++) {
	  TString arg, publisher, selector;

	  // publisher component
	  int ddlno=AliHLTDAQ::DdlIDOffset(3);
	  if (part>1) ddlno+=72+4*slice+(part-2);
	  else ddlno+=2*slice+part;
	  
	  publisher.Form("RCU-DP_%02d_%d", slice, part);
	  arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec %s -silent", ddlno, AliHLTDAQ::HLTSpecificationFromDdlID(ddlno).c_str());
	  handler->CreateConfiguration(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

	  // selector component
	  selector.Form("RCU-chselector_%02d_%d", slice, part);
	  arg="-signal-threshold 1";
	  handler->CreateConfiguration(selector.Data(), "AltroChannelSelector", publisher.Data(), arg.Data());


	  if (sinkChannelSelectors.Length()>0) sinkChannelSelectors+=" ";
	  sinkChannelSelectors+=selector;
	}
      }
      handler->CreateConfiguration("RCU-channelselect", "BlockFilter", sinkChannelSelectors.Data(), "");
    }
  }
  return 0;
}

const char* AliHLTRCUAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
						    AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  return NULL;
}

const char* AliHLTRCUAgent::GetRequiredComponentLibraries() const
{
  // see header file for class documentation
  return NULL;
}

int AliHLTRCUAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
  // see header file for class documentation
  assert(pHandler);
  if (!pHandler) return -EINVAL;
  pHandler->AddComponent(new AliHLTAltroChannelSelectorComponent);
  pHandler->AddComponent(new AliHLTAltroTimebinAverageComponent);
  return 0;
}
