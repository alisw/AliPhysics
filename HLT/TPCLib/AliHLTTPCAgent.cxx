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

/** global instance for agent registration */
AliHLTTPCAgent gAliHLTTPCAgent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCAgent)

AliHLTTPCAgent::AliHLTTPCAgent()
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
					  AliRunLoader* runloader) const
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

const char* AliHLTTPCAgent::GetLocalRecConfigurations(AliRunLoader* runloader) const
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
