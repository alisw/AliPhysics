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

/** @file   AliHLTSampleOfflineSinkComponent.cxx
    @author Matthias Richter
    @date   
    @brief  A sample offline sink component.
*/

#include "AliHLTSampleOfflineSinkComponent.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliLog.h"
#include "TTree.h"

/** global instance for agent registration */
AliHLTSampleOfflineSinkComponent gAliHLTSampleOfflineSinkComponent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTSampleOfflineSinkComponent)

AliHLTSampleOfflineSinkComponent::AliHLTSampleOfflineSinkComponent()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTSampleOfflineSinkComponent::~AliHLTSampleOfflineSinkComponent()
{
  // see header file for class documentation
}

const char* AliHLTSampleOfflineSinkComponent::GetComponentID()
{
  // see header file for class documentation
  return "SampleOfflineDataSink";
}

void AliHLTSampleOfflineSinkComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponent* AliHLTSampleOfflineSinkComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTSampleOfflineSinkComponent;
}

int AliHLTSampleOfflineSinkComponent::DoInit( int /*argc*/, const char** /*argv*/ )
{
  // see header file for class documentation
  int iResult=0;
  return iResult;
}

int AliHLTSampleOfflineSinkComponent::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  return iResult;
}

int AliHLTSampleOfflineSinkComponent::DumpEvent(const AliHLTComponentEventData& /*evtData*/,
						AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation
  int iResult=0;
  AliInfoStream() << "dump event " << GetEventCount() << endl;
  TTree* pTree=(TTree*)GetFirstInputObject(kAliHLTAnyDataType, "TTree");
  while (pTree) {
    AliInfoStream() << " got TTree object with " << pTree->GetEntries() << " entries" << endl;
    pTree=(TTree*)GetNextInputObject();
  }
  return iResult;
}

int AliHLTSampleOfflineSinkComponent::FillESD(int eventNo, AliRunLoader* /*runLoader*/, AliESDEvent* /*esd*/) 
{
  // see header file for class documentation
  int iResult=0;
  AliInfoStream() << "event " << eventNo << endl;
  return iResult;
}

