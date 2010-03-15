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

/** @file   AliHLTGlobalDCSPublisherComponent.cxx
    @author Matthias Richter
    @date   20010-03-10
    @brief  DIM publisher component for global HLT data
*/

#include "AliHLTGlobalDCSPublisherComponent.h"
#include "AliHLTDimServer.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalDCSPublisherComponent)

AliHLTGlobalDCSPublisherComponent::AliHLTGlobalDCSPublisherComponent()
  : AliHLTDataSink()
  , fpServer(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTGlobalDCSPublisherComponent::~AliHLTGlobalDCSPublisherComponent()
{
  // see header file for class documentation

  // file list and file name list are owner of their objects and
  // delete all the objects
}

const char* AliHLTGlobalDCSPublisherComponent::GetComponentID()
{
  // see header file for class documentation
  return "DCSPublisher";
}

void AliHLTGlobalDCSPublisherComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponent* AliHLTGlobalDCSPublisherComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTGlobalDCSPublisherComponent;
}

int AliHLTGlobalDCSPublisherComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;
  iResult=ConfigureFromArgumentString(argc, argv);

  return iResult;
}

int AliHLTGlobalDCSPublisherComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;

  // -servername
  if (argc==0) return 0;
  int i=0;
  const char* serverName=NULL;
  const char* dimdns=NULL;
  TString argument=argv[0];
  if (argument.CompareTo("-servername")==0) {
    if (++i>=argc) return -EPROTO;
    serverName=argv[i];

    // --dimdns
  } else if (argument.BeginsWith("-dimdns")) {
    if (++i>=argc) return -EPROTO;
    dimdns=argv[i];

  } else {
    return -EINVAL;
  }

  fpServer=new AliHLTDimServer(serverName);
  if (!fpServer) return -ENOMEM;
  if ((iResult=fpServer->Init(dimdns))>=0) {
    // add services

    iResult=fpServer->Start();
  }

  return iResult;
}

int AliHLTGlobalDCSPublisherComponent::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  if (!fpServer) return -ENODEV;
  fpServer->Stop();
  delete fpServer;
  fpServer=NULL;

  return iResult;
}

int AliHLTGlobalDCSPublisherComponent::DumpEvent( const AliHLTComponentEventData& /*evtData*/,
						  AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;
  return iResult;
}
