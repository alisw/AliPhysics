// $Id$

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

/** @file   AliHLTOUTComponent.cxx
    @author Matthias Richter
    @date   
    @brief  The HLTOUT data sink component similar to HLTOUT nodes */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTOUTComponent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTComponent)

AliHLTOUTComponent::AliHLTOUTComponent()
  :
  AliHLTOfflineDataSink()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTOUTComponent::~AliHLTOUTComponent()
{
  // see header file for class documentation

  // file list and file name list are owner of their objects and
  // delete all the objects
}

const char* AliHLTOUTComponent::GetComponentID()
{
  // see header file for class documentation
  return "HLTOUT";
}

void AliHLTOUTComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponent* AliHLTOUTComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTOUTComponent;
}

int AliHLTOUTComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    {
      HLTError("unknown argument %s", argument.Data());
      break;
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  if (iResult>=0) {
  }

  return iResult;
}

int AliHLTOUTComponent::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  
  return iResult;
}

int AliHLTOUTComponent::DumpEvent( const AliHLTComponentEventData& evtData,
			 const AliHLTComponentBlockData* blocks, 
			 AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;
  HLTInfo("write %d output blocks", evtData.fBlockCnt);
  for (int n=0; n<(int)evtData.fBlockCnt; n++ ) {
    
  }

  return iResult;
}

int AliHLTOUTComponent::FillESD(int /*eventNo*/, AliRunLoader* /*runLoader*/, AliESDEvent* /*esd*/)
{
  // see header file for class documentation
  int iResult=0;
  return iResult;
}
