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

/** @file   AliHLTTPCDigitDumpComponent.cxx
    @author Matthias Richter
    @date   
    @brief  Special file writer converting TPC digit input to ASCII. */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTPCDigitDumpComponent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCDigitDumpComponent)

AliHLTTPCDigitDumpComponent::AliHLTTPCDigitDumpComponent()
  :
  AliHLTFileWriter()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCDigitDumpComponent::~AliHLTTPCDigitDumpComponent()
{
  // see header file for class documentation
}

const char* AliHLTTPCDigitDumpComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCDigitDump";
}

void AliHLTTPCDigitDumpComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponent* AliHLTTPCDigitDumpComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCDigitDumpComponent;
}

int AliHLTTPCDigitDumpComponent::InitWriter()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCDigitDumpComponent::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    //
    if (argument.CompareTo("-whatsoever")==0) {
    } else {
      iResult=-EINVAL;
    }
  }

  return iResult;
}

int AliHLTTPCDigitDumpComponent::CloseWriter()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCDigitDumpComponent::DumpEvent( const AliHLTComponentEventData& evtData,
			 const AliHLTComponentBlockData* blocks, 
			 AliHLTComponentTriggerData& trigData )
{
  // see header file for class documentation
  int iResult=0;
  return iResult;
}
