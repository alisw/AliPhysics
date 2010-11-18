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

/** @file   AliHLTEsdCollectorComponent.cxx
    @author Matthias Richter
    @date   
    @brief  Base class for writer components to store data in a ROOT file

                                                                          */

#include "AliHLTEsdCollectorComponent.h"
#include "AliHLTEsdManager.h"
#include "TString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTEsdCollectorComponent)

AliHLTEsdCollectorComponent::AliHLTEsdCollectorComponent()
  : AliHLTFileWriter()
  , fpManager(NULL)
  , fTreeName()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  // all event go into the same file, but there are individual files for
  // different data blocks
  SetMode(kConcatenateEvents);
}

AliHLTEsdCollectorComponent::~AliHLTEsdCollectorComponent()
{
  // see header file for class documentation
}

int AliHLTEsdCollectorComponent::InitWriter()
{
  // see header file for class documentation
  int iResult=0;

  // choose .root as default extension
  if (GetExtension().IsNull()) SetExtension("root");

  if ((fpManager=AliHLTEsdManager::New())) {
    TString option="-writelocal";
    if (!GetDirectory().IsNull()) {
      option+=" -directory=";
      option+=GetDirectory();
    }
    if (!fTreeName.IsNull()) {
      option+=" -treename=";
      option+=fTreeName;
    }
    iResult=fpManager->SetOption(option.Data());
  } else {
    HLTError("can not find AliHLTEsdManager class descriptor");
    iResult=-ENODEV;
  }

  return iResult;
}

int AliHLTEsdCollectorComponent::CloseWriter()
{
  // see header file for class documentation
  if (fpManager!=NULL) {
    AliHLTEsdManager::Delete(fpManager);
    fpManager=NULL;
  }
  return 0;
}

int AliHLTEsdCollectorComponent::DumpEvent( const AliHLTComponentEventData& /*evtData*/,
					    const AliHLTComponentBlockData* /*blocks*/, 
					    AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;
  if (!IsDataEvent() && !CheckMode(kWriteAllEvents)) return 0;
  if (!fpManager) return -ENODEV;

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock();
       pBlock && iResult>=0;
       pBlock=GetNextInputBlock()) {
    if (pBlock->fDataType!=kAliHLTDataTypeESDObject &&
	pBlock->fDataType!=kAliHLTDataTypeESDTree) continue;
    HLTInfo("writing ESD, data type %s event %d", (DataType2Text(pBlock->fDataType).c_str()), GetEventCount());
    iResult=fpManager->WriteESD(reinterpret_cast<const AliHLTUInt8_t*>(pBlock->fPtr),
				pBlock->fSize, pBlock->fDataType, NULL, GetEventCount());
  }
  return iResult;
}

int AliHLTEsdCollectorComponent::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation
  if (argc<=0) return 0;
  int iResult=-EINVAL;
  int iArg=0;
  TString argument=argv[iArg];

  // -treename
  if (argument.CompareTo("-treename")==0) {
    if (++iArg==argc) {
      HLTError("expecting parameter for argument '-treename'");
      iResult=-EPROTO;
    } else {
      fTreeName=argv[iArg];
      iResult=++iArg;
    }
  }
  return iResult;
}
