// $Id$

///**************************************************************************
///* This file is property of and copyright by the ALICE HLT Project        * 
///* ALICE Experiment at CERN, All rights reserved.                         *
///*                                                                        *
///* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
///*                  for The ALICE HLT Project.                            *
///*                                                                        *
///* Permission to use, copy, modify and distribute this software and its   *
///* documentation strictly for non-commercial purposes is hereby granted   *
///* without fee, provided that the above copyright notice appears in all   *
///* copies and that both the copyright notice and this permission notice   *
///* appear in the supporting documentation. The authors make no claims     *
///* about the suitability of this software for any purpose. It is          *
///* provided "as is" without express or implied warranty.                  *
///**************************************************************************/

/// @file   AliHLTComponentConfiguration.cxx
/// @author Matthias Richter
/// @date   2010-11-26
/// @brief  HLT configuration description for a single component.
/// @note   The class is used in Offline (AliRoot) context

#include "AliHLTComponentConfiguration.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTConfiguration)

AliHLTComponentConfiguration::AliHLTComponentConfiguration()
  : AliHLTConfiguration()
  , fLibrary()
  , fNodeNames()
  , fOnlineCommand()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTComponentConfiguration::AliHLTComponentConfiguration(const char* id,
							   const char* component,
							   const char* sources,
							   const char* arguments)
  : AliHLTConfiguration(id, component, sources, arguments)
  , fLibrary()
  , fNodeNames()
  , fOnlineCommand()
{
  // constructor
}

AliHLTComponentConfiguration::AliHLTComponentConfiguration(const AliHLTComponentConfiguration& src)
  : AliHLTConfiguration(src)
  , fLibrary(src.fLibrary)
  , fNodeNames(src.fNodeNames)
  , fOnlineCommand(src.fOnlineCommand)
{
  // copy constructor
}

AliHLTComponentConfiguration& AliHLTComponentConfiguration::operator=(const AliHLTComponentConfiguration& src)
{
  // assignment operator
  if (this==&src) return *this;
  AliHLTConfiguration::operator=(src);
  fLibrary=src.fLibrary;
  fNodeNames=src.fNodeNames;
  fOnlineCommand=src.fOnlineCommand;
  return *this;
}

AliHLTComponentConfiguration::~AliHLTComponentConfiguration()
{
  // destructor
}

void AliHLTComponentConfiguration::SetOnlineCommand(const char* cmd)
{
  // set the online command string

  fOnlineCommand=cmd;
}

void AliHLTComponentConfiguration::PrintStatus() const
{
  // see header file for function documentation
  HLTLogKeyword("configuration status");
  if (!fLibrary.IsNull()) HLTMessage("  - component library: \"%s\"",
    fLibrary.Data());
  else HLTMessage("  - component library missing");
  if (!fOnlineCommand.IsNull()) HLTMessage("  - online command: \"%s\"",
    fOnlineCommand.Data());
  else HLTMessage("  - online command missing");
  if (!fNodeNames.IsNull()) HLTMessage("  - online nodes: \"%s\"",
    fNodeNames.Data());
  else HLTMessage("  - no online nodes");
}

void AliHLTComponentConfiguration::Print(const char* option) const
{
  // print information
  AliHLTConfiguration::Print(option);
  if (option && strcmp(option, "status")==0) {
    PrintStatus();
  }
  else {
  HLTLogKeyword("configuration");
  HLTMessage("component library %s, online command %s, online nodes %s",
	     GetComponentLibrary(),
	     GetOnlineCommand(),
	     GetNodeSettings()
	     );
   }
}
