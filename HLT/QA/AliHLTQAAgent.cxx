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

/** @file   AliHLTQAAgent.cxx
    @author Matthias Richter
    @date   2010-03-10
    @brief  Module Agent for HLTqadm library
*/
#include "AliHLTQAAgent.h"

/** global instance for agent registration */
AliHLTQAAgent gAliHLTQAAgent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTQAAgent)

AliHLTQAAgent::AliHLTQAAgent()
  : AliHLTModuleAgent("QA")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTQAAgent::~AliHLTQAAgent()
{
  // see header file for class documentation
}

const char* AliHLTQAAgent::GetQAPlugins() const
{
  // see header file for class documentation
  return "AliHLTTPCQADataMaker";
}
