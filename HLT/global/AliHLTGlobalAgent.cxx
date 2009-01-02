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

/** @file   AliHLTGlobalAgent.cxx
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTGlobal library
*/

#include <cassert>
#include "AliHLTGlobalAgent.h"

// header files of library components
#include "AliHLTGlobalTrackMergerComponent.h"

/** global instance for agent registration */
AliHLTGlobalAgent gAliHLTGlobalAgent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalAgent)

AliHLTGlobalAgent::AliHLTGlobalAgent()
  :
  AliHLTModuleAgent("Global")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTGlobalAgent::~AliHLTGlobalAgent()
{
  // see header file for class documentation
}

int AliHLTGlobalAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
  // see header file for class documentation
  assert(pHandler);
  if (!pHandler) return -EINVAL;
  pHandler->AddComponent(new AliHLTGlobalTrackMergerComponent);
  return 0;
}

int AliHLTGlobalAgent::GetHandlerDescription(AliHLTComponentDataType dt,
					   AliHLTUInt32_t /*spec*/,
					  AliHLTOUTHandlerDesc& desc) const
{
  // see header file for class documentation

  // handler for the HLT readou list and global data data blocks {'HLTRDLST':'HLT '}
  if (dt==AliHLTComponentDataTypeInitializer("HLTRDLST", kAliHLTDataOriginOut) ||
      dt==AliHLTComponentDataTypeInitializer("HLTTRGDT", kAliHLTDataOriginOut)) {
      desc=AliHLTOUTHandlerDesc(kProprietary, dt, GetModuleId());
      return 1;
  }

  return 0;
}

AliHLTOUTHandler* AliHLTGlobalAgent::GetOutputHandler(AliHLTComponentDataType dt,
						   AliHLTUInt32_t /*spec*/)
{
  // see header file for class documentation

  // handler for the HLT readou list and global data data blocks {'HLTRDLST':'HLT '}
  if (dt==AliHLTComponentDataTypeInitializer("HLTRDLST", kAliHLTDataOriginOut) ||
      dt==AliHLTComponentDataTypeInitializer("HLTTRGDT", kAliHLTDataOriginOut)) {
    return NULL;
  }

  return NULL;
}

int AliHLTGlobalAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance)
{
  // see header file for class documentation
  if (pInstance==NULL) return -EINVAL;

  return 0;
}
