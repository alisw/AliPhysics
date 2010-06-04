// @(#) $Id: AliHLTAgentFMD.cxx 25820 2008-05-16 11:47:09Z richterm $

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Hans Hjersing Dalsgaard                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTAgentFMD.cxx
    @author Hans Hjersing Dalsgaard
    @date   
    @brief  Agent of the libAliHLTFMD library
*/

#include <cassert>
#include "AliHLTAgentFMD.h"
#include "TSystem.h"

#include "AliHLTFMDReconstructionComponent.h"

AliHLTAgentFMD gAliHLTAgentFMD;


ClassImp(AliHLTAgentFMD)

AliHLTAgentFMD::AliHLTAgentFMD()
  :
  AliHLTModuleAgent("FMD")
{

}

AliHLTAgentFMD::~AliHLTAgentFMD()
{

}

const char* AliHLTAgentFMD::GetRequiredComponentLibraries() const
{
  // see header file for class documentation
  return "libAliHLTUtil.so libAliHLTFMD.so";
}

int AliHLTAgentFMD::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
  // see header file for class documentation
  assert(pHandler);
  if (!pHandler) return -EINVAL;
  pHandler->AddComponent(new AliHLTFMDReconstructionComponent);
 
  return 0;
}


