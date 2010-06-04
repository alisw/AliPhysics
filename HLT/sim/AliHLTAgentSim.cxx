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

/** @file   AliHLTAgentSim.cxx
    @author Matthias Richter
    @date   
    @brief  Agent of the libHLTsim library
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include <cerrno>
#include "AliHLTAgentSim.h"

// header files of library components
#include "AliHLTOUTComponent.h"

/** global instance for agent registration */
AliHLTAgentSim gAliHLTAgentSim;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAgentSim)

AliHLTAgentSim::AliHLTAgentSim()
  :
  AliHLTModuleAgent("sim")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTAgentSim::~AliHLTAgentSim()
{
  // see header file for class documentation
}

int AliHLTAgentSim::CreateConfigurations(AliHLTConfigurationHandler* /*handler*/,
					  AliRawReader* /*rawReader*/,
					  AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  return 0;
}

const char* AliHLTAgentSim::GetReconstructionChains(AliRawReader* /*rawReader*/,
						     AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  return NULL;
}

const char* AliHLTAgentSim::GetRequiredComponentLibraries() const
{
  // see header file for class documentation
  return NULL;
}

int AliHLTAgentSim::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
  // see header file for class documentation
  assert(pHandler);
  if (!pHandler) return -EINVAL;
  pHandler->AddComponent(new AliHLTOUTComponent);
  return 0;
}
