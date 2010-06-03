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

/** @file   AliHLTCompAgent.cxx
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTComp library
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include <cerrno>
#include "AliHLTCompAgent.h"

// header files of library components
#include "AliHLTCOMPHuffmanAltroComponent.h"
#include "AliHLTCOMPHuffmanAltroCalibComponent.h"

// header file of the module preprocessor
#include "AliHLTCompPreprocessor.h"

/** global instance for agent registration */
AliHLTCompAgent gAliHLTCompAgent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTCompAgent)

AliHLTCompAgent::AliHLTCompAgent()
  :
  AliHLTModuleAgent("Comp")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTCompAgent::~AliHLTCompAgent()
{
  // see header file for class documentation
}

int AliHLTCompAgent::CreateConfigurations(AliHLTConfigurationHandler* /*handler*/,
					  AliRawReader* /*rawReader*/,
					  AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  return 0;
}

const char* AliHLTCompAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
						     AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  return NULL;
}

const char* AliHLTCompAgent::GetRequiredComponentLibraries() const
{
  // see header file for class documentation

  // libAliHLTUtil.so for AliRawReaderPublisher
  //return "libAliHLTUtil.so";
  return NULL;
}

int AliHLTCompAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
  // see header file for class documentation
  assert(pHandler);
  if (!pHandler) return -EINVAL;
  // use fCompressionSwitch = true for decompressed inputtype (i.e. compressed output)
  pHandler->AddComponent(new AliHLTCOMPHuffmanAltroComponent(true));
  // use fCompressionSwitch = false for compressed inputtype (i.e. decompressed output)
  pHandler->AddComponent(new AliHLTCOMPHuffmanAltroComponent(false));
  pHandler->AddComponent(new AliHLTCOMPHuffmanAltroCalibComponent);

  return 0;
}

AliHLTModulePreprocessor* AliHLTCompAgent::GetPreprocessor()
{
  // see header file for class documentation
  return new AliHLTCompPreprocessor;
}
