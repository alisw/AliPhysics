// @(#) $Id$

/** @file   AliHLTAgentUtil.cxx
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTUtil library
*/

#include "AliHLTAgentUtil.h"
#include "AliHLTConfiguration.h"

/** global instance for agent registration */
AliHLTAgentUtil gAliHLTAgentUtil;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAgentUtil)

AliHLTAgentUtil::AliHLTAgentUtil()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTAgentUtil::~AliHLTAgentUtil()
{
  // see header file for class documentation
}

int AliHLTAgentUtil::CreateConfigurations(AliHLTConfigurationHandler* handler,
					  AliRunLoader* runloader) const
{
  // see header file for class documentation
  if (handler) {
  }
  return 0;
}

const char* AliHLTAgentUtil::GetTopConfigurations(AliRunLoader* runloader) const
{
  // see header file for class documentation
  return NULL;
}

const char* AliHLTAgentUtil::GetRequiredComponentLibraries() const
{
  // see header file for class documentation
  return "libAliHLTUtil.so";
}
