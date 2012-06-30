// $Id$

//**************************************************************************
//* This file is property of and copyright by the                          * 
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

/// @file   AliHLTOUTHandler.cxx
/// @author Matthias Richter
/// @date   
/// @brief  Base class implementation of HLTOUT handlers.
///

#include "AliHLTOUTHandler.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTHandler)

AliHLTOUTHandler::AliHLTOUTHandler()
  :
  fState(kHandlerOK)
{ 
  // Base class declaration of HLT output handlers.
  // The library implementation of the AliHLTModuleAgent allows to generate
  // handlers for data blocks of the HLT output. This can be the output of
  // the real HLT coming from the HLTOUT nodes, or simulated HLT output.
  // Note: The created instance of AliHLTOUTHandler is deleted by the framework.
}

AliHLTOUTHandler::~AliHLTOUTHandler()
{
  // destructor
}

int AliHLTOUTHandler::GetProcessedData(const AliHLTUInt8_t* &pData)
{
  // get pointer to processed data
  pData=NULL;
  return 0;
}

int AliHLTOUTHandler::ReleaseProcessedData(const AliHLTUInt8_t* /*pData*/, int /*size*/)
{
  // release the data pointer previously retrieved by GetProcessedData
  return 0;
}

int AliHLTOUTHandler::FinishEvent()
{
  // cleanup the current event processing.
  return 0;
}
