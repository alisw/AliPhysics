// $Id: $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
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

/** @file   AliHLTOUTHandlerIgnore.cxx
    @author Artur Szostak
    @date   7 Jan 2011
    @brief  Implementation of simple HLTOUT handler for ignoring data blocks.
*/

#include "AliHLTOUTHandlerIgnore.h"

ClassImp(AliHLTOUTHandlerIgnore)

int AliHLTOUTHandlerIgnore::ProcessData(AliHLTOUT* /*data*/)
{
  // Ignores the data and returns zero.
  return 0;
}
