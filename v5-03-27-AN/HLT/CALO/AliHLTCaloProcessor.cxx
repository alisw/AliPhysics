// $Id: AliHLTCaloProcessor.cxx 35966 2009-10-26 12:47:19Z odjuvsla $

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Author:  Per Thomas Hille  <perthi@fys.uio.no>                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTCaloProcessor.h"

const AliHLTComponentDataType AliHLTCaloProcessor::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array

AliHLTCaloProcessor::AliHLTCaloProcessor():AliHLTProcessor(), 
					   fCaloEventCount(0)
{
  lineNumber[0] = '\0';
}


AliHLTCaloProcessor::~AliHLTCaloProcessor()
{

}


const char*
AliHLTCaloProcessor::IntToChar(int number)
{
  // comment
  snprintf(lineNumber, 256, "%d", number);
  return lineNumber;
}

