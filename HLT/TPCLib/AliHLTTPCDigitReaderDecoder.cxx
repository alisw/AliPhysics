// $Id: AliHLTTPCDigitReaderDecoder.cxx,v 1.18 2007/11/26 23:19:47 richterm Exp $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
//*                  Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTTPCDigitReaderDecoder.cxx
    @author Kenneth Aamodt, Matthias Richter
    @date   
    @brief  DigitReader implementation for the fast ALTRO Decoder
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCDigitReaderDecoder.h"
#include "AliAltroDecoder.h"
#include "AliAltroData.h"
#include "AliAltroBunch.h"

ClassImp(AliHLTTPCDigitReaderDecoder)

AliHLTTPCDigitReaderDecoder::AliHLTTPCDigitReaderDecoder()
  :
  AliHLTTPCDigitReader()  
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCDigitReaderDecoder::~AliHLTTPCDigitReaderDecoder()
{
  // see header file for class documentation
}

int AliHLTTPCDigitReaderDecoder::InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice)
{
  // see header file for class documentation
  return 0;
}

bool AliHLTTPCDigitReaderDecoder::NextChannel()
{
  // see header file for class documentation
  return false;
}

int AliHLTTPCDigitReaderDecoder::NextBunch()
{
  // see header file for class documentation
  return 0;
}

bool AliHLTTPCDigitReaderDecoder::NextSignal()
{
  // see header file for class documentation
  return false;
}

AliHLTUInt32_t* AliHLTTPCDigitReaderDecoder::GetSignals()
{
  // see header file for class documentation
  return NULL;
}

int AliHLTTPCDigitReaderDecoder::GetRow()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCDigitReaderDecoder::GetPad()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCDigitReaderDecoder::GetSignal()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCDigitReaderDecoder::GetTime()
{
  // see header file for class documentation
  return 0;
}
