// $Id$
///**************************************************************************
///* This file is property of and copyright by the ALICE HLT Project        * 
///* All rights reserved.                                                   *
///*                                                                        *
///* Primary Author: Jenny Wagner  (jwagner@cern.ch)                        *
///*                                                                        *
///* Permission to use, copy, modify and distribute this software and its   *
///* documentation strictly for non-commercial purposes is hereby granted   *
///* without fee, provided that the above copyright notice appears in all   *
///* copies and that both the copyright notice and this permission notice   *
///* appear in the supporting documentation. The authors make no claims     *
///* about the suitability of this software for any purpose. It is          * 
///* provided "as is" without express or implied warranty.                  *
///**************************************************************************

/// @file   ALIHLTCOMPHuffmanCodeData.cxx
/// @author Jenny Wagner
/// @date   29-08-2007
/// @brief  Data class for the Huffman code table of 10-bit-ADC-values
///

#include "AliHLTCOMPHuffmanCodeData.h"
#include "AliHLTStdIncludes.h"

#if __GNUC__ >= 3
using namespace std;
#endif


ClassImp(AliHLTCOMPHuffmanCodeData)

/** construction without any arguments (used for isolated tests) */
AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeData()
  :
  famplitude(0),
  fhuffmancode(0),
  fvalidcodelength(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

/** HuffmanCodeData destructor  */
AliHLTCOMPHuffmanCodeData::~AliHLTCOMPHuffmanCodeData()
{
  /* destructor, see header file for class documentation */
}

void AliHLTCOMPHuffmanCodeData::SetHuffmanCodeData(const AliHLTCOMPHuffmanCodeStruct& codetableentry)
{
  // see header file for class documentation
  famplitude = codetableentry.famplitude;
  fhuffmancode = codetableentry.fhuffmancode;
  fvalidcodelength = codetableentry.fvalidcodelength;
}

AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* AliHLTCOMPHuffmanCodeData::GetHuffmanCodeData(AliHLTCOMPHuffmanCodeStruct* codetableentry) const
{
  // see header file for class documentation
  codetableentry->famplitude = famplitude;
  codetableentry->fhuffmancode = fhuffmancode;
  codetableentry->fvalidcodelength = fvalidcodelength;

  return codetableentry;
}
