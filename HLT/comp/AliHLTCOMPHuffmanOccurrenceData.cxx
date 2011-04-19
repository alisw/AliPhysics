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

/// @file   ALIHLTCOMPHuffmanOccurrenceData.cxx
/// @author Jenny Wagner
/// @date   29-08-2007
/// @brief  Data class for the occurrence table of 10-bit-ADC-values
///

#include "AliHLTCOMPHuffmanOccurrenceData.h"
#include "AliHLTStdIncludes.h"

#if __GNUC__ >= 3
using namespace std;
#endif


ClassImp(AliHLTCOMPHuffmanOccurrenceData)

/** construction without any arguments (used for isolated tests) */
AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanOccurrenceData()
  :
  famplitude(0),
  fabundance(0),
  fcode(2) // has to be initialised to two since reasonable values are 0 and 1 !!!
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

/** HuffmanOccurrenceData destructor  */
AliHLTCOMPHuffmanOccurrenceData::~AliHLTCOMPHuffmanOccurrenceData()
{
  /* destructor, see header file for class documentation */
}

void AliHLTCOMPHuffmanOccurrenceData::SetHuffmanOccurrenceData(const AliHLTCOMPHuffmanDataStruct& occurrencetableentry)
{
  // see header file for class documentation
  famplitude = occurrencetableentry.famplitude;
  fabundance = occurrencetableentry.fabundance;
  fcode = occurrencetableentry.fcode;
}

AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* AliHLTCOMPHuffmanOccurrenceData::GetHuffmanOccurrenceData(AliHLTCOMPHuffmanDataStruct* occurrencetableentry) const
{
  // see header file for class documentation
  occurrencetableentry->famplitude = famplitude;
  occurrencetableentry->fabundance = fabundance;
  occurrencetableentry->fcode = fcode;

  return occurrencetableentry;
}
