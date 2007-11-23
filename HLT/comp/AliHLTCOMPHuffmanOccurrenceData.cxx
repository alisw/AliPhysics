//-*- Mode: C++ -*-
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Author: Jenny Wagner  (jwagner@cern.ch)                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   ALIHLTCOMPHuffmanOccurrenceData.cxx
    @author Jenny Wagner
    @date   29-08-2007
    @brief  Data class for the occurrence table of 10-bit-ADC-values
*/

#include "AliHLTCOMPHuffmanOccurrenceData.h"
#include "AliHLTStdIncludes.h"

#if __GNUC__ >= 3
using namespace std;
#endif


ClassImp(AliHLTCOMPHuffmanOccurrenceData)

/** construction without any arguments (used for isolated tests) */
AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanOccurrenceData()
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

void AliHLTCOMPHuffmanOccurrenceData::SetHuffmanOccurrenceData(AliHLTCOMPHuffmanData_t const& occurrencetableentry)
{
  // see header file for class documentation
  amplitude = occurrencetableentry.amplitude;
  abundance = occurrencetableentry.abundance;
  code = occurrencetableentry.code;
}

AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanData_t* AliHLTCOMPHuffmanOccurrenceData::GetHuffmanOccurrenceData(AliHLTCOMPHuffmanData_t* occurrencetableentry)
{
  // see header file for class documentation
  occurrencetableentry->amplitude = amplitude;
  occurrencetableentry->abundance = abundance;
  occurrencetableentry->code = code;

  return occurrencetableentry;
}
