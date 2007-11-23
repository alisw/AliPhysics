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

/** @file   ALIHLTCOMPHuffmanData.cxx
    @author Jenny Wagner
    @date   29-08-2007
*/

#include "AliHLTCOMPHuffmanData.h"
#include "AliHLTStdIncludes.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTCOMPHuffmanData)

/** construction without any arguments (used for isolated tests) */
AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanData()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

/** HuffmanData destructor  */
AliHLTCOMPHuffmanData::~AliHLTCOMPHuffmanData()
{
  /* destructor, see header file for class documentation */
}

/** get data from OCDB and write them into instance of HuffmanData */
void AliHLTCOMPHuffmanData::InitHuffmanData(AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanData_t* occurrencetable, AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCode_t* codetable)
{
  // see header file for class documentation
  for(Int_t ii = 0; ii < TIMEBINS; ii++)
    {
      fOccurrenceTable[ii].SetHuffmanOccurrenceData(occurrencetable[ii]);
      
      fCodeTable[ii].SetHuffmanCodeData(codetable[ii]);
    }
}

AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanData_t* AliHLTCOMPHuffmanData::GetOccurrenceTable(AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanData_t* occurrencetable) 
{
  // see header file for class documentation
  for (Int_t ii = 0; ii < TIMEBINS; ii++)
    {
      fOccurrenceTable[ii].GetHuffmanOccurrenceData(&(occurrencetable[ii]));
    }

  return occurrencetable;
  
}

AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCode_t* AliHLTCOMPHuffmanData::GetCodeTable(AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCode_t* codetable)
{
  // see header file for class documentation
  for (Int_t ii = 0; ii < TIMEBINS; ii++)
    {
      fCodeTable[ii].GetHuffmanCodeData(&(codetable[ii]));
    }

  return codetable;
}
