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

#include "AliHLTCOMPHuffmanData.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTCOMPHuffmanData)

/** @file   ALIHLTCOMPHuffmanData.cxx
    @author Jenny Wagner
    @date   29-08-2007
            changed on 03-12-2007
    @brief see header file for documentation
*/

/** construction without any arguments (used for isolated tests) */
AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanData()
  :
  fOrigin(kAliHLTVoidDataOrigin),
  fDataSpec(0)
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
void AliHLTCOMPHuffmanData::InitHuffmanData(AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* occurrencetable, AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* codetable)
{
  // see header file for class documentation
  for(Int_t ii = 0; ii < TIMEBINS; ii++)
    {
      fOccurrenceTable[ii].SetHuffmanOccurrenceData(occurrencetable[ii]);
      
      fCodeTable[ii].SetHuffmanCodeData(codetable[ii]);
    }
}

AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* AliHLTCOMPHuffmanData::GetOccurrenceTable(AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* occurrencetable) 
{
  // see header file for class documentation
  for (Int_t ii = 0; ii < TIMEBINS; ii++)
    {
      fOccurrenceTable[ii].GetHuffmanOccurrenceData(&(occurrencetable[ii]));
    }

  return occurrencetable;
  
}

AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* AliHLTCOMPHuffmanData::GetCodeTable(AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* codetable)
{
  // see header file for class documentation
  for (Int_t ii = 0; ii < TIMEBINS; ii++)
    {
      fCodeTable[ii].GetHuffmanCodeData(&(codetable[ii]));
    }

  return codetable;
}

Int_t AliHLTCOMPHuffmanData::SetOCDBSpecifications(TString origin, Int_t dataspec)
{
  // see header file for class documentation
  fOrigin = origin;
  fDataSpec = dataspec;

  return 0;
}

