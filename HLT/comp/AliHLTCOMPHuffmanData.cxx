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

#include "AliHLTCOMPHuffmanData.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTCOMPHuffmanData)

/// @file   ALIHLTCOMPHuffmanData.cxx
/// @author Jenny Wagner
/// @date   29-08-2007
///         changed on 03-12-2007
/// @brief see header file for documentation
///

AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanData()
  : TObject()
  , AliHLTLogging()
  , fOccurrenceTable()
  , fCodeTable()
  , fOrigin(kAliHLTVoidDataOrigin)
  , fDataSpec(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTCOMPHuffmanData::~AliHLTCOMPHuffmanData()
{
  // destructor, see header file for class documentation
}

int AliHLTCOMPHuffmanData::InitHuffmanData(const AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* occurrencetable,
					   const AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* codetable)
{
  // init internal occurrencetable and codetable
  for(Int_t ii = 0; ii < TIMEBINS; ii++) {
    fOccurrenceTable[ii].SetHuffmanOccurrenceData(occurrencetable[ii]);
    fCodeTable[ii].SetHuffmanCodeData(codetable[ii]);
  }

  return 0;
}

AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* AliHLTCOMPHuffmanData::GetOccurrenceTable(AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* occurrencetable) const
{
  // copy occurrencetable to target
  for (Int_t ii = 0; ii < TIMEBINS; ii++) {
    fOccurrenceTable[ii].GetHuffmanOccurrenceData(&(occurrencetable[ii]));
  }

  return 0;
}

AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* AliHLTCOMPHuffmanData::GetCodeTable(AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* codetable) const
{
  // copy code table to target
  for (Int_t ii = 0; ii < TIMEBINS; ii++) {
    fCodeTable[ii].GetHuffmanCodeData(&(codetable[ii]));
  }

  return 0;
}

int AliHLTCOMPHuffmanData::SetOCDBSpecifications(const TString& origin, Int_t dataspec)
{
  // set data block properties
  fOrigin = origin;
  fDataSpec = dataspec;

  return 0;
}
