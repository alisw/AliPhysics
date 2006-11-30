/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ITS SPD digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliITSRawStreamSPD.h"
#include "AliRawReader.h"

ClassImp(AliITSRawStreamSPD)


const Int_t AliITSRawStreamSPD::fgkDDLModuleMap[kDDLsNumber][kModulesPerDDL] = {
  { 0, 1, 4, 5, 80, 81, 84, 85, 88, 89, 92, 93},
  { 8, 9,12,13, 96, 97,100,101,104,105,108,109},
  {16,17,20,21,112,113,116,117,120,121,124,125},
  {24,25,28,29,128,129,132,133,136,137,140,141},
  {32,33,36,37,144,145,148,149,152,153,156,157},
  {40,41,44,45,160,161,164,165,168,169,172,173},
  {48,49,52,53,176,177,180,181,184,185,188,189},
  {56,57,60,61,192,193,196,197,200,201,204,205},
  {64,65,68,69,208,209,212,213,216,217,220,221},
  {72,73,76,77,224,225,228,229,232,233,236,237},
  { 2, 3, 6, 7, 82, 83, 86, 87, 90, 91, 94, 95},
  {10,11,14,15, 98, 99,102,103,106,107,110,111},
  {18,19,22,23,114,115,118,119,122,123,126,127},
  {26,27,30,31,130,131,134,135,138,139,142,143},
  {34,35,38,39,146,147,150,151,154,155,158,159},
  {42,43,46,47,162,163,166,167,170,171,174,175},
  {50,51,54,55,178,179,182,183,186,187,190,191},
  {58,59,62,63,194,195,198,199,202,203,206,207},
  {66,67,70,71,210,211,214,215,218,219,222,223},
  {74,75,78,79,226,227,230,231,234,235,238,239}
};


AliITSRawStreamSPD::AliITSRawStreamSPD(AliRawReader* rawReader) :
  AliITSRawStream(rawReader),
  fData(0),
  fDDLNumber(-1),
  fEventNumber(-1),
  fOffset(0),
  fHitCount(0),
  fDataChar1(0),
  fDataChar2(0),
  fDataChar3(0),
  fDataChar4(0),
  fFirstWord(kTRUE)
{
// create an object to read ITS SPD raw digits

  fRawReader->Select("ITSSPD");
}

Bool_t AliITSRawStreamSPD::ReadNextShort() 
{
  if (fFirstWord) {
    fFirstWord=kFALSE;
    Bool_t b1 = fRawReader->ReadNextChar(fDataChar1);
    if (!b1) return kFALSE;
    Bool_t  b2, b3, b4;
    b2 = fRawReader->ReadNextChar(fDataChar2);
    b3 = fRawReader->ReadNextChar(fDataChar3);
    b4 = fRawReader->ReadNextChar(fDataChar4);
    if (!(b2 && b3 && b4)) {
      return kFALSE;
    }
    fData = fDataChar3+(fDataChar4<<8);
  }
  else {
    fFirstWord=kTRUE;
    fData = fDataChar1+(fDataChar2<<8);
  }

  return kTRUE;
}

void AliITSRawStreamSPD::SkipCalibHeader()
{
  // Checks if there is an extra calibration header 
  // present in the raw data. Reads past this in that case.

  fRawReader->ReadHeader(); // need this to get access to the block attributes
  UChar_t attr = fRawReader->GetBlockAttributes();
  if ((attr & 0x40) == 0x40) { // is the header present?
    Bool_t  b1, b2, b3, b4;
    b1 = fRawReader->ReadNextChar(fDataChar1);
    b2 = fRawReader->ReadNextChar(fDataChar2);
    b3 = fRawReader->ReadNextChar(fDataChar3);
    b4 = fRawReader->ReadNextChar(fDataChar4);
    if (b1 && b2 && b3 && b4) {
      // length of cal header:
      UInt_t calLen = fDataChar1+(fDataChar2<<8)+(fDataChar3<<16)+(fDataChar4<<24);
      // read past the cal header:
      UInt_t tmpData;
      for (UInt_t iword=0; iword<calLen; iword++) {
	fRawReader->ReadNextInt(tmpData);
      }
    }
  }
}

Bool_t AliITSRawStreamSPD::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left

  fPrevModuleID = fModuleID;
  //  while (fRawReader->ReadNextShort(fData)) {
  while (ReadNextShort()) {
    if ((fData & 0xC000) == 0x4000) {         // header
      fHitCount = 0;
      UShort_t eventNumber = (fData >> 4) & 0x007F;
      if (fEventNumber < 0) {
	fEventNumber = eventNumber;
      } 
      else if (eventNumber != fEventNumber) {
	Warning("Next", "mismatching event numbers: %d != %d", 
		eventNumber, fEventNumber);
      }
      UShort_t chipAddr = fData & 0x000F;
      if (chipAddr>9) {
	Error("Next", "overflow chip addr (= %d) , setting it to 9", chipAddr);
	chipAddr=9;
      }
      UShort_t halfStaveNr = (fData & 0x3800)>>11;
      if (halfStaveNr>5 || fRawReader->TestBlockAttribute(halfStaveNr)) {
	Error("Next", "half stave number error(= %d) , setting it to 5", halfStaveNr);
	halfStaveNr=5;
      }
      fDDLNumber = fRawReader->GetDDLID();
      if (fDDLNumber>19 || fDDLNumber<0) {
	Error("Next", "DDL number error (= %d) , setting it to 19", fDDLNumber);
	fDDLNumber=19;
      }
      fModuleID = fgkDDLModuleMap[fDDLNumber][halfStaveNr*2+chipAddr/5];
      fOffset = 32 * (chipAddr % 5);
    } 
    else if ((fData & 0xC000) == 0x0000) {    // trailer
      UShort_t hitCount = fData & 0x1FFF;
      if (hitCount != fHitCount) Error("Next", "wrong number of hits: %d != %d", fHitCount, hitCount);
    } 
    else if ((fData & 0xC000) == 0x8000) {    // pixel hit
      fHitCount++;
      fCoord1 = (fData & 0x001F) + fOffset;
      fCoord2 = (fData >> 5) & 0x00FF;
      return kTRUE;
    } 
    else {                                    // fill word
      if ((fData & 0xC000) != 0xC000) Error("Next", "wrong fill word!");
    }

  }

  return kFALSE;
}
