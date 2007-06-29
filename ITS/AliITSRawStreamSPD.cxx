/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliLog.h"

ClassImp(AliITSRawStreamSPD)


  // this map has to change, waiting for the new geometry
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
  fEventNumber(-1),fChipAddr(0),fHalfStaveNr(0),fCol(0),fRow(0),
  fData(0),fOffset(0),fHitCount(0),
  fDataChar1(0),fDataChar2(0),fDataChar3(0),fDataChar4(0),
  fFirstWord(kTRUE),fPrevEventId(0xffffffff)
{
  // create an object to read ITS SPD raw digits
  fRawReader->Select("ITSSPD");
  // reset calib header words
  for (UInt_t iword=0; iword<kCalHeadLenMax; iword++) {
    fCalHeadWord[iword]=0xffffffff;
  }
  NewEvent();
}

Bool_t AliITSRawStreamSPD::ReadNextShort() 
{
  // read next 16 bit word into fData
  if (fFirstWord) {
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
    if ((*fRawReader->GetEventId())!=fPrevEventId) { // if new event...
      NewEvent();
      fPrevEventId=(*fRawReader->GetEventId());
    }
    fFirstWord=kFALSE;
  }
  else {
    fFirstWord=kTRUE;
    fData = fDataChar1+(fDataChar2<<8);
  }

  return kTRUE;
}

Bool_t AliITSRawStreamSPD::ReadNextInt() 
{
  // reads next 32 bit into fDataChar1..4 
  // (if first 16 bits read already, just completes the present word)
  if (fFirstWord) {
    if (ReadNextShort() && ReadNextShort()) {
      return kTRUE;
    }
  }
  else {
    if (ReadNextShort()) {
      return kTRUE;
    }
  }
  return kFALSE;
}

void AliITSRawStreamSPD::NewEvent()
{
  // call this to reset flags for a new event
  for (UInt_t eqId=0; eqId<20; eqId++) {
    fCalHeadRead[eqId]=kFALSE;
  }
  fEventNumber=-1;
}

Bool_t AliITSRawStreamSPD::ReadCalibHeader()
{
  // read the extra calibration header
  // returns kTRUE if the header is present and has length > 0

  Int_t ddlID = fRawReader->GetDDLID();
  if (ddlID==-1) { // we may need to read one word to get the blockAttr
    if (!ReadNextShort()) return kFALSE;
    ddlID = fRawReader->GetDDLID();
  }
  UChar_t attr = fRawReader->GetBlockAttributes();
  if (ddlID>=0 && ddlID<20) fCalHeadRead[ddlID]=kTRUE;
  if ((attr & 0x40) == 0x40) { // is the header present?
    if (ReadNextInt()) {
      // length of cal header:
      UInt_t calLen = fDataChar1+(fDataChar2<<8)+(fDataChar3<<16)+(fDataChar4<<24);
      if (calLen>kCalHeadLenMax) {
	fRawReader->AddMajorErrorLog(kCalHeaderLengthErr,Form("Header length %d > max = %d",calLen,kCalHeadLenMax));
	AliWarning(Form("Header length problem. %d > %d (max)",calLen,kCalHeadLenMax));
	return kFALSE;
      }
      else if (calLen>0) {
	for (UInt_t iword=0; iword<calLen; iword++) {
	  if (ReadNextInt()) {
	    fCalHeadWord[iword] = fDataChar1+(fDataChar2<<8)+(fDataChar3<<16)+(fDataChar4<<24);
	  }
	  else {
	    fRawReader->AddMajorErrorLog(kCalHeaderLengthErr,"header length problem");
	    AliWarning("header length problem");
	    return kFALSE;
	  }
	}
	return kTRUE;
      }
    }
  }

  return kFALSE;
}

Bool_t AliITSRawStreamSPD::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left

  Int_t ddlID=-1;
  fPrevModuleID = fModuleID;

  while (ReadNextShort()) {

    ddlID = fRawReader->GetDDLID();
    if (ddlID>=0 && ddlID<20) {
      if (!fCalHeadRead[ddlID]) {
	ReadCalibHeader();
      }
    }
    else {
      fRawReader->AddMajorErrorLog(kDDLNumberErr,Form("Wrong DDL number %d",ddlID)); 
      AliWarning(Form("DDL number error (= %d) , setting it to 19", ddlID));
      ddlID=19;
    }

    if ((fData & 0xC000) == 0x4000) {         // header
      fHitCount = 0;
      UShort_t eventNumber = (fData >> 4) & 0x007F;
      if (fEventNumber < 0) {
	fEventNumber = eventNumber;
      } 
      else if (eventNumber != fEventNumber) {
	fRawReader->AddMajorErrorLog(kEventNumberErr,Form("Reading event number %d instead of %d",eventNumber,fEventNumber));
	AliWarning(Form("mismatching event numbers: %d != %d",eventNumber, fEventNumber));
      }
      fChipAddr = fData & 0x000F;
      if (fChipAddr>9) {
	fRawReader->AddMajorErrorLog(kChipAddrErr,Form("Overflow chip address %d - set to 9",fChipAddr));
	AliWarning(Form("overflow chip addr (= %d) , setting it to 9", fChipAddr));
	fChipAddr=9;
      }
      fHalfStaveNr = (fData & 0x3800)>>11;
      if (fHalfStaveNr>5 || fRawReader->TestBlockAttribute(fHalfStaveNr)) {
	fRawReader->AddMajorErrorLog(kStaveNumberErr,Form("Half stave number error %d - set to 5",fHalfStaveNr));
	AliWarning(Form("half stave number error(=%d) , setting it to 5", fHalfStaveNr));
	fHalfStaveNr=5;
      }
      // translate  ("online") ddl, hs, chip nr  to  ("offline") module id :
      fModuleID = fgkDDLModuleMap[ddlID][fHalfStaveNr*2+fChipAddr/5];
      fOffset = 32 * (fChipAddr % 5);
    } 
    else if ((fData & 0xC000) == 0x0000) {    // trailer
      UShort_t hitCount = fData & 0x1FFF;
      if (hitCount != fHitCount){
	fRawReader->AddMajorErrorLog(kNumbHitsErr,Form("Number of hits %d instead of %d",hitCount,fHitCount));
	AliWarning(Form("wrong number of hits: %d != %d", fHitCount, hitCount));
      }
    }
    else if ((fData & 0xC000) == 0x8000) {    // pixel hit
      fHitCount++;
      fCol = (fData & 0x001F);
      fRow = (fData >> 5) & 0x00FF;

      // translate  ("online") chipcol, chiprow  to  ("offline") col (coord1), row (coord2): 
      // This will change, waiting for new geometry!!!
      fCoord1 = fCol;
      //      if      (fModuleID < 80 && ddlID < 10) fCoord1=31-fCoord1;
      //      else if (fModuleID >=80 && ddlID >=10) fCoord1=31-fCoord1;
      fCoord1 += fOffset;
      //      if (ddlID>=10) fCoord1=159-fCoord1;
      fCoord2 = fRow;
      //      if (fModuleID<80) fCoord2=255-fCoord2;

      return kTRUE;
    } 
    else {                                    // fill word
      if ((fData & 0xC000) != 0xC000) {
	fRawReader->AddMajorErrorLog(kWrongWordErr,"Wrong fill word");
	AliWarning("wrong fill word!");
      }
    }

  }

  return kFALSE;
}
Bool_t AliITSRawStreamSPD::GetHalfStavePresent(UInt_t hs) {
  // Reads the half stave present status from the block attributes
  Int_t ddlID = fRawReader->GetDDLID();
  if (ddlID==-1) {
    fRawReader->AddMinorErrorLog(kHalfStaveStatusErr,"DDL ID = -1. Cannot read block attributes.");
    AliWarning("DDL ID = -1. Cannot read block attributes. Return kFALSE.");
    return kFALSE;
  }
  else {
    if (hs>=6) {
      fRawReader->AddMinorErrorLog(kHalfStaveStatusErr,Form( "HS >= 6 requested (%d). Return kFALSE.",hs));
      AliWarning(Form("HS >= 6 requested (%d). Return kFALSE.",hs));
      return kFALSE;
    }
    UChar_t attr = fRawReader->GetBlockAttributes();
    if (((attr>>hs) & 0x01) == 0x01) { // bit set means not present
      return kFALSE;
    }
    else {
      return kTRUE;
    }
  }
}

Bool_t AliITSRawStreamSPD::GetHhalfStaveScanned(UInt_t hs) const {
  if (hs<6) return (Bool_t)((fCalHeadWord[0]>>(6+hs)) & (0x00000001));
  else return kFALSE;
}
Bool_t AliITSRawStreamSPD::GetHchipPresent(UInt_t hs, UInt_t chip) const {
  if (hs<6 && chip<10) return ((( fCalHeadWord[hs/3+3]>>((hs%3)*10+chip)) & 0x00000001) == 1);
  else return kFALSE;
}
UInt_t AliITSRawStreamSPD::GetHdacHigh(UInt_t hs) const {
  if (hs<6) return (fCalHeadWord[hs/2+7]>>(24-16*(hs%2)) & 0x000000ff);
  else return 0;
}
UInt_t AliITSRawStreamSPD::GetHdacLow(UInt_t hs) const {
  if (hs<6) return (fCalHeadWord[hs/2+7]>>(16-16*(hs%2)) & 0x000000ff);
  else return 0;
}
UInt_t AliITSRawStreamSPD::GetHTPAmp(UInt_t hs) const {
  if (hs<6) return fCalHeadWord[hs+10];
  else return 0;
}
Bool_t AliITSRawStreamSPD::GetHminTHchipPresent(UInt_t chip) const {
  if (chip<10) return ((( fCalHeadWord[7]>>(16+chip)) & 0x00000001) == 1);
  else return kFALSE;
}
