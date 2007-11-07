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


const Int_t AliITSRawStreamSPD::fgkDDLModuleMap[kDDLsNumber][kModulesPerDDL] = {
  { 4, 5, 0, 1, 80, 81, 84, 85, 88, 89, 92, 93},
  {12,13, 8, 9, 96, 97,100,101,104,105,108,109},
  {20,21,16,17,112,113,116,117,120,121,124,125},
  {28,29,24,25,128,129,132,133,136,137,140,141},
  {36,37,32,33,144,145,148,149,152,153,156,157},
  {44,45,40,41,160,161,164,165,168,169,172,173},
  {52,53,48,49,176,177,180,181,184,185,188,189},
  {60,61,56,57,192,193,196,197,200,201,204,205},
  {68,69,64,65,208,209,212,213,216,217,220,221},
  {76,77,72,73,224,225,228,229,232,233,236,237},
  { 7, 6, 3, 2, 83, 82, 87, 86, 91, 90, 95, 94},
  {15,14,11,10, 99, 98,103,102,107,106,111,110},
  {23,22,19,18,115,114,119,118,123,122,127,126},
  {31,30,27,26,131,130,135,134,139,138,143,142},
  {39,38,35,34,147,146,151,150,155,154,159,158},
  {47,46,43,42,163,162,167,166,171,170,175,174},
  {55,54,51,50,179,178,183,182,187,186,191,190},
  {63,62,59,58,195,194,199,198,203,202,207,206},
  {71,70,67,66,211,210,215,214,219,218,223,222},
  {79,78,75,74,227,226,231,230,235,234,239,238}
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
      fModuleID = GetOfflineModuleFromOnline(ddlID,fHalfStaveNr,fChipAddr);
      //      fOffset = 32 * (fChipAddr % 5);
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

      fCoord1 = GetOfflineColFromOnline(ddlID,fHalfStaveNr,fChipAddr,fCol);
      fCoord2 = GetOfflineRowFromOnline(ddlID,fHalfStaveNr,fChipAddr,fRow);

      // translate  ("online") chipcol, chiprow  to  ("offline") col (coord1), row (coord2): 
      // This will change, waiting for new geometry!!!
      //      fCoord1 = fCol;
      //      if      (fModuleID < 80 && ddlID < 10) fCoord1=31-fCoord1;
      //      else if (fModuleID >=80 && ddlID >=10) fCoord1=31-fCoord1;
      //      fCoord1 += fOffset;
      //      if (ddlID>=10) fCoord1=159-fCoord1;
      //      fCoord2 = fRow;
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




Int_t AliITSRawStreamSPD::GetModuleNumber(UInt_t iDDL, UInt_t iModule) {
  if (iDDL<20 && iModule<12) return fgkDDLModuleMap[iDDL][iModule];
  else return 240;
}




Bool_t AliITSRawStreamSPD::OfflineToOnline(UInt_t module, UInt_t colM, UInt_t rowM, UInt_t& eq, UInt_t& hs, UInt_t& chip, UInt_t& col, UInt_t& row) {
  // converts offline coordinates to online
  eq = GetOnlineEqIdFromOffline(module);
  hs = GetOnlineHSFromOffline(module);
  chip = GetOnlineChipFromOffline(module,colM);
  col = GetOnlineColFromOffline(module,colM);
  row = GetOnlineRowFromOffline(module,rowM);
  if (eq>=20 || hs>=6 || chip>=10 || col>=32 || row>=256) return kFALSE;
  else return kTRUE;
}


Bool_t AliITSRawStreamSPD::OnlineToOffline(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row, UInt_t& module, UInt_t& colM, UInt_t& rowM) {
  // converts online coordinates to offline
  module = GetOfflineModuleFromOnline(eq,hs,chip);
  colM = GetOfflineColFromOnline(eq,hs,chip,col);
  rowM = GetOfflineRowFromOnline(eq,hs,chip,row);
  if (module>=240 || colM>=160 || rowM>=256) return kFALSE;
  else return kTRUE;
}


UInt_t AliITSRawStreamSPD::GetOnlineEqIdFromOffline(UInt_t module) {
  // offline->online (eq)
  for (UInt_t eqId=0; eqId<20; eqId++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      if (GetModuleNumber(eqId,iModule)==(Int_t)module) return eqId;
    }
  }
  return 20; // error
}

UInt_t AliITSRawStreamSPD::GetOnlineHSFromOffline(UInt_t module) {
  // offline->online (hs)
  for (UInt_t eqId=0; eqId<20; eqId++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      if (GetModuleNumber(eqId,iModule)==(Int_t)module) return iModule/2;
    }
  }
  return 6; // error
}

UInt_t AliITSRawStreamSPD::GetOnlineChipFromOffline(UInt_t module, UInt_t colM) {
  // offline->online (chip)
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      if (GetModuleNumber(eq,iModule)==(Int_t)module) {
	if (eq<10) { // side A
	  return (159-colM)/32 + 5*(iModule%2);
	}
	else { // side C
	  return colM/32 + 5*(iModule%2);
	}
      }
    }
  }
  return 10; // error
}

UInt_t AliITSRawStreamSPD::GetOnlineColFromOffline(UInt_t module, UInt_t colM) {
  // offline->online (col)
  if (module<80) { // inner layer
    return colM%32;
  }
  else if (module<240) { // outer layer
    return (159-colM)%32;
  }
  return 32; // error
}

UInt_t AliITSRawStreamSPD::GetOnlineRowFromOffline(UInt_t module, UInt_t rowM) {
  // offline->online (row)
  if (module<80) { // inner layer
    return rowM;
  }
  else if (module<240) { // outer layer
    return (255-rowM);
  }
  return 256; // error
}





UInt_t AliITSRawStreamSPD::GetOfflineModuleFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip) {
  // online->offline (module)
  if (eqId<20 && hs<6 && chip<10) return fgkDDLModuleMap[eqId][hs*2+chip/5];
  else return 240;
}

UInt_t AliITSRawStreamSPD::GetOfflineColFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col) {
  // online->offline (col)
  if (eqId>=20 || hs>=6 || chip>=10 || col>=32) return 160; // error
  UInt_t offset = 32 * (chip % 5);
  if (hs<2) {
    if (eqId<10) {
      return 159 - (31-col + offset); // inner layer, side A
    }
    else {
      return col + offset; // inner layer, side C
    }
  }
  else {
    if (eqId<10) {
      return 159 - (col + offset); // outer layer, side A
    }
    else {
      return 31-col + offset; // outer layer, side C
    }
  }
}

UInt_t AliITSRawStreamSPD::GetOfflineRowFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t row) {
  // online->offline (row)
  if (eqId>=20 || hs>=6 || chip>=10 || row>=256) return 256; // error
  if (hs<2) {
    return row;
  }
  else {
    return 255-row;
  }
}
