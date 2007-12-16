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
//__________________________________________________________________________
AliITSRawStreamSPD::AliITSRawStreamSPD(AliRawReader* rawReader) :
  AliITSRawStream(rawReader),
  fEventCounter(-1),fChipAddr(0),fHalfStaveNr(0),fCol(0),fRow(0),
  fData(0),fOffset(0),fHitCount(0),
  fDataChar1(0),fDataChar2(0),fDataChar3(0),fDataChar4(0),
  fFirstWord(kTRUE),fPrevEventId(0xffffffff),
  fEqPLBytesRead(0),fEqPLChipHeadersRead(0),fEqPLChipTrailersRead(0),
  fHeaderOrTrailerReadLast(kFALSE),fExpectedHeaderTrailerCount(0),
  fFillOutOfSynch(kFALSE),fDDLID(-1),fLastDDLID(-1),fAdvancedErrorLog(kFALSE),
  fAdvLogger(NULL)
{
  // create an object to read ITS SPD raw digits
  fRawReader->Select("ITSSPD");
  // reset calib header words
  for (UInt_t iword=0; iword<kCalHeadLenMax; iword++) {
    fCalHeadWord[iword]=0xffffffff;
  }
  NewEvent();
}
//__________________________________________________________________________
AliITSRawStreamSPD::AliITSRawStreamSPD(const AliITSRawStreamSPD& rstream) :
  AliITSRawStream(rstream.fRawReader),
  fEventCounter(-1),fChipAddr(0),fHalfStaveNr(0),fCol(0),fRow(0),
  fData(0),fOffset(0),fHitCount(0),
  fDataChar1(0),fDataChar2(0),fDataChar3(0),fDataChar4(0),
  fFirstWord(kTRUE),fPrevEventId(0xffffffff),
  fEqPLBytesRead(0),fEqPLChipHeadersRead(0),fEqPLChipTrailersRead(0),
  fHeaderOrTrailerReadLast(kFALSE),fExpectedHeaderTrailerCount(0),
  fFillOutOfSynch(kFALSE),fDDLID(-1),fLastDDLID(-1),fAdvancedErrorLog(kFALSE),
  fAdvLogger(NULL)
{
  // copy constructor
  AliError("Copy constructor should not be used.");
}
//__________________________________________________________________________
AliITSRawStreamSPD& AliITSRawStreamSPD::operator=(const AliITSRawStreamSPD& rstream) {
  // assignment operator
  if (this!=&rstream) {}
  AliError("Assignment opertator should not be used.");
  return *this;
}
//__________________________________________________________________________
Bool_t AliITSRawStreamSPD::ReadNextShort() {
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
  fEqPLBytesRead+=2;
  
  return kTRUE;
}
//__________________________________________________________________________
Bool_t AliITSRawStreamSPD::ReadNextInt() {
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
//__________________________________________________________________________
void AliITSRawStreamSPD::NewEvent() {
  // call this to reset flags for a new event
  for (UInt_t eqId=0; eqId<20; eqId++) {
    fCalHeadRead[eqId]=kFALSE;
  }
  fEventCounter = -1;
  fDDLID = -1;
  fLastDDLID = -1;
}
//__________________________________________________________________________
Int_t AliITSRawStreamSPD::ReadCalibHeader() {
  // read the extra calibration header
  // returns the length of the header if it is present, -1 otherwise

  Int_t ddlID = fRawReader->GetDDLID();
  if (ddlID==-1) { // we may need to read one word to get the blockAttr
    if (!ReadNextShort()) return -1;
    ddlID = fRawReader->GetDDLID();
  }
  // reset flags and counters
  fEqPLBytesRead = 2;
  fEqPLChipHeadersRead = 0;
  fEqPLChipTrailersRead = 0;
  fHeaderOrTrailerReadLast = kFALSE;
  fFillOutOfSynch = kFALSE;

  // check what number of chip headers/trailers to expect
  fExpectedHeaderTrailerCount = 0;
  for (UInt_t hs=0; hs<6; hs++) {
    if (!fRawReader->TestBlockAttribute(hs)) fExpectedHeaderTrailerCount+=10;
  }

  if (ddlID>=0 && ddlID<20) fCalHeadRead[ddlID]=kTRUE;
  if (fRawReader->TestBlockAttribute(6)) { // is the calib header present?
    if (ReadNextInt()) {
      // length of cal header:
      UInt_t calLen = fDataChar1+(fDataChar2<<8)+(fDataChar3<<16)+(fDataChar4<<24);
      if (calLen>kCalHeadLenMax) {
	TString errMess = Form("Header length %d > max = %d",calLen,kCalHeadLenMax);
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kCalHeaderLengthErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kCalHeaderLengthErr,ddlID,-1,-1,errMess.Data());
	return -1;
      }
      else {
	for (UInt_t iword=0; iword<calLen; iword++) {
	  if (ReadNextInt()) {
	    fCalHeadWord[iword] = fDataChar1+(fDataChar2<<8)+(fDataChar3<<16)+(fDataChar4<<24);
	  }
	  else {
	    TString errMess = "Header length problem";
	    AliError(errMess.Data());
	    fRawReader->AddMajorErrorLog(kCalHeaderLengthErr,errMess.Data());
	    if (fAdvancedErrorLog) fAdvLogger->ProcessError(kCalHeaderLengthErr,ddlID,-1,-1,errMess.Data());
	    return -1;
	  }
	}
	return calLen;
      }
    }
  }

  return -1;
}
//__________________________________________________________________________
Bool_t AliITSRawStreamSPD::Next() {
// read the next raw digit
// returns kFALSE if there is no digit left

  fPrevModuleID = fModuleID;

  while (ReadNextShort()) {

    fLastDDLID = fDDLID;
    fDDLID = fRawReader->GetDDLID();
    if (fDDLID>=0 && fDDLID<20) {
      if (!fCalHeadRead[fDDLID]) {
	if (fLastDDLID!=-1) {  // if not the first equipment for this event
	  fEqPLBytesRead -= 2;
	  CheckHeaderAndTrailerCount(fLastDDLID);
	}
	if (ReadCalibHeader()>=0) continue;
      }
    }
    else {
      TString errMess = Form("Error in DDL number (=%d) , setting it to 19",fDDLID);
      AliError(errMess.Data());
      fRawReader->AddMajorErrorLog(kDDLNumberErr,errMess.Data()); 
      if (fAdvancedErrorLog) fAdvLogger->AddMessage(errMess.Data());
      fDDLID=19;
    }

    if ((fData & 0xC000) == 0x4000) {         // header
      if (fHeaderOrTrailerReadLast) {
	TString errMess = "Chip trailer missing";
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kTrailerMissingErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kTrailerMissingErr,fLastDDLID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
      }
      fHeaderOrTrailerReadLast = kTRUE;
      fEqPLChipHeadersRead++;
      fHitCount = 0;
      UShort_t eventCounter = (fData >> 4) & 0x007F;
      if (fEventCounter < 0) {
	fEventCounter = eventCounter;
      } 
      else if (eventCounter != fEventCounter) {
	TString errMess;
	if (fEqPLChipHeadersRead==1) {
	  errMess = Form("Mismatching event counters between this equipment and the previous: %d != %d",
			 eventCounter,fEventCounter);
	}
	else {
	  errMess = Form("Mismatching event counters between this chip header and the previous: %d != %d",
			 eventCounter,fEventCounter);
	}
	fEventCounter = eventCounter;
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kEventCounterErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kEventCounterErr,fDDLID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
      }
      fChipAddr = fData & 0x000F;
      if (fChipAddr>9) {
	TString errMess = Form("Overflow chip address %d - set to 0",fChipAddr);
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kChipAddrErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kChipAddrErr,fDDLID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
	fChipAddr=0;
      }
      fHalfStaveNr = (fData & 0x3800)>>11;
      if (fHalfStaveNr>5 || fRawReader->TestBlockAttribute(fHalfStaveNr)) {
	TString errMess = Form("Half stave number error: %d - set to 0",fHalfStaveNr);
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kHSNumberErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kHSNumberErr,fDDLID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
	fHalfStaveNr=0;
      }
      // translate  ("online") ddl, hs, chip nr  to  ("offline") module id :
      fModuleID = GetOfflineModuleFromOnline(fDDLID,fHalfStaveNr,fChipAddr);
    } 
    else if ((fData & 0xC000) == 0x0000) {    // trailer
      if ( (fEqPLBytesRead+fFillOutOfSynch*2)%4 != 0 ) {
	TString errMess = "Fill word is missing";
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kFillMissingErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kFillMissingErr,fDDLID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
	if (fFillOutOfSynch) fFillOutOfSynch = kFALSE;
	else fFillOutOfSynch = kTRUE;
      }
      if (!fHeaderOrTrailerReadLast) {
	TString errMess = "Trailer without previous header";
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kTrailerWithoutHeaderErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kTrailerWithoutHeaderErr,fLastDDLID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
      }
      fHeaderOrTrailerReadLast = kFALSE;
      fEqPLChipTrailersRead++;
      UShort_t hitCount = fData & 0x1FFF;
      if (hitCount != fHitCount){
	TString errMess = Form("Number of hits %d, while %d expected",fHitCount,hitCount);
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kNumberHitsErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kNumberHitsErr,fDDLID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
      }
    }
    else if ((fData & 0xC000) == 0x8000) {    // pixel hit
      if (!fHeaderOrTrailerReadLast) {
	TString errMess = "Chip header missing";
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kHeaderMissingErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kHeaderMissingErr,fDDLID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
	fHeaderOrTrailerReadLast = kTRUE;
      }
      fHitCount++;
      fCol = (fData & 0x001F);
      fRow = (fData >> 5) & 0x00FF;

      fCoord1 = GetOfflineColFromOnline(fDDLID,fHalfStaveNr,fChipAddr,fCol);
      fCoord2 = GetOfflineRowFromOnline(fDDLID,fHalfStaveNr,fChipAddr,fRow);

      return kTRUE;
    } 
    else {                                    // fill word
      if ((fData & 0xC000) != 0xC000) {
	TString errMess = "Wrong fill word!";
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kWrongFillWordErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kWrongFillWordErr,fDDLID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
      }
      if ( (fEqPLBytesRead+fFillOutOfSynch*2)%4 != 2 ) {
	TString errMess = "Fill word is unexpected";
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kFillUnexpectErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kFillUnexpectErr,fDDLID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
	if (fFillOutOfSynch) fFillOutOfSynch = kFALSE;
	else fFillOutOfSynch = kTRUE;
      }
    }

  }
  if (fDDLID>=0 && fDDLID<20) {
    CheckHeaderAndTrailerCount(fDDLID);
  }
  return kFALSE;
}
//__________________________________________________________________________
void AliITSRawStreamSPD::CheckHeaderAndTrailerCount(Int_t ddlID) {
  // Checks that the number of header and trailers found for the ddl are as expected
  if (fEqPLChipHeadersRead != fExpectedHeaderTrailerCount) {
    TString errMess = Form("Chip header count inconsistent %d != %d (expected) for ddl %d",
			   fEqPLChipHeadersRead,fExpectedHeaderTrailerCount,ddlID);
    AliError(errMess.Data());
    fRawReader->AddMajorErrorLog(kHeaderCountErr,errMess.Data());
    if (fAdvancedErrorLog) fAdvLogger->ProcessError(kHeaderCountErr,ddlID,-1,-1,errMess.Data());
  }
  if (fEqPLChipTrailersRead != fExpectedHeaderTrailerCount) {
    TString errMess = Form("Chip trailer count inconsistent %d != %d (expected) for ddl %d",
			   fEqPLChipTrailersRead,fExpectedHeaderTrailerCount,ddlID);
    AliError(errMess.Data());
    fRawReader->AddMajorErrorLog(kHeaderCountErr,errMess.Data());
    if (fAdvancedErrorLog) fAdvLogger->ProcessError(kHeaderCountErr,ddlID,-1,-1,errMess.Data());
  }
}
//__________________________________________________________________________
void AliITSRawStreamSPD::ActivateAdvancedErrorLog(Bool_t activate, AliITSRawStreamSPDErrorLog* advLogger) {
  // Activate the advanced error logging.
  // Logger object has to be created elsewhere and a pointer sent here.
  fAdvancedErrorLog=activate;
  if (activate && advLogger!=NULL) {
    fAdvLogger = advLogger;
  }
  if (fAdvLogger==NULL) {
    fAdvancedErrorLog=kFALSE;
  }
}
//__________________________________________________________________________
const Char_t* AliITSRawStreamSPD::GetErrorName(UInt_t errorCode) {
  // Returns a string for each error code
  if      (errorCode==kTotal)                   return "All Errors";
  else if (errorCode==kHeaderMissingErr)        return "Header Missing";
  else if (errorCode==kTrailerMissingErr)       return "Trailer Missing";
  else if (errorCode==kTrailerWithoutHeaderErr) return "Trailer Unexpected";
  else if (errorCode==kHeaderCountErr)          return "Header Count Wrong";
  else if (errorCode==kTrailerCountErr)         return "Trailer Count Wrong";
  else if (errorCode==kFillUnexpectErr)         return "Fill Unexpected";
  else if (errorCode==kFillMissingErr)          return "Fill Missing";
  else if (errorCode==kWrongFillWordErr)        return "Fill Word Wrong";
  else if (errorCode==kNumberHitsErr)           return "Hit Count Wrong";
  else if (errorCode==kEventCounterErr)         return "Event Counter Error";
  else if (errorCode==kDDLNumberErr)            return "DDL Number Error";
  else if (errorCode==kHSNumberErr)             return "HS Number Error";
  else if (errorCode==kChipAddrErr)             return "Chip Address Error";
  else if (errorCode==kCalHeaderLengthErr)      return "Calib Header Length Error";
  else return "";
}
//__________________________________________________________________________
Bool_t AliITSRawStreamSPD::GetHalfStavePresent(UInt_t hs) {
  // Reads the half stave present status from the block attributes
  Int_t ddlID = fRawReader->GetDDLID();
  if (ddlID==-1) {
    AliWarning("DDL ID = -1. Cannot read block attributes. Return kFALSE.");
    return kFALSE;
  }
  else {
    if (hs>=6) {
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
//__________________________________________________________________________
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
//__________________________________________________________________________
Int_t AliITSRawStreamSPD::GetModuleNumber(UInt_t iDDL, UInt_t iModule) {
  if (iDDL<20 && iModule<12) return fgkDDLModuleMap[iDDL][iModule];
  else return 240;
}
//__________________________________________________________________________
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
//__________________________________________________________________________
Bool_t AliITSRawStreamSPD::OnlineToOffline(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row, UInt_t& module, UInt_t& colM, UInt_t& rowM) {
  // converts online coordinates to offline
  module = GetOfflineModuleFromOnline(eq,hs,chip);
  colM = GetOfflineColFromOnline(eq,hs,chip,col);
  rowM = GetOfflineRowFromOnline(eq,hs,chip,row);
  if (module>=240 || colM>=160 || rowM>=256) return kFALSE;
  else return kTRUE;
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOnlineEqIdFromOffline(UInt_t module) {
  // offline->online (eq)
  for (UInt_t eqId=0; eqId<20; eqId++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      if (GetModuleNumber(eqId,iModule)==(Int_t)module) return eqId;
    }
  }
  return 20; // error
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOnlineHSFromOffline(UInt_t module) {
  // offline->online (hs)
  for (UInt_t eqId=0; eqId<20; eqId++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      if (GetModuleNumber(eqId,iModule)==(Int_t)module) return iModule/2;
    }
  }
  return 6; // error
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOnlineChipFromOffline(UInt_t module, UInt_t colM) {
  // offline->online (chip)
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      if (GetModuleNumber(eq,iModule)==(Int_t)module) {
	return colM/32 + 5*(iModule%2);
      }
    }
  }
  return 10; // error
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOnlineColFromOffline(UInt_t module, UInt_t colM) {
  // offline->online (col)
  if (colM<160) return colM%32;
  else return 32; // error
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOnlineRowFromOffline(UInt_t module, UInt_t rowM) {
  // offline->online (row)
  if (rowM<256) return rowM;
  else return 256; // error
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOfflineModuleFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip) {
  // online->offline (module)
  if (eqId<20 && hs<6 && chip<10) return fgkDDLModuleMap[eqId][hs*2+chip/5];
  else return 240;
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOfflineColFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col) {
  // online->offline (col)
  if (eqId>=20 || hs>=6 || chip>=10 || col>=32) return 160; // error
  UInt_t offset = 32 * (chip % 5);
  return col+offset;
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOfflineRowFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t row) {
  // online->offline (row)
  if (eqId>=20 || hs>=6 || chip>=10 || row>=256) return 256; // error
  return row;
}




