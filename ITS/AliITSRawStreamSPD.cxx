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
//__________________________________________________________________________
AliITSRawStreamSPD::AliITSRawStreamSPD(AliRawReader* rawReader) :
  AliITSRawStream(rawReader),
  fEventCounter(-1),fChipAddr(0),fHalfStaveNr(0),fCol(0),fRow(0),fCalHeadLen(0),
  fData(0),fOffset(0),fHitCount(0),
  fDataChar1(0),fDataChar2(0),fDataChar3(0),fDataChar4(0),
  fFirstWord(kTRUE),fPrevEventId(0xffffffff),
  fEqPLBytesRead(0),fEqPLChipHeadersRead(0),fEqPLChipTrailersRead(0),
  fHeaderOrTrailerReadLast(kFALSE),fExpectedHeaderTrailerCount(0),
  fFillOutOfSynch(kFALSE),fDDLID(-1),fLastDDLID(-1),fAdvancedErrorLog(kFALSE),
  fAdvLogger(NULL)
{
  // create an object to read ITS SPD raw digits
  fRawReader->Reset();
  fRawReader->Select("ITSSPD");
  // reset calib header words
  for (UInt_t iword=0; iword<kCalHeadLenMax; iword++) {
    fCalHeadWord[iword]=0xffffffff;
  }
  for (UInt_t eq=0; eq<20; eq++) {
    fActiveEq[eq]=kFALSE;
    for (UInt_t hs=0; hs<6; hs++) {
      fActiveHS[eq][hs]=kFALSE;
      for (UInt_t chip=0; chip<10; chip++) {
	fActiveChip[eq][hs][chip]=kFALSE;
	fEventCounterFull[eq][hs][chip] = -1;
      }
    }
  }
  NewEvent();
}
//__________________________________________________________________________
AliITSRawStreamSPD::AliITSRawStreamSPD(const AliITSRawStreamSPD& rstream) :
  AliITSRawStream(rstream.fRawReader),
  fEventCounter(-1),fChipAddr(0),fHalfStaveNr(0),fCol(0),fRow(0),fCalHeadLen(0),
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
    const UInt_t *intPtr = fRawReader->GetEventId();
    if (intPtr!=0) {
      if (*intPtr!=fPrevEventId) { // if new event...
	NewEvent();
	fPrevEventId=*intPtr;
      }
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
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	fFastOrSignal[eq][hs][chip] = kFALSE;
      }
    }
  }
}
//__________________________________________________________________________
Int_t AliITSRawStreamSPD::ReadCalibHeader() {
  // needs to be called in the beginning of each equipment data
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

  // check router error bits:
  UInt_t statusBits = fRawReader->GetStatusBits();
  if ((statusBits >> 5) & 1) { // linkrx/detector fatal error bit
    TString errMess = "LinkRx Error Bit Set";
    AliError(errMess.Data());
    fRawReader->AddMajorErrorLog(kLinkRxDetectorFatalErr,errMess.Data());
    if (fAdvancedErrorLog) fAdvLogger->ProcessError(kLinkRxDetectorFatalErr,ddlID,-1,-1,errMess.Data());
  }
  if ((statusBits >> 12) & 1) { // trigger sequence monitor error bit
    TString errMess = "TSM Trigger Error Bit Set";
    AliError(errMess.Data());
    fRawReader->AddMajorErrorLog(kTSMtriggerErr,errMess.Data());
    if (fAdvancedErrorLog) fAdvLogger->ProcessError(kTSMtriggerErr,ddlID,-1,-1,errMess.Data());
  }
  if (fRawReader->TestBlockAttribute(7)) { // bunch crossing difference error bit
    TString errMess = "High Multiplicity Event Flag Set";
    AliError(errMess.Data());
    fRawReader->AddMajorErrorLog(kHighMultiplicityFlag,errMess.Data());
    if (fAdvancedErrorLog) fAdvLogger->ProcessError(kHighMultiplicityFlag,ddlID,-1,-1,errMess.Data());
  }

  // set eq active and the participating half-staves (for access from outside this class), 
  // set what number of chip headers/trailers to expect (for the moment not used)
  fActiveEq[ddlID]=kTRUE;
  fExpectedHeaderTrailerCount = 0;
  for (UInt_t hs=0; hs<6; hs++) {
    if (!fRawReader->TestBlockAttribute(hs)) {
      fActiveHS[ddlID][hs]=kTRUE;
      fExpectedHeaderTrailerCount+=10;
    }
    else {
      fActiveHS[ddlID][hs]=kFALSE;
    }
  }

  if (ddlID>=0 && ddlID<20) fCalHeadRead[ddlID]=kTRUE;
  if (fRawReader->TestBlockAttribute(6)) { // is the calib header present?
    if (ReadNextInt()) {
      // length of cal header:
      fCalHeadLen = fDataChar1+(fDataChar2<<8)+(fDataChar3<<16)+(fDataChar4<<24);
      if (fCalHeadLen>kCalHeadLenMax) {
	TString errMess = Form("Header length %d > max = %d",fCalHeadLen,kCalHeadLenMax);
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kCalHeaderLengthErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kCalHeaderLengthErr,ddlID,-1,-1,errMess.Data());
	return -1;
      }
      else {
	for (UInt_t iword=0; iword<fCalHeadLen; iword++) {
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
	return fCalHeadLen;
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
	  // The next line is commented, since we do not have the information how many chips to expect
	  // CheckHeaderAndTrailerCount(fLastDDLID);
	}
	if (ReadCalibHeader()>=0) continue; // skip to next word if we found a calib header, otherwise parse this word as regular data (go on with this while loop)
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
      fActiveChip[fDDLID][fHalfStaveNr][fChipAddr]=kTRUE;
      fEventCounterFull[fDDLID][fHalfStaveNr][fChipAddr] = eventCounter;
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
      UShort_t hitCount = fData & 0x0FFF;
      if (hitCount != fHitCount){
	TString errMess = Form("Number of hits %d, while %d expected",fHitCount,hitCount);
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kNumberHitsErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kNumberHitsErr,fDDLID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
      }
      Bool_t errorBit = fData & 0x1000;
      if (errorBit) {
	TString errMess = Form("Trailer error bit set for chip %d,%d,%d",fDDLID,fHalfStaveNr,fChipAddr);
	AliError(errMess.Data());
	fRawReader->AddMajorErrorLog(kTrailerErrorBitErr,errMess.Data());
	if (fAdvancedErrorLog) fAdvLogger->ProcessError(kTrailerErrorBitErr,fDDLID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
      }
      Bool_t fastorBit = fData & 0x2000;
      fFastOrSignal[fDDLID][fHalfStaveNr][fChipAddr] = fastorBit;
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
    // The next line is commented, since we do not have the information how many chips to expect
    // CheckHeaderAndTrailerCount(fDDLID);
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
    if (fAdvancedErrorLog) fAdvLogger->ProcessError(kHeaderCountErr,ddlID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
  }
  if (fEqPLChipTrailersRead != fExpectedHeaderTrailerCount) {
    TString errMess = Form("Chip trailer count inconsistent %d != %d (expected) for ddl %d",
			   fEqPLChipTrailersRead,fExpectedHeaderTrailerCount,ddlID);
    AliError(errMess.Data());
    fRawReader->AddMajorErrorLog(kHeaderCountErr,errMess.Data());
    if (fAdvancedErrorLog) fAdvLogger->ProcessError(kHeaderCountErr,ddlID,fEqPLBytesRead,fEqPLChipHeadersRead,errMess.Data());
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
  else if (errorCode==kAdvEventCounterErr)      return "Event Counter Error (Adv)";
  else if (errorCode==kAdvEventCounterOrderErr) return "Event Counter Jump Error (Adv)";
  else if (errorCode==kTrailerErrorBitErr)      return "Trailer Error Bit Set";
  else if (errorCode==kLinkRxDetectorFatalErr)  return "LinkRx/Detector Fatal Error Bit Set";
  else if (errorCode==kTSMtriggerErr)           return "TSM Trigger Error Bit Set";
  else if (errorCode==kHighMultiplicityFlag)    return "High Multiplicity Event Flag Set";
  else return "";
}
//__________________________________________________________________________
Bool_t AliITSRawStreamSPD::GetFastOrSignal(UInt_t eq, UInt_t hs, UInt_t chip) {
  // returns if there was a fastor signal from this chip
  if (eq>=20 || hs>=6 || chip>=10) {
    TString errMess = Form("eq,hs,chip = %d,%d,%d out of bounds. Return kFALSE.",eq,hs,chip);
    AliError(errMess.Data());
    return kFALSE;
  }
  return fFastOrSignal[eq][hs][chip];
}
//__________________________________________________________________________
Bool_t AliITSRawStreamSPD::IsActiveEq(UInt_t eq) const {
  // returns if the eq is active (seen in data)
  if (eq>=20) {
    TString errMess = Form("eq = %d out of bounds. Return kFALSE.",eq);
    AliError(errMess.Data());
    return kFALSE;
  }
  return fActiveEq[eq];
}
//__________________________________________________________________________
Bool_t AliITSRawStreamSPD::IsActiveHS(UInt_t eq, UInt_t hs) const {
  // returns if the hs is active (info from block attr)
  if (eq>=20 || hs>=6) {
    TString errMess = Form("eq,hs = %d,%d out of bounds. Return kFALSE.",eq,hs);
    AliError(errMess.Data());
    return kFALSE;
  }
  return fActiveHS[eq][hs];
}
//__________________________________________________________________________
Bool_t AliITSRawStreamSPD::IsActiveChip(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // returns if the chip is active (seen in data)
  if (eq>=20 || hs>=6 || chip>=10) {
    TString errMess = Form("eq,hs,chip = %d,%d,%d out of bounds. Return kFALSE.",eq,hs,chip);
    AliError(errMess.Data());
    return kFALSE;
  }
  return fActiveChip[eq][hs][chip];
}
//__________________________________________________________________________
Bool_t AliITSRawStreamSPD::GetHalfStavePresent(UInt_t hs) {
  // Reads the half stave present status from the block attributes
  // This is not really needed anymore (kept for now in case it is still used somewhere).
  // The same information can be reached through the "IsActiveHS" method instead.
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
Bool_t AliITSRawStreamSPD::IsEventCounterFullConsistent() const {
  // checks if the event counter values are consistent within the event
  Short_t reference = -1;
  for (UInt_t eq=0; eq<20; eq++) {
    if (IsActiveEq(eq)) {
      for (UInt_t hs=0; hs<6; hs++) {
	if (IsActiveHS(eq,hs)) {
	  for (UInt_t chip=0; chip<10; chip++) {
	    if (fEventCounterFull[eq][hs][chip]!=-1) {
	      if (reference==-1) reference = fEventCounterFull[eq][hs][chip];
	      if (fEventCounterFull[eq][hs][chip] != reference) return kFALSE;
	    }
	  }
	}
      }
    }
  }
  return kTRUE;
}
//__________________________________________________________________________
Short_t AliITSRawStreamSPD::GetEventCounterFullEq(UInt_t eq) const {
  // if the eq is active; returns the event counter value
  if (eq>=20) {
    TString errMess = Form("eq (%d) out of bounds",eq);
    AliError(errMess.Data());
    return -1;
  }
  if (IsActiveEq(eq)) {
    for (UInt_t hs=0; hs<6; hs++) {
      if (IsActiveHS(eq,hs)) {
	for (UInt_t chip=0; chip<10; chip++) {
	  if (fEventCounterFull[eq][hs][chip]!=-1) {
	    return fEventCounterFull[eq][hs][chip];
	  }
	}
      }
    }
  }
  return -1;
}
//__________________________________________________________________________
Short_t AliITSRawStreamSPD::GetEventCounterFullHS(UInt_t eq, UInt_t hs) const {
  // if the eq,hs is active; returns the event counter value
  if (eq>=20 || hs>=6) {
    TString errMess = Form("eq,hs (%d,%d) out of bounds",eq,hs);
    AliError(errMess.Data());
    return -1;
  }
  if (IsActiveEq(eq)) {
    if (IsActiveHS(eq,hs)) {
      for (UInt_t chip=0; chip<10; chip++) {
	if (fEventCounterFull[eq][hs][chip]!=-1) {
	  return fEventCounterFull[eq][hs][chip];
	}
      }
    }
  }
  return -1;
}
//__________________________________________________________________________
Short_t AliITSRawStreamSPD::GetEventCounterFullChip(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // if the eq,hs,chip is active; returns the event counter value
  if (eq>=20 || hs>=6 || chip>=10) {
    TString errMess = Form("eq,hs,chip (%d,%d,%d) out of bounds",eq,hs,chip);
    AliError(errMess.Data());
    return -1;
  }
  if (IsActiveEq(eq)) {
    if (IsActiveHS(eq,hs)) {
      if (IsActiveChip(eq,hs,chip)) {
	return fEventCounterFull[eq][hs][chip];
      }
    }
  }
  return -1;
}
//__________________________________________________________________________
Int_t AliITSRawStreamSPD::GetHword(UInt_t index) {
  if (index<kCalHeadLenMax) return fCalHeadWord[index];
  else return 0;
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
UInt_t AliITSRawStreamSPD::GetFOHnumDacs() const {
  return (fCalHeadLen-37)/2;
}
UInt_t AliITSRawStreamSPD::GetFOHdacIndex(UInt_t index) const {
  if (index>=GetFOHnumDacs()) {
    TString errMess = Form("Only %d DACs in this run, returning 0",GetFOHnumDacs());
    AliError(errMess.Data());
    return 0;
  }
  return fCalHeadWord[7+index*2];
}
UInt_t AliITSRawStreamSPD::GetFOHdacValue(UInt_t index) const {
  if (index>=GetFOHnumDacs()) {
    TString errMess = Form("Only %d DACs in this run, returning 0",GetFOHnumDacs());
    AliError(errMess.Data());
    return 0;
  }
  return fCalHeadWord[7+1+index*2];
}
UInt_t AliITSRawStreamSPD::GetFOHchipCount(UInt_t hs, UInt_t chip) const {
  if (hs<6 && chip<10) {
    if (chip%2==0) {
      return ((fCalHeadWord[7 + GetFOHnumDacs()*2 + (hs*10 + chip)/2] >> 16) & 0x0000ffff);
    }
    else {
      return (fCalHeadWord[7 + GetFOHnumDacs()*2 + (hs*10 + chip)/2] & 0x0000ffff);
    }
  }
  else return 0;
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
	if (module<80) {
	  if (eq<10) { // side A
	    return (159-colM)/32 + 5*(iModule%2);
	  }
	  else { // side C
	    return colM/32 + 5*(iModule%2);
	  }
	}
	else if (module<240) {
	  if (eq<10) { // side A
	    return colM/32 + 5*(iModule%2);
	  }
	  else { // side C
	    return (159-colM)/32 + 5*(iModule%2);
	  }
	}
      }
    }
  }
  return 10; // error
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOnlineColFromOffline(UInt_t module, UInt_t colM) {
  // offline->online (col)
  if (module<80) { // inner layer
    return colM%32;
  }
  else if (module<240) { // outer layer
    return colM%32;
  }
  return 32; // error
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOnlineRowFromOffline(UInt_t module, UInt_t rowM) {
  // offline->online (row)
  if (module<80) { // inner layer
    return (255-rowM);
  }
  else if (module<240) { // outer layer
    return (255-rowM);
  }
  return 256; // error
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOfflineModuleFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip) {
  // online->offline (module)
  if (eqId<20 && hs<6 && chip<10) return fgkDDLModuleMap[eqId][hs*2+chip/5];
  else return 240;
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOfflineChipKeyFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip) {
  // online->offline (chip key: 0-1199)
  if (eqId<20 && hs<6 && chip<10) {
    UInt_t module = GetOfflineModuleFromOnline(eqId,hs,chip);
    UInt_t chipInModule = ( chip>4 ? chip-5 : chip ); 
    if(eqId>9) chipInModule = 4 - chipInModule;  // side C only
    return (module*5 + chipInModule);
  } else return 1200;
}
//__________________________________________________________________________
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
      return (col + offset); // outer layer, side A
    }
    else {
      return 159 - (31-col + offset); // outer layer, side C
    }
  }
}
//__________________________________________________________________________
UInt_t AliITSRawStreamSPD::GetOfflineRowFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t row) {
  // online->offline (row)
  if (eqId>=20 || hs>=6 || chip>=10 || row>=256) return 256; // error
  return 255-row;
}

