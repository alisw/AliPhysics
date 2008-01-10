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

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading the HMPID raw data
/// The format of the raw data corresponds to the one
/// which was documented by Paolo Martinengo.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliHMPIDRawStream.h"
#include "AliRawReader.h"
#include "AliLog.h"

Int_t fPos[170000];
Int_t iPos = 0;
static Bool_t stDeb = kTRUE;

ClassImp(AliHMPIDRawStream)

//_____________________________________________________________________________
AliHMPIDRawStream::AliHMPIDRawStream(AliRawReader* rawReader) :
  fDDLNumber(-1),
  fRawReader(rawReader),
  fData(NULL),
  fPosition(-1)
{
  // Constructor
  fNPads = 0;
  fCharge = 0x0;
  fPad =0x0;

  for(Int_t l=1; l < kSumErr; l++) fNumOfErr[l]=0;               //reset errors
    
  fRawReader->Reset();
  fRawReader->Select("HMPID");
}
//-----------------------------------------------------------------------------
AliHMPIDRawStream::AliHMPIDRawStream() :
  fDDLNumber(-1),
  fRawReader(0x0),
  fData(NULL),
  fPosition(-1)
{
  // Constructor
  for(Int_t l=1; l < kSumErr; l++) fNumOfErr[l]=0;               //reset errors
}
//_____________________________________________________________________________
AliHMPIDRawStream::~AliHMPIDRawStream()
{
  // destructor
    DelVars();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::Reset()
{
  // reset raw stream params
  // Reinitalize the containers
  fDDLNumber = -1;
  fPosition = -1;
  fData = NULL;
  if (fRawReader) fRawReader->Reset();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::Next()
{
  // read next DDL raw data from the HMPID raw data stream
  // return kFALSE in case of error or no data left
  AliDebug(1,"Start.");
  do {
    if (!fRawReader->ReadNextData(fData)) return kFALSE;
  } while (fRawReader->GetDataSize() == 0);
  
     
  /*
  Event type is selected as in $ALICE_ROOT/RAW/event.h  
  #define START_OF_RUN                    ((eventTypeType) 1)
  #define END_OF_RUN                      ((eventTypeType) 2)
  #define START_OF_RUN_FILES              ((eventTypeType) 3)
  #define END_OF_RUN_FILES                ((eventTypeType) 4)
  #define START_OF_BURST                  ((eventTypeType) 5)
  #define END_OF_BURST                    ((eventTypeType) 6)
  #define PHYSICS_EVENT                   ((eventTypeType) 7) <<---------------  
  #define CALIBRATION_EVENT               ((eventTypeType) 8)
  #define EVENT_FORMAT_ERROR              ((eventTypeType) 9)
  #define START_OF_DATA                   ((eventTypeType)10)
  #define END_OF_DATA                     ((eventTypeType)11)
  #define SYSTEM_SOFTWARE_TRIGGER_EVENT   ((eventTypeType)12)
  #define DETECTOR_SOFTWARE_TRIGGER_EVENT ((eventTypeType)13)
  #define EVENT_TYPE_MIN                  1
  #define EVENT_TYPE_MAX                  13 
  */

  fPosition = 0;
  Bool_t status;
  
  if(fRawReader->GetType() == 7)  {                             //New: Select Physics events, Old: Raw data size is not 0 and not 47148 (pedestal)
    fDDLNumber = fRawReader->GetDDLID();
    Printf("DDL %i started to be decoded!.",fDDLNumber);
    InitVars(fRawReader->GetDataSize()/4);
    status = ReadHMPIDRawData();
    if(status) Printf("Event DDL %i successfully decoded!.",fDDLNumber);
    else Printf("Event DDL %i ERROR in decoding!.",fDDLNumber);
    if(stDeb) DumpData(fRawReader->GetDataSize());
//    stDeb=kFALSE;
  }
//  return status;
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::InitVars(Int_t n)
{
  fNPads = 0;
  fCharge = new Int_t[n];
  fPad = new Int_t[n];
  //for debug purpose
  for(Int_t i=0; i < 170000; i++) fPos[i]=0;                     //reset debug
  iPos = 0;
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::DelVars()
{
  fNPads = 0;
  delete fCharge;
  delete fPad;
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::ReadHMPIDRawData()
{
  Int_t cntGlob = fRawReader->GetDataSize()/4;
  Int_t cnt = cntGlob;
  Int_t nwSeg;
  Int_t cntSegment;

  fWord = GetWord(cnt);cnt--;

  
  while (cnt>0) {
    
    nwSeg = (fWord >> kbit8) & 0xfff;
    if(!CheckSegment()) return kFALSE;
    if(!ReadSegment(cntSegment)) return kFALSE;

    if(nwSeg != cntSegment) {Printf("Error in Segment counters: %i different wrt %i",nwSeg,cntSegment);return kFALSE;}
//    Printf(" cnt %i cntSegment %i",cnt,cntSegment);
    fWord = GetWord(cntSegment+1,kBwd);
    cnt-=cntSegment+1;
  }
  
  return kTRUE;
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::ReadSegment(Int_t &cntSegment)
{
  cntSegment = (fWord >> kbit8) & 0xfff;
  Int_t cnt = cntSegment;
  Int_t cntRow;
  Int_t nwRow;

  fWord = GetWord(cnt,kBwd);
  
  while (cnt>0) {

    cntRow  = (fWord >> kbit16) & 0xfff;
    if(!CheckRowMarker()) return kFALSE;
    if(!ReadRow(nwRow)) return kFALSE;

    if(nwRow != cntRow) {Printf("Error in Row counters: %i different wrt %i",nwRow,cntRow);return kFALSE;}
    fWord = GetWord(cntRow+1);
    cnt -= cntRow+1;
    
  }

  cntSegment -= cnt;
  
  return kTRUE;
    
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::ReadRow(Int_t &cntRow)
{
  Int_t cnt;
  Int_t cntDilogic;
  Int_t nwDil;
  
  cntRow  = (fWord >> kbit16) & 0xfff;
  cnt = cntRow;  
  
  fWord = GetWord(cntRow);
  
  while (cnt>0) {
    
    if(!CheckEoE(nwDil)) return kFALSE;
    if(!ReadDilogic(cntDilogic)) return kFALSE;
    
    if(nwDil != cntDilogic) {Printf("Error in Dilogic counters: %i different wrt %i",nwDil,cntDilogic);return kFALSE;}
    cnt -= cntDilogic;
    fWord = GetWord(1,kBwd); // go to next Dilogic bank...
    cnt--;
//    Printf(" cnt %i cntDilogic %i ",cnt,cntDilogic);
  }
  
  cntRow -= cnt;  
  
  return kTRUE;
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::ReadDilogic(Int_t &cntDilogic)
{
  cntDilogic = fWord & 0x7f;
  
  Int_t cnt = cntDilogic;
  
//  Printf(" cnt DILOGIC %i at %i word %08X",cnt,fPosition,fWord);

  for(Int_t iDil=0;iDil<cntDilogic;iDil++) {
    UInt_t dilogic = 0, row = 0;
    fWord = GetWord(1,kBwd);
//check on row number      
    cnt--;
    row = (fWord >> kbit22) & 0x1f;
    if(!CheckRow(row)) continue;
//check dilogic number     
    dilogic = (fWord >> kbit18) & 0xf;                                              //dilogic info in raw word is between bits: 18...21
    if(!CheckDilogic(dilogic)) continue;
//check pad number
    UInt_t pad = (fWord >> kbit12) & 0x3f;                                          //pad info in raw word is between bits: 12...17
    if(!CheckPad(pad)) continue;
    Int_t charge = fWord & 0xfff;
    fPad[fNPads] = GetPad(fDDLNumber,row,dilogic,pad);
    fCharge[fNPads] = charge;
    fNPads++;
  }

  cntDilogic -= cnt;  
  return kTRUE;  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckSegment()
{
  UInt_t markSegment = 0xAB0;
  /*
  if (iRow%8 == 0) {
    UInt_t segWord = GetWord();          
    if ((segWord >> 20) != statusSegWord) {
      fRawReader->AddMajorErrorLog(kBadSegWordErr);
      AliWarning(Form("Wrong segment word signature: %x, expected 0xab0!",(segWord >> 20)));
      fNumOfErr[kBadSegWordErr]++;
      return kTRUE;
    }
*/
  UInt_t segMarker = (fWord >> kbit20) & 0xfff;
  if (segMarker != markSegment ) {
    fRawReader->AddMajorErrorLog(kWrongSegErr,Form("Segment marker %0X wrong (expected %0X) at %i in word %0X!",segMarker,markSegment,fPosition,fWord));
    AliWarning(Form("Segment marker %X wrong (expected %0X)! at %i in word %0X!",segMarker,markSegment,fPosition,fWord));
    fNumOfErr[kWrongSegErr]++;
    return kFALSE;
  }
  
  UInt_t segAddress = fWord & 0xff;
  if (segAddress<1 ||segAddress>3) {
    fRawReader->AddMajorErrorLog(kWrongSegErr,Form("Segment address %d not in the valid range [1-3] at %i in word %0X",segAddress,fPosition,fWord));
    AliWarning(Form("Segment address %d not in the valid range [1-3]",segAddress));
    fNumOfErr[kWrongSegErr]++;
    return kFALSE;
  }
//  Printf("Segment Marker found at %i! Number of segment is %i",fPosition,segAddress);
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckRow(UInt_t row)
{
//check on row number      
//  Printf("ROW %i word %0X",row,fWord);
  if(row>=1 && row <=kNRows) return kTRUE;
  
  fRawReader->AddMajorErrorLog(kWrongRowErr,Form("row %d",row));
  AliWarning(Form("Wrong row index: %d, expected (1 -> %d) word %0X at %i...",row,kNRows,fWord,fPosition));
  fNumOfErr[kWrongRowErr]++;
  return kFALSE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckDilogic(UInt_t dilogic)
{
//check dilogic number     
  if (dilogic>= 1 && dilogic <=kNDILOGICAdd) return kTRUE;

  fRawReader->AddMajorErrorLog(kWrongDilogicErr,Form("dil %d",dilogic));
  AliWarning(Form("Wrong DILOGIC index: %d, expected (1 -> %d)!",dilogic,kNDILOGICAdd));
  fNumOfErr[kWrongDilogicErr]++;
  //dilogic = iDILOGIC;
  return kFALSE;
}      
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckPad(UInt_t pad)
{
//check pad number     
  if (pad < kNPadAdd) return kTRUE;
  
  fRawReader->AddMajorErrorLog(kWrongPadErr,Form("pad %d",pad));
  AliWarning(Form("Wrong pad index: %d, expected (0 -> %d)!",pad,kNPadAdd));
  fNumOfErr[kWrongPadErr]++;
  return kFALSE;
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckEoE(Int_t &nDil)
{
  if (!((fWord >> kbit27) & 0x1)) {                                                //check 27th bit in EoE. It must be 1!
    fRawReader->AddMajorErrorLog(kEoEFlagErr);
    AliWarning(Form("Missing end-of-event flag! (%08X) at %i",fWord,fPosition));
    fNumOfErr[kEoEFlagErr]++;
    return kFALSE;
  }
  nDil = fWord & 0x7f;
  if(nDil < 0 || nDil > 48 ) {
    fRawReader->AddMajorErrorLog(kEoESizeErr,Form("EoE size=%d",nDil));
    AliWarning(Form("Wrong end-of-event word-count: %08X",fWord));
    fNumOfErr[kEoESizeErr]++;
    return kFALSE;
  }
//  UInt_t da = (eOfEvent >> 18) & 0xf;
//  if (cntData!=0 && da != dilogic) {
//    fRawReader->AddMajorErrorLog(kEoEDILOGICErr,Form("eoe dil %d != %d",da,dilogic));
//    AliWarning(Form("Wrong DILOGIC address found in end-of-event: %d, expected %d!",da,dilogic));
//    fNumOfErr[kEoEDILOGICErr]++;
//    return kFALSE;  AliQAChecker::Instance()->Run(AliQA::kHMPID, task, obj) ;  

//  }
//  UInt_t ca = (eOfEvent >> 22) & 0x1f;
//  if (cntData!=0 &&  ca != row) {
//    fRawReader->AddMajorErrorLog(kEoERowErr,Form("eoe row %d != %d",ca,row));
//    AliWarning(Form("Wrong row index found in end-of-event: %d, expected %d!",ca,row));
//    fNumOfErr[kEoERowErr]++;
//    return kFALSE;
//  }
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckRowMarker()
{
  UInt_t nMAXwordsInRow = 0x1EA;
  UInt_t statusControlRow = 0x32a8; // 0x36a8 for zero suppression
//First check on row marker    
  UInt_t rowControlWord = fWord >> kbit0 & 0xfbff;
  
  if(rowControlWord != statusControlRow) {
    fRawReader->AddMajorErrorLog(kRowMarkerErr);
    AliWarning(Form("Wrong row marker %x expected 0x32a8!",rowControlWord));
    fNumOfErr[kRowMarkerErr]++;
    return kFALSE;
  }
//Second check on row marker    
   UInt_t wordsInRow = fWord >> kbit16 & 0x0fff;                // Number of words after the row marker
  
  if (wordsInRow > nMAXwordsInRow) {
    fRawReader->AddMajorErrorLog(kRowMarkerSizeErr);
    AliWarning(Form(" FATAL: Number of words %x in a row exceeds the expected value: 0x1EA !",wordsInRow));
    fNumOfErr[kRowMarkerSizeErr]++;
    return kFALSE;
  }
  
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
UInt_t AliHMPIDRawStream::GetWord(Int_t n,EDirection dir)
{
  // This method returns the n-th 32 bit word
  // inside the raw data payload.
  // The method is supposed to be endian (platform)
  // independent.
  if (!fData || fPosition < 0) AliFatal("Raw data payload buffer is not yet initialized !");
  
  UInt_t word = 0;
  
  if(dir==kBwd) n = -n; 
  fPosition+=4*n-4;
  
  StorePosition();
  
  word |= fData[fPosition++];
  word |= fData[fPosition++] << 8;
  word |= fData[fPosition++] << 16;
  word |= fData[fPosition++] << 24;

  return word;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::DumpData(Int_t nw)
{
   for(Int_t i=0;i<nw;i+=4) {
     if(!(i%16)) printf(" \n %8i) ",i);
     printf("%02X%02X%02X%02X [ %06i ] ",fData[i+3],fData[i+2],fData[i+1],fData[i+0],fPos[i]);
   }
   Printf(" \n -----end of dump ----------- ");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::StorePosition()
{
  if(fPos[fPosition]!=0) {
//    Printf("Position already stored!!! Value %i at address %i",fPos[fPosition],fPosition); 
    return;
  }
  iPos++;
  fPos[fPosition] = iPos;
//  if(stDeb)Printf("%i - Actual position %i",iPos,fPosition); 
}
