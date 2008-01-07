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
//static Bool_t stDeb = kTRUE;

ClassImp(AliHMPIDRawStream)

//_____________________________________________________________________________
AliHMPIDRawStream::AliHMPIDRawStream(AliRawReader* rawReader) :
  fDDLNumber(-1),
  fRawReader(rawReader),
  fData(NULL),
  fPosition(-1)
{
  // Constructor
  Init();

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
  Init();
}
//_____________________________________________________________________________
AliHMPIDRawStream::~AliHMPIDRawStream()
{
  // destructor
}

//_____________________________________________________________________________
void AliHMPIDRawStream::Init()
{
  // Initalize the container
  // with the pad charges
  Int_t n=0;
  iPos = 0;
  for(Int_t h = 0; h < kNDDL; h++) {  
    for(Int_t i = 0; i < kNRows; i++){
      for(Int_t j = 0; j < kNDILOGICAdd; j++){
	for(Int_t k = 0; k < kNPadAdd; k++){
	  fCharge[h][i][j][k] = -1;
	  fPad[h][i][j][k]=-1;
          fPos[++n] = 0;
	}
      }
    }
  }
  fZeroSup=kTRUE;
  
  for(Int_t l=1; l < kSumErr; l++) fNumOfErr[l]=0;               //reset errors
}//Init()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::Reset()
{
  // reset raw stream params
  // Reinitalize the containers
  Init();
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
    
    Init();
    
    status = ReadHMPIDRawData();
//    if(status) Printf("Event DDL %i successfully decoded!.",fDDLNumber);
// Just for test...
//    for(Int_t i=0;i<fRawReader->GetDataSize()/4;i++) {
//      GetWord();
//    }
//...
    
//    if(stDeb) DumpData(fRawReader->GetDataSize());
//    stDeb=kFALSE;
  }
  return status;
}
    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::ReadHMPIDRawData()
{
  Int_t cntGlob = fRawReader->GetDataSize()/4;
  Int_t cnt = cntGlob;
  UInt_t word32;
  Int_t nwSeg;
  Int_t cntSegment;

  word32 = GetWord(cnt);cnt--;

  
  while (cnt>0) {
    
    nwSeg = (word32 >> kbit8) & 0xfff;
    if(!CheckSegment(word32)) return kFALSE;
    if(!ReadSegment(word32,cntSegment)) return kFALSE;

    if(nwSeg != cntSegment) {Printf("Error in Segment counters: %i different wrt %i",nwSeg,cntSegment);return kFALSE;}
    Printf(" cnt %i cntSegment %i",cnt,cntSegment);
    word32 = GetWord(cntSegment+1,kBwd);
    cnt-=cntSegment+1;
  }
  
  return kTRUE;
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::ReadSegment(UInt_t word32,Int_t &cntSegment)
{
  cntSegment = (word32 >> kbit8) & 0xfff;
  Int_t cnt = cntSegment;
  Int_t cntRow;
  Int_t nwRow;

  word32 = GetWord(cnt,kBwd);
  cntRow  = (word32 >> kbit16) & 0xfff;
  
  while (cnt>0) {

    if(!CheckRowMarker(word32)) return kFALSE;
    if(!ReadRow(word32,nwRow)) return kFALSE;

    if(nwRow != cntRow) {Printf("Error in Row counters: %i different wrt %i",nwRow,cntRow);return kFALSE;}
    word32 = GetWord(cntRow+1);
    cnt -= cntRow+1;
    
  }

  cntSegment -= cnt;
  
  return kTRUE;
    
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::ReadRow(UInt_t word32,Int_t &cntRow)
{
  Int_t cnt;
  Int_t cntDilogic;
  Int_t nwDil;
  
  cntRow  = (word32 >> kbit16) & 0xfff;
  cnt = cntRow;  
  
  word32 = GetWord(cntRow);
  
  while (cnt>0) {
    
    if(!CheckEoE(word32,nwDil)) return kFALSE;
    if(!ReadDilogic(word32,cntDilogic)) return kFALSE;
    
    if(nwDil != cntDilogic) {Printf("Error in Dilogic counters: %i different wrt %i",nwDil,cntDilogic);return kFALSE;}
    cnt -= cntDilogic;
    word32 = GetWord(1,kBwd); // go to next Dilogic bank...
    cnt--;
//    Printf(" cnt %i cntDilogic %i ",cnt,cntDilogic);
  }
  
  cntRow -= cnt;  
  
  return kTRUE;
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::ReadDilogic(UInt_t word32,Int_t &cntDilogic)
{
  cntDilogic = word32 & 0x7f;
  
  Int_t cnt = cntDilogic;

  for(Int_t iDil=0;iDil<cntDilogic;iDil++) {
    UInt_t dilogic = 0, row = 0;
    word32 = GetWord(1,kBwd);
//check on row number      
    cnt--;
    row = (word32 >> kbit22) & 0xf;
    if(!CheckRow(row)) continue;
//check dilogic number     
    dilogic = (word32 >> kbit18) & 0xf;                                              //dilogic info in raw word is between bits: 18...21
    if(!CheckDilogic(dilogic)) continue;
//check pad number
    UInt_t pad = (word32 >> kbit12) & 0x3f;                                          //pad info in raw word is between bits: 12...17
    if(!CheckPad(pad)) continue;
    fCharge[fDDLNumber][row][dilogic][pad] = word32 & 0xfff;
//            Printf(" (word %08X) DDL %i row %i dil %i pad %i ",word32,fDDLNumber,row,dilogic,pad);
  }

  cntDilogic -= cnt;  
  return kTRUE;  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckSegment(UInt_t word)
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
  UInt_t segMarker = (word >> kbit20) & 0xfff;
  if (segMarker != markSegment ) {
    fRawReader->AddMajorErrorLog(kWrongSegErr,Form("Segment marker %0X wrong (expected %0X) at address %i in word %0X!",segMarker,markSegment,fPosition/4,word));
    AliWarning(Form("Segment marker %X wrong (expected %0X)! at address %i in word %0X!",segMarker,markSegment,fPosition/4,word));
    fNumOfErr[kWrongSegErr]++;
    return kFALSE;
  }
  
  UInt_t segAddress = word & 0xff;
  if (segAddress<1 ||segAddress>3) {
    fRawReader->AddMajorErrorLog(kWrongSegErr,Form("Segment address %d not in the valid range [1-3] at address %i in word %0X",segAddress,fPosition/4,word));
    AliWarning(Form("Segment address %d not in the valid range [1-3]",segAddress));
    fNumOfErr[kWrongSegErr]++;
    return kFALSE;
  }
  Printf("Segment Marker found! Number of segment is %i",segAddress);
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckRow(UInt_t row)
{
//check on row number      
  if(row>=1 && row <=kNRows) return kTRUE;
  
  fRawReader->AddMajorErrorLog(kWrongRowErr,Form("row %d",row));
  AliWarning(Form("Wrong row index: %d, expected (1 -> %d)!",row,kNRows));
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
Bool_t AliHMPIDRawStream::CheckEoE(UInt_t word,Int_t &nDil)
{
  if (!((word >> kbit27) & 0x1)) {                                                //check 27th bit in EoE. It must be 1!
    fRawReader->AddMajorErrorLog(kEoEFlagErr);
    AliWarning(Form("Missing end-of-event flag! (%08X)",word));
    fNumOfErr[kEoEFlagErr]++;
    return kFALSE;
  }
  nDil = word & 0x7f;
  if(nDil < 1 || nDil > 48 ) {
    fRawReader->AddMajorErrorLog(kEoESizeErr,Form("EoE size=%d",nDil));
    AliWarning(Form("Wrong end-of-event word-count: %08X",word));
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
Bool_t AliHMPIDRawStream::CheckRowMarker(UInt_t word)
{
  UInt_t nMAXwordsInRow = 0x1EA;
  UInt_t statusControlRow = 0x32a8; // 0x36a8 for zero suppression
//First check on row marker    
  UInt_t rowControlWord = word >> kbit0 & 0xfbff;
  
  if(rowControlWord != statusControlRow) {
    fRawReader->AddMajorErrorLog(kRowMarkerErr);
    AliWarning(Form("Wrong row marker %x expected 0x32a8!",rowControlWord));
    fNumOfErr[kRowMarkerErr]++;
    return kFALSE;
  }
//Second check on row marker    
   UInt_t wordsInRow = word >> kbit16 & 0x0fff;                // Number of words after the row marker
  
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
//  if(fPos[fPosition]!=0) {Printf("Position already stored!!! Value %i at address %i",fPos[fPosition],fPosition); return;}
  iPos++;
  fPos[fPosition] = iPos;
//  if(stDeb)Printf("%i - Actual position %i",iPos,fPosition); 
}
