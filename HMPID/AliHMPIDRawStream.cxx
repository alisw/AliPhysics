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

ClassImp(AliHMPIDRawStream)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRawStream::AliHMPIDRawStream(AliRawReader* rawReader) :
  fNPads(0),
  fCharge(0x0),
  fPad(0x0),
  fDDLNumber(-1),
  fnDDLInStream(0x0),
  fnDDLOutStream(0x0),
  fLDCNumber( 0),
  fTimeStamp( 0),
  fRawReader(rawReader),
  fData(0x0),
  fNumOfErr(0x0),
  fPosition(-1),
  fWord(0),
  fZeroSup(kTRUE),
  fPos(0x0),
  fiPos(0),
  fTurbo(kFALSE),
  fRawDataSize(0)
{
  //
  // Constructor
  //
  fNumOfErr = new Int_t*[kNDDL];                                 // Store the numner of errors for a given error type and a given DD
  for(Int_t i=0;i<kNDDL;i++) {
    fNumOfErr[i] = new Int_t [kSumErr];
  }
  
  fnDDLInStream=new Int_t[kNDDL];
  fnDDLOutStream=new Int_t[kNDDL];
  for(Int_t iddl=0;iddl<kNDDL;iddl++) { fnDDLInStream[iddl]=-1;fnDDLOutStream[iddl]=-1;}
 
  for(Int_t iddl=0;iddl<kNDDL;iddl++) 
    for(Int_t ierr=0; ierr < kSumErr; ierr++) fNumOfErr[iddl][ierr]=0;               //reset errors  
  fRawReader->Reset();
  fRawReader->Select("HMPID");
  
  
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRawStream::AliHMPIDRawStream() :
  fNPads(0),
  fCharge(0x0),
  fPad(0x0),
  fDDLNumber(-1),
  fnDDLInStream(0x0),
  fnDDLOutStream(0x0),
  fLDCNumber( 0),
  fTimeStamp( 0),
  fRawReader(0x0),
  fData(0x0),
  fNumOfErr(0x0),  
  fPosition(-1),
  fWord(0),
  fZeroSup(kTRUE),
  fPos(0x0),
  fiPos(0),
  fTurbo(kFALSE) ,
  fRawDataSize(0)
{
  //
  // Constructor
  //
  fNumOfErr = new Int_t*[kNDDL];                                 // Store the numner of errors for a given error type and a given DD
  for(Int_t i=0;i<kNDDL;i++) {
    fNumOfErr[i] = new Int_t [kSumErr];
  }
  fnDDLInStream=new Int_t[kNDDL];
  fnDDLOutStream=new Int_t[kNDDL];
  for(Int_t iddl=0;iddl<kNDDL;iddl++) { fnDDLInStream[iddl]=-1;fnDDLOutStream[iddl]=-1;}

  
  for(Int_t iddl=0;iddl<kNDDL;iddl++) 
    for(Int_t ierr=0; ierr < kSumErr; ierr++) fNumOfErr[iddl][ierr]=0;               //reset errors  

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRawStream::~AliHMPIDRawStream()
{
  //
  // destructor
  //
  DelVars();
  
  
  fDDLNumber=0;
  fLDCNumber=0;
  fTimeStamp=0;
  fPosition=0;
  fWord=0;
  fZeroSup=0;
  fTurbo=0;
  fRawDataSize=0;
  for(Int_t i=0;i<kNDDL;i++) delete [] fNumOfErr[i]; 
  delete [] fNumOfErr; 

  delete [] fnDDLInStream;
  delete [] fnDDLOutStream;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::Reset()
{
  // reset raw stream params
  // Reinitalize the containers
  fDDLNumber = -1;
  fLDCNumber =  0;
  fTimeStamp =  0;
  fPosition = -1;
  fData = NULL;
  if (fRawReader) fRawReader->Reset();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::Turbo()
{
  
  Int_t row,dilogic;UInt_t pad;
  Int_t cntGlob = fRawReader->GetDataSize()/4;
  fPosition=0;
  fNPads=0;
//  Int_t gw=0;
  for(Int_t i=1;i<cntGlob;i++) {
    if(!GetWord(1)) return kFALSE;
    if (((fWord >> kbit27) & 1)) continue;
    UInt_t statusControlRow = 0x32a8; 
    UInt_t rowControlWord = fWord >> kbit0 & 0xfbff;
    if(rowControlWord == statusControlRow) continue;

    row = (fWord >> kbit22) & 0x1f;
    dilogic = (fWord >> kbit18) & 0xf;                                              //dilogic info in raw word is between bits: 18...21
    
    pad = (fWord >> kbit12) & 0x3f;                                          //pad info in raw word is between bits: 12...17
    if(!CheckPad(pad)) continue;
    Int_t charge = fWord & 0xfff;
    if(GetPad(fDDLNumber,row,dilogic,pad)<0) continue;
    fPad[fNPads] = GetPad(fDDLNumber,row,dilogic,pad);
    fCharge[fNPads] = charge; 
    fNPads++;
    if(charge==0) fNumOfErr[fDDLNumber][kPedQZero]++;
  }//word loop
  //Printf("Size: %i  DDL %i row %i dilogic %i pad %i fPos %i fNPads: %i Charge: %d Word %4.4x GoodW: %i",cntGlob,fDDLNumber,row,dilogic,pad,fPosition,fNPads,fCharge[fNPads-1],fWord,gw++);
  return kTRUE;
}//Turbo()
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
  Bool_t status=kFALSE;
  fRawDataSize=0;        
  fDDLNumber = fRawReader->GetDDLID();
  if(fRawReader->GetType() == 7 || fRawReader->GetType() == 8 )  {           //New: Select Physics events, Old: Raw data size is not 0 and not 47148 (pedestal)
    fnDDLInStream[fDDLNumber]=1; fnDDLOutStream[fDDLNumber]=0;
    
    fLDCNumber = fRawReader->GetLDCId();
    fTimeStamp = fRawReader->GetTimestamp();
    
    AliDebug(1,Form("DDL %i started to be decoded!",fDDLNumber));
    fRawDataSize=fRawReader->GetDataSize()/4;
    DelVars();                                         //We have to delete the variables initialized in the InitVars before recall IntiVars!!!
    InitVars(fRawDataSize);                            //To read the charge and pads we cannot delete before the status return
    
    if(fTurbo==kTRUE) status=Turbo();
    else status = ReadHMPIDRawData();
   
    if(status) AliDebug(1,Form("Event DDL %i successfully decoded!.",fDDLNumber));
    else AliDebug(1,Form("Event DDL %i ERROR in decoding!.",fDDLNumber));
//      DumpData(fRawReader->GetDataSize());
  
   
    }
    if(status==kTRUE) {fnDDLOutStream[fDDLNumber]++; }//Printf("fnDDLOutStream[%d]=%d",fDDLNumber,fnDDLOutStream[fDDLNumber]); } //Count the number of events when the DDL was succesfully decoded
   
//    return status;  // temporary solution...
   return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::InitVars(Int_t n)
{
  //
  //
  //
  fNPads = 0;
  fCharge = new Int_t[n]; 
  fPad = new Int_t[n];
  //for debug purpose
  fPos = new Int_t[4*n+4];                     //reset debug
  fiPos = 0;
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::DelVars()
{
  //
  //
  // 
  //Clean the initvars!!!!!!!!
  fNPads = 0; 
  if (fCharge)     { delete [] fCharge;    fCharge = 0x0; }
  if (fPad)        { delete [] fPad;       fPad = 0x0;       }   
  if (fPos)        { delete [] fPos;       fPos = 0x0;       }     
 
  fiPos=0;
  
     
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::ReadHMPIDRawData()
{
  //Here the loop on the decoding the raw bank 
  //for one ddl starts.
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK
  Int_t cntGlob = fRawReader->GetDataSize()/4;
  if(cntGlob==0) {fNumOfErr[fDDLNumber][kRawDataSizeErr]++; return kFALSE; }
  
  Int_t cnt = cntGlob;
  Int_t nwSeg;
  Int_t cntSegment;

  if(!GetWord(cnt)) return kFALSE;
  cnt--;

  
  while (cnt>0) {
    
    nwSeg = (fWord >> kbit8) & 0xfff;
    if(!CheckSegment()) return kFALSE;
    if(!ReadSegment(cntSegment)) return kFALSE;

    if(nwSeg != cntSegment) {AliDebug(1,Form("Error in Segment counters: %i different wrt %i",nwSeg,cntSegment)); return kFALSE;}
    if(!GetWord(cntSegment+1,kBwd)) return kFALSE;
    cnt-=cntSegment+1;
  }
  
  return kTRUE;
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::ReadSegment(Int_t &cntSegment)
{
  //Read the segment
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK
  cntSegment = (fWord >> kbit8) & 0xfff;
  Int_t cnt = cntSegment;
  Int_t cntRow;
  Int_t nwRow;

  if(!GetWord(cnt,kBwd)) return kFALSE;
  
  while (cnt>0) {

    cntRow  = (fWord >> kbit16) & 0xfff;
    if(!CheckRowMarker()) return kFALSE;
    if(!ReadRow(nwRow)) return kFALSE;

    if(nwRow != cntRow) {AliDebug(1,Form("Error in Row counters: %i different wrt %i",nwRow,cntRow)); return kFALSE;}
    if(!GetWord(cntRow+1)) return kFALSE;
    cnt -= cntRow+1;
    
  }

  cntSegment -= cnt;
  
  return kTRUE;
    
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::ReadRow(Int_t &cntRow)
{
  // Read the row
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK

  Int_t cnt;
  Int_t cntDilogic;
  Int_t nwDil;
  
  cntRow  = (fWord >> kbit16) & 0xfff;
  cnt = cntRow;  
  
  if(!GetWord(cntRow)) return kFALSE;
  
  while (cnt>0) {
    
    if(!CheckEoE(nwDil)) return kFALSE;
    if(!ReadDilogic(cntDilogic)) return kFALSE;
    
    if(nwDil != cntDilogic) {AliDebug(1,Form("Error in Dilogic counters: %i different wrt %i",nwDil,cntDilogic));return kFALSE;}
    cnt -= cntDilogic;
    if(!GetWord(1,kBwd)) return kFALSE; // go to next Dilogic bank...
    cnt--;
//    Printf(" cnt %i cntDilogic %i ",cnt,cntDilogic);
  }
  
  cntRow -= cnt;  
  
  return kTRUE;
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::ReadDilogic(Int_t &cntDilogic)
{
  // Read the dilogic bank
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK

  cntDilogic = fWord & 0x7f;
  
  Int_t cnt = cntDilogic;
  
//  Printf(" cnt DILOGIC %i at %i word %08X",cnt,fPosition,fWord);

  for(Int_t iDil=0;iDil<cntDilogic;iDil++) {
    UInt_t dilogic = 0, row = 0;
    if(!GetWord(1,kBwd)) return kFALSE;
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
    if(GetPad(fDDLNumber,row,dilogic,pad)<0) continue;
    fPad[fNPads] = GetPad(fDDLNumber,row,dilogic,pad);
    fCharge[fNPads] = charge; 
    fNPads++;
    
    if(charge==0) 
    {
      AliDebug(1,Form("If PEDESTAL run -> WARNING: ZERO charge is read from DDL: %d row: %d dil: %d pad: %d",fDDLNumber,row,dilogic,pad));
      fNumOfErr[fDDLNumber][kPedQZero]++;
    }
    
  }//iDil

  cntDilogic -= cnt;  
  return kTRUE;  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckSegment()
{
  // Check the segment marker
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK

  UInt_t markSegment = 0xAB0;
  /*
  if (iRow%8 == 0) {
    UInt_t segWord = GetWord();          
    if ((segWord >> 20) != statusSegWord) {
      fRawReader->AddMajorErrorLog(kBadSegWordErr);
      AliDebug(1,Form("Wrong segment word signature: %x, expected 0xab0!",(segWord >> 20)));
      fNumOfErr[kBadSegWordErr]++;
      return kTRUE;
    }
*/
  UInt_t segMarker = (fWord >> kbit20) & 0xfff;
  if (segMarker != markSegment ) {
    //fRawReader->AddMajorErrorLog(kWrongSegErr,Form("Segment marker %0X wrong (expected %0X) at %i in word %0X!",segMarker,markSegment,fPosition,fWord));
    AliDebug(1,Form("Segment marker %X wrong (expected %0X)! at %i in word %0X!",segMarker,markSegment,fPosition,fWord));
    fNumOfErr[fDDLNumber][kWrongSegErr]++;
    return kFALSE;
  }
  
  UInt_t segAddress = fWord & 0xff;
  if (segAddress<1 ||segAddress>3) {
    //fRawReader->AddMajorErrorLog(kWrongSegErr,Form("Segment address %d not in the valid range [1-3] at %i in word %0X",segAddress,fPosition,fWord));
    AliDebug(1,Form("Segment address %d not in the valid range [1-3]",segAddress));
    fNumOfErr[fDDLNumber][kWrongSegErr]++;
    return kFALSE;
  }
//  Printf("Segment Marker found at %i! Number of segment is %i",fPosition,segAddress);
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckRow(UInt_t row)
{
  //check on row number      
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK
//  Printf("ROW %i word %0X",row,fWord);
  if(row>=1 && row <=kNRows) return kTRUE;
  
  //fRawReader->AddMajorErrorLog(kWrongRowErr,Form("row %d",row));
  AliDebug(1,Form("Wrong row index: %d, expected (1 -> %d) word %0X at %i...",row,kNRows,fWord,fPosition));
  fNumOfErr[fDDLNumber][kWrongRowErr]++;
  return kFALSE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckDilogic(UInt_t dilogic)
{
//check dilogic number     
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK
  if (dilogic>= 1 && dilogic <=kNDILOGICAdd) return kTRUE;

  //fRawReader->AddMajorErrorLog(kWrongDilogicErr,Form("dil %d",dilogic));
  AliDebug(1,Form("Wrong DILOGIC index: %d, expected (1 -> %d)!",dilogic,kNDILOGICAdd));
  fNumOfErr[fDDLNumber][kWrongDilogicErr]++;
  //dilogic = iDILOGIC;
  return kFALSE;
}      
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckPad(UInt_t pad)
{
//check pad number     
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK
  if (pad < kNPadAdd) return kTRUE;
  
  //fRawReader->AddMajorErrorLog(kWrongPadErr,Form("pad %d",pad));
  AliDebug(1,Form("Wrong pad index: %d, expected (0 -> %d)!",pad,kNPadAdd));
  fNumOfErr[fDDLNumber][kWrongPadErr]++;
  return kFALSE;
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckEoE(Int_t &nDil)
{
  //check the End of Event
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK
  if (!((fWord >> kbit27) & 0x1)) {                                                //check 27th bit in EoE. It must be 1!
    //fRawReader->AddMajorErrorLog(kEoEFlagErr);
    AliDebug(1,Form("Missing end-of-event flag! (%08X) at %i",fWord,fPosition));
    fNumOfErr[fDDLNumber][kEoEFlagErr]++;
    return kFALSE;
  }
  nDil = fWord & 0x7f;                                               //nDil=EoE word count
  if(nDil < 0 || nDil > 48 ) { 

    //fRawReader->AddMajorErrorLog(kEoESizeErr,Form("EoE size=%d",nDil));
    AliDebug(1,Form("Wrong end-of-event word-count: %08X",fWord));
    fNumOfErr[fDDLNumber][kEoESizeErr]++;
    return kFALSE;
  }
//  UInt_t da = (eOfEvent >> 18) & 0xf;
//  if (cntData!=0 && da != dilogic) {
//    fRawReader->AddMajorErrorLog(kEoEDILOGICErr,Form("eoe dil %d != %d",da,dilogic));
//    AliDebug(1,Form("Wrong DILOGIC address found in end-of-event: %d, expected %d!",da,dilogic));
//    fNumOfErr[kEoEDILOGICErr]++;
//    return kFALSE;  AliQAChecker::Instance()->Run(AliQAv1::kHMPID, task, obj) ;  

//  }
//  UInt_t ca = (eOfEvent >> 22) & 0x1f;
//  if (cntData!=0 &&  ca != row) {
//    fRawReader->AddMajorErrorLog(kEoERowErr,Form("eoe row %d != %d",ca,row));
//    AliDebug(1,Form("Wrong row index found in end-of-event: %d, expected %d!",ca,row));
//    fNumOfErr[kEoERowErr]++;
//    return kFALSE;
//  }
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::CheckRowMarker()
{
  //check the row marker
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK
  UInt_t nMAXwordsInRow = 0x1EA;
  UInt_t statusControlRow = 0x32a8; // 0x36a8 for zero suppression
//First check on row marker    
  UInt_t rowControlWord = fWord >> kbit0 & 0xfbff;
  
  if(rowControlWord != statusControlRow) {
    //fRawReader->AddMajorErrorLog(kRowMarkerErr);
    AliDebug(1,Form("Wrong row marker %x expected 0x32a8!",rowControlWord));
              fNumOfErr[fDDLNumber][kRowMarkerErr]++;
              return kFALSE;
   }
//Second check on row marker    
   UInt_t wordsInRow = fWord >> kbit16 & 0x0fff;                // Number of words after the row marker, bit 10 is skipped in this check
  
  if (wordsInRow > nMAXwordsInRow) {
    //fRawReader->AddMajorErrorLog(kRowMarkerSizeErr);
    AliDebug(1,Form(" FATAL: Number of words %x in a row exceeds the expected value: 0x1EA !",wordsInRow));
    fNumOfErr[fDDLNumber][kRowMarkerSizeErr]++;
    return kFALSE;
  }
  
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::GetWord(Int_t n,EDirection dir)
{
  // This method returns the n-th 32 bit word
  // inside the raw data payload.
  // The method is supposed to be endian (platform)
  // independent.
  
  fWord = 0;
  if (fPosition < 0) {
    AliError("fPosition < 0 !!! Event skipped.");
    fRawReader->AddMajorErrorLog(kRawDataSizeErr,"fPosition<0");
    return kFALSE;
  }

  if(dir==kBwd) n = -n; 
  fPosition+=4*n-4;

  if(fPosition==-4) return kTRUE;
  
  if(fPosition<0 || fPosition > fRawReader->GetDataSize()) {
    AliWarning(Form("fPosition out of boundaries %i",fPosition));
    return kFALSE;
  }
    
  StorePosition();
  
  fWord |= fData[fPosition++];
  fWord |= fData[fPosition++] << 8;
  fWord |= fData[fPosition++] << 16;
  fWord |= fData[fPosition++] << 24;

  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::DumpData(Int_t nw)
{
  //just a simple raw data dump
  // in () is the position in bytes
  //--
   for(Int_t i=0;i<nw;i+=4) {
     if(!(i%16)) printf(" \n %8i) ",i);
     printf("%02X%02X%02X%02X [ %06i ] ",fData[i+3],fData[i+2],fData[i+1],fData[i+0],fPos[i]);
   }
   Printf(" \n -----end of dump ----------- ");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::StorePosition()
{
  //just for debug purpose
  // it stores the position
  //read for the first time
//  Printf("@@@@@@@@@ fPos: %x fPosition: %d",fPos,fPosition);
  if(fPos[fPosition]!=0) {
//    Printf("Position already stored!!! Value %i at address %i",fPos[fPosition],fPosition); 
    return;
  }
  fiPos++;
  fPos[fPosition] = fiPos;
//  if(stDeb)Printf("%i - Actual position %i",iPos,fPosition); 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
