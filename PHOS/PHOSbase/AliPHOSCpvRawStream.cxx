/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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

#include "AliPHOSCpvRawStream.h"
#include "AliRawReader.h"
#include <stdio.h>

ClassImp(AliPHOSCpvRawStream);
using namespace std;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliPHOSCpvRawStream::AliPHOSCpvRawStream(AliRawReader* rawReader) :
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
  fPosition(0),
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
  Int_t kNDDL = AliPHOSCpvParam::kNDDL;
  fNumOfErr = new Int_t*[kNDDL];     // Store the number of errors for a given error type and a given DD
  for(Int_t iddl=0; iddl<kNDDL; iddl++) {
    fNumOfErr[iddl] = new Int_t [kSumErr];
    for(Int_t ierr=0; ierr < kSumErr; ierr++) {
      fNumOfErr[iddl][ierr]=0;       // reset errors
    }
  }
  
  fnDDLInStream  = new Int_t[kNDDL];
  fnDDLOutStream = new Int_t[kNDDL];
  for(Int_t i=0; i<kNDDL; i++) {
    fnDDLInStream [i] = -1;
    fnDDLOutStream[i] = -1;
  }
 
  fRawReader->Reset();
  fRawReader->Select("CPV");
} // Constructor

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliPHOSCpvRawStream::AliPHOSCpvRawStream() :
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
  fPosition(0),
  fWord(0),
  fZeroSup(kTRUE),
  fPos(0x0),
  fiPos(0),
  fTurbo(kFALSE) ,
  fRawDataSize(0)
{
  //
  // Default constructor
  //
  Int_t kNDDL = AliPHOSCpvParam::kNDDL;
  fNumOfErr = new Int_t*[kNDDL]; // Store the number of errors for a given error type and a given DDL
  for(Int_t iddl=0; iddl<kNDDL; iddl++) {
    fNumOfErr[iddl] = new Int_t [kSumErr];
    for(Int_t ierr=0; ierr < kSumErr; ierr++) {
      fNumOfErr[iddl][ierr]=0;   // reset errors
    }
  }

  fnDDLInStream  = new Int_t[kNDDL];
  fnDDLOutStream = new Int_t[kNDDL];
  for(Int_t i=0; i<kNDDL; i++) {
    fnDDLInStream [i] = -1;
    fnDDLOutStream[i] = -1;
  }
} // Default constructor
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliPHOSCpvRawStream::~AliPHOSCpvRawStream()
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
  for(Int_t i=0; i<AliPHOSCpvParam::kNDDL; i++) delete [] fNumOfErr[i]; 
  delete [] fNumOfErr; 
  delete [] fnDDLInStream;
  delete [] fnDDLOutStream;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliPHOSCpvRawStream::Reset()
{
  // reset raw stream params
  // Reinitalize the containers
  fDDLNumber = -1;
  fLDCNumber =  0;
  fTimeStamp =  0;
  fPosition = 0;
  fData = NULL;
  if (fRawReader) fRawReader->Reset();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::Turbo()
{
  
  UInt_t row,_3G,pad;
  Int_t cntGlob = fRawReader->GetDataSize()/4;
  fPosition = 0;
  fNPads    = 0;

  //std::cout<<"words to be read: "<<cntGlob<<std::endl;

  for(Int_t i=1; i<cntGlob; i++) {
    if (!GetWord(1)) return kFALSE;
    //std::cout<<"i've passed getWord(1)"<<std::endl;
    if (fPosition/4 <= 5) continue; // Skip first 5 words 
    
    row = ((fWord >> kbit22) & 0x1f) - 1;
    _3G = ((fWord >> kbit18) & 0xf) - 1; // 3GASSIPLEX info in raw word is between bits: 18...21    
    pad = (fWord >> kbit12) & 0x3f;      // pad info in raw word is between bits: 12...17
    //std::cout<<"row = "<<row<<", 3Gassiplex = "<<_3G<<", pad info = "<<pad<<std::endl;
    
    Int_t charge, abs, eType;
    if(!AliPHOSCpvParam::DecodeRawWord(fDDLNumber,fWord,abs,charge,eType)) {
      if (eType > 0) {
	fNumOfErr[fDDLNumber][eType]++;
	//std::cout<<"AliPHOSCpvRawStream::Turbo(): I cannot decode word!"<<std::endl;
	//cout<<"DDL = "<<  fDDLNumber << "; word = "<< fWord <<"; abs = " << abs 
	//    <<"; charge = "<< charge <<"; etype = "<< eType << endl;
      }
      continue;
    }
    else {
      fPad   [fNPads] = abs;
      fCharge[fNPads] = charge;
      fNPads++;
      if(charge==0) fNumOfErr[fDDLNumber][kPedQZero]++;
      //Printf("Size: %i  DDL %i row %i 3G %i pad %i fPos %i fNPads: %i Charge: %d fWord %4.4x",cntGlob,fDDLNumber,row,_3G,pad,fPosition,fNPads,fCharge[fNPads-1],fWord);
    }
  }//word loop
  return kTRUE;
}//Turbo()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::Next()
{
  // read next DDL raw data from the CPV raw data stream
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
  fDDLNumber = -1;
  if (fRawReader->GetDDLID()>=0)
    fDDLNumber = fRawReader->GetDDLID();

  //debug!!!
  //fDDLNumber=0;

  if(fDDLNumber<0) {
    AliWarning(Form("fDDLNumber not a acceptable value %i",fDDLNumber));
    return kFALSE;
  }
    
  if(fRawReader->GetType() == 7 || fRawReader->GetType() == 8 )  {  //New: Select Physics events, Old: Raw data size is not 0 and not 47148 (pedestal)
    fnDDLInStream[fDDLNumber]=1; 
    fnDDLOutStream[fDDLNumber]=0;
    
    fLDCNumber = fRawReader->GetLDCId();
    fTimeStamp = fRawReader->GetTimestamp();
    
    fRawDataSize=fRawReader->GetDataSize()/4;
    DelVars();                    //We have to delete the variables initialized in the InitVars before recall InitVars!!!
    InitVars(fRawDataSize);       //To read the charge and pads we cannot delete before the status return
  
    if(fTurbo==kTRUE) status=Turbo();
    else status = ReadCPVRawData(); // If there is any error in raw event, the whole event is thrown away
   
    if(status) AliDebug(1,Form("Event DDL %i successfully decoded!.",fDDLNumber));
    else       AliDebug(1,Form("Event DDL %i ERROR in decoding!.",fDDLNumber));
    //DumpData(fRawReader->GetDataSize());
  } // if(Select Physics events)

  if(status==kTRUE) fnDDLOutStream[fDDLNumber]++; //Count the number of events when the DDL was succesfully decoded
   
   return status;
} // Next()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliPHOSCpvRawStream::InitVars(Int_t n)
{
  fNPads = 0;
  fCharge = new Int_t[n]; 
  fPad = new Int_t[n];
  fPos = new Int_t[4*n];                 //reset debug 
  for(Int_t ie = 0 ; ie < 4*n; ie++) fPos[ie] = 0; //initialize for 0, otherwise the position is considered filled and will not be updated for the dump
  fiPos = 0;
} // InitVars()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliPHOSCpvRawStream::DelVars()
{
  //Clean the initvars!!!!!!!!
  fNPads = 0; 
  if (fCharge)     { delete [] fCharge;    fCharge = 0x0; }
  if (fPad)        { delete [] fPad;       fPad = 0x0;       }   
  if (fPos)        { delete [] fPos;       fPos = 0x0;       }     
 
  fiPos=0;
} // DelVars()
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::ReadCPVRawData()
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

  while (cnt>20) { //counter limit is changed from 0 to 20 to take into account (skip) the 5 extra words in the equipment header
    nwSeg = (fWord >> kbit8) & 0xfff;
    if(!CheckSegment()) return kFALSE; 
    if(!ReadSegment(cntSegment)) return kFALSE;
    AliDebug(1,"ReadCPVRawData(): ReadSegment() passed");
    if(nwSeg != cntSegment) return kFALSE;
    if(!GetWord(cntSegment+1,kBwd)) return kFALSE;
    cnt -= cntSegment+1;
  }
  return kTRUE;
} // ReadCPVRawData()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::ReadSegment(Int_t &cntSegment)
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
    //    cntRow  = (fWord >> kbit16) & 0xfff;
    cntRow = 490;
    if(!CheckRowMarker()) return kFALSE;
    if(!ReadRow(nwRow))   return kFALSE;
    if(nwRow != cntRow) {AliDebug(1,Form("Error in Row counters: %i different wrt %i",nwRow,cntRow)); return kFALSE;}
    if(!GetWord(cntRow+1)) return kFALSE;
    cnt -= cntRow+1;
  }
  cntSegment -= cnt;
  return kTRUE;
} // ReadSegment()    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::ReadRow(Int_t &cntRow)
{
  // Read the row
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK

  // 3G = 3-gassiplex
  Int_t cnt;
  Int_t cnt3G;
  Int_t nw3G;
  
  //  cntRow  = (fWord >> kbit16) & 0xfff;
  cntRow = 490;
  cnt = cntRow;  
  
  if(!GetWord(cntRow)) return kFALSE;

  while (cnt>480) {
    if(!CheckEoE(nw3G)) return kFALSE;
    if(!Read3G(cnt3G)) return kFALSE;
    if(nw3G != cnt3G) {AliDebug(1,Form("Error in 3gassiplex counters: %i different wrt %i",nw3G,cnt3G));return kFALSE;}
    if(!GetWord(1,kBwd)) return kFALSE; // go to next 3gassiplex...
    cnt--;
  }
  
  cntRow -= cnt;  
  cntRow += 480;
  return kTRUE;
} // ReadRow()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::Read3G(Int_t &cnt3G)
{
  // Read the 3gassiplex
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK

  cnt3G = fWord & 0x7f;
  
  Int_t cnt = cnt3G;

  for(Int_t i3G=0; i3G<cnt3G; i3G++) {
    UInt_t _3G = 0, row = 0;
    if(!GetWord(1,kBwd)) return kFALSE;
    //check on row number      
    cnt--;
    row = ((fWord >> kbit22) & 0x1f) - 1;
    if(!CheckRow(row)) continue;
    //check 3gassiplex number     
    _3G = ((fWord >> kbit18) & 0xf) - 1;     //3gassiplex info in raw word is between bits: 18...21
    if(!Check3G(_3G)) continue;
    //check pad number
    UInt_t pad = (fWord >> kbit12) & 0x3f; //pad info in raw word is between bits: 12...17
    if(!CheckPad(pad)) continue;
    Int_t charge = fWord & 0xfff;
    if(AliPHOSCpvParam::Abs(fDDLNumber,row,_3G,pad)<0) continue;

    fPad[fNPads] = AliPHOSCpvParam::Abs(fDDLNumber,row,_3G,pad);
    fCharge[fNPads] = charge; 
    fNPads++;
    
    if(charge==0) 
    {
      AliDebug(1,Form("If PEDESTAL run -> WARNING: ZERO charge is read from DDL: %d row: %d 3G: %d pad: %d",fDDLNumber,row,_3G,pad));
      fNumOfErr[fDDLNumber][kPedQZero]++;
    }
    
  }//i3G

  cnt3G -= cnt;  
  return kTRUE;  
} // ReadDilogic
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::CheckSegment()
{
  // Check the segment marker
  // It returns: kFALSE if any error occurs
  //             kTRUE  if all OK

  UInt_t markSegment = 0xAB0;
  UInt_t segMarker = (fWord >> kbit20) & 0xfff;

  if (segMarker != markSegment ) {
    AliDebug(1,Form("Segment marker %X wrong (expected %0X)! at %i in word %0X!",segMarker,markSegment,fPosition,fWord));
    fNumOfErr[fDDLNumber][kWrongSegErr]++;
    return kFALSE;
  }
  
  UInt_t segAddress = fWord & 0xff;
  if (segAddress<1 ||segAddress>3) {
    AliDebug(1,Form("Segment address %d not in the valid range [1-3]",segAddress));
    fNumOfErr[fDDLNumber][kWrongSegErr]++;
    return kFALSE;
  }
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::CheckRow(UInt_t row)
{
  //check on row number      
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK
  if(row>=0 && row <AliPHOSCpvParam::kNRows) return kTRUE;
  AliDebug(1,Form("Wrong row index: %d, expected (0 -> %d) word %0X at %i...",row,AliPHOSCpvParam::kNRows-1,fWord,fPosition));
  fNumOfErr[fDDLNumber][kWrongRowErr]++;
  return kFALSE;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::Check3G(UInt_t _3G)
{
  //check 3gassiplex number     
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK
  if (_3G>= 0 && _3G <AliPHOSCpvParam::kN3GAdd) return kTRUE;
  AliDebug(1,Form("Wrong 3GASSIPLEX index: %d, expected (0 -> %d)!",_3G,AliPHOSCpvParam::kN3GAdd-1));
  fNumOfErr[fDDLNumber][kWrong3GErr]++;
  return kFALSE;
}      
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::CheckPad(UInt_t pad)
{
  //check pad number     
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK
  if (pad < AliPHOSCpvParam::kNPadAdd) return kTRUE;
  AliDebug(1,Form("Wrong pad index: %d, expected (0 -> %d)!",pad,AliPHOSCpvParam::kNPadAdd - 1));
  fNumOfErr[fDDLNumber][kWrongPadErr]++;
  return kFALSE;
}    
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::CheckEoE(Int_t &nDil)
{
  //check the End of Event
  // "End of event" is not end of event. It checks subordinate words between 3gassiplex blocks
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK

  if (!((fWord >> kbit27) & 0x1)) {      //check 27th bit in EoE. It must be 1!
    AliDebug(1,Form("Missing end-of-event flag! (%08X) at %i",fWord,fPosition/4));
    fNumOfErr[fDDLNumber][kEoEFlagErr]++;
    return kFALSE;
  }
  nDil = fWord & 0x7f;           //nDil=EoE word count
  if(nDil < 0 || nDil > 48 ) { 
    AliDebug(1,Form("Wrong end-of-event word-count: %08X",fWord));
    fNumOfErr[fDDLNumber][kEoESizeErr]++;
    return kFALSE;
  }
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::CheckRowMarker()
{
  //check the row marker
  //It returns: kFALSE if any error occurs
  //            kTRUE  if all OK
  UInt_t nMAXwordsInRow = 0x1EA;
  UInt_t statusControlRow = 0x32A8; // 0x36A8 for zero suppression
//First check on row marker    
  UInt_t rowControlWord = fWord >> kbit0 & 0xffff;

  if((fWord >> kbit27) & 0x1) return kTRUE;

  if(rowControlWord != statusControlRow) {
    AliDebug(1,Form("Wrong row marker %x expected 0x32a8!",rowControlWord));
    fNumOfErr[fDDLNumber][kRowMarkerErr]++;
    return kFALSE;
  }

  //Second check on row marker    
  UInt_t wordsInRow = fWord >> kbit16 & 0x0fff;    // Number of words after the row marker, bit 10 is skipped in this check
  
  if (wordsInRow > nMAXwordsInRow) {
    AliDebug(1,Form(" FATAL: Number of words %x in a row exceeds the expected value: 0x1EA !",wordsInRow));
    fNumOfErr[fDDLNumber][kRowMarkerSizeErr]++;
    return kFALSE;
  }
  
  return kTRUE;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvRawStream::GetWord(Int_t n,EDirection dir)
{
  // This method returns the n-th 32 bit word
  // inside the raw data payload.
  // The method is supposed to be platform independent.
  
  fWord = 0;
  if (fPosition < 0) {
    AliError("fPosition < 0 !!! Event skipped.");
    fRawReader->AddMajorErrorLog(kRawDataSizeErr,"fPosition<0");
    return kFALSE;
  }

  if(dir==kBwd) n = -n; 
  fPosition += 4*n-4;

  if(fPosition == -4) return kTRUE;
  
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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliPHOSCpvRawStream::DumpData(Int_t nw)
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliPHOSCpvRawStream::StorePosition()
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
