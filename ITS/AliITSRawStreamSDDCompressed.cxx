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

/* $Id $*/

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to decode compressed SDD Raw Data format                //
// The 32 bits for a data word are defined as follows:           //
//   31 control bit (0=data word, 1= control word)               //
//   30 -                                                        //
//   29  |                                                       //
//   28  |-> 4 bits to identify the Carlos (0-11) inside the DDL //
//   27 -                                                        //
//   26 detecor side (0= left, =right)                           //
//   25 -                                                        //
//   24  |                                                       //
//   23  |                                                       //
//   22  |                                                       //
//   21  |-> 8 bits to identify the anode number (0-255)         //
//   20  |                                                       //
//   19  |                                                       //
//   18 -                                                        //
//   17 -                                                        //
//   16  |                                                       //
//   15  |                                                       //
//   14  |                                                       //
//   13  |-> 8 bits to identify the time bin (0-255)             //
//   12  |                                                       //
//   11  |                                                       //
//   10 -                                                        //
//    9 -                                                        //
//    8  |                                                       //
//    7  |                                                       //
//    6  |                                                       //
//    5  |                                                       //
//    4  |-> 10 bit for the ADC counts                           //
//    3  |                                                       //
//    2  |                                                       //
//    1  |                                                       //
//    0 -                                                        //
//                                                               //
// Plus 3 types of control words                                 //
// identified by the the 4 most significant bits (31-28)         //
// 1) Jitter word     = 1000                                     //
// 2) JTAG answer     = 1100                                     //
// 3) End of module data (needed by the Cluster Finder)   = 1111 //
//                                                               //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////



#include "AliITSRawStreamSDDCompressed.h"
#include "AliRawReader.h"
#include "AliLog.h"

ClassImp(AliITSRawStreamSDDCompressed)
  


//______________________________________________________________________
AliITSRawStreamSDDCompressed::AliITSRawStreamSDDCompressed(AliRawReader* rawReader) :
  AliITSRawStream(rawReader),
  fDDLModuleMap(0),
  fData(0),
  fCarlosId(-1),
  fChannel(0),
  fJitter(0),
  fDDL(0),
  fADCEncoded(0)
{
// create an object to read ITS SDD raw digits
  for(Int_t im=0;im<kSDDModules;im++){
    fLowThresholdArray[im][0]=0;
    fLowThresholdArray[im][1]=0;
  }
  fRawReader->Reset();
  fRawReader->Select("ITSSDD");


}

//______________________________________________________________________
AliITSRawStreamSDDCompressed::AliITSRawStreamSDDCompressed(const AliITSRawStreamSDDCompressed& rs) :
  AliITSRawStream(rs.fRawReader),
  fDDLModuleMap(rs.fDDLModuleMap),
  fData(0),
  fCarlosId(-1),
  fChannel(0),
  fJitter(0),
  fDDL(0),
  fADCEncoded(0)
{
  // copy constructor
  AliError("Copy constructor should not be used.");
}
//__________________________________________________________________________
AliITSRawStreamSDDCompressed& AliITSRawStreamSDDCompressed::operator=(const AliITSRawStreamSDDCompressed& rs) {
  // assignment operator
  if (this!=&rs) {}
  AliError("Assignment opertator should not be used.");
  return *this;
}

//______________________________________________________________________
AliITSRawStreamSDDCompressed::~AliITSRawStreamSDDCompressed(){
  if(fDDLModuleMap) delete fDDLModuleMap;
}


//______________________________________________________________________
Int_t AliITSRawStreamSDDCompressed::DecompAmbra(Int_t value) const
{
  // AMBRA decompression (from 8 to 10 bit)
  
  if ((value & 0x80) == 0) {
    return value & 0x7f;
  } else if ((value & 0x40) == 0) {
    if(value&1) return 0x080 + ((value & 0x3f) << 1);
    return 0x081 + ((value & 0x3f) << 1);
  } else if ((value & 0x20) == 0) {
    if(value&1) return 0x103 + ((value & 0x1f) << 3);
    return 0x104 + ((value & 0x1f) << 3);
  } else {
    if(value&1) return 0x207 + ((value & 0x1f) << 4);
    return 0x208 + ((value & 0x1f) << 4);
  }
  
}
//______________________________________________________________________
Bool_t AliITSRawStreamSDDCompressed::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left
// returns kTRUE and fCompletedModule=kFALSE and fCompletedDDL=kFALSE when a digit is found
// returns kTRUE and fCompletedModule=kTRUE  and fCompletedDDL=kFALSE when a module is completed (=3x3FFFFFFF footer words)
// returns kTRUE and fCompletedModule=kFALSE and fCompletedDDL=kTRUE  when a DDL is completed (=jitter word)


  UInt_t idJit=8;    // Jitter word has the 4 most significant bits = 1000
  UInt_t idEom=15;   // end of module has the 4 most significant bits = 1111
  UInt_t idJtag=12;  // end of module has the 4 most significant bits = 1100

  UInt_t maskmod=15;   // last 4 bits for module number in end of module word
  UInt_t maskCarlos=15<<27; // 4 bits  (27-30) for CarlosId in data word
  UInt_t maskSide=1<<26;    // 1 bit   (26)    for side     in data word
  UInt_t maskAnode=255<<18; // 8 bits  (18-25) for Nanode   in data word
  UInt_t maskTb=255<<10;    // 8 bits  (10-27) for Ntimebin in data word
  UInt_t maskADC=1023;      // 10 bits (0-9)   for ADC      in data word
  UInt_t maskCode=7;        // 3 bits (0-2)    for ADC range in encoded-ADC case
    
  while(kTRUE){
    if (!fRawReader->ReadNextInt(fData)) return kFALSE;  // read next word
    UInt_t mostsigbits=fData>>28; 
    if(fData==0xFFFFFFFF){ 
      // CarlosRX header do nothing
    } else if(mostsigbits==idEom){
      // end of module word
      fCarlosId=fData&maskmod;
      fDDL=fRawReader->GetDDLID();
      fModuleID = GetModuleNumber(fDDL,fCarlosId);
      fCompletedDDL=kFALSE;
      fCompletedModule=kTRUE;
      return kTRUE;
    } else if(mostsigbits==idJit){
      // jitter word
      fJitter = fData&0x000000ff;      
      fCompletedModule=kFALSE;
      fCompletedDDL=kTRUE;
      return kTRUE;
    } else if(mostsigbits==idJtag){
      // jtag word -> skipped
      continue;
    }else if(mostsigbits<8){
      // data word
      fCarlosId=(fData&maskCarlos)>>27;
      fDDL=fRawReader->GetDDLID();
      fModuleID = GetModuleNumber(fDDL,fCarlosId);
      fChannel=(fData&maskSide)>>26;
      fCoord1=(fData&maskAnode)>>18;
      fCoord2=(fData&maskTb)>>10;
      Int_t sig8bit;
      if(fADCEncoded){
	UInt_t code=fData&maskCode;
	if (code < 2 || code > 7){ 
	  AliError(Form("Wrong ADC code value %d",code));
	  continue;
	}
	UInt_t adcmask=(1<<code)-1;
	sig8bit=((fData&(adcmask<<3))>>3) + (1<<code);
      }else{      
	sig8bit=fData&maskADC;
      }
      sig8bit+=fLowThresholdArray[fModuleID-kSPDModules][fChannel];
      fSignal=DecompAmbra(sig8bit);
      fCompletedModule=kFALSE;
      fCompletedDDL=kFALSE;
      return kTRUE;
    }
  }
  return kFALSE;
}


