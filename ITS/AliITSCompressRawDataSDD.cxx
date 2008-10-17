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

/* $Id$*/

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to decode the SDD Raw Data from the CarlosRX format to  //
// a compressed format consisting in a word of 32 bit per cell   //
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
// Plus 2 typs of control words:                                 //
// - End of module data (needed by the Cluster Finder)           //
//       first 4 most significant bits                   = 1111  //
// - Jitter word                                                 //
//       first 4 most significant bits                   = 1000  //
//                                                               //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include"AliLog.h"
#include "AliITSCompressRawDataSDD.h"
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliITSRawStreamSDD.h"


ClassImp(AliITSCompressRawDataSDD)

AliITSCompressRawDataSDD::AliITSCompressRawDataSDD():
TObject(),
fRawReader(0),
fPointerToData(0),
fSizeInMemory(0)
{
  // default constructor
}
//______________________________________________________________________
AliITSCompressRawDataSDD::~AliITSCompressRawDataSDD(){
  // raw reader is passed from outdside, don't delete it
}
//______________________________________________________________________
UInt_t AliITSCompressRawDataSDD::CompressEvent(UChar_t* inputPtr){
  // Method to be used in HLT
  UInt_t siz=0;
  memcpy(fPointerToData,inputPtr,32); // event header, 8 words
  fPointerToData+=32;
  siz+=32;
  UInt_t word=0;
  AliITSRawStreamSDD s(fRawReader);
  s.SetDecompressAmbra(kFALSE);
  Int_t mask1=0xFF000000;
  Int_t mask2=0x00FF0000;
  Int_t mask3=0x0000FF00;
  Int_t mask4=0x000000FF;
  while(s.Next()){
    if(s.IsCompletedModule()==kTRUE){
      word=15<<28;
      word+=s.GetCarlosId();
      if(siz+4<fSizeInMemory){
	*(fPointerToData)=(word&mask4);
	++fPointerToData;
	*(fPointerToData)=(word&mask3)>>8;
	++fPointerToData;
	*(fPointerToData)=(word&mask2)>>16;
	++fPointerToData;
	*(fPointerToData)=(word&mask1)>>24;
	++fPointerToData;
	siz+=4;
      }
    }else if(s.IsCompletedDDL()==kTRUE){
      word=8<<28;
      word+=s.GetJitter();
      if(siz+4<fSizeInMemory){
	*(fPointerToData)=(word&mask4);
	++fPointerToData;
	*(fPointerToData)=(word&mask3)>>8;
	++fPointerToData;
	*(fPointerToData)=(word&mask2)>>16;
	++fPointerToData;
	*(fPointerToData)=(word&mask1)>>24;
	++fPointerToData;
	siz+=4;
      }
    }else{
      word=s.GetCarlosId()<<27;
      word+=s.GetChannel()<<26;
      word+=s.GetCoord1()<<18;
      word+=s.GetCoord2()<<10;
      word+=s.GetEightBitSignal();
      if(siz+4<fSizeInMemory){
	*(fPointerToData)=(word&mask4);
	++fPointerToData;
	*(fPointerToData)=(word&mask3)>>8;
	++fPointerToData;
	*(fPointerToData)=(word&mask2)>>16;
	++fPointerToData;
	*(fPointerToData)=(word&mask1)>>24;
	++fPointerToData;
	siz+=4;
      }
    }
  }
  return siz;
}
