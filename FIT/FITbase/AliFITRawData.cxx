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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  FIT raw data conversion class                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//#include <Riostream.h>
//#include <TTree.h>
#include <TMap.h>
#include "AliFIT.h"
#include "AliFITRawData.h"
#include "AliFITDigit.h"
#include "AliBitPacking.h"
#include "AliRawDataHeader.h"
#include "AliRawDataHeaderSim.h"
#include "AliBitPacking.h"
#include "AliFstream.h"
#include "AliRunLoader.h"
#include "AliDAQ.h"
#include "AliLog.h"
#include <iostream>

using std::cout;
using std::endl;
using std::ios_base;

ClassImp(AliFITRawData)

//_____________________________________________________________________________
  AliFITRawData::AliFITRawData():TObject(),
				 fFITdigitArray(NULL),
				 fVerbose(0),      
				 fIndex(-1) ,     
				 fEventNumber(0), 
				 fDataHeaderPos(0),
				 fFile(0x0)   
    
{
  /*
    -  160 channels (2 words each as in TOF DDL) for :
    word 1 :0-5bit number of PMT; word 2: 0-7 error sign, 8-31 TDC
    and the same but for amplified signal. Now I wrote the same time because
    CDF are not ready and differences didn't measured yet.
    
    -  160x2 channel for amplitude: very preliminary, QTC features are not
    known now, preliminary i put as T1 time signal for this PMT in first
    channel and T1+A in second, where A=Log(Amplitude);
    and the same for amplified but A=Log(10*Amplitude).
    
    - Or FIT-A and Or FIT-C 2 channels
    - FITA-FITC vertex information
    - Time Meaner (FITA+FITC)/2
    - 4 MPD  multiplicity signals sum amp both sides
  */
  //open the output file
  // char fileName[15];
  TString fileName = Form("%s",AliDAQ::DdlFileName("FIT",0));
  fFile = new AliFstream(fileName.Data());
  memset(fBuffer,0,512*sizeof(UInt_t));
  
  //get event number 
  AliRunLoader *runloader = AliRunLoader::Instance();
  if (runloader) {
    fEventNumber = runloader->GetEventNumber();
  }
   for ( Int_t k=0; k<1000; k++)   fAllData[k] = -1;

}

//_____________________________________________________________________________

AliFITRawData::AliFITRawData(const AliFITRawData &r):TObject(),
						     fFITdigitArray(NULL),
						     fVerbose(0),      
						     fIndex(-1) ,     
						     fEventNumber(0), 
						     fDataHeaderPos(0),
						     fFile(0x0)  
{
  //
  // AliFITrawData copy constructor
  //
  
  ((AliFITRawData &) r).Copy(*this);
  for ( Int_t k=0; k<1000; k++)   fAllData[k] = -1;

  
}

//_____________________________________________________________________________
AliFITRawData::~AliFITRawData()
{
  //
  // Destructor
  //
}

//_____________________________________________________________________________
AliFITRawData &AliFITRawData::operator=(const AliFITRawData &r)
{
  //
  // Assignment operator
  //
  
  if (this != &r) ((AliFITRawData &) r).Copy(*this);
  return *this;
  
}

//_____________________________________________________________________________
void AliFITRawData::GetDigits()
{
  
  //This method packs the FIT digits in a proper 32 bits structure

  //read FIT digits and fill TDC and ADC arrays
  
  Int_t digit = -1;
  Int_t ndigits = fFITdigitArray->GetEntries();
  AliDebug(2, Form(" Number of read digits = %d",ndigits));
  AliFITDigit *digs;
  for(Int_t i=0; i<1000; i++) fAllData[i]=-1;
  // loop on FIT digits
  for (digit=0; digit<ndigits; digit++) {
    digs = (AliFITDigit*)fFITdigitArray->UncheckedAt(digit);
    Int_t pmt = digs->NPMT();
    fAllData[pmt] = digs->TimeCFD();
    fAllData[pmt+240] = digs->TimeQT0();
    fAllData[pmt+480] = digs->TimeQT1();
  }
  //  Int_t error=0;
  
}
//---------------------------------------------------------------------------------------

Int_t AliFITRawData::RawDataFIT(TBranch* branch)
{
  //This method creates the Raw data files for FIT detector
  
  
  //max number of digits per DDL file times 2
  //  UInt_t fBuffer[kSize];
  //  UInt_t baseWord;
  // UInt_t word;
  cout.setf( ios_base::hex, ios_base::basefield );
  
  fIndex=-1;
  Int_t ch[4] = {0,2,4,6};
  
  fFITdigitArray = * (TClonesArray**) branch->GetAddress();

  AliRawDataHeaderSim header;
  WriteDataHeader(kTRUE, kFALSE);
  
  branch->GetEvent();
  
  GetDigits();
  
  Int_t time,  positionOfTRMHeader, iTDC, channel;
  //space for DRM header
  fIndex += 6;
  Int_t trm1words=0;
  
  for (Int_t itrm=0; itrm <4 ; itrm++) {
    Int_t startTRM=fIndex;
    //space for 1st TRM header
    fIndex ++;
    positionOfTRMHeader= fIndex;
    
    for (Int_t chain=0; chain <2 ; chain++) {
      // space for chain  header
      fIndex ++;
      WriteChainDataHeader(chain+1, 1); // 
      //data TRM 1 chain 1 
      for (Int_t det = 0; det < 60; det++) {
	time = fAllData[det + itrm*120 + chain*60];
	if (time >0 && time !=999999) {
	  fIndex++;
	  iTDC = det / 4;
	  channel = ch[det - iTDC*4];
	  AliDebug(2, Form("det %i  alldata %i trm %i chain %i tdc %i channel %i \n",det, det + itrm*120 + chain*60, itrm, chain, iTDC, det - iTDC*4) );
	  FillTime(channel,iTDC, time);
	}
      } 
      fIndex++;
      WriteChainDataTrailer (chain+1);
    }
    fIndex++;
    WriteTrailer(15, 0,fEventNumber,5); // 1st TRM trailer
    
    trm1words = fIndex - startTRM;
    WriteTRMDataHeader(itrm, trm1words , positionOfTRMHeader);
  }
  
  //DRM trailer
  fIndex++;
  WriteTrailer(1,0,fEventNumber,5);
  WriteDRMDataHeader();
  
  //write packing digits    
  fFile->WriteBuffer((char*) fBuffer,((fIndex+1)*sizeof(UInt_t)));
  //write real data header on its place
  WriteDataHeader(kFALSE, kFALSE);
  
  
  //end for
  
  return 0;  
  
}

//_____________________________________________________________________________

void  AliFITRawData::WriteDRMDataHeader()
{
  //Write a (dummy or real) DDL DRM  data header, 
  //set the compression bit if compressed
  //  UInt_t drmheader[4];  
  cout.setf( ios_base::hex, ios_base::basefield );
  UInt_t word;
  UInt_t baseWord=0;
  //fill DRM headers
  //DRM Global Header
  word = 1;
  AliBitPacking::PackWord(word,baseWord,0, 3); // 0001 
  word = fIndex ;
  AliBitPacking::PackWord(word,baseWord,4, 20); // event words 
  word=124;
  AliBitPacking::PackWord(word,baseWord, 21, 27); // event words 
  word=4;
  AliBitPacking::PackWord(word,baseWord,28, 31);// 0100 marks header
  fBuffer[0]=  baseWord;
  //DRM status header 1
  word = 1;
  AliBitPacking::PackWord(word,baseWord,0, 3); // 0001 
  word = 1;
  AliBitPacking::PackWord(word,baseWord,4, 14); // slotID now 0000000001
  word = 1;
  AliBitPacking::PackWord(word,baseWord,15, 15); //if 1  LHC clock is coorectly recieved from CPDM 
  word=0;
  AliBitPacking::PackWord(word,baseWord,16,27); // reserve for future use
  word=4;
  AliBitPacking::PackWord(word,baseWord,28,31); // 0100 marks header
  fBuffer[1] = baseWord;
  
  word=0;
  baseWord=0;
  
  //DRM status header 2
  word = 1;
  AliBitPacking::PackWord(word,baseWord, 0, 3); // 0001 
  word = 3;
   AliBitPacking::PackWord(word,baseWord, 4, 14); //enable slotID now 00000000011
   word = 0;
   AliBitPacking::PackWord(word,baseWord, 15, 15); // something
   word=0;
   AliBitPacking::PackWord(word,baseWord, 16, 27); // fault ID for simulation 0
   word=4;
   AliBitPacking::PackWord(word,baseWord,28,31); // 0100 marks header
   fBuffer[2]=  baseWord;
   
   word=0;
   baseWord=0;
   //DRM status header 3
   word = 1;
   AliBitPacking::PackWord(word,baseWord,0, 3); // 0001 
   word = 0;
   AliBitPacking::PackWord(word,baseWord,4, 27); // TTC event counter
   word=4;
   AliBitPacking::PackWord(word,baseWord,28,31); // 0100 marks header
   fBuffer[3]=  baseWord;
    
   // new DRM format
   fBuffer[4]=  baseWord;
   fBuffer[5]=  baseWord;
   
   word=0;
   baseWord=0;
   
}
  
//_____________________________________________________________________________

void  AliFITRawData::WriteTRMDataHeader(UInt_t slotID, Int_t nWordsInTRM,
					Int_t  positionOfTRMHeader)
{
  //Write a (dummy or real) DDL TRM  data header, 
  //set the compression bit if compressed
  //  UInt_t trmheader;  
  UInt_t word;
  UInt_t baseWord=0;
  //fill TRM headers
  //TRM Global Header
  word = slotID;
  AliBitPacking::PackWord(word,baseWord,0, 3); // slotID
  word = nWordsInTRM;
 //+this word - DRM header 
  
  AliBitPacking::PackWord(word,baseWord,4, 16); // event words 
  word=0;
  AliBitPacking::PackWord(word,baseWord,17,18); // ACQ
  word=0;
  AliBitPacking::PackWord(word,baseWord,19,19); //  L SEY inside LUT
  word=0;
  AliBitPacking::PackWord(word,baseWord,20,27); //  MBZ
  word=4;
  AliBitPacking::PackWord(word,baseWord,28,31); // 0100 marks header
  fBuffer[positionOfTRMHeader] =  baseWord;
  word=0; 
  baseWord=0;
}

//_____________________________________________________________________________

void  AliFITRawData::WriteChainDataHeader(UInt_t chainNumber,UInt_t slotID)
{
  //Write a (dummy or real) DDL Chain  data header, 
  //set the compression bit if compressed
  //  chainNumber 00 or 10
  UInt_t word;
  UInt_t baseWord=0;
  //fill TRM headers
  //TRM Global Header
  word = slotID; // ask Tatiana 7 or 9 
  AliBitPacking::PackWord(word,baseWord,0, 3); // slotID
  word = 0;
  AliBitPacking::PackWord(word,baseWord,4, 15); // bunchID
  word=0;
  AliBitPacking::PackWord(word,baseWord,16,23); // PB24 temperature
  word=0;
  AliBitPacking::PackWord(word,baseWord,24,26); //  PB24 ID
  word=0;
  AliBitPacking::PackWord(word,baseWord,27,27); //  TS
  word=chainNumber;
  AliBitPacking::PackWord(word,baseWord,28,31); // 0100 marks header
  fBuffer[fIndex] =  baseWord;
  word=0;
  baseWord=0;     
}
//_____________________________________________________________________________

void  AliFITRawData::WriteChainDataTrailer(UInt_t chainNumber )
{
  //Write a (dummy or real) DDL Chain  data trailer 
  //set the compression bit if compressed
  //  chainNumber 00 or 10
  UInt_t word;
  UInt_t baseWord=0;
  word = 0; // ask Tatiana 7 or 9 
  AliBitPacking::PackWord(word,baseWord,0, 3); // status
  word = 0;
  AliBitPacking::PackWord(word,baseWord,4, 15); // MBZ
  word=fEventNumber;
  AliBitPacking::PackWord(word,baseWord,16,27); // event counter
  word=chainNumber;
  AliBitPacking::PackWord(word,baseWord,28,31); // chain number
 // fIndex++;
  fBuffer[fIndex] =  baseWord;

  word=0;
  baseWord=0;       
}
//_____________________________________________________________________________

void  AliFITRawData::WriteDataHeader(Bool_t dummy, Bool_t compressed)
{
  //Write a (dummy or real) DDL data header, 
  //set the compression bit if compressed
  
  AliRawDataHeaderSim header;
  
  if (dummy) {
    //if size=0 it means that this data header is a dummy data header
    fDataHeaderPos = fFile->Tellp();
    fFile->WriteBuffer((char*)(&header), sizeof(header));
  } else {
    UInt_t currentFilePos = fFile->Tellp();
    fFile->Seekp(fDataHeaderPos);
    header.fSize = currentFilePos-fDataHeaderPos;
    header.SetAttribute(0);  // valid data
    if (compressed) header.SetAttribute(1); 
    fFile->WriteBuffer((char*)(&header), sizeof(header));
    fFile->Seekp(currentFilePos);
  }
  
}

//___ __________________________________________________________________________
void  AliFITRawData::WriteTrailer(UInt_t slot, Int_t word1, UInt_t word2, UInt_t word3)
{
  //Write a (dummy or real) DDL Chain  data trailer 
  
  UInt_t word;
  UInt_t baseWord=0;
  word = slot;
  AliBitPacking::PackWord(word,baseWord,0, 3); // 0001 
  word=word1;
  AliBitPacking::PackWord(word,baseWord,4, 15); // CRC ?
  word = word2;
  AliBitPacking::PackWord(word,baseWord,16,27); // event counter
  word=word3;
  AliBitPacking::PackWord(word,baseWord,28,31); //  marks trailer
  fBuffer[fIndex] =  baseWord;
  word=0;
  baseWord=0;
  
}

//---------------------------------------------------------------------------------------
void  AliFITRawData::FillTime(Int_t ch, Int_t iTDC, Int_t time)
{
  //  put all time fields thother in 1 word
  
  UInt_t word;
  UInt_t baseWord=0;
  
  word=time;
  AliBitPacking::PackWord(word,baseWord, 0, 20); // Time 
  
  word=ch;
  AliBitPacking::PackWord(word,baseWord, 21, 23); // number of channel 
  word=iTDC;
  AliBitPacking::PackWord(word,baseWord, 24, 27); // TDC ID
  
  word=0;
  AliBitPacking::PackWord(word,baseWord, 28, 28); // E = 0 in simulation
  word=0;
  AliBitPacking::PackWord(word,baseWord, 29, 30); // PS bit data 00
  word=1;
  AliBitPacking::PackWord(word,baseWord, 31, 31); // 1
  fBuffer[fIndex]=baseWord;
  word=0;
  baseWord=0;
}
