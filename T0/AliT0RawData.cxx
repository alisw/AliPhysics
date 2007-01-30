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
//  T0 raw data conversion class                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TTree.h>

#include "AliT0.h"
#include "AliT0RawData.h"
#include "AliT0digit.h"
#include "AliBitPacking.h"
#include "AliRawDataHeader.h"
#include "AliRawDataHeaderSim.h"
#include "AliBitPacking.h"
#include "AliFstream.h"
#include "AliRunLoader.h"

ClassImp(AliT0RawData)

//_____________________________________________________________________________
AliT0RawData::AliT0RawData():TObject()
{
  /*
-  48 channels (2 words each as in TOF DDL) for :
word 1 :0-5bit number of PMT; word 2: 0-7 error sign, 8-31 TDC
and the same but for amplified signal. Now I wrote the same time because
CDF are not ready and differences didn't measured yet.

-  48 channel for amplitude: very preliminary, QTC features are not
known now, preliminary i put as T1 time signal for this PMT in first
channel and T1+A in second, where A=Log(Amplitude);
and the same for amplified but A=Log(10*Amplitude).

- T0-A and T0-C 2 channels
- T0A-T0C vertex information
- Time Meaner where T0C TOF increase to the T0A TOF distance
- 6 multiplicity signals the same way as amplitude and with the same
uncertances
  */

  fIndex=-1;
  fDigits = NULL;

  fTimeCFD = new TArrayI(24);
  fADC1 = new TArrayI(24);
  fTimeLED = new TArrayI(24);
  fADC0 = new TArrayI(24);
  fFile = NULL;
  fDataHeaderPos = 0;
  fDRMDataHeaderPos = 0; 
  memset(fBuffer,0,512*sizeof(UInt_t));

  //open the output file
  char fileName[15];
  sprintf(fileName,"T0_%d.ddl", 0xd00);
  fFile = new AliFstream(fileName);
  //get event number 
  AliRunLoader *runloader = AliRunLoader::GetRunLoader();
  if (runloader) {
    fEventNumber = runloader->GetEventNumber();
  }
}

//_____________________________________________________________________________
AliT0RawData::AliT0RawData(const AliT0RawData &r):TObject()
{
  //
  // AliT0rawData copy constructor
  //

  ((AliT0RawData &) r).Copy(*this);

}

//_____________________________________________________________________________
AliT0RawData::~AliT0RawData()
{
  //
  // Destructor
  //
  if (fDigits) {
    delete fDigits;
    fDigits = NULL;
  }
  delete fTimeCFD;
  delete fADC1;
  delete fTimeLED;
  delete fADC0;
}

//_____________________________________________________________________________
AliT0RawData &AliT0RawData::operator=(const AliT0RawData &r)
{
  //
  // Assignment operator
  //

  if (this != &r) ((AliT0RawData &) r).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliT0RawData::GetDigits(AliT0digit *fDigits)
{
 
  //This method packs the T0 digits in a proper 32 bits structure

  //read T0 digits and fill TDC and ADC arrays


  //  Int_t error=0;
  Int_t time,  positionOfTRMHeader;
  
  // Get the digits array
  
  fDigits->GetTime(*fTimeCFD);
  fDigits->GetADC(*fADC1);
  fDigits->GetTimeAmp(*fTimeLED);
  fDigits->GetADCAmp(*fADC0);
  Int_t meantime = fDigits->MeanTime(); 
  Int_t timediff = fDigits->TimeDiff(); 
  Int_t mult0=fDigits->SumMult();
  Int_t mult1=fDigits->SumMult();
  Int_t timeA = fDigits->BestTimeLeft();
  Int_t timeC = fDigits->BestTimeRight();
  
  
  TArrayI  *allData = new TArrayI(110);
  Int_t i=0;
  allData->AddAt(0,0);
  for (i=1; i<25; i++) {
    allData->AddAt(fTimeLED->At(i-1),i);
    allData->AddAt(fTimeCFD->At(i-1),i+24);
    allData->AddAt(fADC0->At(i-1),i+54);
    allData->AddAt(fADC1->At(i-1),i+78);
    //    cout<<i<<" led "<<fTimeLED->At(i-1)<<" cfd "<<fTimeCFD->At(i-1)<<" qt0 "<<fADC0->At(i-1)<<" qt1 "<<fADC1->At(i-1)<<endl;
  }
  allData->AddAt(meantime,49);
  allData->AddAt(timediff,50);
  allData->AddAt(timediff,103); //trigger vertex
  allData->AddAt(timeA,51);
  allData->AddAt(timeA,104); //trigger T0A
  allData->AddAt(timeC,52);
  allData->AddAt(timeC,105); //trigger T0C
  allData->AddAt(mult0,53);
  allData->AddAt(mult1,106); //trigger central
  allData->AddAt(mult1,54);
  allData->AddAt(mult1,107); //trigger semi-central

  // cout<<" new Event "<<endl;
  //   for (Int_t ia=0; ia<110; ia++) cout<<ia<<" "<<allData->At(ia)<<endl;
  //space for DRM header
  fIndex += 4;

  //space for 1st TRM header
  fIndex ++;
  positionOfTRMHeader= fIndex;

  //space for chain  header
  fIndex ++;

  // Loop through all PMT
  Int_t chain=0; 
  Int_t iTDC = 0;
  Int_t channel=0;
  Int_t trm1words=0;
  Int_t fWordsIn1stTRM=0;
  //LED
  for (Int_t det = 0; det < 55; det++) {
    time = allData->At(det);

    if (time >0) { 
      FillTime(channel,  iTDC,  time); 
      trm1words++;   
     }
    if (channel < 6) channel +=2;
    else {
      channel = 0; 
      iTDC++;
      if (iTDC>15) { chain++; iTDC=0;}
    }
    //  cout<<det<<" "<<time<<" "<<iTDC<<" "<<channel<<" "<<endl;
  }
  
  WriteTrailer(0,0,fEventNumber,1); // 1st chain trailer
  WriteTrailer(15,0,fEventNumber,5); // 1st TRM trailer
  fWordsIn1stTRM = trm1words + 4;
  //  WriteTRMDataHeader(3, trm1words , positionOfTRMHeader);
  WriteTRMDataHeader(0, trm1words , positionOfTRMHeader);


  //space for 2st TRM header
  fIndex ++;
  positionOfTRMHeader= fIndex;

  //space for chain  header
  fIndex ++;


  chain=0; 
  iTDC = 0;
  channel=0;
  Int_t trm2words=0;
  for (Int_t det = 55; det < 108; det++) {
    time = allData->At(det);

    if (time >0) { 
      FillTime(channel,  iTDC,  time); 
      trm2words++;}
    if (channel < 6) channel +=2;
    else {
      channel = 0; 
      iTDC++;
      if (iTDC>15) { chain++; iTDC=0;}
    }
    //    cout<<det<<" "<<time<<" "<<channel<<" "<<iTDC<<endl;
 }

  WriteTrailer(0,0,fEventNumber,1); // 1st chain trailer
  WriteTrailer(15,0,fEventNumber,5); // 1st TRM trailer
  //  WriteTRMDataHeader(5,trm2words,positionOfTRMHeader);
  WriteTRMDataHeader(1,trm2words,positionOfTRMHeader);

  WriteTrailer(1,fEventNumber,0,5); //DRM trailer
  WriteDRMDataHeader();

}
//------------------------------------------------------------------------------
void AliT0RawData::PackWord(UInt_t &BaseWord, UInt_t Word, Int_t StartBit, Int_t StopBit)
{
  
  // Build mask
  Int_t len=StopBit-StartBit+1;
  UInt_t mask=0;
  for(Int_t jb=0; jb<len; mask|=1<<jb++);
  // Check consistency
  if(Word > mask){
    Error("PackWord", "Word to be filled is not within desired length\n"
          "Word:%d Start bit:%d Stop Bit:%d",Word,StartBit,StopBit);
    return;
  }
  BaseWord=(BaseWord&~(mask<<StartBit))|Word<<StartBit;

}


//_____________________________________________________________________________

void  AliT0RawData::WriteDRMDataHeader()
{
//Write a (dummy or real) DDL DRM  data header, 
//set the compression bit if compressed
//  UInt_t drmheader[4];  
  UInt_t word;
  UInt_t baseWord=0;
  //fill DRM headers
  //DRM Global Header
  word = 1;
  PackWord(baseWord,word,0, 3); // 0001 
  word = fIndex ;

  PackWord(baseWord,word,4, 20); // event words 
  word=124;
  PackWord(baseWord,word,21,27); // DRM ID for T0 - 124
  word=4;
  PackWord(baseWord,word,28,31); // 0100 marks header
  fBuffer[0]=  baseWord;

  //DRM status header 1
  word = 1;
  PackWord(baseWord,word,0, 3); // 0001 
  word = 1;
  PackWord(baseWord,word,4, 14); // slotID now 0000000001
  word = 1;
  PackWord(baseWord,word,15, 15); //if 1  LHC clock is coorectly recieved from CPDM 
  word=0;
  PackWord(baseWord,word,16,27); // reserve for future use
  word=4;
  PackWord(baseWord,word,28,31); // 0100 marks header
   fBuffer[1] = baseWord;

  word=0;
  baseWord=0;

    //DRM status header 2
    word = 1;
    PackWord(baseWord,word, 0, 3); // 0001 
    word = 3;
    PackWord(baseWord,word, 4, 14); //enable slotID now 00000000011
    word = 0;
    PackWord(baseWord,word, 15, 15); // something
    word=0;
    PackWord(baseWord,word, 16, 27); // fault ID for simulation 0
    word=4;
    PackWord(baseWord,word,28,31); // 0100 marks header
    fBuffer[2]=  baseWord;

    
    word=0;
    baseWord=0;
    //DRM status header 3
    word = 1;
    PackWord(baseWord,word,0, 3); // 0001 
    word = 0;
    PackWord(baseWord,word,4, 27); // TTC event counter
    word=4;
    PackWord(baseWord,word,28,31); // 0100 marks header
    fBuffer[3]=  baseWord;


    word=0;
    baseWord=0;
    
}
  
//_____________________________________________________________________________

void  AliT0RawData::WriteTRMDataHeader(UInt_t slotID, Int_t nWordsInTRM,
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
  PackWord(baseWord,word,0, 3); // slotID
  word = nWordsInTRM;
 //+this word - DRM header 

  PackWord(baseWord,word,4, 16); // event words 
  word=0;
  PackWord(baseWord,word,17,18); // ACQ
  word=0;
  PackWord(baseWord,word,19,19); //  L SEY inside LUT
  word=0;
  PackWord(baseWord,word,20,27); //  MBZ
  word=4;
  PackWord(baseWord,word,28,31); // 0100 marks header
  fBuffer[positionOfTRMHeader] =  baseWord;
 cout<<" TRM header "<<baseWord<<endl;   
  word=0; 
  baseWord=0;
     
}

//_____________________________________________________________________________

void  AliT0RawData::WriteChainDataHeader(UInt_t chainNumber,UInt_t slotID)
{
//Write a (dummy or real) DDL Chain  data header, 
//set the compression bit if compressed
//  chainNumber 00 or 10
  UInt_t word;
  UInt_t baseWord=0;
  //fill TRM headers
  //TRM Global Header
  word = slotID; // ask Tatiana 7 or 9 
  PackWord(baseWord,word,0, 3); // slotID
  word = 0;
  PackWord(baseWord,word,4, 15); // bunchID
  word=0;
  PackWord(baseWord,word,16,23); // PB24 temperature
  word=0;
  PackWord(baseWord,word,24,26); //  PB24 ID
  word=0;
  PackWord(baseWord,word,27,27); //  TS
  word=chainNumber;
  PackWord(baseWord,word,28,31); // 0100 marks header
  fBuffer[4] =  baseWord;

  word=0;
  baseWord=0;     
  
}
//_____________________________________________________________________________

void  AliT0RawData::WriteDataHeader(Bool_t dummy, Bool_t compressed)
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


void  AliT0RawData::WriteTrailer(UInt_t slot, Int_t word1, UInt_t word2, UInt_t word3)
{

  UInt_t word;
  UInt_t baseWord=0;
  word = slot;
  PackWord(baseWord,word,0, 3); // 0001 
  word=word1;
  PackWord(baseWord,word,4, 15); // CRC ?
  word = word2;
  PackWord(baseWord,word,16,27); // event counter
  word=word3;
  PackWord(baseWord,word,28,31); //  marks trailer
  fIndex++;
  fBuffer[fIndex] =  baseWord;
  
  word=0;
  baseWord=0;

}

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
void  AliT0RawData::FillTime(Int_t ch, Int_t iTDC, Int_t time)
{
  UInt_t word;
  UInt_t baseWord=0;

  word=time;
  PackWord(baseWord,word, 0, 20); // Time 

  word=ch;
  PackWord(baseWord,word, 21, 23); // number of channel 
  word=iTDC;
  PackWord(baseWord,word, 24, 27); // TDC ID

  word=0;
  PackWord(baseWord,word, 28, 28); // E = 0 in simulation
  word=0;
  PackWord(baseWord,word, 29, 30); // PS bit data 00
  word=1;
  PackWord(baseWord,word, 31, 31); // 1
  fIndex++;
  fBuffer[fIndex]=baseWord;
  word=0;
  baseWord=0;
	
}
//---------------------------------------------------------------------------------------

Int_t AliT0RawData::RawDataT0(AliT0digit *fDigits)
{
   //This method creates the Raw data files for T0 detector


  // const Int_t kSize=512; //2*AliTOFGeometry::NpadXSector() 
                          //max number of digits per DDL file times 2
  //  UInt_t fBuffer[kSize];
  //  UInt_t baseWord;
  // UInt_t word;

  fIndex=-1;
 
  AliRawDataHeaderSim header;
  //loop over TOF DDL files
    //write Dummy DATA HEADER
   WriteDataHeader(kTRUE, kFALSE);
  GetDigits(fDigits);
  //write packing digits
  fFile->WriteBuffer((char*) fBuffer,((fIndex+1)*sizeof(UInt_t)));
  //write real data header on its place
   WriteDataHeader(kFALSE, kFALSE);
  
  
  //end for
  
  return 0;  
  
}
