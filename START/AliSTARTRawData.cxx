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
//  START raw data conversion class                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TTree.h>

#include "AliSTART.h"
#include "AliSTARTRawData.h"
#include "AliSTARTdigit.h"
#include "AliBitPacking.h"
#include "AliRawDataHeader.h"
#include "AliBitPacking.h"

ClassImp(AliSTARTRawData)

//_____________________________________________________________________________
AliSTARTRawData::AliSTARTRawData():TObject()
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
  fADC = new TArrayI(24);
  fTimeLED = new TArrayI(24);
  fADC0 = new TArrayI(24);
   //   this->Dump();
  
}

//_____________________________________________________________________________
AliSTARTRawData::AliSTARTRawData(const AliSTARTRawData &r):TObject()
{
  //
  // AliSTARTrawData copy constructor
  //

  ((AliSTARTRawData &) r).Copy(*this);

}

//_____________________________________________________________________________
AliSTARTRawData::~AliSTARTRawData()
{
  //
  // Destructor
  //
  if (fDigits) {
    delete fDigits;
    fDigits = NULL;
  }
  delete fTimeCFD;
  delete fADC;
  delete fTimeLED;
  delete fADC0;
}

//_____________________________________________________________________________
AliSTARTRawData &AliSTARTRawData::operator=(const AliSTARTRawData &r)
{
  //
  // Assignment operator
  //

  if (this != &r) ((AliSTARTRawData &) r).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliSTARTRawData::GetDigits(AliSTARTdigit *fDigits, UInt_t *buf)
{
 
  //This method packs the START digits in a proper 32 bits structure

  //read START digits and fill TDC and ADC arrays


  UInt_t word;
  UInt_t baseWord=0;
  Int_t error=0;

    // Get the digits array

  fDigits->GetTime(*fTimeCFD);
  fDigits->GetADC(*fADC);
  fDigits->GetTimeAmp(*fTimeLED);
  fDigits->GetADCAmp(*fADC0);
     
  // Loop through all PMT
 
  for (Int_t det = 0; det < 24; det++) {
    Int_t timeLED=fTimeLED->At(det);
    //   if( timeLED > 0){
      // DDL 1 0-5 -#PMT, 6-31 - empty
      //LED
      word=det;;
      PackWord(baseWord,word, 0, 8); 
      fIndex++;
      buf[fIndex]=baseWord;
      
      word=0;
      baseWord=0;
      word=error;
      PackWord(baseWord,word,0, 7); // Error flag
      word=timeLED;
      PackWord(baseWord,word,8,31); // time-of-flight
      fIndex++;
      buf[fIndex]=baseWord;
      word=0;
      baseWord=0;
      //  }
  }
   for (Int_t det = 0; det < 24; det++) {
    //CDF
     Int_t timeCFD=fTimeCFD->At(det);
     //   if ( timeCFD >0 ) {
       // DDL2 2 0-5 -#PMT, 6-31 - empty
       word=0;
       baseWord=0;
       word=det+24;
       PackWord(baseWord,word, 0, 8); // number of PMT on the right side
       fIndex++;
       buf[fIndex]=baseWord;
       word=0;
       baseWord=0;
       word=error;
       PackWord(baseWord,word,0, 7); // Error flag
       word=timeCFD;
       PackWord(baseWord,word,8,31); // time-of-flight
       fIndex++;
       buf[fIndex]=baseWord;
        word=0;
       baseWord=0;
       //  }
   }
  for (Int_t det = 0; det < 24; det++) {
    //conver ADC to time (preliminary algorithm)
    Int_t qtc=fADC->At(det);
    //    if ( qtc > 0 )
    //  {
	word=det+48;
	PackWord(baseWord,word, 0, 8); // number of PMT on the right side
	fIndex++;
	buf[fIndex]=baseWord;
	baseWord=0;
	word=error;
	PackWord(baseWord,word,0, 7); // Error flag
	word=qtc;
	PackWord(baseWord,word,8,31); // Q->Time
	fIndex++;
	buf[fIndex]=baseWord;
 	word=0;
	baseWord=0;
	//     }
  }
  
  for (Int_t det = 0; det < 24; det++) {
    Int_t qtcAmp=fADC0->At(det);
    
    // DDL 4 amplified QTC charge * 10
    
     //Amplified  ADC -> TDC 
    
    word=det+72;
    PackWord(baseWord,word, 0, 8); // number of PMT on the right side
    fIndex++;
    buf[fIndex]=baseWord;
    baseWord=0;
    word=error;
    PackWord(baseWord,word,0, 7); // Error flag
    word=qtcAmp;
    PackWord(baseWord,word,8,31); // Q->T amplified
    fIndex++;
    buf[fIndex]=baseWord;
    
    word=0;
    baseWord=0;
  }


  word=0;
  baseWord=0;
  fIndex++;
  word=97;
  PackWord(baseWord,word, 0, 8); // ?????????????????
  buf[fIndex]=baseWord;
  baseWord=0;
  word=error;
  PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->MeanTime();
  PackWord(baseWord,word,8,31); // MEANER
 
  fIndex++;
  buf[fIndex]=baseWord;



  // besttime right & left
  word=98;
  PackWord(baseWord,word, 0, 8); // T0-A sign
  fIndex++;
  buf[fIndex]=baseWord;

  baseWord=0;
  word=error;
  PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->BestTimeRight();
  PackWord(baseWord,word,8,31); // time-of-flight T0-A
  fIndex++;
  buf[fIndex]=baseWord;

  word=99;
  PackWord(baseWord,word, 0, 8); // T0-C sign
  fIndex++;
  buf[fIndex]=baseWord;

  baseWord=0;
  word=error;
  PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->BestTimeLeft();
  PackWord(baseWord,word,8,31); // time-of-flight T0-C 
  fIndex++;
  buf[fIndex]=baseWord;
 
  // time difference
  word=100;
  PackWord(baseWord,word, 0, 8); // TVDS sign
  fIndex++;
  buf[fIndex]=baseWord;

  baseWord=0;

  word=error;
  PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->TimeDiff();
  PackWord(baseWord,word,8,31); // T0verex
  fIndex++;
  buf[fIndex]=baseWord;


  // multiplicity 
  
  Int_t mult=fDigits->SumMult();
  word=101;
  PackWord(baseWord,word, 0, 8); 
  fIndex++;
  buf[fIndex]=baseWord;
  
  baseWord=0;
  word=error;
  PackWord(baseWord,word,0, 7); // Error flag
  word=mult;
  PackWord(baseWord,word,8,31); // time amplitude
  fIndex++;
  buf[fIndex]=baseWord;
  


  // trigger channels
   // besttime right & left
  word=102;
  PackWord(baseWord,word, 0, 8); // T0-A sign
  fIndex++;
  buf[fIndex]=baseWord;

  baseWord=0;
  word=error;
  PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->BestTimeRight();
  PackWord(baseWord,word,8,31); // time-of-flight T0-A
  fIndex++;
  buf[fIndex]=baseWord;
 
  word=103;
  PackWord(baseWord,word, 0, 8); // T0-C sign
  fIndex++;
  buf[fIndex]=baseWord;

  baseWord=0;

  word=error;
  PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->BestTimeLeft();
  PackWord(baseWord,word,8,31); // time-of-flight T0-C 
  fIndex++;
  buf[fIndex]=baseWord;

  // time difference
  word=104;
  PackWord(baseWord,word, 0, 8); // TVDS sign
  word=fDigits->TimeDiff();
   PackWord(baseWord,word, 6, 31); // T0 vertex 
  fIndex++;
  buf[fIndex]=baseWord;

  baseWord=0;

  word=error;
  PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->TimeDiff();
  PackWord(baseWord,word,8,31); // T0verex
  fIndex++;
  buf[fIndex]=baseWord;

  // multiplicity 
  
  mult=fDigits->SumMult();
  word=105;
  PackWord(baseWord,word, 0, 8); 
  word=mult;
  PackWord(baseWord,word, 6, 31); // sum amplitude
  fIndex++;
  buf[fIndex]=baseWord;
  
  baseWord=0;
  word=error;
 PackWord(baseWord,word,0, 7); // Error flag
  word=mult;
 PackWord(baseWord,word,8,31); // time amplitude
  fIndex++;
  buf[fIndex]=baseWord;
  
  // multiplicity 
  
   mult=fDigits->SumMult();
  word=106;
  PackWord(baseWord,word, 0, 8); 
  fIndex++;
  buf[fIndex]=baseWord;
  
  baseWord=0;
  word=error;
  PackWord(baseWord,word,0, 7); // Error flag
  word=mult;
  PackWord(baseWord,word,8,31); // time amplitude
  fIndex++;
  buf[fIndex]=baseWord;
  cout<<endl;
}

//-----------------------------------------------------------------------------------

void AliSTARTRawData::PackWord(UInt_t &BaseWord, UInt_t Word, Int_t StartBit, Int_t StopBit)
{
  //This method packs a word into the Baseword buffer starting form the "StartBit"
  //and tacking StopBit-StartBit+1 bits
  UInt_t dummyWord,offSet;
  Int_t  length;
  UInt_t sum;
  //The BaseWord is being filled with 1 from StartBit to StopBit
  length=StopBit-StartBit+1;
  sum=(UInt_t)TMath::Power(2,length)-1;

  if(Word > sum){
    Error("PackWord", "Word to be filled is not within desired length\n"
          "Word:%d Start bit:%d Stop Bit:%d",Word,StartBit,StopBit);
    return;
  }
  offSet=sum;
  offSet<<=StartBit;
  BaseWord=BaseWord|offSet;

  //The Word to be filled is shifted to the position StartBit
  //and the remaining  Left and Right bits are filled with 1
  sum=(UInt_t)TMath::Power(2,StartBit)-1;
  dummyWord=0xFFFFFFFF<<length;
  dummyWord +=Word;
  dummyWord<<=StartBit;
  dummyWord+=sum;
  BaseWord=BaseWord&dummyWord;
  return;
}

//---------------------------------------------------------------------------------------

Int_t AliSTARTRawData::RawDataSTART(AliSTARTdigit *fDigits)
{
   //This method creates the Raw data files for START detector


  const Int_t kSize=512; //2*AliTOFGeometry::NpadXSector() 
                          //max number of digits per DDL file times 2
  UInt_t buf[kSize];
  UInt_t baseWord;
  UInt_t word;

  fIndex=-1;

  char fileName[15];
  ofstream outfile;         // logical name of the output file 
  AliRawDataHeader header;
  AliBitPacking *pack ;
  //loop over TOF DDL files
  sprintf(fileName,"START_%d.ddl", 0xd00);
  //   sprintf(fileName,"START_0xd00.ddl"); //The name of the output file
#ifndef __DECCXX
    outfile.open(fileName,ios::binary);
#else
    outfile.open(fileName);
#endif
    //write Dummy DATA HEADER
    UInt_t dataHeaderPosition=outfile.tellp();
    outfile.write((char*)(&header),sizeof(header));

    baseWord=0;
    word=0;
    pack-> PackWord(baseWord,word,0, 31); // Number of DDL file

    fIndex++;
    buf[fIndex]=baseWord;
    GetDigits(fDigits,buf);

    outfile.write((char *)buf,((fIndex+1)*sizeof(UInt_t)));
    for(Int_t ii=0;ii<(fIndex+1);ii++) buf[ii]=0;
    fIndex=-1;
    
    //Write REAL DATA HEADER
    UInt_t currentFilePosition=outfile.tellp();
    outfile.seekp(dataHeaderPosition);
    header.fSize=currentFilePosition-dataHeaderPosition;
    header.SetAttribute(0);  // valid data
    outfile.write((char*)(&header),sizeof(header));
    outfile.close();

 //end for
   
  return 0;  
  
}
