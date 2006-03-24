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
  AliBitPacking *pack ;
    // Get the digits array

  fDigits->GetTime(*fTimeCFD);
  fDigits->GetADC(*fADC);
  fDigits->GetTimeAmp(*fTimeLED);
  fDigits->GetADCAmp(*fADC0);
     
  // Loop through all PMT
 
  for (Int_t det = 0; det < 24; det++) {
    Int_t timeLED=fTimeLED->At(det);

     // DDL 1 0-5 -#PMT, 6-31 - empty
    //LED
     word=det;;
     pack->PackWord(baseWord,word, 0, 5); 
     fIndex++;
     buf[fIndex]=baseWord;

     word=0;
     baseWord=0;
     word=error;
     pack->PackWord(baseWord,word,0, 7); // Error flag
     word=timeLED;
     cout<<det<<" led "<<timeLED<<" "<<word<<" ";
    pack->PackWord(baseWord,word,8,31); // time-of-flight
     fIndex++;
     buf[fIndex]=baseWord;
     word=0;
     baseWord=0;
  }
  cout<<endl;
   for (Int_t det = 0; det < 24; det++) {
    //CDF
    Int_t timeCFD=fTimeCFD->At(det);
    // DDL2 2 0-5 -#PMT, 6-31 - empty
     word=0;
     baseWord=0;
     word=det+24;
     pack->PackWord(baseWord,word, 0, 5); // number of PMT on the right side
     fIndex++;
     buf[fIndex]=baseWord;
     word=0;
     baseWord=0;
     word=error;
    pack-> PackWord(baseWord,word,0, 7); // Error flag
     word=timeCFD;
     cout<<det<<" CFD "<<timeCFD<<" "<<word<<" ";
     pack->PackWord(baseWord,word,8,31); // time-of-flight
     fIndex++;
     buf[fIndex]=baseWord;
     word=0;
     baseWord=0;
   }
   cout<<endl;
  for (Int_t det = 0; det < 24; det++) {
    //conver ADC to time (preliminary algorithm)
    Int_t qtc=fADC->At(det);
    word=det+48;
    pack->PackWord(baseWord,word, 0, 5); // number of PMT on the right side
    fIndex++;
    buf[fIndex]=baseWord;
    baseWord=0;
    word=error;
    pack->PackWord(baseWord,word,0, 7); // Error flag
    word=qtc;
    cout<<det<<" "<<fIndex<<" adc1 "<<word<<" ";
    pack->PackWord(baseWord,word,8,31); // Q->Time
    fIndex++;
    buf[fIndex]=baseWord;
    word=0;
    baseWord=0;
  }
  cout<<endl;
  
  for (Int_t det = 0; det < 24; det++) {
    Int_t qtcAmp=fADC0->At(det);
    
    // DDL 4 amplified QTC charge * 10
    
     //Amplified  ADC -> TDC 
    
    word=det+72;
    pack->PackWord(baseWord,word, 0, 5); // number of PMT on the right side
    fIndex++;
    buf[fIndex]=baseWord;
    baseWord=0;
    word=error;
    pack->PackWord(baseWord,word,0, 7); // Error flag
    word=qtcAmp;
    cout<<det<<" "<<fIndex<<" adc2 "<<word<<" ";
    pack->PackWord(baseWord,word,8,31); // Q->T amplified
    fIndex++;
    buf[fIndex]=baseWord;
    
    word=0;
    baseWord=0;
  }


  word=0;
  baseWord=0;
  fIndex++;
  word=97;
  pack->PackWord(baseWord,word, 0, 5); // ?????????????????
  buf[fIndex]=baseWord;
  baseWord=0;
  word=error;
  pack->PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->MeanTime();
  cout<<fIndex<<" "<<" mean "<< word<<endl;
  pack->PackWord(baseWord,word,8,31); // MEANER
  fIndex++;
  buf[fIndex]=baseWord;
    cout<<fIndex<<" meaner "<<word<<" ";



  // besttime right & left
  word=98;
 pack-> PackWord(baseWord,word, 0, 5); // T0-A sign
  fIndex++;
  buf[fIndex]=baseWord;

  baseWord=0;
  word=error;
  pack->PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->BestTimeRight();
  pack-> PackWord(baseWord,word,8,31); // time-of-flight T0-A
  fIndex++;
  buf[fIndex]=baseWord;
   cout<<fIndex<<" T0-A "<<word<<" ";

  word=99;
   pack->PackWord(baseWord,word, 0, 5); // T0-C sign
  fIndex++;
  buf[fIndex]=baseWord;

  baseWord=0;

  word=error;
  pack->PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->BestTimeLeft();
 pack-> PackWord(baseWord,word,8,31); // time-of-flight T0-C 
  fIndex++;
  buf[fIndex]=baseWord;

   cout<<fIndex<<" T0-C "<<word<<" ";
  // time difference
  word=100;
 pack-> PackWord(baseWord,word, 0, 5); // TVDS sign
  word=fDigits->TimeDiff();
 pack-> PackWord(baseWord,word, 6, 31); // T0 vertex 
  fIndex++;
  buf[fIndex]=baseWord;

  baseWord=0;

  word=error;
  pack->PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->TimeDiff();
  pack->PackWord(baseWord,word,8,31); // T0verex
  fIndex++;
  buf[fIndex]=baseWord;
  cout<<fIndex<<" diff "<<word<<endl;

  // multiplicity 
  
   Int_t mult=fDigits->SumMult();
  word=100;
 pack-> PackWord(baseWord,word, 0, 5); 
  word=mult;
 pack-> PackWord(baseWord,word, 6, 31); // sum amplitude
  fIndex++;
  buf[fIndex]=baseWord;
  
  baseWord=0;
  word=error;
 pack-> PackWord(baseWord,word,0, 7); // Error flag
  word=mult;
 pack-> PackWord(baseWord,word,8,31); // time amplitude
  fIndex++;
  buf[fIndex]=baseWord;
  
  cout<<endl;

  // trigger channels
   // besttime right & left
  word=101;
  pack->PackWord(baseWord,word, 0, 5); // T0-A sign
  fIndex++;
  buf[fIndex]=baseWord;

  baseWord=0;
  word=error;
 pack-> PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->BestTimeRight();
 pack-> PackWord(baseWord,word,8,31); // time-of-flight T0-A
  fIndex++;
  buf[fIndex]=baseWord;
   cout<<fIndex<<" T0-A "<<word<<" ";

  word=102;
 pack-> PackWord(baseWord,word, 0, 5); // T0-C sign
  fIndex++;
  buf[fIndex]=baseWord;

  baseWord=0;

  word=error;
 pack-> PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->BestTimeLeft();
 pack-> PackWord(baseWord,word,8,31); // time-of-flight T0-C 
  fIndex++;
  buf[fIndex]=baseWord;

   cout<<fIndex<<" T0-C "<<word<<" ";
  // time difference
  word=103;
 pack-> PackWord(baseWord,word, 0, 5); // TVDS sign
  word=fDigits->TimeDiff();
 pack-> PackWord(baseWord,word, 6, 31); // T0 vertex 
  fIndex++;
  buf[fIndex]=baseWord;

  baseWord=0;

  word=error;
 pack-> PackWord(baseWord,word,0, 7); // Error flag
  word=fDigits->TimeDiff();
 pack-> PackWord(baseWord,word,8,31); // T0verex
  fIndex++;
  buf[fIndex]=baseWord;
  cout<<fIndex<<" diff "<<word<<endl;

  // multiplicity 
  
   mult=fDigits->SumMult();
  word=104;
 pack-> PackWord(baseWord,word, 0, 5); 
  word=mult;
pack->  PackWord(baseWord,word, 6, 31); // sum amplitude
  fIndex++;
  buf[fIndex]=baseWord;
  
  baseWord=0;
  word=error;
 pack-> PackWord(baseWord,word,0, 7); // Error flag
  word=mult;
 pack-> PackWord(baseWord,word,8,31); // time amplitude
  fIndex++;
  buf[fIndex]=baseWord;
  
  // multiplicity 
  
   mult=fDigits->SumMult();
  word=105;
 pack-> PackWord(baseWord,word, 0, 5); 
  word=mult;
 pack-> PackWord(baseWord,word, 6, 31); // sum amplitude
  fIndex++;
  buf[fIndex]=baseWord;
  
  baseWord=0;
  word=error;
 pack-> PackWord(baseWord,word,0, 7); // Error flag
  word=mult;
 pack-> PackWord(baseWord,word,8,31); // time amplitude
  fIndex++;
  buf[fIndex]=baseWord;
  cout<<endl;
 

}

//-------------------------------------------------------------------------------------
/*
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
*/
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
