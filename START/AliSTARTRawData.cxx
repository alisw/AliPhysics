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

#include "AliSTART.h"
#include "AliSTARTRawData.h"
#include "AliSTARTdigit.h"
#include "AliSTARTLoader.h"

#include <AliLoader.h>
#include <AliRunLoader.h>
 
#include "/home/alla/AliRoot/verynew/RAW/AliRawDataHeader.h"
#include "/home/alla/AliRoot/verynew/RAW/AliRawData.h"

ClassImp(AliSTARTRawData)

//_____________________________________________________________________________
AliSTARTRawData::AliSTARTRawData():TObject()
{


  fIndex=-1;
  fDigits = NULL;
  
  ftimeTDC = new TArrayI(24); 
   fADC = new TArrayI(24); 
   
  this->Dump();
  
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
    delete ftimeTDC;
    delete fADC;
  }

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

  cout<<"GetDigits(AliSTARTdigit *fDigits, UInt_t *buf) "<<endl;
  UInt_t word;
  UInt_t baseWord=0;
  Int_t error=0;

    fDigits->Print();

   // Get the digits array

      fDigits->GetTime(*ftimeTDC);
      fDigits->GetADC(*fADC);
      fBestTimeRight=fDigits->GetBestTimeRight();
      fBestTimeLeft=fDigits->GetBestTimeLeft();
      fMeanTime = fDigits-> GetMeanTime();

     
  // Loop through all PMT
 

  for (Int_t det = 0; det < 24; det++) {
    Int_t time=ftimeTDC->At(det);
    Int_t ADC=fADC->At(det);
    printf(" right det %x time %x ADC %x \n",det,time,ADC);
    //conver ADC to time (preliminary algorithm)


     // DDL 1 0-5 -#PMT, 6-31 - empty

     word=det;;
     PackWord(baseWord,word, 0, 5); 
     fIndex++;
     buf[fIndex]=baseWord;

     word=0;
     baseWord=0;

     //TDC   
     word=error;
     PackWord(baseWord,word,0, 7); // Error flag
     word=time;
     PackWord(baseWord,word,8,31); // time-of-flight
     fIndex++;
     buf[fIndex]=baseWord;
     word=0;
     baseWord=0;
    
    // DDL2 2 0-5 -#PMT, 6-31 - empty
     word=det;;
     PackWord(baseWord,word, 0, 5); // number of PMT on the right side
     fIndex++;
     buf[fIndex]=baseWord;
     word=0;
     baseWord=0;

     //  amplified TDC    
     word=error;
     PackWord(baseWord,word,0, 7); // Error flag
     word=time;
     PackWord(baseWord,word,8,31); // time-of-flight
     fIndex++;
     buf[fIndex]=baseWord;
     word=0;
     baseWord=0;

     // DDL 3
     word=det;;
     PackWord(baseWord,word, 0, 5); // number of PMT on the right side
     fIndex++;
     buf[fIndex]=baseWord;
     word=0;
     baseWord=0;

     // ADC -> TDC     
     word=error;
     PackWord(baseWord,word,0, 7); // Error flag
     word=ADC;
     PackWord(baseWord,word,8,31); // time-of-flight
     fIndex++;
     buf[fIndex]=baseWord;

     //Amplified  ADC -> TDC 
     word=0;
     baseWord=0;

     word=det;;
     PackWord(baseWord,word, 0, 5); // number of PMT on the right side
     fIndex++;
     buf[fIndex]=baseWord;
     baseWord=0;
     word=error;
     PackWord(baseWord,word,0, 7); // Error flag
     word=ADC;
     PackWord(baseWord,word,8,31); // time-of-flight
     fIndex++;
     buf[fIndex]=baseWord;
     word=0;
     baseWord=0;


 }
  /*
  //timemean
  fIndex++;
  buf[fIndex]=baseWord;
  word=25;
  PackWord(baseWord,word, 0, 5); // number of PMT on the right side
  word=fMeanTime;
  PackWord(baseWord,word, 6, 31); // TDC on the right side from Marin
  printf("meantime buf[%i]=%x\n",fIndex,buf[fIndex]);
  
  fIndex++;
  buf[fIndex]=baseWord;
  
  baseWord=0;
  
  word=error;
  PackWord(baseWord,word,0, 7); // Error flag
  word=fMeanTime;
  PackWord(baseWord,word,8,31); // time-of-flight
  
  fIndex++;
  buf[fIndex]=baseWord;
  
  printf("meantime buf[%i]=%x\n",fIndex,buf[fIndex]);

     
  // besttime right & left
  fIndex++;
  cout<<" left "<<fBestTimeLeft<<" right "<<fBestTimeRight<<endl;
  buf[fIndex]=baseWord;
  word=26;
  PackWord(baseWord,word, 0, 5); // number of PMT on the right side
  word=fBestTimeRight;
  PackWord(baseWord,word, 6, 31); // TDC on the right side from Marin
  printf("best buf[%i]=%x\n",fIndex,buf[fIndex]);
  
  fIndex++;
  buf[fIndex]=baseWord;
  
  baseWord=0;
  
  word=error;
  PackWord(baseWord,word,0, 7); // Error flag
  word=fBestTimeRight;
  PackWord(baseWord,word,8,31); // time-of-flight
  
  fIndex++;
  buf[fIndex]=baseWord;
  
  printf("4 right buf[%i]=%x\n",fIndex,buf[fIndex]);
  
  word=27;
  PackWord(baseWord,word, 0, 5); // number of PMT on the right side
  word=fBestTimeLeft;
  PackWord(baseWord,word, 6, 31); // TDC on the right side from Marin
  printf("5 left buf[%i]=%x\n",fIndex,buf[fIndex]);
  
  fIndex++;
  buf[fIndex]=baseWord;
  
  baseWord=0;
  
  word=error;
  PackWord(baseWord,word,0, 7); // Error flag
  word=fBestTimeLeft;
  PackWord(baseWord,word,8,31); // time-of-flight
  
  fIndex++;
  buf[fIndex]=baseWord;
  
  printf("5 left buf[%i]=%x\n",fIndex,buf[fIndex]);
  */   
  word=0;
  baseWord=0;
   
}

//-------------------------------------------------------------------------------------

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

Int_t AliSTARTRawData::RawDataSTART(AliSTARTdigit *fDigits){
  
  //This method creates the Raw data files for TOF detector
  const Int_t kSize=512; //2*AliTOFGeometry::NpadXSector() 
                          //max number of digits per DDL file times 2
  UInt_t buf[kSize];
  UInt_t baseWord;
  UInt_t word;

  fIndex=-1;

  char fileName[15];
  ofstream outfile;         // logical name of the output file 
  AliRawDataHeader header;
  cout<<" AliRawDataHeader header; start "<<endl;
  //loop over TOF DDL files
     sprintf(fileName,"START_0xd00.ddl"); //The name of the output file
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
    PackWord(baseWord,word,0, 31); // Number of DDL file

    fIndex++;
    buf[fIndex]=baseWord;

    //   branch->GetEvent();

    //For each DDL file, buf contains the array of data words in Binary format
    //fIndex gives the number of 32 bits words in the buffer for each DDL file
    cout<<" AliSTARTRawData::RawDataSTART "<<fDigits<<endl;

    GetDigits(fDigits,buf);
    cout<<"REAL DATA "<<fIndex<<endl;
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
