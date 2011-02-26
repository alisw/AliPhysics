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

//#include <Riostream.h>
//#include <TTree.h>
#include <TMap.h>
#include "AliT0.h"
#include "AliT0RawData.h"
#include "AliT0digit.h"
#include "AliBitPacking.h"
#include "AliRawDataHeader.h"
#include "AliRawDataHeaderSim.h"
#include "AliBitPacking.h"
#include "AliFstream.h"
#include "AliRunLoader.h"
#include "AliDAQ.h"
#include "AliT0LookUpValue.h"
#include "AliT0LookUpKey.h"

ClassImp(AliT0RawData)

//_____________________________________________________________________________
  AliT0RawData::AliT0RawData():TObject(),
			       fVerbose(0),      
			       fIndex(-1) ,     
			       fEventNumber(0), 
			       fTimeCFD(new TArrayI(24)),    
			       fADC1( new TArrayI(24)),     
			       fTimeLED( new TArrayI(24)), 
			       fADC0( new TArrayI(24)),     
			       fFile(0x0),   
			       fDataHeaderPos(0),
			       fDRMDataHeaderPos(0),
			       fTRMDataHeaderPos(0),
			       fParam(0),
			       fLookUp(0)
  

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

  //open the output file
  char fileName[15];
  strcpy(fileName,AliDAQ::DdlFileName("T0",0)); //The name of the output file
  fFile = new AliFstream(fileName);
  memset(fBuffer,0,512*sizeof(UInt_t));

  //get event number 
  AliRunLoader *runloader = AliRunLoader::Instance();
  if (runloader) {
    fEventNumber = runloader->GetEventNumber();
  }

  // Inverse lookup table for simulation

  fParam = AliT0Parameters::Instance();
  fParam->Init();
  AliT0LookUpKey* lookkey;//= new AliT0LookUpKey();
  AliT0LookUpValue*  lookvalue;//= new AliT0LookUpValue();
  TMap *lookup = fParam->GetMapLookup();
  TMapIter iter(lookup);

  for( Int_t iline=0; iline<106; iline++)
    {
      lookvalue = ( AliT0LookUpValue*) iter.Next();
      lookkey = (AliT0LookUpKey*) lookup->GetValue(lookvalue);
      fLookUp.Add(lookkey, lookvalue);
      //lookkey= new AliT0LookUpKey();
      //lookvalue= new AliT0LookUpValue();
    }
    
}

//_____________________________________________________________________________

AliT0RawData::AliT0RawData(const AliT0RawData &r):TObject(),
						  fVerbose(0),      
						  fIndex(-1) ,     
						  fEventNumber(0), 
						  fTimeCFD(new TArrayI(24)),    
						  fADC1( new TArrayI(24)),     
						  fTimeLED( new TArrayI(24)), 
						  fADC0( new TArrayI(24)),     
						  fFile(0x0),   
						  fDataHeaderPos(0),
						  fDRMDataHeaderPos(0),
						  fTRMDataHeaderPos(0),
						  fParam(0),
						  fLookUp(0)

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
  delete fTimeCFD;
  delete fADC1;
  delete fTimeLED;
  delete fADC0;
  delete fFile;
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
//void AliT0RawData::GetDigits(fDigits)
{
 
  //This method packs the T0 digits in a proper 32 bits structure

  //read T0 digits and fill TDC and ADC arrays


  //  Int_t error=0;
  Int_t time,  positionOfTRMHeader;
  
  // Get the digits array
  
  fDigits->GetTimeCFD(*fTimeCFD);
  fDigits->GetQT0(*fADC1);
  fDigits->GetTimeLED(*fTimeLED);
  fDigits->GetQT1(*fADC0);
  Int_t meantime = fDigits->MeanTime(); 
  Int_t timediff = fDigits->TimeDiff(); 
  Int_t mult0=fDigits->SumMult();
  Int_t mult1=fDigits->SumMult();
  Int_t timeA = fDigits->BestTimeC();
  Int_t timeC = fDigits->BestTimeA();
  
  
  //  TArrayI  *allData = new TArrayI(110);
  Int_t allData[110][1];
  for (Int_t i=0; i<110; i++) allData[i][0] = 0;

  allData[0][0]=0;
  for (Int_t i=1; i<13; i++) {
    allData[i][0]    = fTimeCFD->At(i-1);
    allData[i+12][0] = fTimeLED->At(i-1);
    allData[i+56][0] = fTimeCFD->At(i-1+12);
    allData[i+68][0] = fTimeLED->At(i-1+12);
  }
  
  for (Int_t iii=0; iii<12; iii++) {
    allData[2*iii+25][0] = fADC1->At(iii);
    allData[2*iii+26][0] = fADC0->At(iii);
  }
  for (Int_t ii=12; ii<24; ii++) {
    allData[2*ii+57][0] = fADC1->At(ii);
    allData[2*ii+58][0] = fADC0->At(ii);
  }
  
  allData[49][0] = meantime;
  allData[50][0] = timediff;
  allData[51][0] = timeA;
  allData[52][0] = timeC;
  allData[53][0] = mult0;
  allData[54][0] = mult1;
  allData[55][0] = mult0;
  allData[56][0] = mult1;

  //    cout.setf( ios_base::hex, ios_base::basefield );
  //space for DRM header
  fIndex += 6;


  Int_t startTRM=fIndex;
  //space for 1st TRM header
  fIndex ++;
  positionOfTRMHeader= fIndex;
  //space for chain  header
  fIndex ++;
  WriteChainDataHeader(1, 1); // 

  //  fIndex++;
  // Loop through all PMT
  Int_t chain=0; 
  Int_t iTDC = 0;
  Int_t channel=0;
  Int_t trm1words=0;
  Int_t itrm=7;
  Int_t inside =0;
  Int_t isData = 0;
  AliT0LookUpKey * lookkey  = new AliT0LookUpKey();
  AliT0LookUpValue * lookvalue ;//= new AliT0LookUpValue(trm,tdc,chain,channel);
  for (Int_t det = 0; det < 105; det++) {
    time = allData[det][0];
    if (time >0 && time !=99999) {
      lookkey->SetKey(det);
      lookvalue = (AliT0LookUpValue*) fLookUp.GetValue((TObject*)lookkey);     
      if (lookvalue ) 
	{
	  isData++;
	  itrm= lookvalue->GetTRM();
	  if (det >56 &&inside == 0)  {
	    WriteChainDataTrailer(1); // 1st chain trailer
	    fIndex++;
	    WriteChainDataHeader(2, 1);
	    //	    fIndex++;
	    inside++;
	  }	    
	  chain = lookvalue->GetChain();
	  iTDC = lookvalue->GetTDC();
	  channel = lookvalue->GetChannel();
	  FillTime(channel,  iTDC,  time);
	  AliDebug(1,Form("look %i  itrm %i ,  chain %i , iTDC %i, channel %i time %i", det,itrm,chain,iTDC,channel, time));
	}
      else
	{
	  cout<<" no lookup value for key "<<det<<endl;
	  //  break;
	}
    }
    
  }
  if (inside==0) {
    WriteChainDataTrailer(1); // 1st chain trailer
    fIndex++;
    WriteChainDataHeader(2, 1);
  }
    //  WriteChainDataHeader(2, 1); // 
  WriteChainDataTrailer(3); // 2st chain trailer
  WriteTrailer(15,0,fEventNumber,5); // 1st TRM trailer
  
  
  trm1words = fIndex - startTRM;
  //space for 2st TRM header
  
  WriteTRMDataHeader(itrm, trm1words , positionOfTRMHeader);
  
  //DRM trailer
  WriteTrailer(1,0,fEventNumber,5);
    
    WriteDRMDataHeader();
    
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
  //cout<<" chain header "<<baseWord<<" number "<<chainNumber<<endl;
  word=0;
  baseWord=0;     
  
}
//_____________________________________________________________________________

void  AliT0RawData::WriteChainDataTrailer(UInt_t chainNumber )
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
  fIndex++;
  fBuffer[fIndex] =  baseWord;

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
  fIndex++;
  fBuffer[fIndex] =  baseWord;

  word=0;
  baseWord=0;

}

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
void  AliT0RawData::FillTime(Int_t ch, Int_t iTDC, Int_t time)
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
  fIndex++;
  fBuffer[fIndex]=baseWord;

  word=0;
  baseWord=0;
}
//---------------------------------------------------------------------------------------

Int_t AliT0RawData::RawDataT0(AliT0digit *fDigits)
  //Int_t AliT0RawData::RawDataT0(*fDigits)
{
   //This method creates the Raw data files for T0 detector


  // const Int_t kSize=512; //2*AliTOFGeometry::NpadXSector() 
                          //max number of digits per DDL file times 2
  //  UInt_t fBuffer[kSize];
  //  UInt_t baseWord;
  // UInt_t word;

  fIndex=-1;
 

   WriteDataHeader(kTRUE, kFALSE);
  GetDigits(fDigits);
  //write packing digits
  
  
  fFile->WriteBuffer((char*) fBuffer,((fIndex+1)*sizeof(UInt_t)));
  //write real data header on its place
   WriteDataHeader(kFALSE, kFALSE);
  
  
  //end for
  
  return 0;  
  
}
